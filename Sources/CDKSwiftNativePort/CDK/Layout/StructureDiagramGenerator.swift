import Foundation

#if canImport(CoreGraphics)
    import CoreGraphics
#endif

enum CDKStructureDiagramGenerator {
    // Phase-2 tuned constants (benchmark set: norbornane, decalin, naphthalene, steroid-like fused core).
    enum Tuning {
        static let anchorDriftWeight: CGFloat = 260
        static let hardOverlapRatio: CGFloat = 0.76
        static let softOverlapRatio: CGFloat = 1.12
        static let hardOverlapPenalty: CGFloat = 300
        static let softOverlapPenalty: CGFloat = 70
        static let intraRingHardRatio: CGFloat = 0.86
        static let intraRingPenalty: CGFloat = 115
        static let edgeCrossPenalty: CGFloat = 280
        static let edgeNearRatio: CGFloat = 0.42
        static let edgeNearPenalty: CGFloat = 34

        static let relaxIterations: Int = 150
        static let bondSpringFactor: CGFloat = 0.40
        static let nonBondedMinDistanceRatio: CGFloat = 1.14
        static let nonBondedPushFactor: CGFloat = 0.28
        static let crossingPushRatio: CGFloat = 0.15

        static let chainPreferredAngle: CGFloat = 120.0 * .pi / 180.0
        static let branchFanSpread: CGFloat = .pi / 2.2
        static let branchOpenSpread: CGFloat = .pi * 1.0
        static let chainPassLimit: Int = 32
        static let bondFlipGainThreshold: CGFloat = 0.97
    }

    static func apply(to molecule: Molecule) -> Molecule {
        guard molecule.atoms.count >= 2 else { return molecule }

        let standardBond: CGFloat = 1.4
        let graph = CDKMolecularGraph(molecule: molecule)
        let connectedComponents = graph.connectedComponents()
        let isMarkush =
            connectedComponents.count > 1
            && connectedComponents.contains { componentRGroupMembership($0, molecule: molecule) != nil }
        let components = connectedComponents.sorted {
            compareComponents($0, $1, molecule: molecule, markush: isMarkush)
        }
        let sssr = CDKRingSearch.sssrLikeRings(in: molecule, graph: graph)

        var positions: [Int: CGPoint] = [:]
        var offsetX: CGFloat = 0

        for component in components {
            let ringsInComponent =
                sssr
                .filter { !Set($0).isDisjoint(with: component) }
                .sorted { lhs, rhs in
                    if lhs.count != rhs.count { return lhs.count < rhs.count }
                    return lhs.lexicographicallyPrecedes(rhs)
                }
            let sharedEdgeMultiplicity = ringEdgeMultiplicity(ringsInComponent)
            let ringSystems = connectedRingSystems(ringsInComponent)
                .sorted { lhs, rhs in
                    let lFused = fusedEdgeCount(lhs, multiplicity: sharedEdgeMultiplicity)
                    let rFused = fusedEdgeCount(rhs, multiplicity: sharedEdgeMultiplicity)
                    if lFused != rFused { return lFused > rFused }
                    let lCount = Set(lhs.flatMap { $0 }).count
                    let rCount = Set(rhs.flatMap { $0 }).count
                    if lCount != rCount { return lCount > rCount }
                    return lhs.count > rhs.count
                }
            let rigidRingAtoms = simpleRingLockAtoms(ringSystems)

            if let primaryRingSystem = ringSystems.first {
                placeRingSystem(
                    primaryRingSystem,
                    graph: graph,
                    component: component,
                    positions: &positions,
                    center: CGPoint(x: offsetX, y: 0),
                    bondLength: standardBond)
            }

            let ringAtoms = Set(ringsInComponent.flatMap { $0 })
            let ringEdges = Set(ringsInComponent.flatMap { CDKRingSearch.edgeKeys(for: $0) })

            if !component.contains(where: { positions[$0] != nil }) {
                let seed = chooseSeed(
                    in: component,
                    graph: graph,
                    ringSet: Set(ringsInComponent.flatMap { $0 }),
                    preferAcyclicTerminalSeed: !isMarkush)
                positions[seed] = CGPoint(x: offsetX, y: 0)
            }

            var rounds = 0
            var progressed = true
            while progressed && rounds < Tuning.chainPassLimit {
                progressed = false
                if placeLongestUnplacedChains(
                    component: component,
                    graph: graph,
                    molecule: molecule,
                    ringAtoms: ringAtoms,
                    ringEdges: ringEdges,
                    positions: &positions,
                    bondLength: standardBond)
                {
                    progressed = true
                }
                if placeDistributedPartners(
                    component: component,
                    graph: graph,
                    molecule: molecule,
                    ringsInComponent: ringsInComponent,
                    positions: &positions,
                    bondLength: standardBond,
                    fallbackCenter: CGPoint(x: offsetX, y: 0))
                {
                    progressed = true
                }
                rounds += 1
            }

            // Second ring-placement pass after chain growth, useful for distal/fused systems.
            for system in ringSystems
            where system.contains(where: { ring in ring.contains(where: { positions[$0] == nil }) }) {
                let anchorCandidates = Set(system.flatMap { $0 }).compactMap { positions[$0] }
                let anchor =
                    anchorCandidates.isEmpty
                    ? CGPoint(x: offsetX, y: 0)
                    : CGPoint(
                        x: anchorCandidates.reduce(0) { $0 + $1.x } / CGFloat(anchorCandidates.count),
                        y: anchorCandidates.reduce(0) { $0 + $1.y } / CGFloat(anchorCandidates.count))
                placeRingSystem(
                    system,
                    graph: graph,
                    component: component,
                    positions: &positions,
                    center: anchor,
                    bondLength: standardBond)
            }

            _ = placeDistributedPartners(
                component: component,
                graph: graph,
                molecule: molecule,
                ringsInComponent: ringsInComponent,
                positions: &positions,
                bondLength: standardBond,
                fallbackCenter: CGPoint(x: offsetX, y: 0))

            // Handle any leftover atoms in this component.
            for atomID in component where positions[atomID] == nil {
                let neighbors = graph.neighbors(of: atomID).filter { component.contains($0) }
                let placed = neighbors.filter { positions[$0] != nil }
                if let anchor = placed.first, let anchorPos = positions[anchor] {
                    let angle = CGFloat((atomID * 37) % 360) * .pi / 180.0
                    positions[atomID] = CGPoint(
                        x: anchorPos.x + cos(angle) * standardBond,
                        y: anchorPos.y + sin(angle) * standardBond
                    )
                } else {
                    positions[atomID] = CGPoint(x: offsetX, y: 0)
                }
            }

            // Keep aromatic systems rigid while allowing aliphatic fused/bridged systems
            // to relax and untangle if they start with a rough placement.
            let lockedAtoms = Set(component.filter { isAromaticAtom($0, molecule: molecule) })
                .union(rigidRingAtoms)
            optimizeByBondFlips(
                component: component,
                graph: graph,
                positions: &positions,
                locked: lockedAtoms,
                ringEdges: ringEdges,
                bondLength: standardBond)
            relax(
                component: component,
                graph: graph,
                positions: &positions,
                bondLength: standardBond,
                locked: lockedAtoms,
                iterations: Tuning.relaxIterations)
            if isMarkush {
                let usedAttachmentOrientation = orientMarkushAttachmentComponent(
                    component: component,
                    graph: graph,
                    molecule: molecule,
                    positions: &positions)
                if !usedAttachmentOrientation {
                    selectOrientation(
                        component: component,
                        graph: graph,
                        molecule: molecule,
                        positions: &positions,
                        widthDiff: standardBond,
                        alignDiff: 1)
                }
                orientMarkushRootComponent(
                    component: component,
                    molecule: molecule,
                    positions: &positions)
            } else if ringsInComponent.isEmpty {
                selectOrientation(
                    component: component,
                    graph: graph,
                    molecule: molecule,
                    positions: &positions,
                    widthDiff: standardBond,
                    alignDiff: 1)
            }

            if let box = boundingBox(component: component, positions: positions) {
                let dx = offsetX - box.minX
                let dy = -box.midY
                for atomID in component {
                    if let p = positions[atomID] {
                        positions[atomID] = CGPoint(x: p.x + dx, y: p.y + dy)
                    }
                }
                if let placedBox = boundingBox(component: component, positions: positions) {
                    offsetX += max(6.0, placedBox.width + 4.0)
                } else {
                    offsetX += 6.0
                }
            } else {
                offsetX += 6.0
            }
        }

        arrangeCoordinateBondedSandwichRings(
            molecule: molecule,
            rings: sssr,
            positions: &positions,
            bondLength: standardBond)

        if isMarkush {
            arrangeMarkushComponents(
                components: components,
                molecule: molecule,
                positions: &positions,
                spacing: standardBond)
        }

        var out = molecule
        out.atoms = molecule.atoms.map { atom in
            var a = atom
            if let p = positions[a.id] { a.position = p }
            return a
        }
        return out
    }
}
