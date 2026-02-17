import Foundation
import CoreGraphics

// CDK-inspired SDG port entry point and layout support types.
private struct CDKEdgeKey: Hashable {
    let a: Int
    let b: Int

    init(_ u: Int, _ v: Int) {
        a = min(u, v)
        b = max(u, v)
    }
}

private struct CDKMolecularGraph {
    let atomIDs: [Int]
    let neighborsByAtom: [Int: [Int]]
    let edgeBonds: [CDKEdgeKey: Bond]

    init(molecule: Molecule) {
        atomIDs = molecule.atoms.map(\.id).sorted()

        var neighborSets: [Int: Set<Int>] = [:]
        var edgeBonds: [CDKEdgeKey: Bond] = [:]
        for id in atomIDs {
            neighborSets[id] = []
        }
        for bond in molecule.bonds {
            neighborSets[bond.a1, default: []].insert(bond.a2)
            neighborSets[bond.a2, default: []].insert(bond.a1)
            edgeBonds[CDKEdgeKey(bond.a1, bond.a2)] = bond
        }

        self.neighborsByAtom = neighborSets.mapValues { Array($0).sorted() }
        self.edgeBonds = edgeBonds
    }

    var edgeCount: Int { edgeBonds.count }

    func neighbors(of atomID: Int) -> [Int] {
        neighborsByAtom[atomID] ?? []
    }

    func connectedComponents() -> [Set<Int>] {
        var out: [Set<Int>] = []
        var seen: Set<Int> = []

        for seed in atomIDs where !seen.contains(seed) {
            var component: Set<Int> = [seed]
            var stack: [Int] = [seed]
            seen.insert(seed)

            while let cur = stack.popLast() {
                for nxt in neighbors(of: cur) where !seen.contains(nxt) {
                    seen.insert(nxt)
                    component.insert(nxt)
                    stack.append(nxt)
                }
            }
            out.append(component)
        }

        return out
    }
}

private enum CDKRingSearch {
    static func sssrLikeRings(in molecule: Molecule, graph: CDKMolecularGraph) -> [[Int]] {
        let components = graph.connectedComponents().count
        let rank = max(0, graph.edgeCount - graph.atomIDs.count + components)
        guard rank > 0 else { return [] }

        // Approximate CDK ring-search by first generating shortest edge-based cycle candidates
        // (basis-like), then filling with deterministic simple-cycle fallbacks if needed.
        let basisCandidates = shortestCycleCandidates(graph: graph, maxCycleSize: 12)
        let fallbackCandidates = molecule.simpleCycles(maxSize: 12)
        let candidates = deduplicateCycles(basisCandidates + fallbackCandidates)
            .sorted { lhs, rhs in
                if lhs.count != rhs.count { return lhs.count < rhs.count }
                return lhs.lexicographicallyPrecedes(rhs)
            }
        guard !candidates.isEmpty else { return [] }

        var selected: [[Int]] = []
        var coveredEdges: Set<CDKEdgeKey> = []

        for ring in candidates {
            let edges = edgeKeys(for: ring)
            if edges.contains(where: { !coveredEdges.contains($0) }) {
                selected.append(ring)
                for e in edges { coveredEdges.insert(e) }
            }
            if selected.count >= rank { break }
        }

        // If rank wasn't reached due edge-overlap pruning, fill deterministically.
        if selected.count < rank {
            for ring in candidates where !selected.contains(where: { $0 == ring }) {
                selected.append(ring)
                if selected.count >= rank { break }
            }
        }

        return Array(selected.prefix(rank))
    }

    static func edgeKeys(for ring: [Int]) -> [CDKEdgeKey] {
        guard ring.count >= 2 else { return [] }
        return (0..<ring.count).map { i in
            CDKEdgeKey(ring[i], ring[(i + 1) % ring.count])
        }
    }

    private static func shortestCycleCandidates(graph: CDKMolecularGraph, maxCycleSize: Int) -> [[Int]] {
        let edges = Array(graph.edgeBonds.keys).sorted { lhs, rhs in
            if lhs.a != rhs.a { return lhs.a < rhs.a }
            return lhs.b < rhs.b
        }

        var out: [[Int]] = []
        for edge in edges {
            guard let path = shortestPathExcludingEdge(graph: graph,
                                                       start: edge.a,
                                                       end: edge.b,
                                                       excluding: edge,
                                                       maxDepth: maxCycleSize - 1),
                  path.count >= 3,
                  path.count <= maxCycleSize else { continue }
            out.append(canonicalCycle(path))
        }
        return out
    }

    private static func shortestPathExcludingEdge(graph: CDKMolecularGraph,
                                                  start: Int,
                                                  end: Int,
                                                  excluding: CDKEdgeKey,
                                                  maxDepth: Int) -> [Int]? {
        var queue: [Int] = [start]
        var head = 0
        var seen: Set<Int> = [start]
        var parent: [Int: Int] = [:]
        var depth: [Int: Int] = [start: 0]

        while head < queue.count {
            let cur = queue[head]
            head += 1
            let curDepth = depth[cur] ?? 0
            if cur == end { break }
            if curDepth >= maxDepth { continue }

            for nxt in graph.neighbors(of: cur) {
                if CDKEdgeKey(cur, nxt) == excluding { continue }
                if seen.contains(nxt) { continue }
                seen.insert(nxt)
                parent[nxt] = cur
                depth[nxt] = curDepth + 1
                queue.append(nxt)
            }
        }

        guard seen.contains(end) else { return nil }
        var path: [Int] = [end]
        var cur = end
        while cur != start {
            guard let p = parent[cur] else { return nil }
            cur = p
            path.append(cur)
        }
        return Array(path.reversed())
    }

    private static func deduplicateCycles(_ cycles: [[Int]]) -> [[Int]] {
        var seen: Set<[Int]> = []
        var out: [[Int]] = []
        for cycle in cycles {
            let canonical = canonicalCycle(cycle)
            guard seen.insert(canonical).inserted else { continue }
            out.append(canonical)
        }
        return out
    }

    private static func canonicalCycle(_ cycle: [Int]) -> [Int] {
        guard !cycle.isEmpty else { return [] }
        let n = cycle.count

        func rotations(of seq: [Int]) -> [[Int]] {
            (0..<n).map { i in
                Array(seq[i..<n]) + Array(seq[0..<i])
            }
        }

        let forward = rotations(of: cycle)
        let reverse = rotations(of: Array(cycle.reversed()))
        return (forward + reverse).min(by: isLexicographicallySmaller) ?? cycle
    }

    private static func isLexicographicallySmaller(_ lhs: [Int], _ rhs: [Int]) -> Bool {
        for (a, b) in zip(lhs, rhs) {
            if a != b { return a < b }
        }
        return lhs.count < rhs.count
    }
}

enum CDKStructureDiagramGenerator {
    // Phase-2 tuned constants (benchmark set: norbornane, decalin, naphthalene, steroid-like fused core).
    private enum Tuning {
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
        let components = graph.connectedComponents()
            .sorted { ($0.min() ?? 0) < ($1.min() ?? 0) }
        let sssr = CDKRingSearch.sssrLikeRings(in: molecule, graph: graph)

        var positions: [Int: CGPoint] = [:]
        var offsetX: CGFloat = 0

        for component in components {
            let ringsInComponent = sssr
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

            if let primaryRingSystem = ringSystems.first {
                placeRingSystem(primaryRingSystem,
                                graph: graph,
                                component: component,
                                positions: &positions,
                                center: CGPoint(x: offsetX, y: 0),
                                bondLength: standardBond)
            }

            let ringAtoms = Set(ringsInComponent.flatMap { $0 })
            let ringEdges = Set(ringsInComponent.flatMap { CDKRingSearch.edgeKeys(for: $0) })

            if !component.contains(where: { positions[$0] != nil }) {
                let seed = chooseSeed(in: component, graph: graph, ringSet: Set(ringsInComponent.flatMap { $0 }))
                positions[seed] = CGPoint(x: offsetX, y: 0)
            }

            var rounds = 0
            var progressed = true
            while progressed && rounds < Tuning.chainPassLimit {
                progressed = false
                if placeLongestUnplacedChains(component: component,
                                              graph: graph,
                                              molecule: molecule,
                                              ringAtoms: ringAtoms,
                                              ringEdges: ringEdges,
                                              positions: &positions,
                                              bondLength: standardBond) {
                    progressed = true
                }
                if placeDistributedPartners(component: component,
                                           graph: graph,
                                           molecule: molecule,
                                           ringsInComponent: ringsInComponent,
                                           positions: &positions,
                                           bondLength: standardBond,
                                           fallbackCenter: CGPoint(x: offsetX, y: 0)) {
                    progressed = true
                }
                rounds += 1
            }

            // Second ring-placement pass after chain growth, useful for distal/fused systems.
            for system in ringSystems where system.contains(where: { ring in ring.contains(where: { positions[$0] == nil }) }) {
                let anchorCandidates = Set(system.flatMap { $0 }).compactMap { positions[$0] }
                let anchor = anchorCandidates.isEmpty
                    ? CGPoint(x: offsetX, y: 0)
                    : CGPoint(x: anchorCandidates.reduce(0) { $0 + $1.x } / CGFloat(anchorCandidates.count),
                              y: anchorCandidates.reduce(0) { $0 + $1.y } / CGFloat(anchorCandidates.count))
                placeRingSystem(system,
                                graph: graph,
                                component: component,
                                positions: &positions,
                                center: anchor,
                                bondLength: standardBond)
            }

            _ = placeDistributedPartners(component: component,
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
            optimizeByBondFlips(component: component,
                                graph: graph,
                                positions: &positions,
                                locked: lockedAtoms,
                                ringEdges: ringEdges,
                                bondLength: standardBond)
            relax(component: component,
                  graph: graph,
                  positions: &positions,
                  bondLength: standardBond,
                  locked: lockedAtoms,
                  iterations: Tuning.relaxIterations)

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

        var out = molecule
        out.atoms = molecule.atoms.map { atom in
            var a = atom
            if let p = positions[a.id] { a.position = p }
            return a
        }
        return out
    }

    private static func connectedRingSystems(_ rings: [[Int]]) -> [[[Int]]] {
        guard !rings.isEmpty else { return [] }
        let ringSets = rings.map(Set.init)
        var seen: Set<Int> = []
        var out: [[[Int]]] = []

        for i in rings.indices where !seen.contains(i) {
            seen.insert(i)
            var stack = [i]
            var group: [Int] = [i]

            while let cur = stack.popLast() {
                for j in rings.indices where !seen.contains(j) {
                    if !ringSets[cur].isDisjoint(with: ringSets[j]) {
                        seen.insert(j)
                        stack.append(j)
                        group.append(j)
                    }
                }
            }

            out.append(group.map { rings[$0] })
        }
        return out
    }

    private static func ringEdgeMultiplicity(_ rings: [[Int]]) -> [CDKEdgeKey: Int] {
        var out: [CDKEdgeKey: Int] = [:]
        for ring in rings {
            for edge in CDKRingSearch.edgeKeys(for: ring) {
                out[edge, default: 0] += 1
            }
        }
        return out
    }

    private static func fusedEdgeCount(_ ringSystem: [[Int]],
                                       multiplicity: [CDKEdgeKey: Int]) -> Int {
        let edges = Set(ringSystem.flatMap { CDKRingSearch.edgeKeys(for: $0) })
        return edges.reduce(into: 0) { acc, edge in
            if (multiplicity[edge] ?? 0) > 1 { acc += 1 }
        }
    }

    // CDK RingPlacer-like attachment sequencing:
    // fused -> spiro -> bridged, traversed from already placed ring order.
    private enum RingAttachmentMode: Int {
        case fused = 0
        case spiro = 1
        case bridged = 2
        case isolated = 3
    }

    private static func ringAttachmentMode(_ ring: [Int], to placedRing: [Int]) -> RingAttachmentMode {
        let shared = Set(ring).intersection(placedRing)
        if shared.count >= 2 {
            let ringEdges = Set(CDKRingSearch.edgeKeys(for: ring))
            let placedEdges = Set(CDKRingSearch.edgeKeys(for: placedRing))
            return ringEdges.isDisjoint(with: placedEdges) ? .bridged : .fused
        }
        if shared.count == 1 {
            return .spiro
        }
        return .isolated
    }

    private static func placeRingSystem(_ rings: [[Int]],
                                        graph: CDKMolecularGraph,
                                        component: Set<Int>,
                                        positions: inout [Int: CGPoint],
                                        center: CGPoint,
                                        bondLength: CGFloat) {
        guard !rings.isEmpty else { return }
        let ordered = rings.sorted { lhs, rhs in
            if lhs.count != rhs.count { return lhs.count < rhs.count }
            return lhs.lexicographicallyPrecedes(rhs)
        }

        // Seed with a ring that already has most anchors, otherwise place the first ring regularly.
        let seedIndex = ordered.indices.max { a, b in
            let ca = ordered[a].filter { positions[$0] != nil }.count
            let cb = ordered[b].filter { positions[$0] != nil }.count
            if ca != cb { return ca < cb }
            return ordered[a].count > ordered[b].count
        } ?? 0

        if ordered[seedIndex].allSatisfy({ positions[$0] == nil }) {
            placeRegularRing(ordered[seedIndex],
                             center: center,
                             bondLength: bondLength,
                             positions: &positions)
        }

        var placedRings: Set<Int> = []
        var placementOrder: [Int] = []
        if ordered[seedIndex].contains(where: { positions[$0] != nil }) {
            placedRings.insert(seedIndex)
            placementOrder.append(seedIndex)
        }

        var progress = true
        while progress {
            progress = false
            var candidates: [(idx: Int, mode: RingAttachmentMode, parentRank: Int, shared: Int)] = []

            for idx in ordered.indices where !placedRings.contains(idx) {
                let ring = ordered[idx]
                let sharedAnchors = ring.filter { positions[$0] != nil }
                guard !sharedAnchors.isEmpty else { continue }

                var bestMode: RingAttachmentMode = .isolated
                var bestParentRank = Int.max
                for (rank, parentIdx) in placementOrder.enumerated() {
                    guard placedRings.contains(parentIdx) else { continue }
                    let mode = ringAttachmentMode(ring, to: ordered[parentIdx])
                    if mode.rawValue < bestMode.rawValue || (mode == bestMode && rank < bestParentRank) {
                        bestMode = mode
                        bestParentRank = rank
                    }
                }
                if bestMode == .isolated {
                    bestMode = .bridged
                    bestParentRank = Int.max / 2
                }

                candidates.append((idx, bestMode, bestParentRank, sharedAnchors.count))
            }

            candidates.sort { lhs, rhs in
                if lhs.mode != rhs.mode { return lhs.mode.rawValue < rhs.mode.rawValue }
                if lhs.parentRank != rhs.parentRank { return lhs.parentRank < rhs.parentRank }
                if lhs.shared != rhs.shared { return lhs.shared > rhs.shared }
                let lRing = ordered[lhs.idx]
                let rRing = ordered[rhs.idx]
                if lRing.count != rRing.count { return lRing.count < rRing.count }
                return lRing.lexicographicallyPrecedes(rRing)
            }

            for candidateDesc in candidates {
                let idx = candidateDesc.idx
                if placedRings.contains(idx) { continue }
                let ring = ordered[idx]
                let shared = ring.filter { positions[$0] != nil }
                guard !shared.isEmpty else { continue }

                if let candidate = bestRingPlacementCandidate(ring: ring,
                                                              sharedAnchors: shared,
                                                              graph: graph,
                                                              component: component,
                                                              positions: positions,
                                                              bondLength: bondLength) {
                    var didPlace = false
                    for atomID in ring where positions[atomID] == nil {
                        if let p = candidate[atomID] {
                            positions[atomID] = p
                            didPlace = true
                        }
                    }
                    if didPlace || ring.allSatisfy({ positions[$0] != nil }) {
                        placedRings.insert(idx)
                        placementOrder.append(idx)
                        progress = true
                    }
                }
            }
        }
    }

    private static func placeDistributedPartners(component: Set<Int>,
                                                 graph: CDKMolecularGraph,
                                                 molecule: Molecule,
                                                 ringsInComponent: [[Int]],
                                                 positions: inout [Int: CGPoint],
                                                 bondLength: CGFloat,
                                                 fallbackCenter: CGPoint) -> Bool {
        var progressedAny = false
        var localProgress = true
        var guardPasses = 0

        while localProgress && guardPasses < max(4, component.count * 2) {
            localProgress = false
            guardPasses += 1

            for center in component.sorted() where positions[center] != nil {
                let neighbors = graph.neighbors(of: center).filter { component.contains($0) }
                let placed = neighbors.filter { positions[$0] != nil }
                var unplaced = neighbors.filter { positions[$0] == nil }
                guard !unplaced.isEmpty else { continue }

                // Try ring placement before branching to preserve ring geometry.
                let candidateRings = ringsInComponent.filter { ring in
                    ring.contains(center) && ring.contains(where: { positions[$0] == nil })
                }
                for ring in candidateRings {
                    placeRingSystem([ring],
                                    graph: graph,
                                    component: component,
                                    positions: &positions,
                                    center: positions[center] ?? fallbackCenter,
                                    bondLength: bondLength)
                }

                unplaced = neighbors.filter { positions[$0] == nil }
                guard !unplaced.isEmpty else { continue }

                let dirs = proposedDirections(center: center,
                                              placedNeighbors: placed,
                                              unplacedCount: unplaced.count,
                                              positions: positions,
                                              molecule: molecule,
                                              graph: graph)
                for (idxUnplaced, atomID) in unplaced.enumerated() {
                    guard let centerPos = positions[center] else { continue }
                    let dir = dirs[min(idxUnplaced, dirs.count - 1)]
                    positions[atomID] = CGPoint(
                        x: centerPos.x + dir.dx * bondLength,
                        y: centerPos.y + dir.dy * bondLength
                    )
                    localProgress = true
                    progressedAny = true
                }
            }
        }

        return progressedAny
    }

    private static func placeLongestUnplacedChains(component: Set<Int>,
                                                   graph: CDKMolecularGraph,
                                                   molecule: Molecule,
                                                   ringAtoms: Set<Int>,
                                                   ringEdges: Set<CDKEdgeKey>,
                                                   positions: inout [Int: CGPoint],
                                                   bondLength: CGFloat) -> Bool {
        var anyPlaced = false
        var passes = 0

        while passes < Tuning.chainPassLimit {
            guard let chain = bestUnplacedChain(component: component,
                                                graph: graph,
                                                molecule: molecule,
                                                ringAtoms: ringAtoms,
                                                ringEdges: ringEdges,
                                                positions: positions),
                  chain.count >= 2 else { break }

            let anchor = chain[0]
            let first = chain[1]
            let initial = initialChainVector(anchor: anchor,
                                             firstUnplaced: first,
                                             component: component,
                                             graph: graph,
                                             positions: positions)
            placeLinearChain(chain,
                             initialVector: initial,
                             component: component,
                             positions: &positions,
                             bondLength: bondLength)
            anyPlaced = true
            passes += 1
        }

        return anyPlaced
    }

    private static func bestUnplacedChain(component: Set<Int>,
                                          graph: CDKMolecularGraph,
                                          molecule: Molecule,
                                          ringAtoms: Set<Int>,
                                          ringEdges: Set<CDKEdgeKey>,
                                          positions: [Int: CGPoint]) -> [Int]? {
        var best: [Int] = []

        for anchor in component.sorted() where positions[anchor] != nil {
            let starts = graph.neighbors(of: anchor)
                .filter { component.contains($0) && positions[$0] == nil }
                .sorted()

            for start in starts {
                let chain = longestUnplacedChain(anchor: anchor,
                                                 start: start,
                                                 component: component,
                                                 graph: graph,
                                                 molecule: molecule,
                                                 ringAtoms: ringAtoms,
                                                 ringEdges: ringEdges,
                                                 positions: positions)
                if chain.count > best.count ||
                    (chain.count == best.count && pathLexicographicallyPrecedes(chain, best)) {
                    best = chain
                }
            }
        }

        return best.count >= 2 ? best : nil
    }

    private static func longestUnplacedChain(anchor: Int,
                                             start: Int,
                                             component: Set<Int>,
                                             graph: CDKMolecularGraph,
                                             molecule: Molecule,
                                             ringAtoms: Set<Int>,
                                             ringEdges: Set<CDKEdgeKey>,
                                             positions: [Int: CGPoint]) -> [Int] {
        guard positions[start] == nil else { return [anchor] }
        var best: [Int] = [anchor, start]
        var visited: Set<Int> = [anchor, start]

        func updateBest(_ path: [Int]) {
            if path.count > best.count ||
                (path.count == best.count && pathLexicographicallyPrecedes(path, best)) {
                best = path
            }
        }

        func dfs(prev: Int, cur: Int, path: [Int]) {
            if ringAtoms.contains(cur) {
                updateBest(path)
                return
            }
            guard chainEligibleAtom(cur, molecule: molecule, ringAtoms: ringAtoms) else {
                updateBest(path)
                return
            }

            var extended = false
            let neighbors = graph.neighbors(of: cur)
                .filter { component.contains($0) && $0 != prev && positions[$0] == nil && !visited.contains($0) }
                .sorted()

            for nxt in neighbors {
                if ringAtoms.contains(nxt) {
                    updateBest(path + [nxt])
                    continue
                }
                guard chainEligibleAtom(nxt, molecule: molecule, ringAtoms: ringAtoms) else { continue }
                guard chainEligibleEdge(cur, nxt, graph: graph, ringEdges: ringEdges) else { continue }
                visited.insert(nxt)
                dfs(prev: cur, cur: nxt, path: path + [nxt])
                visited.remove(nxt)
                extended = true
            }

            if !extended { updateBest(path) }
        }

        if ringAtoms.contains(start) {
            return best
        }
        dfs(prev: anchor, cur: start, path: [anchor, start])
        return best
    }

    private static func chainEligibleAtom(_ atomID: Int,
                                          molecule: Molecule,
                                          ringAtoms: Set<Int>) -> Bool {
        guard !ringAtoms.contains(atomID), let atom = molecule.atom(id: atomID) else { return false }
        return atom.element.uppercased() != "H"
    }

    private static func chainEligibleEdge(_ a: Int,
                                          _ b: Int,
                                          graph: CDKMolecularGraph,
                                          ringEdges: Set<CDKEdgeKey>) -> Bool {
        let key = CDKEdgeKey(a, b)
        guard !ringEdges.contains(key) else { return false }
        guard let bond = graph.edgeBonds[key] else { return false }
        return bond.order == .single
    }

    private static func placeLinearChain(_ chain: [Int],
                                         initialVector: CGVector,
                                         component: Set<Int>,
                                         positions: inout [Int: CGPoint],
                                         bondLength: CGFloat) {
        guard chain.count >= 2, let anchorPos = positions[chain[0]] else { return }
        let startDir = normalize(initialVector) ?? unitVector(angle: 0)
        if positions[chain[1]] == nil {
            positions[chain[1]] = CGPoint(x: anchorPos.x + startDir.dx * bondLength,
                                          y: anchorPos.y + startDir.dy * bondLength)
        }
        if chain.count < 3 { return }

        let centroid = placedCentroid(component: component, positions: positions, excluding: Set(chain.dropFirst()))
        var lastSign: CGFloat = 0

        for i in 2..<chain.count {
            let a = chain[i - 2]
            let b = chain[i - 1]
            let c = chain[i]
            guard let pa = positions[a], let pb = positions[b] else { continue }

            let back = normalize(CGVector(dx: pa.x - pb.x, dy: pa.y - pb.y)) ?? unitVector(angle: 0)
            let dir1 = rotate(back, by: Tuning.chainPreferredAngle)
            let dir2 = rotate(back, by: -Tuning.chainPreferredAngle)
            let p1 = CGPoint(x: pb.x + dir1.dx * bondLength, y: pb.y + dir1.dy * bondLength)
            let p2 = CGPoint(x: pb.x + dir2.dx * bondLength, y: pb.y + dir2.dy * bondLength)
            let s1 = chainPointScore(point: p1, centroid: centroid, positions: positions, bondLength: bondLength)
            let s2 = chainPointScore(point: p2, centroid: centroid, positions: positions, bondLength: bondLength)

            let sign1: CGFloat = 1
            let sign2: CGFloat = -1
            let pickFirst: Bool
            if lastSign == 0 {
                pickFirst = s1 <= s2
            } else {
                let preferredSign = -lastSign
                let preferredScore = preferredSign > 0 ? s1 : s2
                let alternateScore = preferredSign > 0 ? s2 : s1
                if preferredScore <= alternateScore * 1.08 {
                    pickFirst = preferredSign > 0
                } else {
                    pickFirst = s1 <= s2
                }
            }

            positions[c] = pickFirst ? p1 : p2
            lastSign = pickFirst ? sign1 : sign2
        }
    }

    private static func initialChainVector(anchor: Int,
                                           firstUnplaced: Int,
                                           component: Set<Int>,
                                           graph: CDKMolecularGraph,
                                           positions: [Int: CGPoint]) -> CGVector {
        guard let center = positions[anchor] else {
            let angle = CGFloat((anchor * 41 + firstUnplaced * 17) % 360) * .pi / 180.0
            return unitVector(angle: angle)
        }

        var sum = CGVector.zero
        for n in graph.neighbors(of: anchor) where component.contains(n) {
            guard let p = positions[n], let u = normalize(CGVector(dx: p.x - center.x, dy: p.y - center.y)) else { continue }
            sum.dx += u.dx
            sum.dy += u.dy
        }

        if let open = normalize(CGVector(dx: -sum.dx, dy: -sum.dy)) {
            return open
        }
        let angle = CGFloat((anchor * 41 + firstUnplaced * 17) % 360) * .pi / 180.0
        return unitVector(angle: angle)
    }

    private static func chainPointScore(point: CGPoint,
                                        centroid: CGPoint?,
                                        positions: [Int: CGPoint],
                                        bondLength: CGFloat) -> CGFloat {
        var score: CGFloat = 0
        let hard = bondLength * 0.95
        let soft = bondLength * 1.20
        for q in positions.values {
            let d = point.distance(to: q)
            if d < hard {
                let x = hard - d
                score += x * x * 180
            } else if d < soft {
                let x = soft - d
                score += x * x * 24
            }
        }
        if let centroid {
            score -= point.distance(to: centroid) * 0.22
        }
        return score
    }

    private static func placedCentroid(component: Set<Int>,
                                       positions: [Int: CGPoint],
                                       excluding: Set<Int> = []) -> CGPoint? {
        let points = component
            .filter { !excluding.contains($0) }
            .compactMap { positions[$0] }
        guard !points.isEmpty else { return nil }
        let sx = points.reduce(CGFloat(0)) { $0 + $1.x }
        let sy = points.reduce(CGFloat(0)) { $0 + $1.y }
        return CGPoint(x: sx / CGFloat(points.count), y: sy / CGFloat(points.count))
    }

    private static func pathLexicographicallyPrecedes(_ lhs: [Int], _ rhs: [Int]) -> Bool {
        guard !rhs.isEmpty else { return true }
        for (a, b) in zip(lhs, rhs) {
            if a != b { return a < b }
        }
        return lhs.count < rhs.count
    }

    private static func placeRegularRing(_ ring: [Int],
                                         center: CGPoint,
                                         bondLength: CGFloat,
                                         positions: inout [Int: CGPoint]) {
        guard ring.count >= 3 else { return }
        let n = CGFloat(ring.count)
        let radius = bondLength / (2 * sin(.pi / n))
        let base = CGFloat((ring[0] * 11) % 360) * .pi / 180.0
        let step = (2 * .pi) / n

        for (idx, atomID) in ring.enumerated() {
            let angle = base + CGFloat(idx) * step
            positions[atomID] = CGPoint(
                x: center.x + cos(angle) * radius,
                y: center.y + sin(angle) * radius
            )
        }
    }

    private static func regularRingLocalCoordinates(_ ring: [Int], bondLength: CGFloat) -> [Int: CGPoint] {
        guard ring.count >= 3 else { return [:] }
        let n = CGFloat(ring.count)
        let radius = bondLength / (2 * sin(.pi / n))
        let step = (2 * .pi) / n
        var out: [Int: CGPoint] = [:]
        for (idx, atomID) in ring.enumerated() {
            let angle = CGFloat(idx) * step
            out[atomID] = CGPoint(x: cos(angle) * radius, y: sin(angle) * radius)
        }
        return out
    }

    private static func bestRingPlacementCandidate(ring: [Int],
                                                   sharedAnchors: [Int],
                                                   graph: CDKMolecularGraph,
                                                   component: Set<Int>,
                                                   positions: [Int: CGPoint],
                                                   bondLength: CGFloat) -> [Int: CGPoint]? {
        var candidates: [[Int: CGPoint]] = []
        let local = regularRingLocalCoordinates(ring, bondLength: bondLength)
        let localMirror = local.mapValues { CGPoint(x: $0.x, y: -$0.y) }

        if sharedAnchors.count >= 2 {
            for pair in anchorPairs(in: ring, shared: sharedAnchors, maxPairs: 4) {
                guard let t1 = positions[pair.0], let t2 = positions[pair.1] else { continue }
                if let c1 = transformedRing(local: local, ring: ring, anchor1: pair.0, anchor2: pair.1, target1: t1, target2: t2) {
                    candidates.append(c1)
                }
                if let c2 = transformedRing(local: localMirror, ring: ring, anchor1: pair.0, anchor2: pair.1, target1: t1, target2: t2) {
                    candidates.append(c2)
                }
            }
        } else if let s = sharedAnchors.first, let target = positions[s] {
            candidates += singleAnchorCandidates(ring: ring,
                                                 local: local,
                                                 localMirror: localMirror,
                                                 sharedAtom: s,
                                                 target: target,
                                                 graph: graph,
                                                 positions: positions)
        }

        guard !candidates.isEmpty else { return nil }
        return candidates.min { lhs, rhs in
            let sl = ringPlacementScore(candidate: lhs,
                                        ring: ring,
                                        sharedAnchors: sharedAnchors,
                                        graph: graph,
                                        component: component,
                                        positions: positions,
                                        bondLength: bondLength)
            let sr = ringPlacementScore(candidate: rhs,
                                        ring: ring,
                                        sharedAnchors: sharedAnchors,
                                        graph: graph,
                                        component: component,
                                        positions: positions,
                                        bondLength: bondLength)
            return sl < sr
        }
    }

    private static func anchorPairs(in ring: [Int], shared: [Int], maxPairs: Int) -> [(Int, Int)] {
        guard shared.count >= 2 else { return [] }
        struct AnchorPairKey: Hashable {
            let a: Int
            let b: Int
            init(_ x: Int, _ y: Int) {
                if x <= y {
                    a = x
                    b = y
                } else {
                    a = y
                    b = x
                }
            }
        }

        let sharedSet = Set(shared)
        var pairs: [(Int, Int)] = []
        var seen: Set<AnchorPairKey> = []

        func appendPair(_ a: Int, _ b: Int) {
            let key = AnchorPairKey(a, b)
            guard seen.insert(key).inserted else { return }
            pairs.append((a, b))
        }

        // Prefer adjacent shared anchors (fused rings), then broader-separation pairs.
        for i in 0..<ring.count {
            let a = ring[i]
            let b = ring[(i + 1) % ring.count]
            if sharedSet.contains(a) && sharedSet.contains(b) {
                appendPair(a, b)
            }
        }

        var byGap: [(gap: Int, a: Int, b: Int)] = []
        for i in 0..<shared.count {
            for j in (i + 1)..<shared.count {
                let a = shared[i]
                let b = shared[j]
                guard let ia = ring.firstIndex(of: a), let ib = ring.firstIndex(of: b) else { continue }
                let diff = abs(ia - ib)
                let gap = min(diff, ring.count - diff)
                byGap.append((gap, a, b))
            }
        }
        byGap.sort { lhs, rhs in
            if lhs.gap != rhs.gap { return lhs.gap > rhs.gap }
            if lhs.a != rhs.a { return lhs.a < rhs.a }
            return lhs.b < rhs.b
        }
        for item in byGap {
            appendPair(item.a, item.b)
            if pairs.count >= maxPairs { break }
        }
        return Array(pairs.prefix(max(1, maxPairs)))
    }

    private static func transformedRing(local: [Int: CGPoint],
                                        ring: [Int],
                                        anchor1: Int,
                                        anchor2: Int,
                                        target1: CGPoint,
                                        target2: CGPoint) -> [Int: CGPoint]? {
        guard let l1 = local[anchor1], let l2 = local[anchor2] else { return nil }
        let lv = CGVector(dx: l2.x - l1.x, dy: l2.y - l1.y)
        let tv = CGVector(dx: target2.x - target1.x, dy: target2.y - target1.y)
        let ll = hypot(lv.dx, lv.dy)
        let tl = hypot(tv.dx, tv.dy)
        guard ll > 0.0001, tl > 0.0001 else { return nil }

        let scale = tl / ll
        let angle = atan2(tv.dy, tv.dx) - atan2(lv.dy, lv.dx)
        let c = cos(angle)
        let s = sin(angle)

        var out: [Int: CGPoint] = [:]
        for atomID in ring {
            guard let p = local[atomID] else { continue }
            let sx = (p.x - l1.x) * scale
            let sy = (p.y - l1.y) * scale
            let rx = sx * c - sy * s
            let ry = sx * s + sy * c
            out[atomID] = CGPoint(x: target1.x + rx, y: target1.y + ry)
        }
        return out
    }

    private static func singleAnchorCandidates(ring: [Int],
                                               local: [Int: CGPoint],
                                               localMirror: [Int: CGPoint],
                                               sharedAtom: Int,
                                               target: CGPoint,
                                               graph: CDKMolecularGraph,
                                               positions: [Int: CGPoint]) -> [[Int: CGPoint]] {
        guard let idx = ring.firstIndex(of: sharedAtom) else { return [] }
        let n1 = ring[(idx + 1) % ring.count]
        let n2 = ring[(idx - 1 + ring.count) % ring.count]
        let preferred = preferredExpansionDirection(for: sharedAtom, graph: graph, positions: positions)

        func makeCandidate(from base: [Int: CGPoint], neighbor: Int) -> [Int: CGPoint]? {
            guard let ls = base[sharedAtom], let ln = base[neighbor] else { return nil }
            let localVec = CGVector(dx: ln.x - ls.x, dy: ln.y - ls.y)
            guard let uLocal = normalize(localVec) else { return nil }
            guard let uPref = normalize(preferred) else { return nil }
            let angle = atan2(uPref.dy, uPref.dx) - atan2(uLocal.dy, uLocal.dx)
            let c = cos(angle)
            let s = sin(angle)

            var out: [Int: CGPoint] = [:]
            for atomID in ring {
                guard let p = base[atomID] else { continue }
                let x = p.x - ls.x
                let y = p.y - ls.y
                out[atomID] = CGPoint(
                    x: target.x + x * c - y * s,
                    y: target.y + x * s + y * c
                )
            }
            return out
        }

        var out: [[Int: CGPoint]] = []
        if let c1 = makeCandidate(from: local, neighbor: n1) { out.append(c1) }
        if let c2 = makeCandidate(from: localMirror, neighbor: n2) { out.append(c2) }
        return out
    }

    private static func preferredExpansionDirection(for atomID: Int,
                                                    graph: CDKMolecularGraph,
                                                    positions: [Int: CGPoint]) -> CGVector {
        guard let center = positions[atomID] else { return CGVector(dx: 1, dy: 0) }
        let placedNeighbors = graph.neighbors(of: atomID).compactMap { id -> CGPoint? in
            positions[id]
        }
        guard !placedNeighbors.isEmpty else {
            let angle = CGFloat((atomID * 53) % 360) * .pi / 180.0
            return CGVector(dx: cos(angle), dy: sin(angle))
        }

        var sum = CGVector.zero
        for p in placedNeighbors {
            if let u = normalize(CGVector(dx: p.x - center.x, dy: p.y - center.y)) {
                sum.dx += u.dx
                sum.dy += u.dy
            }
        }
        if let open = normalize(CGVector(dx: -sum.dx, dy: -sum.dy)) {
            return open
        }
        if let first = placedNeighbors.first {
            let v = CGVector(dx: first.x - center.x, dy: first.y - center.y)
            return normalize(CGVector(dx: -v.dy, dy: v.dx)) ?? CGVector(dx: 1, dy: 0)
        }
        return CGVector(dx: 1, dy: 0)
    }

    private static func ringPlacementScore(candidate: [Int: CGPoint],
                                           ring: [Int],
                                           sharedAnchors: [Int],
                                           graph: CDKMolecularGraph,
                                           component: Set<Int>,
                                           positions: [Int: CGPoint],
                                           bondLength: CGFloat) -> CGFloat {
        var score: CGFloat = 0
        let sharedSet = Set(sharedAnchors)
        let ringSet = Set(ring)

        for atomID in sharedAnchors {
            if let p = candidate[atomID], let fixed = positions[atomID] {
                score += p.distance(to: fixed) * Tuning.anchorDriftWeight
            }
        }

        let existingAtoms = component.filter { positions[$0] != nil && !ringSet.contains($0) }
        let newAtoms = ring.filter { !sharedSet.contains($0) }

        for atomID in newAtoms {
            guard let p = candidate[atomID] else { continue }
            for other in existingAtoms {
                guard let q = positions[other] else { continue }
                let d = p.distance(to: q)
                let hard = bondLength * Tuning.hardOverlapRatio
                let soft = bondLength * Tuning.softOverlapRatio
                if d < hard {
                    let x = hard - d
                    score += x * x * Tuning.hardOverlapPenalty
                } else if d < soft {
                    let x = soft - d
                    score += x * x * Tuning.softOverlapPenalty
                }
            }
        }

        for i in 0..<newAtoms.count {
            for j in (i + 1)..<newAtoms.count {
                guard let p1 = candidate[newAtoms[i]], let p2 = candidate[newAtoms[j]] else { continue }
                let d = p1.distance(to: p2)
                let hard = bondLength * Tuning.intraRingHardRatio
                if d < hard {
                    let x = hard - d
                    score += x * x * Tuning.intraRingPenalty
                }
            }
        }

        let ringEdges = CDKRingSearch.edgeKeys(for: ring)
        let existingEdges: [(Int, Int)] = graph.edgeBonds.values.compactMap { b in
            guard component.contains(b.a1), component.contains(b.a2) else { return nil }
            guard positions[b.a1] != nil, positions[b.a2] != nil else { return nil }
            let key = CDKEdgeKey(b.a1, b.a2)
            if ringEdges.contains(key) { return nil }
            return (b.a1, b.a2)
        }

        for edge in ringEdges {
            let a = edge.a
            let b = edge.b
            guard let pa = candidate[a], let pb = candidate[b] else { continue }
            for (u, v) in existingEdges {
                if a == u || a == v || b == u || b == v { continue }
                guard let pu = positions[u], let pv = positions[v] else { continue }
                if segmentsIntersect(pa, pb, pu, pv) {
                    score += Tuning.edgeCrossPenalty
                } else {
                    let near = bondLength * Tuning.edgeNearRatio
                    let d = segmentDistance(pa, pb, pu, pv)
                    if d < near {
                        let x = near - d
                        score += x * x * Tuning.edgeNearPenalty
                    }
                }
            }
        }

        return score
    }

    private static func chooseSeed(in component: Set<Int>,
                                   graph: CDKMolecularGraph,
                                   ringSet: Set<Int>) -> Int {
        component.max { lhs, rhs in
            let lScore = (ringSet.contains(lhs) ? 100 : 0) + graph.neighbors(of: lhs).count
            let rScore = (ringSet.contains(rhs) ? 100 : 0) + graph.neighbors(of: rhs).count
            if lScore != rScore { return lScore < rScore }
            return lhs < rhs
        } ?? (component.min() ?? 0)
    }

    private static func isAromaticAtom(_ atomID: Int, molecule: Molecule) -> Bool {
        if molecule.atom(id: atomID)?.aromatic == true { return true }
        return molecule.bonds(forAtom: atomID).contains { $0.order == .aromatic }
    }

    private static func proposedDirections(center: Int,
                                           placedNeighbors: [Int],
                                           unplacedCount: Int,
                                           positions: [Int: CGPoint],
                                           molecule: Molecule,
                                           graph: CDKMolecularGraph) -> [CGVector] {
        guard unplacedCount > 0 else { return [] }
        guard let centerPos = positions[center] else {
            return (0..<unplacedCount).map { i in
                unitVector(angle: CGFloat(i) * (2 * .pi / CGFloat(max(1, unplacedCount))))
            }
        }

        if placedNeighbors.isEmpty {
            let base = CGFloat((center * 47) % 360) * .pi / 180.0
            return fanDirections(count: unplacedCount, baseAngle: base, totalSpread: 2 * .pi)
        }

        if placedNeighbors.count == 1, let parentPos = positions[placedNeighbors[0]] {
            let toParent = CGVector(dx: parentPos.x - centerPos.x, dy: parentPos.y - centerPos.y)
            let uParent = normalize(toParent) ?? CGVector(dx: -1, dy: 0)
            let target = preferredAngle(at: center, molecule: molecule, graph: graph)

            if unplacedCount == 1 {
                let sign: CGFloat = ((center + placedNeighbors[0]) % 2 == 0) ? 1.0 : -1.0
                return [rotate(uParent, by: sign * target)]
            }

            let opposite = rotate(uParent, by: .pi)
            let base = atan2(opposite.dy, opposite.dx)
            let spread = min(Tuning.branchOpenSpread, target + 0.5)
            return fanDirections(count: unplacedCount, baseAngle: base, totalSpread: spread * 2)
        }

        let vecs = placedNeighbors.compactMap { id -> CGVector? in
            guard let p = positions[id] else { return nil }
            return normalize(CGVector(dx: p.x - centerPos.x, dy: p.y - centerPos.y))
        }
        let sum = vecs.reduce(CGVector.zero) { acc, v in
            CGVector(dx: acc.dx + v.dx, dy: acc.dy + v.dy)
        }
        let open = normalize(CGVector(dx: -sum.dx, dy: -sum.dy))
            ?? normalize(CGVector(dx: -vecs.first!.dy, dy: vecs.first!.dx))
            ?? CGVector(dx: 1, dy: 0)

        let base = atan2(open.dy, open.dx)
        let spread = Tuning.branchFanSpread
        return fanDirections(count: unplacedCount, baseAngle: base, totalSpread: spread)
    }

    private static func preferredAngle(at atomID: Int, molecule: Molecule, graph: CDKMolecularGraph) -> CGFloat {
        let neighbors = graph.neighbors(of: atomID)
        guard let atom = molecule.atom(id: atomID) else { return 109.5 * .pi / 180.0 }

        let hasPiLike = molecule.bonds(forAtom: atomID).contains { $0.order == .double || $0.order == .triple || $0.order == .aromatic }
            || atom.aromatic
            || neighbors.contains(where: { n in
                molecule.bonds(forAtom: n).contains { b in
                    let other = (b.a1 == n) ? b.a2 : b.a1
                    return other != atomID && (b.order == .double || b.order == .triple || b.order == .aromatic)
                }
            })

        switch atom.element.uppercased() {
        case "C", "N":
            return hasPiLike ? (120.0 * .pi / 180.0) : (109.5 * .pi / 180.0)
        case "O", "S":
            return hasPiLike ? (120.0 * .pi / 180.0) : (109.5 * .pi / 180.0)
        default:
            return hasPiLike ? (120.0 * .pi / 180.0) : (109.5 * .pi / 180.0)
        }
    }

    private static func optimizeByBondFlips(component: Set<Int>,
                                            graph: CDKMolecularGraph,
                                            positions: inout [Int: CGPoint],
                                            locked: Set<Int>,
                                            ringEdges: Set<CDKEdgeKey>,
                                            bondLength: CGFloat) {
        let candidateBonds = graph.edgeBonds.values
            .filter { component.contains($0.a1) && component.contains($0.a2) && $0.order == .single }
            .sorted { $0.id < $1.id }

        for bond in candidateBonds {
            let edge = CDKEdgeKey(bond.a1, bond.a2)
            if ringEdges.contains(edge) { continue }

            let left = sideComponent(start: bond.a1,
                                     blockedA: bond.a1,
                                     blockedB: bond.a2,
                                     graph: graph,
                                     component: component)
            if left.contains(bond.a2) { continue } // not a bridge in this component
            let right = component.subtracting(left)

            let sideToFlip = (left.count <= right.count) ? left : right
            guard !sideToFlip.isEmpty else { continue }
            if sideToFlip.contains(where: { locked.contains($0) }) { continue }
            guard let p1 = positions[bond.a1], let p2 = positions[bond.a2] else { continue }

            let before = layoutPenalty(component: component,
                                       graph: graph,
                                       positions: positions,
                                       bondLength: bondLength)
            var trial = positions
            for atomID in sideToFlip {
                guard let p = trial[atomID] else { continue }
                trial[atomID] = reflect(point: p, acrossLineFrom: p1, to: p2)
            }
            let after = layoutPenalty(component: component,
                                      graph: graph,
                                      positions: trial,
                                      bondLength: bondLength)
            if after < before * Tuning.bondFlipGainThreshold {
                positions = trial
            }
        }
    }

    private static func sideComponent(start: Int,
                                      blockedA: Int,
                                      blockedB: Int,
                                      graph: CDKMolecularGraph,
                                      component: Set<Int>) -> Set<Int> {
        var seen: Set<Int> = [start]
        var stack: [Int] = [start]
        while let cur = stack.popLast() {
            for nxt in graph.neighbors(of: cur) where component.contains(nxt) {
                if (cur == blockedA && nxt == blockedB) || (cur == blockedB && nxt == blockedA) {
                    continue
                }
                if seen.insert(nxt).inserted {
                    stack.append(nxt)
                }
            }
        }
        return seen
    }

    private static func layoutPenalty(component: Set<Int>,
                                      graph: CDKMolecularGraph,
                                      positions: [Int: CGPoint],
                                      bondLength: CGFloat) -> CGFloat {
        var score: CGFloat = 0
        let atoms = component.sorted()
        let bonds = graph.edgeBonds.values.filter { component.contains($0.a1) && component.contains($0.a2) }
        let bonded = Set(bonds.map { CDKEdgeKey($0.a1, $0.a2) })

        for i in 0..<atoms.count {
            for j in (i + 1)..<atoms.count {
                let a = atoms[i]
                let b = atoms[j]
                if bonded.contains(CDKEdgeKey(a, b)) { continue }
                guard let p = positions[a], let q = positions[b] else { continue }
                let d = p.distance(to: q)
                let hard = bondLength * 0.95
                let soft = bondLength * 1.20
                if d < hard {
                    let x = hard - d
                    score += x * x * 120
                } else if d < soft {
                    let x = soft - d
                    score += x * x * 16
                }
            }
        }

        for i in 0..<bonds.count {
            for j in (i + 1)..<bonds.count {
                let b1 = bonds[i]
                let b2 = bonds[j]
                if b1.a1 == b2.a1 || b1.a1 == b2.a2 || b1.a2 == b2.a1 || b1.a2 == b2.a2 {
                    continue
                }
                guard let p1 = positions[b1.a1], let p2 = positions[b1.a2],
                      let q1 = positions[b2.a1], let q2 = positions[b2.a2] else { continue }
                if segmentsIntersect(p1, p2, q1, q2) {
                    score += 160
                } else {
                    let d = segmentDistance(p1, p2, q1, q2)
                    if d < bondLength * 0.35 {
                        let x = bondLength * 0.35 - d
                        score += x * x * 40
                    }
                }
            }
        }

        return score
    }

    private static func relax(component: Set<Int>,
                              graph: CDKMolecularGraph,
                              positions: inout [Int: CGPoint],
                              bondLength: CGFloat,
                              locked: Set<Int>,
                              iterations: Int) {
        let atoms = component.sorted()
        guard atoms.count >= 2 else { return }
        let bonded = Set(graph.edgeBonds.keys)
        let bonds = graph.edgeBonds.values.filter { component.contains($0.a1) && component.contains($0.a2) }

        for _ in 0..<max(1, iterations) {
            for bond in bonds {
                guard let p1 = positions[bond.a1], let p2 = positions[bond.a2] else { continue }
                let dx = p2.x - p1.x
                let dy = p2.y - p1.y
                let d = max(0.0001, hypot(dx, dy))
                let ux = dx / d
                let uy = dy / d
                let delta = (d - bondLength) * Tuning.bondSpringFactor

                if !locked.contains(bond.a1) {
                    positions[bond.a1] = CGPoint(x: p1.x + ux * delta, y: p1.y + uy * delta)
                }
                if !locked.contains(bond.a2) {
                    positions[bond.a2] = CGPoint(x: p2.x - ux * delta, y: p2.y - uy * delta)
                }
            }

            // Non-bonded overlap push.
            for i in 0..<atoms.count {
                for j in (i + 1)..<atoms.count {
                    let a = atoms[i]
                    let b = atoms[j]
                    if bonded.contains(CDKEdgeKey(a, b)) { continue }
                    guard let pa = positions[a], let pb = positions[b] else { continue }
                    let dx = pb.x - pa.x
                    let dy = pb.y - pa.y
                    let d = max(0.0001, hypot(dx, dy))
                    let minD = bondLength * Tuning.nonBondedMinDistanceRatio
                    guard d < minD else { continue }
                    let push = (minD - d) * Tuning.nonBondedPushFactor
                    let ux = dx / d
                    let uy = dy / d

                    if !locked.contains(a) {
                        positions[a] = CGPoint(x: pa.x - ux * push, y: pa.y - uy * push)
                    }
                    if !locked.contains(b) {
                        positions[b] = CGPoint(x: pb.x + ux * push, y: pb.y + uy * push)
                    }
                }
            }

            // Untangle edge crossings.
            for i in 0..<bonds.count {
                for j in (i + 1)..<bonds.count {
                    let b1 = bonds[i]
                    let b2 = bonds[j]
                    if b1.a1 == b2.a1 || b1.a1 == b2.a2 || b1.a2 == b2.a1 || b1.a2 == b2.a2 {
                        continue
                    }
                    guard let p1 = positions[b1.a1], let p2 = positions[b1.a2],
                          let q1 = positions[b2.a1], let q2 = positions[b2.a2] else { continue }
                    guard segmentsIntersect(p1, p2, q1, q2) else { continue }

                    let v1 = normalize(CGVector(dx: p2.x - p1.x, dy: p2.y - p1.y)) ?? CGVector(dx: 1, dy: 0)
                    let v2 = normalize(CGVector(dx: q2.x - q1.x, dy: q2.y - q1.y)) ?? CGVector(dx: 0, dy: 1)
                    let sign: CGFloat = (v1.dx * v2.dy - v1.dy * v2.dx) >= 0 ? 1 : -1
                    let perp1 = CGVector(dx: -v1.dy * sign, dy: v1.dx * sign)
                    let perp2 = CGVector(dx: v2.dy * sign, dy: -v2.dx * sign)
                    let push = bondLength * Tuning.crossingPushRatio

                    if !locked.contains(b1.a1), let p = positions[b1.a1] {
                        positions[b1.a1] = CGPoint(x: p.x + perp1.dx * push, y: p.y + perp1.dy * push)
                    }
                    if !locked.contains(b1.a2), let p = positions[b1.a2] {
                        positions[b1.a2] = CGPoint(x: p.x + perp1.dx * push, y: p.y + perp1.dy * push)
                    }
                    if !locked.contains(b2.a1), let p = positions[b2.a1] {
                        positions[b2.a1] = CGPoint(x: p.x + perp2.dx * push, y: p.y + perp2.dy * push)
                    }
                    if !locked.contains(b2.a2), let p = positions[b2.a2] {
                        positions[b2.a2] = CGPoint(x: p.x + perp2.dx * push, y: p.y + perp2.dy * push)
                    }
                }
            }
        }
    }

    private static func boundingBox(component: Set<Int>, positions: [Int: CGPoint]) -> CGRect? {
        let points = component.compactMap { positions[$0] }
        guard let first = points.first else { return nil }
        var minX = first.x
        var minY = first.y
        var maxX = first.x
        var maxY = first.y
        for p in points {
            minX = min(minX, p.x)
            minY = min(minY, p.y)
            maxX = max(maxX, p.x)
            maxY = max(maxY, p.y)
        }
        return CGRect(x: minX, y: minY, width: max(0.0001, maxX - minX), height: max(0.0001, maxY - minY))
    }

    private static func fanDirections(count: Int, baseAngle: CGFloat, totalSpread: CGFloat) -> [CGVector] {
        guard count > 0 else { return [] }
        if count == 1 { return [unitVector(angle: baseAngle)] }

        let start = baseAngle - totalSpread * 0.5
        let step = totalSpread / CGFloat(count - 1)
        return (0..<count).map { i in
            unitVector(angle: start + CGFloat(i) * step)
        }
    }

    private static func unitVector(angle: CGFloat) -> CGVector {
        CGVector(dx: cos(angle), dy: sin(angle))
    }

    private static func normalize(_ v: CGVector) -> CGVector? {
        let len = hypot(v.dx, v.dy)
        guard len > 0.0001 else { return nil }
        return CGVector(dx: v.dx / len, dy: v.dy / len)
    }

    private static func rotate(_ v: CGVector, by angle: CGFloat) -> CGVector {
        let c = cos(angle)
        let s = sin(angle)
        return CGVector(
            dx: v.dx * c - v.dy * s,
            dy: v.dx * s + v.dy * c
        )
    }

    private static func segmentsIntersect(_ p1: CGPoint, _ p2: CGPoint, _ p3: CGPoint, _ p4: CGPoint) -> Bool {
        let o1 = orientation(p1, p2, p3)
        let o2 = orientation(p1, p2, p4)
        let o3 = orientation(p3, p4, p1)
        let o4 = orientation(p3, p4, p2)
        return o1 * o2 < 0 && o3 * o4 < 0
    }

    private static func orientation(_ a: CGPoint, _ b: CGPoint, _ c: CGPoint) -> CGFloat {
        (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
    }

    private static func reflect(point p: CGPoint, acrossLineFrom a: CGPoint, to b: CGPoint) -> CGPoint {
        let vx = b.x - a.x
        let vy = b.y - a.y
        let len2 = vx * vx + vy * vy
        guard len2 > 0.0001 else { return p }
        let t = ((p.x - a.x) * vx + (p.y - a.y) * vy) / len2
        let proj = CGPoint(x: a.x + t * vx, y: a.y + t * vy)
        return CGPoint(x: 2 * proj.x - p.x, y: 2 * proj.y - p.y)
    }

    private static func segmentDistance(_ p1: CGPoint, _ p2: CGPoint, _ p3: CGPoint, _ p4: CGPoint) -> CGFloat {
        if segmentsIntersect(p1, p2, p3, p4) { return 0 }
        return min(
            pointSegmentDistance(p1, a: p3, b: p4),
            pointSegmentDistance(p2, a: p3, b: p4),
            pointSegmentDistance(p3, a: p1, b: p2),
            pointSegmentDistance(p4, a: p1, b: p2)
        )
    }

    private static func pointSegmentDistance(_ p: CGPoint, a: CGPoint, b: CGPoint) -> CGFloat {
        let vx = b.x - a.x
        let vy = b.y - a.y
        let len2 = vx * vx + vy * vy
        guard len2 > 0.0001 else { return p.distance(to: a) }
        let t = max(0, min(1, ((p.x - a.x) * vx + (p.y - a.y) * vy) / len2))
        let proj = CGPoint(x: a.x + t * vx, y: a.y + t * vy)
        return p.distance(to: proj)
    }
}
