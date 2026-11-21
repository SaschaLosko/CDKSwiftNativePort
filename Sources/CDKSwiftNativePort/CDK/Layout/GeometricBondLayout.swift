import Foundation

#if canImport(CoreGraphics)
    import CoreGraphics
#endif

enum CDKGeometricBondLayout {
    static func applying(to molecule: Molecule) -> Molecule {
        let geometricBonds = molecule.bonds.filter { $0.order == .double && $0.doubleBondStereo != nil }
        guard !geometricBonds.isEmpty else { return molecule }

        let graph = CDKMolecularGraph(molecule: molecule)
        let directionCenters = directionCenters(in: molecule)
        let fixedStereoCenters = Set(
            molecule.bonds.compactMap { bond -> Int? in
                guard bond.order == .single, !isDirectionMarker(bond, centers: directionCenters) else { return nil }
                switch bond.stereo {
                case .up, .down: return bond.a1
                case .upReversed, .downReversed: return bond.a2
                case .none, .either: return nil
                }
            })
        var positions = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0.position) })

        for bond in geometricBonds {
            guard let configuration = bond.doubleBondStereo,
                let refs = bond.stereoReferenceAtomIDs, refs.count == 4,
                CDKEdgeKey(refs[1], refs[2]) == CDKEdgeKey(bond.a1, bond.a2),
                let first = positions[refs[0]], let start = positions[refs[1]],
                let end = positions[refs[2]], let last = positions[refs[3]]
            else { continue }

            let dx = end.x - start.x
            let dy = end.y - start.y
            let lengthSquared = dx * dx + dy * dy
            let firstSide = dx * (first.y - start.y) - dy * (first.x - start.x)
            let lastSide = dx * (last.y - start.y) - dy * (last.x - start.x)
            guard lengthSquared > 1e-10, abs(firstSide) > 1e-8, abs(lastSide) > 1e-8,
                (firstSide * lastSide > 0) != (configuration == .cis),
                let left = fragment(graph: graph, from: refs[1], excluding: refs[2]),
                let right = fragment(graph: graph, from: refs[2], excluding: refs[1])
            else { continue }

            // Reflect one connected fragment across the double bond. A ring has
            // no separable fragment. Preserve any previously assigned wedges.
            let candidates = left.count <= right.count ? [left, right] : [right, left]
            guard let movable = candidates.first(where: { $0.isDisjoint(with: fixedStereoCenters) }) else {
                continue
            }
            for id in movable {
                guard let point = positions[id] else { continue }
                let fraction = ((point.x - start.x) * dx + (point.y - start.y) * dy) / lengthSquared
                positions[id] = CGPoint(
                    x: 2 * (start.x + dx * fraction) - point.x,
                    y: 2 * (start.y + dy * fraction) - point.y)
            }
        }

        var result = molecule
        for index in result.atoms.indices {
            if let position = positions[result.atoms[index].id] {
                result.atoms[index].position = position
            }
        }
        return result
    }

    /// The SMILES parser stores slash directions as single-bond stereo markers
    /// and records their shared configuration on the double bond.
    static func directionCenters(in molecule: Molecule) -> Set<Int> {
        Set(
            molecule.bonds.lazy.filter { $0.order == .double && $0.stereo == .either && $0.doubleBondStereo != nil }
                .flatMap { [$0.a1, $0.a2] })
    }

    static func isDirectionMarker(_ bond: Bond, centers: Set<Int>) -> Bool {
        guard bond.order == .single, centers.contains(bond.a1) || centers.contains(bond.a2) else { return false }
        return bond.stereo == .up || bond.stereo == .down
    }

    private static func fragment(graph: CDKMolecularGraph, from start: Int, excluding other: Int) -> Set<Int>? {
        let cut = CDKEdgeKey(start, other)
        var visited: Set<Int> = [start]
        var pending = [start]
        while let atom = pending.popLast() {
            for neighbor in graph.neighbors(of: atom) where CDKEdgeKey(atom, neighbor) != cut {
                guard neighbor != other else { return nil }
                if visited.insert(neighbor).inserted { pending.append(neighbor) }
            }
        }
        return visited
    }
}
