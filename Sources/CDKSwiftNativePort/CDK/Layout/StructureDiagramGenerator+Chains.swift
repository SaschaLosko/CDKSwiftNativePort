import Foundation

#if canImport(CoreGraphics)
    import CoreGraphics
#endif

extension CDKStructureDiagramGenerator {
    static func placeDistributedPartners(
        component: Set<Int>,
        graph: CDKMolecularGraph,
        molecule: Molecule,
        ringsInComponent: [[Int]],
        positions: inout [Int: CGPoint],
        bondLength: CGFloat,
        fallbackCenter: CGPoint
    ) -> Bool {
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
                    placeRingSystem(
                        [ring],
                        graph: graph,
                        component: component,
                        positions: &positions,
                        center: positions[center] ?? fallbackCenter,
                        bondLength: bondLength)
                }

                unplaced = neighbors.filter { positions[$0] == nil }
                guard !unplaced.isEmpty else { continue }

                let dirs = proposedDirections(
                    center: center,
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

    static func placeLongestUnplacedChains(
        component: Set<Int>,
        graph: CDKMolecularGraph,
        molecule: Molecule,
        ringAtoms: Set<Int>,
        ringEdges: Set<CDKEdgeKey>,
        positions: inout [Int: CGPoint],
        bondLength: CGFloat
    ) -> Bool {
        var anyPlaced = false
        var passes = 0

        while passes < Tuning.chainPassLimit {
            guard
                let chain = bestUnplacedChain(
                    component: component,
                    graph: graph,
                    molecule: molecule,
                    ringAtoms: ringAtoms,
                    ringEdges: ringEdges,
                    positions: positions),
                chain.count >= 3
            else { break }

            let anchor = chain[0]
            let first = chain[1]
            let initial = initialChainVector(
                anchor: anchor,
                firstUnplaced: first,
                component: component,
                graph: graph,
                positions: positions)
            placeLinearChain(
                chain,
                initialVector: initial,
                component: component,
                positions: &positions,
                bondLength: bondLength)
            anyPlaced = true
            passes += 1
        }

        return anyPlaced
    }

    static func bestUnplacedChain(
        component: Set<Int>,
        graph: CDKMolecularGraph,
        molecule: Molecule,
        ringAtoms: Set<Int>,
        ringEdges: Set<CDKEdgeKey>,
        positions: [Int: CGPoint]
    ) -> [Int]? {
        var best: [Int] = []

        for anchor in component.sorted() where positions[anchor] != nil {
            let starts = graph.neighbors(of: anchor)
                .filter {
                    component.contains($0)
                        && positions[$0] == nil
                        && !ringEdges.contains(CDKEdgeKey(anchor, $0))
                }
                .sorted()

            for start in starts {
                let chain = longestUnplacedChain(
                    anchor: anchor,
                    start: start,
                    component: component,
                    graph: graph,
                    molecule: molecule,
                    ringAtoms: ringAtoms,
                    ringEdges: ringEdges,
                    positions: positions)
                if chain.count > best.count
                    || (chain.count == best.count && pathLexicographicallyPrecedes(chain, best))
                {
                    best = chain
                }
            }
        }

        return best.count >= 3 ? best : nil
    }

    static func longestUnplacedChain(
        anchor: Int,
        start: Int,
        component: Set<Int>,
        graph: CDKMolecularGraph,
        molecule: Molecule,
        ringAtoms: Set<Int>,
        ringEdges: Set<CDKEdgeKey>,
        positions: [Int: CGPoint]
    ) -> [Int] {
        guard positions[start] == nil else { return [anchor] }
        var best: [Int] = [anchor, start]
        var visited: Set<Int> = [anchor, start]

        func updateBest(_ path: [Int]) {
            if path.count > best.count
                || (path.count == best.count && pathLexicographicallyPrecedes(path, best))
            {
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
                .filter {
                    component.contains($0) && $0 != prev && positions[$0] == nil && !visited.contains($0)
                }
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

    static func chainEligibleAtom(
        _ atomID: Int,
        molecule: Molecule,
        ringAtoms: Set<Int>
    ) -> Bool {
        guard !ringAtoms.contains(atomID), let atom = molecule.atom(id: atomID) else { return false }
        return atom.element.uppercased() != "H"
    }

    static func chainEligibleEdge(
        _ a: Int,
        _ b: Int,
        graph: CDKMolecularGraph,
        ringEdges: Set<CDKEdgeKey>
    ) -> Bool {
        let key = CDKEdgeKey(a, b)
        guard !ringEdges.contains(key) else { return false }
        guard let bond = graph.edgeBonds[key] else { return false }
        return bond.order == .single || bond.order == .double
    }

    static func placeLinearChain(
        _ chain: [Int],
        initialVector: CGVector,
        component: Set<Int>,
        positions: inout [Int: CGPoint],
        bondLength: CGFloat
    ) {
        guard chain.count >= 2, let anchorPos = positions[chain[0]] else { return }
        let startDir = normalize(initialVector) ?? unitVector(angle: 0)
        if positions[chain[1]] == nil {
            positions[chain[1]] = CGPoint(
                x: anchorPos.x + startDir.dx * bondLength,
                y: anchorPos.y + startDir.dy * bondLength)
        }
        if chain.count < 3 { return }

        let centroid = placedCentroid(
            component: component, positions: positions, excluding: Set(chain.dropFirst()))
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
            let s1 = chainPointScore(
                point: p1, centroid: centroid, positions: positions, bondLength: bondLength)
            let s2 = chainPointScore(
                point: p2, centroid: centroid, positions: positions, bondLength: bondLength)

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

    static func initialChainVector(
        anchor: Int,
        firstUnplaced: Int,
        component: Set<Int>,
        graph: CDKMolecularGraph,
        positions: [Int: CGPoint]
    ) -> CGVector {
        guard let center = positions[anchor] else {
            let angle = CGFloat((anchor * 41 + firstUnplaced * 17) % 360) * .pi / 180.0
            return unitVector(angle: angle)
        }

        var sum = CGVector.zero
        for n in graph.neighbors(of: anchor) where component.contains(n) {
            guard let p = positions[n],
                let u = normalize(CGVector(dx: p.x - center.x, dy: p.y - center.y))
            else { continue }
            sum.dx += u.dx
            sum.dy += u.dy
        }

        if let open = normalize(CGVector(dx: -sum.dx, dy: -sum.dy)) {
            return open
        }
        let angle = CGFloat((anchor * 41 + firstUnplaced * 17) % 360) * .pi / 180.0
        return unitVector(angle: angle)
    }

    static func chainPointScore(
        point: CGPoint,
        centroid: CGPoint?,
        positions: [Int: CGPoint],
        bondLength: CGFloat
    ) -> CGFloat {
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

    static func placedCentroid(
        component: Set<Int>,
        positions: [Int: CGPoint],
        excluding: Set<Int> = []
    ) -> CGPoint? {
        let points =
            component
            .filter { !excluding.contains($0) }
            .compactMap { positions[$0] }
        guard !points.isEmpty else { return nil }
        let sx = points.reduce(CGFloat(0)) { $0 + $1.x }
        let sy = points.reduce(CGFloat(0)) { $0 + $1.y }
        return CGPoint(x: sx / CGFloat(points.count), y: sy / CGFloat(points.count))
    }

    static func pathLexicographicallyPrecedes(_ lhs: [Int], _ rhs: [Int]) -> Bool {
        guard !rhs.isEmpty else { return true }
        for (a, b) in zip(lhs, rhs) {
            if a != b { return a < b }
        }
        return lhs.count < rhs.count
    }
}
