import Foundation

#if canImport(CoreGraphics)
    import CoreGraphics
#endif

extension CDKStructureDiagramGenerator {
    static func chooseSeed(
        in component: Set<Int>,
        graph: CDKMolecularGraph,
        ringSet: Set<Int>,
        preferAcyclicTerminalSeed: Bool
    ) -> Int {
        if preferAcyclicTerminalSeed,
            ringSet.isEmpty,
            let acyclicSeed = acyclicTerminalSeed(in: component, graph: graph)
        {
            return acyclicSeed
        }

        return component.max { lhs, rhs in
            let lScore = (ringSet.contains(lhs) ? 100 : 0) + graph.neighbors(of: lhs).count
            let rScore = (ringSet.contains(rhs) ? 100 : 0) + graph.neighbors(of: rhs).count
            if lScore != rScore { return lScore < rScore }
            return lhs < rhs
        } ?? (component.min() ?? 0)
    }

    static func acyclicTerminalSeed(
        in component: Set<Int>,
        graph: CDKMolecularGraph
    ) -> Int? {
        let leaves =
            component
            .filter { graph.neighbors(of: $0).filter(component.contains).count <= 1 }
            .sorted()
        guard !leaves.isEmpty else { return component.min() }

        return leaves.max { lhs, rhs in
            let lhsDistance = maximumGraphDistance(from: lhs, in: component, graph: graph)
            let rhsDistance = maximumGraphDistance(from: rhs, in: component, graph: graph)
            if lhsDistance != rhsDistance { return lhsDistance < rhsDistance }
            return lhs > rhs
        }
    }

    static func maximumGraphDistance(
        from seed: Int,
        in component: Set<Int>,
        graph: CDKMolecularGraph
    ) -> Int {
        var queue: [Int] = [seed]
        var head = 0
        var distance: [Int: Int] = [seed: 0]
        var maxDistance = 0

        while head < queue.count {
            let current = queue[head]
            head += 1
            let currentDistance = distance[current] ?? 0
            maxDistance = max(maxDistance, currentDistance)

            for neighbor in graph.neighbors(of: current) where component.contains(neighbor) {
                guard distance[neighbor] == nil else { continue }
                distance[neighbor] = currentDistance + 1
                queue.append(neighbor)
            }
        }

        return maxDistance
    }

    static func isAromaticAtom(_ atomID: Int, molecule: Molecule) -> Bool {
        if molecule.atom(id: atomID)?.aromatic == true { return true }
        return molecule.bonds(forAtom: atomID).contains { $0.order == .aromatic }
    }

    static func proposedDirections(
        center: Int,
        placedNeighbors: [Int],
        unplacedCount: Int,
        positions: [Int: CGPoint],
        molecule: Molecule,
        graph: CDKMolecularGraph
    ) -> [CGVector] {
        guard unplacedCount > 0 else { return [] }
        guard let centerPos = positions[center] else {
            return (0..<unplacedCount).map { i in
                unitVector(angle: CGFloat(i) * (2 * .pi / CGFloat(max(1, unplacedCount))))
            }
        }

        if placedNeighbors.isEmpty {
            let base = CGFloat((center * 47) % 360) * .pi / 180.0
            return radialDirections(count: unplacedCount, baseAngle: base)
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
        let open =
            normalize(CGVector(dx: -sum.dx, dy: -sum.dy))
            ?? normalize(CGVector(dx: -vecs.first!.dy, dy: vecs.first!.dx))
            ?? CGVector(dx: 1, dy: 0)

        let base = atan2(open.dy, open.dx)
        let spread = Tuning.branchFanSpread
        return fanDirections(count: unplacedCount, baseAngle: base, totalSpread: spread)
    }

    static func preferredAngle(at atomID: Int, molecule: Molecule, graph: CDKMolecularGraph)
        -> CGFloat
    {
        let neighbors = graph.neighbors(of: atomID)
        guard let atom = molecule.atom(id: atomID) else { return 109.5 * .pi / 180.0 }

        let hasPiLike =
            molecule.bonds(forAtom: atomID).contains {
                $0.order == .double || $0.order == .triple || $0.order == .aromatic
            }
            || atom.aromatic
            || neighbors.contains(where: { n in
                molecule.bonds(forAtom: n).contains { b in
                    let other = (b.a1 == n) ? b.a2 : b.a1
                    return other != atomID
                        && (b.order == .double || b.order == .triple || b.order == .aromatic)
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

    static func optimizeByBondFlips(
        component: Set<Int>,
        graph: CDKMolecularGraph,
        positions: inout [Int: CGPoint],
        locked: Set<Int>,
        ringEdges: Set<CDKEdgeKey>,
        bondLength: CGFloat
    ) {
        let candidateBonds = graph.edgeBonds.values
            .filter { component.contains($0.a1) && component.contains($0.a2) && $0.order == .single }
            .sorted { $0.id < $1.id }

        for bond in candidateBonds {
            let edge = CDKEdgeKey(bond.a1, bond.a2)
            if ringEdges.contains(edge) { continue }

            let left = sideComponent(
                start: bond.a1,
                blockedA: bond.a1,
                blockedB: bond.a2,
                graph: graph,
                component: component)
            if left.contains(bond.a2) { continue }  // not a bridge in this component
            let right = component.subtracting(left)

            let sideToFlip = (left.count <= right.count) ? left : right
            guard !sideToFlip.isEmpty else { continue }
            if sideToFlip.contains(where: { locked.contains($0) }) { continue }
            guard let p1 = positions[bond.a1], let p2 = positions[bond.a2] else { continue }

            let before = layoutPenalty(
                component: component,
                graph: graph,
                positions: positions,
                bondLength: bondLength)
            var trial = positions
            for atomID in sideToFlip {
                guard let p = trial[atomID] else { continue }
                trial[atomID] = reflect(point: p, acrossLineFrom: p1, to: p2)
            }
            let after = layoutPenalty(
                component: component,
                graph: graph,
                positions: trial,
                bondLength: bondLength)
            if after < before * Tuning.bondFlipGainThreshold {
                positions = trial
            }
        }
    }

    static func sideComponent(
        start: Int,
        blockedA: Int,
        blockedB: Int,
        graph: CDKMolecularGraph,
        component: Set<Int>
    ) -> Set<Int> {
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

    static func layoutPenalty(
        component: Set<Int>,
        graph: CDKMolecularGraph,
        positions: [Int: CGPoint],
        bondLength: CGFloat
    ) -> CGFloat {
        var score: CGFloat = 0
        let atoms = component.sorted()
        let bonds = graph.edgeBonds.values.filter {
            component.contains($0.a1) && component.contains($0.a2)
        }
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
                    let q1 = positions[b2.a1], let q2 = positions[b2.a2]
                else { continue }
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

    static func relax(
        component: Set<Int>,
        graph: CDKMolecularGraph,
        positions: inout [Int: CGPoint],
        bondLength: CGFloat,
        locked: Set<Int>,
        iterations: Int
    ) {
        let atoms = component.sorted()
        guard atoms.count >= 2 else { return }
        let bonded = Set(graph.edgeBonds.keys)
        let bonds = graph.edgeBonds.values.filter {
            component.contains($0.a1) && component.contains($0.a2)
        }

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
                        let q1 = positions[b2.a1], let q2 = positions[b2.a2]
                    else { continue }
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
}
