import Foundation

#if canImport(CoreGraphics)
    import CoreGraphics
#endif

extension CDKStructureDiagramGenerator {
    static func placeRegularRing(
        _ ring: [Int],
        center: CGPoint,
        bondLength: CGFloat,
        positions: inout [Int: CGPoint]
    ) {
        guard ring.count >= 3 else { return }
        let n = CGFloat(ring.count)
        let radius = bondLength / (2 * sin(.pi / n))
        let base: CGFloat
        if ring.count.isMultiple(of: 2) {
            // Canonical even-membered rings with near-vertical side edges read
            // closer to the conventional "benzene-like" hexagon than an
            // arbitrary atom-id-derived rotation.
            base = -(.pi / n)
        } else {
            base = .pi / 2
        }
        let step = (2 * .pi) / n

        for (idx, atomID) in ring.enumerated() {
            let angle = base + CGFloat(idx) * step
            positions[atomID] = CGPoint(
                x: center.x + cos(angle) * radius,
                y: center.y + sin(angle) * radius
            )
        }
    }

    static func regularRingLocalCoordinates(_ ring: [Int], bondLength: CGFloat) -> [Int:
        CGPoint]
    {
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

    static func bestRingPlacementCandidate(
        ring: [Int],
        sharedAnchors: [Int],
        graph: CDKMolecularGraph,
        component: Set<Int>,
        positions: [Int: CGPoint],
        bondLength: CGFloat
    ) -> [Int: CGPoint]? {
        var candidates: [[Int: CGPoint]] = []
        let local = regularRingLocalCoordinates(ring, bondLength: bondLength)
        let localMirror = local.mapValues { CGPoint(x: $0.x, y: -$0.y) }

        if sharedAnchors.count >= 2 {
            for pair in anchorPairs(in: ring, shared: sharedAnchors, maxPairs: 4) {
                guard let t1 = positions[pair.0], let t2 = positions[pair.1] else { continue }
                if let c1 = transformedRing(
                    local: local, ring: ring, anchor1: pair.0, anchor2: pair.1, target1: t1, target2: t2)
                {
                    candidates.append(c1)
                }
                if let c2 = transformedRing(
                    local: localMirror, ring: ring, anchor1: pair.0, anchor2: pair.1, target1: t1, target2: t2
                ) {
                    candidates.append(c2)
                }
            }
        } else if let s = sharedAnchors.first, let target = positions[s] {
            candidates += singleAnchorCandidates(
                ring: ring,
                local: local,
                localMirror: localMirror,
                sharedAtom: s,
                target: target,
                graph: graph,
                positions: positions)
        }

        guard !candidates.isEmpty else { return nil }
        return candidates.min { lhs, rhs in
            let sl = ringPlacementScore(
                candidate: lhs,
                ring: ring,
                sharedAnchors: sharedAnchors,
                graph: graph,
                component: component,
                positions: positions,
                bondLength: bondLength)
            let sr = ringPlacementScore(
                candidate: rhs,
                ring: ring,
                sharedAnchors: sharedAnchors,
                graph: graph,
                component: component,
                positions: positions,
                bondLength: bondLength)
            return sl < sr
        }
    }

    static func anchorPairs(in ring: [Int], shared: [Int], maxPairs: Int) -> [(Int, Int)] {
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

    static func transformedRing(
        local: [Int: CGPoint],
        ring: [Int],
        anchor1: Int,
        anchor2: Int,
        target1: CGPoint,
        target2: CGPoint
    ) -> [Int: CGPoint]? {
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

    static func singleAnchorCandidates(
        ring: [Int],
        local: [Int: CGPoint],
        localMirror: [Int: CGPoint],
        sharedAtom: Int,
        target: CGPoint,
        graph: CDKMolecularGraph,
        positions: [Int: CGPoint]
    ) -> [[Int: CGPoint]] {
        guard ring.contains(sharedAtom) else { return [] }
        let preferred = preferredExpansionDirection(for: sharedAtom, graph: graph, positions: positions)

        func makeCandidate(from base: [Int: CGPoint]) -> [Int: CGPoint]? {
            guard let ls = base[sharedAtom] else { return nil }
            let localCenterVector = CGVector(dx: -ls.x, dy: -ls.y)
            guard let uLocal = normalize(localCenterVector) else { return nil }
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
        if let c1 = makeCandidate(from: local) { out.append(c1) }
        if let c2 = makeCandidate(from: localMirror) { out.append(c2) }
        return out
    }

    static func preferredExpansionDirection(
        for atomID: Int,
        graph: CDKMolecularGraph,
        positions: [Int: CGPoint]
    ) -> CGVector {
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

    static func ringPlacementScore(
        candidate: [Int: CGPoint],
        ring: [Int],
        sharedAnchors: [Int],
        graph: CDKMolecularGraph,
        component: Set<Int>,
        positions: [Int: CGPoint],
        bondLength: CGFloat
    ) -> CGFloat {
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
}
