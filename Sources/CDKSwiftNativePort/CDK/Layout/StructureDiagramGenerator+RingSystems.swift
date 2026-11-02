import Foundation

#if canImport(CoreGraphics)
    import CoreGraphics
#endif

extension CDKStructureDiagramGenerator {
    static func connectedRingSystems(_ rings: [[Int]]) -> [[[Int]]] {
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

    static func simpleRingLockAtoms(_ ringSystems: [[[Int]]]) -> Set<Int> {
        Set(
            ringSystems.compactMap { system -> [Int]? in
                guard system.count == 1 else { return nil }
                return system[0]
            }.flatMap { $0 })
    }

    static func ringEdgeMultiplicity(_ rings: [[Int]]) -> [CDKEdgeKey: Int] {
        var out: [CDKEdgeKey: Int] = [:]
        for ring in rings {
            for edge in CDKRingSearch.edgeKeys(for: ring) {
                out[edge, default: 0] += 1
            }
        }
        return out
    }

    static func fusedEdgeCount(
        _ ringSystem: [[Int]],
        multiplicity: [CDKEdgeKey: Int]
    ) -> Int {
        let edges = Set(ringSystem.flatMap { CDKRingSearch.edgeKeys(for: $0) })
        return edges.reduce(into: 0) { acc, edge in
            if (multiplicity[edge] ?? 0) > 1 { acc += 1 }
        }
    }

    // CDK RingPlacer-like attachment sequencing:
    // fused -> spiro -> bridged, traversed from already placed ring order.
    enum RingAttachmentMode: Int {
        case fused = 0
        case spiro = 1
        case bridged = 2
        case isolated = 3
    }

    static func ringAttachmentMode(_ ring: [Int], to placedRing: [Int]) -> RingAttachmentMode {
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

    static func placeRingSystem(
        _ rings: [[Int]],
        graph: CDKMolecularGraph,
        component: Set<Int>,
        positions: inout [Int: CGPoint],
        center: CGPoint,
        bondLength: CGFloat
    ) {
        guard !rings.isEmpty else { return }
        let ordered = rings.sorted { lhs, rhs in
            if lhs.count != rhs.count { return lhs.count < rhs.count }
            return lhs.lexicographicallyPrecedes(rhs)
        }

        // Seed with a ring that already has most anchors, otherwise place the first ring regularly.
        let seedIndex =
            ordered.indices.max { a, b in
                let ca = ordered[a].filter { positions[$0] != nil }.count
                let cb = ordered[b].filter { positions[$0] != nil }.count
                if ca != cb { return ca < cb }
                return ordered[a].count > ordered[b].count
            } ?? 0

        let seedRing = ordered[seedIndex]
        if seedRing.allSatisfy({ positions[$0] == nil }) {
            placeRegularRing(
                seedRing,
                center: center,
                bondLength: bondLength,
                positions: &positions)
        } else if seedRing.contains(where: { positions[$0] == nil }) {
            let shared = seedRing.filter { positions[$0] != nil }
            if let candidate = bestRingPlacementCandidate(
                ring: seedRing,
                sharedAnchors: shared,
                graph: graph,
                component: component,
                positions: positions,
                bondLength: bondLength)
            {
                for atomID in seedRing where positions[atomID] == nil {
                    positions[atomID] = candidate[atomID]
                }
            }
        }

        var placedRings: Set<Int> = []
        var placementOrder: [Int] = []
        if ordered[seedIndex].allSatisfy({ positions[$0] != nil }) {
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

                if let candidate = bestRingPlacementCandidate(
                    ring: ring,
                    sharedAnchors: shared,
                    graph: graph,
                    component: component,
                    positions: positions,
                    bondLength: bondLength)
                {
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
}
