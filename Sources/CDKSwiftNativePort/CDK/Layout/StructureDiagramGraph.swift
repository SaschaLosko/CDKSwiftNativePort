import Foundation

#if canImport(CoreGraphics)
    import CoreGraphics
#endif

// CDK-inspired SDG port entry point and layout support types.
struct CDKEdgeKey: Hashable {
    let a: Int
    let b: Int

    init(_ u: Int, _ v: Int) {
        a = min(u, v)
        b = max(u, v)
    }
}

struct CDKMolecularGraph {
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

enum CDKRingSearch {
    static func sssrLikeRings(in molecule: Molecule, graph: CDKMolecularGraph) -> [[Int]] {
        _ = molecule
        let components = graph.connectedComponents().count
        let rank = max(0, graph.edgeCount - graph.atomIDs.count + components)
        guard rank > 0 else { return [] }

        // CDK-like ring basis pass: deterministic candidate generation + GF(2) cycle-space selection.
        let maxCycle = max(12, min(40, graph.atomIDs.count))
        let basisCandidates = fundamentalCycleCandidates(graph: graph)
        let shortestCandidates = shortestCycleCandidates(graph: graph, maxCycleSize: maxCycle)
        let candidates = deduplicateCycles(basisCandidates + shortestCandidates)
            .sorted { lhs, rhs in
                if lhs.count != rhs.count { return lhs.count < rhs.count }
                return lhs.lexicographicallyPrecedes(rhs)
            }
        guard !candidates.isEmpty else { return [] }

        let edgeUniverse = Array(graph.edgeBonds.keys).sorted(by: edgeLess)
        let edgeIndex = Dictionary(
            uniqueKeysWithValues: edgeUniverse.enumerated().map { ($0.element, $0.offset) })
        return selectIndependentCycles(candidates, rank: rank, edgeIndex: edgeIndex)
    }

    static func edgeKeys(for ring: [Int]) -> [CDKEdgeKey] {
        guard ring.count >= 2 else { return [] }
        return (0..<ring.count).map { i in
            CDKEdgeKey(ring[i], ring[(i + 1) % ring.count])
        }
    }

    private static func shortestCycleCandidates(graph: CDKMolecularGraph, maxCycleSize: Int)
        -> [[Int]]
    {
        let edges = Array(graph.edgeBonds.keys).sorted { lhs, rhs in
            if lhs.a != rhs.a { return lhs.a < rhs.a }
            return lhs.b < rhs.b
        }

        var out: [[Int]] = []
        for edge in edges {
            guard
                let path = shortestPathExcludingEdge(
                    graph: graph,
                    start: edge.a,
                    end: edge.b,
                    excluding: edge,
                    maxDepth: maxCycleSize - 1),
                path.count >= 3,
                path.count <= maxCycleSize
            else { continue }
            out.append(canonicalCycle(path))
        }
        return out
    }

    private static func fundamentalCycleCandidates(graph: CDKMolecularGraph) -> [[Int]] {
        var visited: Set<Int> = []
        var parent: [Int: Int] = [:]
        var depth: [Int: Int] = [:]
        var out: [[Int]] = []
        let atomIDs = graph.atomIDs.sorted()

        func pathToAncestor(from start: Int, ancestor: Int) -> [Int]? {
            var path: [Int] = [start]
            var current = start
            while current != ancestor {
                guard let p = parent[current] else { return nil }
                current = p
                path.append(current)
            }
            return path
        }

        func dfs(_ current: Int) {
            visited.insert(current)
            let neighbors = graph.neighbors(of: current).sorted()
            for next in neighbors {
                if parent[current] == next { continue }
                if !visited.contains(next) {
                    parent[next] = current
                    depth[next] = (depth[current] ?? 0) + 1
                    dfs(next)
                    continue
                }

                let curDepth = depth[current] ?? 0
                let nextDepth = depth[next] ?? 0
                guard nextDepth < curDepth else { continue }
                guard let chain = pathToAncestor(from: current, ancestor: next), chain.count >= 3 else {
                    continue
                }
                out.append(canonicalCycle(chain))
            }
        }

        for seed in atomIDs where !visited.contains(seed) {
            parent[seed] = seed
            depth[seed] = 0
            dfs(seed)
            parent[seed] = nil
        }

        return out
    }

    private static func shortestPathExcludingEdge(
        graph: CDKMolecularGraph,
        start: Int,
        end: Int,
        excluding: CDKEdgeKey,
        maxDepth: Int
    ) -> [Int]? {
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

    private static func selectIndependentCycles(
        _ cycles: [[Int]],
        rank: Int,
        edgeIndex: [CDKEdgeKey: Int]
    ) -> [[Int]] {
        guard rank > 0 else { return [] }

        typealias Row = (pivot: Int, bits: Set<Int>)
        var basis: [Row] = []
        var selected: [[Int]] = []
        selected.reserveCapacity(rank)

        func reduce(_ bitset: Set<Int>) -> Set<Int> {
            var out = bitset
            for row in basis {
                if out.contains(row.pivot) {
                    out.formSymmetricDifference(row.bits)
                }
            }
            return out
        }

        for cycle in cycles {
            let edges = edgeKeys(for: cycle)
            let bits = Set(edges.compactMap { edgeIndex[$0] })
            if bits.isEmpty { continue }

            let reduced = reduce(bits)
            guard let pivot = reduced.max() else { continue }

            basis.append((pivot: pivot, bits: reduced))
            basis.sort { $0.pivot > $1.pivot }
            selected.append(cycle)

            if selected.count >= rank {
                break
            }
        }

        return selected
    }

    private static func edgeLess(_ lhs: CDKEdgeKey, _ rhs: CDKEdgeKey) -> Bool {
        if lhs.a != rhs.a { return lhs.a < rhs.a }
        return lhs.b < rhs.b
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
