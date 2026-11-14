import Foundation

struct MCSGraph {
    let molecule: Molecule
    let atoms: [Int: Atom]
    let edges: [Int: [Int: Bond]]
    var ringBonds: Set<Int> = []
    var ringAtoms: Set<Int> = []
    var relevantRings: [Set<Int>] = []

    init(_ molecule: Molecule) throws {
        if !molecule.atoms.isEmpty {
            try QueryGraphFacts.validate(molecule)
        } else if !molecule.bonds.isEmpty {
            throw CDKMCSError.invalidInput("Bonds without atoms.")
        }
        self.molecule = molecule
        atoms = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0) })
        var edges: [Int: [Int: Bond]] = [:]
        for bond in molecule.bonds {
            edges[bond.a1, default: [:]][bond.a2] = bond
            edges[bond.a2, default: [:]][bond.a1] = bond
        }
        self.edges = edges
    }

    mutating func prepareRings(complete: Bool, budget: MarkushMatchBudget) throws {
        // An edge is cyclic exactly when its endpoints remain connected after removing it.
        for edge in molecule.bonds {
            var visited: Set<Int> = [edge.a1]
            var queue = [edge.a1]; var cursor = 0
            while cursor < queue.count && !visited.contains(edge.a2) {
                try budget.check()
                for (neighbor, bond) in edges[queue[cursor], default: [:]] where bond.id != edge.id {
                    if visited.insert(neighbor).inserted { queue.append(neighbor) }
                }
                cursor += 1
            }
            if visited.contains(edge.a2) { ringBonds.insert(edge.id); ringAtoms.formUnion([edge.a1, edge.a2]) }
        }
        guard complete else { return }
        // Enumerate simple cycles under the same search budget. A cycle is relevant iff
        // it is independent of all strictly shorter cycles (union of minimum cycle bases).
        var cycles: Set<Set<Int>> = []
        for root in ringAtoms.sorted() {
            var path: Set<Int> = [root]
            func walk(_ at: Int, _ bondIDs: [Int]) throws {
                try budget.check()
                for neighbor in edges[at, default: [:]].keys.sorted() {
                    let bond = edges[at]![neighbor]!
                    guard ringBonds.contains(bond.id) else { continue }
                    if neighbor == root && path.count >= 3 {
                        cycles.insert(Set(bondIDs + [bond.id]))
                    } else if neighbor > root && !path.contains(neighbor) {
                        path.insert(neighbor)
                        try walk(neighbor, bondIDs + [bond.id])
                        path.remove(neighbor)
                    }
                }
            }
            try walk(root, [])
        }
        let ordered = cycles.sorted { a, b in
            a.count == b.count ? a.sorted().lexicographicallyPrecedes(b.sorted()) : a.count < b.count
        }
        var basis: [Int: Set<Int>] = [:]
        func reduce(_ cycle: Set<Int>) -> Set<Int> {
            var vector = cycle
            while let pivot = vector.max(), let existing = basis[pivot] { vector.formSymmetricDifference(existing) }
            return vector
        }
        var cursor = 0
        while cursor < ordered.count {
            let count = ordered[cursor].count
            var end = cursor
            while end < ordered.count && ordered[end].count == count { end += 1 }
            for cycle in ordered[cursor..<end] {
                try budget.check()
                if !reduce(cycle).isEmpty { relevantRings.append(cycle) }
            }
            for cycle in ordered[cursor..<end] {
                let vector = reduce(cycle)
                if let pivot = vector.max() { basis[pivot] = vector }
            }
            cursor = end
        }
    }

    func ringsComplete(atomIDs: Set<Int>, bondIDs: Set<Int>) -> Bool {
        let whole = relevantRings.filter { $0.isSubset(of: bondIDs) }.reduce(into: Set<Int>()) { $0.formUnion($1) }
        guard bondIDs.intersection(ringBonds).isSubset(of: whole) else { return false }
        let wholeAtoms = Set(molecule.bonds.filter { whole.contains($0.id) }.flatMap { [$0.a1, $0.a2] })
        return atomIDs.intersection(ringAtoms).isSubset(of: wholeAtoms)
    }
}
