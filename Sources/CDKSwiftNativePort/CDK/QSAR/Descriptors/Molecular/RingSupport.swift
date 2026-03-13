import Foundation

enum CDKDescriptorRingSupport {
    struct Cycle: Hashable {
        let atoms: [Int]
        let bondIDs: Set<Int>

        var size: Int {
            bondIDs.count
        }

        var sortedBondIDs: [Int] {
            bondIDs.sorted()
        }
    }

    static func aromaticBondIDs(in molecule: Molecule, perceiveAromaticity: Bool) -> Set<Int> {
        var bondIDs = Set(molecule.bonds.filter { $0.order == .aromatic }.map(\.id))
        guard perceiveAromaticity else { return bondIDs }

        for ring in molecule.aromaticDisplayRings() {
            for index in ring.indices {
                let atomA = ring[index]
                let atomB = ring[(index + 1) % ring.count]
                if let bond = molecule.bond(between: atomA, and: atomB) {
                    bondIDs.insert(bond.id)
                }
            }
        }
        return bondIDs
    }

    static func aromaticAtomIDs(in molecule: Molecule, perceiveAromaticity: Bool) -> Set<Int> {
        var atomIDs = Set(molecule.atoms.filter(\.aromatic).map(\.id))
        guard perceiveAromaticity else { return atomIDs }

        for ring in molecule.aromaticDisplayRings() {
            atomIDs.formUnion(ring)
        }
        return atomIDs
    }

    static func hasThreeMemberedRing(containing atomID: Int, in molecule: Molecule) -> Bool {
        let neighbors = molecule.neighbors(of: atomID)
        guard neighbors.count >= 2 else { return false }

        for firstIndex in 0..<(neighbors.count - 1) {
            for secondIndex in (firstIndex + 1)..<neighbors.count {
                if molecule.bond(between: neighbors[firstIndex], and: neighbors[secondIndex]) != nil {
                    return true
                }
            }
        }
        return false
    }

    static func smallestRingBasis(in molecule: Molecule, maximumCycleSize: Int = 24) -> [Cycle] {
        let heavyAtomIDs = CDKDescriptorSupport.heavyAtomIDs(in: molecule)
        guard !heavyAtomIDs.isEmpty else { return [] }

        let heavyBonds = molecule.bonds.filter { heavyAtomIDs.contains($0.a1) && heavyAtomIDs.contains($0.a2) }
        let cycleRank = max(0, heavyBonds.count - heavyAtomIDs.count + CDKDescriptorSupport.heavyAtomConnectedComponentCount(in: molecule))
        guard cycleRank > 0 else { return [] }

        let adjacency = heavyAdjacency(in: molecule, heavyAtomIDs: heavyAtomIDs)
        var candidatesByBondSet: [Set<Int>: Cycle] = [:]

        for bond in heavyBonds {
            guard let path = shortestPath(from: bond.a1,
                                          to: bond.a2,
                                          excludingBondID: bond.id,
                                          adjacency: adjacency) else {
                continue
            }

            let pathBondIDs = bondIDs(for: path, in: molecule)
            guard !pathBondIDs.isEmpty else { continue }
            let cycleBondIDs = pathBondIDs.union([bond.id])
            candidatesByBondSet[cycleBondIDs] = Cycle(atoms: path, bondIDs: cycleBondIDs)
        }

        if candidatesByBondSet.count < cycleRank {
            let cappedCycleSize = min(maximumCycleSize, max(3, heavyAtomIDs.count))
            for ring in molecule.simpleCycles(maxSize: cappedCycleSize) {
                guard ring.count >= 3 else { continue }
                guard ring.allSatisfy({ heavyAtomIDs.contains($0) }) else { continue }
                let ringBondIDs = bondIDs(in: ring, molecule: molecule)
                guard ringBondIDs.count == ring.count else { continue }
                candidatesByBondSet[ringBondIDs] = Cycle(atoms: ring, bondIDs: ringBondIDs)
            }
        }

        let sortedCandidates = candidatesByBondSet.values.sorted { lhs, rhs in
            if lhs.size != rhs.size { return lhs.size < rhs.size }
            if lhs.atoms.count != rhs.atoms.count { return lhs.atoms.count < rhs.atoms.count }
            return lhs.sortedBondIDs.lexicographicallyPrecedes(rhs.sortedBondIDs)
        }

        var selected: [Cycle] = []
        var basisRows: [Int: Set<Int>] = [:]

        for candidate in sortedCandidates {
            let reduced = reduce(candidate.bondIDs, by: basisRows)
            guard !reduced.isEmpty else { continue }
            guard let pivot = reduced.max() else { continue }

            selected.append(candidate)
            basisRows[pivot] = reduced

            if selected.count == cycleRank {
                break
            }
        }

        return selected
    }

    private static func heavyAdjacency(in molecule: Molecule,
                                       heavyAtomIDs: Set<Int>) -> [Int: [(neighbor: Int, bondID: Int)]] {
        var adjacency: [Int: [(neighbor: Int, bondID: Int)]] = [:]
        for atomID in heavyAtomIDs {
            adjacency[atomID] = []
        }
        for bond in molecule.bonds where heavyAtomIDs.contains(bond.a1) && heavyAtomIDs.contains(bond.a2) {
            adjacency[bond.a1, default: []].append((bond.a2, bond.id))
            adjacency[bond.a2, default: []].append((bond.a1, bond.id))
        }
        return adjacency.mapValues { entries in
            entries.sorted {
                if $0.neighbor != $1.neighbor { return $0.neighbor < $1.neighbor }
                return $0.bondID < $1.bondID
            }
        }
    }

    private static func shortestPath(from start: Int,
                                     to goal: Int,
                                     excludingBondID: Int,
                                     adjacency: [Int: [(neighbor: Int, bondID: Int)]]) -> [Int]? {
        guard start != goal else { return nil }

        var queue = [start]
        var visited: Set<Int> = [start]
        var previous: [Int: Int] = [:]
        var cursor = 0

        while cursor < queue.count {
            let current = queue[cursor]
            cursor += 1

            for entry in adjacency[current, default: []] where entry.bondID != excludingBondID {
                if !visited.insert(entry.neighbor).inserted {
                    continue
                }
                previous[entry.neighbor] = current
                if entry.neighbor == goal {
                    var path = [goal]
                    var step = goal
                    while let predecessor = previous[step] {
                        path.append(predecessor)
                        if predecessor == start {
                            return path.reversed()
                        }
                        step = predecessor
                    }
                    return nil
                }
                queue.append(entry.neighbor)
            }
        }

        return nil
    }

    private static func bondIDs(for atomPath: [Int], in molecule: Molecule) -> Set<Int> {
        guard atomPath.count >= 2 else { return [] }
        var bondIDs = Set<Int>()
        for index in 0..<(atomPath.count - 1) {
            guard let bond = molecule.bond(between: atomPath[index], and: atomPath[index + 1]) else {
                return []
            }
            bondIDs.insert(bond.id)
        }
        return bondIDs
    }

    private static func bondIDs(in cycle: [Int], molecule: Molecule) -> Set<Int> {
        guard cycle.count >= 3 else { return [] }
        var bondIDs = Set<Int>()
        for index in cycle.indices {
            let atomA = cycle[index]
            let atomB = cycle[(index + 1) % cycle.count]
            guard let bond = molecule.bond(between: atomA, and: atomB) else {
                return []
            }
            bondIDs.insert(bond.id)
        }
        return bondIDs
    }

    private static func reduce(_ vector: Set<Int>, by basisRows: [Int: Set<Int>]) -> Set<Int> {
        var remainder = vector
        for pivot in basisRows.keys.sorted(by: >) {
            guard remainder.contains(pivot), let basisVector = basisRows[pivot] else { continue }
            remainder.formSymmetricDifference(basisVector)
        }
        return remainder
    }
}
