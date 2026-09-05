import XCTest
@testable import CDKSwiftNativePort

final class DescriptorRingSupportTests: XCTestCase {
    // The former edge-shortest-path candidate set has only four cycles for
    // this rank-five graph, which triggered exhaustive 24-atom enumeration.
    private let deficientEdges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5),
                                  (0, 2), (2, 5), (1, 4), (0, 4), (3, 5)]

    func testDeficientShortestEdgeCandidatesProduceCompleteMinimumBasis() {
        let molecule = makeMolecule(atomCount: 6, edges: deficientEdges)
        let basis = CDKDescriptorRingSupport.smallestRingBasis(in: molecule)
        XCTAssertEqual(basis.count, 5)
        XCTAssertEqual(basis.map(\.size).sorted(), [3, 3, 3, 3, 4])
        assertValidIndependentBasis(basis, in: molecule)
    }

    func testHighlyCyclicDisconnectedGraphDoesNotEnumerateAllSimpleCycles() {
        var edges = deficientEdges
        for y in 0..<8 {
            for x in 0..<8 {
                let atom = 6 + y * 8 + x
                if x < 7 { edges.append((atom, atom + 1)) }
                if y < 7 { edges.append((atom, atom + 8)) }
            }
        }
        let molecule = makeMolecule(atomCount: 70, edges: edges)
        let basis = CDKDescriptorRingSupport.smallestRingBasis(in: molecule)
        XCTAssertEqual(basis.count, 54)
        XCTAssertEqual(basis.map(\.size).reduce(0, +), 49 * 4 + 16)
        assertValidIndependentBasis(basis, in: molecule)
    }

    func testMinimumBasisMatchesExhaustiveOracleForEveryFiveAtomGraph() {
        let possibleEdges = (0..<5).flatMap { a in ((a + 1)..<5).map { (a, $0) } }
        for mask in 0..<(1 << possibleEdges.count) {
            let edges = possibleEdges.enumerated().compactMap { index, edge in
                mask & (1 << index) == 0 ? nil : edge
            }
            let molecule = makeMolecule(atomCount: 5, edges: edges)
            let basis = CDKDescriptorRingSupport.smallestRingBasis(in: molecule)
            // Exhaustive enumeration is deliberately confined to five vertices.
            let cycles = molecule.simpleCycles(maxSize: 5).map { atoms in
                Set(atoms.indices.compactMap { index in
                    molecule.bond(between: atoms[index], and: atoms[(index + 1) % atoms.count])?.id
                })
            }
            let oracle = independentCycles(cycles.sorted { $0.count < $1.count })
            XCTAssertEqual(basis.count, oracle.count, "graph mask \(mask)")
            XCTAssertEqual(basis.map(\.size).reduce(0, +), oracle.map(\.count).reduce(0, +), "graph mask \(mask)")
            assertValidIndependentBasis(basis, in: molecule)
        }
    }

    func testMacrocycleAndHydrogenDoNotTruncateTheBasis() {
        let edges = (0..<30).map { ($0, ($0 + 1) % 30) }
        var molecule = makeMolecule(atomCount: 31, edges: edges + [(0, 30)])
        molecule.atoms[30].element = "H"
        let basis = CDKDescriptorRingSupport.smallestRingBasis(in: molecule)
        XCTAssertEqual(basis.map(\.size), [30])
        XCTAssertFalse(basis.flatMap(\.atoms).contains(30))
    }

    private func assertValidIndependentBasis(_ basis: [CDKDescriptorRingSupport.Cycle], in molecule: Molecule,
                                            file: StaticString = #filePath, line: UInt = #line) {
        XCTAssertEqual(independentCycles(basis.map(\.bondIDs)).count, basis.count, file: file, line: line)
        for cycle in basis {
            XCTAssertEqual(Set(cycle.atoms).count, cycle.size, file: file, line: line)
            for index in cycle.atoms.indices {
                let bond = molecule.bond(between: cycle.atoms[index], and: cycle.atoms[(index + 1) % cycle.atoms.count])
                XCTAssertTrue(bond.map { cycle.bondIDs.contains($0.id) } ?? false, file: file, line: line)
            }
        }
    }

    private func independentCycles(_ cycles: [Set<Int>]) -> [Set<Int>] {
        var rows: [Int: Set<Int>] = [:]
        var chosen: [Set<Int>] = []
        for cycle in cycles {
            var remainder = cycle
            while let pivot = remainder.max(), let row = rows[pivot] {
                remainder.formSymmetricDifference(row)
            }
            if let pivot = remainder.max() {
                rows[pivot] = remainder
                chosen.append(cycle)
            }
        }
        return chosen
    }

    private func makeMolecule(atomCount: Int, edges: [(Int, Int)]) -> Molecule {
        Molecule(atoms: (0..<atomCount).map { Atom(id: $0, element: "C", position: .zero) },
                 bonds: edges.enumerated().map { id, edge in
                     Bond(id: id, a1: edge.0, a2: edge.1, order: .single)
                 })
    }
}
