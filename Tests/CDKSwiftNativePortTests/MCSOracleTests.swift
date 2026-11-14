import XCTest
@testable import CDKSwiftNativePort

final class MCSOracleTests: XCTestCase {
    // Independent exhaustive oracle: enumerate vertex subsets and target permutations,
    // then build the common-edge matrix and connected components from scratch.
    func testAllFourVertexTopologiesAgainstIndependentSubsetOracle() throws {
        let edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        func graph(_ mask: Int, nitrogen: Bool) -> Molecule {
            let atoms = (0..<4).map {
                Atom(id: $0 * 11 - 23, element: nitrogen && $0 == 3 ? "N" : "C", position: .zero)
            }
            let bonds = edges.enumerated().compactMap { index, edge -> Bond? in
                guard mask & (1 << index) != 0 else { return nil }
                return Bond(id: index * 7 + 1, a1: atoms[edge.0].id, a2: atoms[edge.1].id, order: .single)
            }
            return Molecule(atoms: atoms, bonds: bonds)
        }
        for sourceMask in 0..<64 {
            for targetMask in [0, 13, 45, 63] {
                let query = graph(sourceMask, nitrogen: false)
                let target = graph(targetMask, nitrogen: true)
                for disconnected in [false, true] {
                    let options = CDKMCSOptions(minimumBondsPerFragment: disconnected ? 1 : nil)
                    let actual = try CDKMCS(query: query, options: options).search(in: target)
                    let expected = oracle(query, target, disconnected: disconnected)
                    XCTAssertTrue(actual.isMaximumProven)
                    XCTAssertEqual(
                        [actual.mappedAtomCount, actual.mappedBondCount], expected,
                        "query mask \(sourceMask), target mask \(targetMask), disconnected \(disconnected)")
                    for witness in actual.mappings {
                        XCTAssertEqual(witness.atomIDByQueryAtomID.count, expected[0])
                        XCTAssertEqual(witness.bondIDByQueryBondID.count, expected[1])
                        XCTAssertEqual(Set(witness.atomIDByQueryAtomID.values).count, expected[0])
                    }
                }
            }
        }
    }

    private func oracle(_ query: Molecule, _ target: Molecule, disconnected: Bool) -> [Int] {
        var best = [0, 0]
        for mask in 1..<(1 << query.atoms.count) {
            let selected = query.atoms.indices.filter { mask & (1 << $0) != 0 }
            func assign(_ targetIndices: [Int]) {
                if targetIndices.count < selected.count {
                    for index in target.atoms.indices where !targetIndices.contains(index) {
                        assign(targetIndices + [index])
                    }
                    return
                }
                guard zip(selected, targetIndices).allSatisfy({ query.atoms[$0].element == target.atoms[$1].element })
                else { return }
                let n = selected.count
                var labels = Array(0..<n)
                var bonds = 0
                var degrees = Array(repeating: 0, count: n)
                for a in 0..<n {
                    for b in (a + 1)..<n {
                        if let left = query.bond(
                            between: query.atoms[selected[a]].id, and: query.atoms[selected[b]].id),
                            let right = target.bond(
                                between: target.atoms[targetIndices[a]].id, and: target.atoms[targetIndices[b]].id),
                            left.order == right.order
                        {
                            bonds += 1; degrees[a] += 1; degrees[b] += 1
                            let old = labels[b], replacement = labels[a]
                            labels = labels.map { $0 == old ? replacement : $0 }
                        }
                    }
                }
                if disconnected {
                    guard degrees.allSatisfy({ $0 > 0 }) else { return }
                    if bonds > best[1] || (bonds == best[1] && n > best[0]) { best = [n, bonds] }
                } else if Set(labels).count == 1 && (n > best[0] || (n == best[0] && bonds > best[1])) {
                    best = [n, bonds]
                }
            }
            assign([])
        }
        return best
    }

    func testAtomPermutationsSourceNonmutationAndSymmetricScores() throws {
        let parser = CDKSmilesParser(generateCoordinates: false)
        let source = try parser.parseSmiles("CC(=O)NC")
        var target = try parser.parseSmiles("CCC(=O)N")
        target.atoms.reverse(); target.bonds.reverse()
        let saved = target
        let forward = try CDKMCS(query: source).search(in: target)
        let reverse = try CDKMCS(query: target).search(in: source)
        XCTAssertEqual(
            [forward.mappedAtomCount, forward.mappedBondCount], [reverse.mappedAtomCount, reverse.mappedBondCount])
        XCTAssertEqual(target, saved)
        let encoder = JSONEncoder()
        let decoded = try JSONDecoder().decode(
            CDKSubstructureMatch.self, from: encoder.encode(XCTUnwrap(forward.mappings.first)))
        XCTAssertEqual(decoded, forward.mappings.first)
    }

    func testLargerBoundedSearchPreservesWitnessAndActiveCancellation() async throws {
        let parser = CDKSmilesParser(generateCoordinates: false)
        let query = try parser.parseSmiles(String(repeating: "C", count: 18))
        let target = try parser.parseSmiles("C1" + String(repeating: "C", count: 19) + "1")
        let result = try CDKMCS(query: query, options: .init(maximumMappings: 2, maximumSearchSteps: 5_000)).search(
            in: target)
        XCTAssertEqual(result.termination, .searchLimit)
        XCTAssertFalse(result.isMaximumProven)
        XCTAssertFalse(result.mappings.isEmpty)
        let task = Task.detached {
            try CDKMCS(query: query, options: .init(maximumSearchSteps: 10_000_000, timeoutSeconds: 5)).search(
                in: target)
        }
        try await Task.sleep(for: .milliseconds(10))
        task.cancel()
        let cancelled = try await task.value
        XCTAssertEqual(cancelled.termination, .cancelled)
        XCTAssertFalse(cancelled.isMaximumProven)
    }
}
