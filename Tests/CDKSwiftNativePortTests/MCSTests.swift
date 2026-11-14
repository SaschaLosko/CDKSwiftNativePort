import XCTest
@testable import CDKSwiftNativePort

final class MCSTests: XCTestCase {
    private func molecule(_ text: String) throws -> Molecule {
        try CDKSmilesParser(generateCoordinates: false).parseSmiles(text)
    }
    private func search(_ query: String, _ target: String, _ options: CDKMCSOptions = .init()) throws -> CDKMCSResult {
        try CDKMCS(query: molecule(query), options: options).search(in: molecule(target))
    }

    func testUpstreamAmideExampleAndNonInducedCommonEdges() throws {
        // CDK 2.14.1 MCS.java API example; no upstream runtime is required.
        let query = try molecule("CC(=O)NC"), target = try molecule("CCC(=O)N")
        let result = try CDKMCS(query: query).search(in: target)
        XCTAssertTrue(result.isMaximumProven)
        XCTAssertEqual(result.mappedAtomCount, 4)
        XCTAssertEqual(result.mappedBondCount, 3)
        XCTAssertEqual(result.mappings.count, 1)
        let expected = Dictionary(uniqueKeysWithValues: (0..<4).map { (query.atoms[$0].id, target.atoms[$0 + 1].id) })
        XCTAssertEqual(result.mappings.first?.atomIDByQueryAtomID, expected)
        let broken = try search("C1CC1", "CCC")
        XCTAssertTrue(broken.isMaximumProven)
        XCTAssertEqual(broken.mappedAtomCount, 3)
        XCTAssertEqual(broken.mappedBondCount, 2)
    }

    func testMappingLimitDoesNotPreventMaximumProof() throws {
        let all = try search("C1CC1", "C1CC1")
        XCTAssertEqual(all.mappings.count, 6)
        XCTAssertTrue(all.allMaximumMappingsReturned)
        let limited = try search("C1CC1", "C1CC1", .init(maximumMappings: 2))
        XCTAssertTrue(limited.isMaximumProven)
        XCTAssertTrue(limited.mappingsTruncated)
        XCTAssertFalse(limited.allMaximumMappingsReturned)
        XCTAssertEqual(limited.mappings.count, 2)
        XCTAssertEqual(Set(limited.mappings).count, 2)
    }

    func testDisconnectedUsesCommonBondMinimumPerFragment() throws {
        let connected = try search("CCNCC", "CCOCC")
        XCTAssertEqual(connected.mappedAtomCount, 2)
        let disconnected = try search("CCNCC", "CCOCC", .init(minimumBondsPerFragment: 1))
        XCTAssertTrue(disconnected.isMaximumProven)
        XCTAssertEqual(disconnected.mappedAtomCount, 4)
        XCTAssertEqual(disconnected.mappedBondCount, 2)
        let tooSmall = try search("CCNCC", "CCOCC", .init(minimumBondsPerFragment: 2))
        XCTAssertTrue(tooSmall.isMaximumProven)
        XCTAssertTrue(tooSmall.mappings.isEmpty)
        let isolated = try search("C.C", "C.C", .init(minimumBondsPerFragment: 1))
        XCTAssertTrue(isolated.mappings.isEmpty)
    }

    func testCompleteRelevantRingsAndChainControls() throws {
        let options = CDKMCSOptions(completeRings: true, maximumSearchSteps: 5_000_000)
        let broken = try search("C1CC1", "CCC", options)
        XCTAssertTrue(broken.isMaximumProven)
        XCTAssertEqual(broken.mappedAtomCount, 0)
        let fused = try search("c1ccccc1", "c1ccc2ccccc2c1", options)
        XCTAssertTrue(fused.isMaximumProven)
        XCTAssertEqual(fused.mappedAtomCount, 6)
        XCTAssertEqual(fused.mappedBondCount, 6)
        let bridged = try search("C1CCCCC1", "C1CC2CCC1C2", options)
        XCTAssertTrue(bridged.isMaximumProven)
        XCTAssertEqual(bridged.mappedAtomCount, 0)
        let substituent = try search("CC1CC1", "CCC", options)
        XCTAssertEqual(substituent.mappedAtomCount, 1)
        let spiro = try search("C1CC1", "C1CC12CC2", options)
        XCTAssertEqual(spiro.mappedAtomCount, 3)
        // The six-membered perimeter of norbornane is not a relevant ring.
        let budget = MarkushMatchBudget(.init(maximumSearchSteps: 1_000_000))
        var graph = try MCSGraph(molecule("C1CC2CCC1C2"))
        try graph.prepareRings(complete: true, budget: budget)
        XCTAssertEqual(graph.relevantRings.map(\.count).sorted(), [5, 5])
        // Equal paths in bicyclo[1.1.1]pentane yield three relevant rings, although rank is two.
        var theta = try MCSGraph(molecule("C1C2CC1C2"))
        try theta.prepareRings(complete: true, budget: budget)
        XCTAssertEqual(theta.relevantRings.map(\.count).sorted(), [4, 4, 4])
    }

    func testTetrahedralStereoIsPartialMappingConstraintAndUnspecifiedIsCompatible() throws {
        let options = CDKMCSOptions(matchStereochemistry: true)
        XCTAssertEqual(try search("F[C@](Cl)(Br)I", "F[C@](Cl)(Br)I", options).mappedAtomCount, 5)
        XCTAssertEqual(try search("F[C@](Cl)(Br)I", "F[C@@](Cl)(Br)I", options).mappedAtomCount, 3)
        XCTAssertEqual(try search("F[C@](Cl)(Br)I", "FC(Cl)(Br)I", options).mappedAtomCount, 5)
        XCTAssertEqual(try search("FC(Cl)(Br)I", "F[C@@](Cl)(Br)I", options).mappedAtomCount, 5)
        XCTAssertEqual(try search("F[C@](Cl)(Br)I", "F[C@@](Cl)(Br)I").mappedAtomCount, 5)
        var permuted = try molecule("F[C@](Cl)(Br)I")
        permuted.atoms.reverse(); permuted.bonds.reverse()
        XCTAssertEqual(
            try CDKMCS(query: molecule("F[C@](Cl)(Br)I"), options: options).search(in: permuted).mappedAtomCount, 5)
    }

    func testAlkeneStereoUnspecifiedAndAlternativeCarriers() throws {
        let options = CDKMCSOptions(matchStereochemistry: true)
        XCTAssertEqual(try search("F/C=C/Cl", "F/C=C/Cl", options).mappedAtomCount, 4)
        XCTAssertEqual(try search("F/C=C/Cl", "F/C=C\\Cl", options).mappedAtomCount, 3)
        XCTAssertEqual(try search("F/C=C/Cl", "FC=CCl", options).mappedAtomCount, 4)
        XCTAssertEqual(try search("FC=CCl", "F/C=C/Cl", options).mappedAtomCount, 4)
        XCTAssertEqual(try search("F/C=C/Cl", "Cl/C=C/F", options).mappedAtomCount, 4)
    }

    func testRelativeGroupsInvertTogether() throws {
        let options = CDKMCSOptions(matchStereochemistry: true, maximumSearchSteps: 2_000_000)
        let left = "F[C@H](Cl)[C@H](Br)I |&1:1,3|"
        let inverted = "F[C@@H](Cl)[C@@H](Br)I"
        let mixed = "F[C@@H](Cl)[C@H](Br)I"
        XCTAssertEqual(try search(left, inverted, options).mappedAtomCount, 6)
        XCTAssertLessThan(try search(left, mixed, options).mappedAtomCount, 6)
        XCTAssertEqual(try search("F[C@H](Cl)[C@H](Br)I |r|", inverted, options).mappedAtomCount, 6)
    }

    func testExplicitHydrogensCountAndIsotopeChargePolicies() throws {
        XCTAssertEqual(try search("[H]C([H])([H])[H]", "[H]C([H])([H])[H]").mappedAtomCount, 5)
        XCTAssertEqual(try search("[13CH3]O", "CO").mappedAtomCount, 2)
        XCTAssertEqual(try search("[13CH3]O", "CO", .init(matchIsotopes: true)).mappedAtomCount, 1)
        XCTAssertEqual(try search("C[NH3+]", "CN", .init(matchFormalCharge: true)).mappedAtomCount, 1)
        XCTAssertEqual(try search("*O", "CO").mappedAtomCount, 1)
        XCTAssertEqual(try search("*O", "CO", .init(atomMatching: .element)).mappedAtomCount, 2)
        XCTAssertEqual(try search("CO", "C=O", .init(ignoreBondOrder: true)).mappedAtomCount, 2)
    }

    func testQueryPredicatesAndRecursiveAnchoring() throws {
        let query = try CDKSmartsParser.parse("C(=O)[N,O]")
        let result = try CDKMCS(query: query).search(in: molecule("CC(=O)O"))
        XCTAssertTrue(result.isMaximumProven)
        XCTAssertEqual(result.mappedAtomCount, 3)
        let stereoPattern = try CDKSmartsParser.parse("N[C@@H](C)C(=O)O")
        XCTAssertEqual(try CDKMCS(query: stereoPattern).search(in: molecule("N[C@H](C)C(=O)O")).mappedAtomCount, 6)
        XCTAssertThrowsError(
            try CDKMCS(query: query, options: .init(ignoreBondOrder: true)).search(in: molecule("CC=O")))
        let recursive = try CDKSmartsParser.parse("[$(C=O)]")
        XCTAssertEqual(try CDKMCS(query: recursive).search(in: molecule("CC=O")).mappedAtomCount, 1)
        XCTAssertEqual(try CDKMCS(query: recursive).search(in: molecule("CCO")).mappedAtomCount, 0)
        XCTAssertThrowsError(
            try CDKMCS(query: query, options: .init(matchStereochemistry: true)).search(in: molecule("CC=O")))
    }

    func testBudgetTimeoutCancellationAndEmptyGraphs() async throws {
        let query = try molecule("C1CCCCC1")
        let target = try molecule("C1CCC2CCCCC2C1")
        let limited = try CDKMCS(query: query, options: .init(maximumSearchSteps: 100)).search(in: target)
        XCTAssertEqual(limited.termination, .searchLimit)
        XCTAssertFalse(limited.isMaximumProven)
        XCTAssertThrowsError(try CDKMCS(query: query, options: .init(maximumSearchSteps: 1)).match(in: target))
        let timed = try CDKMCS(query: query, options: .init(timeoutSeconds: 1e-12)).search(in: target)
        XCTAssertEqual(timed.termination, .timedOut)
        let task = Task {
            withUnsafeCurrentTask { $0?.cancel() }
            return try CDKMCS(query: query).search(in: target)
        }
        let cancelled = try await task.value
        XCTAssertEqual(cancelled.termination, .cancelled)
        let empty = try CDKMCS(query: Molecule()).search(in: target)
        XCTAssertTrue(empty.isMaximumProven)
        XCTAssertTrue(empty.mappings.isEmpty)
        XCTAssertThrowsError(try CDKMCS(query: query, options: .init(minimumBondsPerFragment: 0)).search(in: target))
    }
}
