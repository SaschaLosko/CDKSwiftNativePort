import XCTest
@testable import CDKSwiftNativePort

final class QueryVariableAttachmentTests: XCTestCase {
    private func molecule(_ text: String) throws -> Molecule {
        try CDKSmilesParser(generateCoordinates: false).parseSmiles(text)
    }

    func testDirectPositionVariationAndQueryBondOrder() throws {
        let family = try molecule("O*.c1ccccc1 |m:1:2.3|")
        let query = try CDKQueryPattern(molecule: family)
        let positive = CDKQueryMatcher.match(query, in: try molecule("Oc1ccccc1"))
        XCTAssertEqual(positive.verdict, .match, positive.detail ?? "")
        XCTAssertEqual(positive.atomMapping?.count, family.atoms.count)
        XCTAssertEqual(CDKQueryMatcher.match(query, in: try molecule("Nc1ccccc1")).verdict, .noMatch)
        var carbonyl = try CDKQueryPattern(molecule: molecule("O=*.CC |m:1:2.3|"))
        XCTAssertEqual(CDKQueryMatcher.match(carbonyl, in: try molecule("CC=O")).verdict, .match)
        XCTAssertEqual(CDKQueryMatcher.match(carbonyl, in: try molecule("CCO")).verdict, .noMatch)
        let variableBond = try XCTUnwrap(carbonyl.molecule.bonds.first)
        carbonyl.bonds[variableBond.id] = .or(.order(.single), .order(.double))
        XCTAssertEqual(CDKQueryMatcher.match(carbonyl, in: try molecule("CCO")).verdict, .match)
    }

    func testBothVariableEndsAndBacktrackingOverlappingChoices() throws {
        var graph = Molecule()
        graph.atoms =
            (0..<4).map { Atom(id: $0, element: "*", position: .zero) }
            + (4..<7).map { Atom(id: $0, element: "C", position: .zero, isotopeMassNumber: $0 + 9) }
        graph.bonds = [Bond(id: 0, a1: 0, a2: 1, order: .single), Bond(id: 1, a1: 2, a2: 3, order: .single)]
        graph.sgroups = [
            .init(kind: .extMulticenter, keyword: "m", atomIDs: [0, 4, 6]),
            .init(kind: .extMulticenter, keyword: "m", atomIDs: [1, 5]),
            .init(kind: .extMulticenter, keyword: "m", atomIDs: [2, 4]),
            .init(kind: .extMulticenter, keyword: "m", atomIDs: [3, 5]),
        ]
        let query = try CDKQueryPattern(molecule: graph)
        let target = try molecule("[13CH3][14CH2][15CH3]")
        let result = CDKQueryMatcher.match(query, in: target, mode: .wholeGraph)
        XCTAssertEqual(result.verdict, .match, result.detail ?? "")
        XCTAssertEqual(result.atomMapping?[0], target.atoms[2].id)
        XCTAssertEqual(result.atomMapping?[2], target.atoms[0].id)
        graph.sgroups[0].atomIDs = [0, 4]
        XCTAssertEqual(CDKQueryMatcher.match(try CDKQueryPattern(molecule: graph), in: target).verdict, .noMatch)
    }

    func testFixedAndVariableEdgesCannotConsumeSameBond() throws {
        var graph = try molecule("**.CN |m:0:2,m:1:3|")
        let query = try CDKQueryPattern(molecule: graph)
        XCTAssertEqual(CDKQueryMatcher.match(query, in: try molecule("CN")).verdict, .noMatch)
        graph.bonds.removeAll { $0.a1 == graph.atoms[2].id || $0.a2 == graph.atoms[2].id }
        XCTAssertEqual(
            CDKQueryMatcher.match(try CDKQueryPattern(molecule: graph), in: try molecule("CN")).verdict, .match)
    }

    func testInvalidMetadataAndBudgetsRemainUnresolved() throws {
        var graph = try molecule("O*.CC |m:1:2.3|")
        graph.cxState?.positionalVariations[1] = [999]
        let result = CDKQueryMatcher.match(try CDKQueryPattern(molecule: graph), in: try molecule("CCO"))
        XCTAssertEqual(result.verdict, .unresolved); XCTAssertEqual(result.reason, .invalidInput)
        let query = try CDKQueryPattern(molecule: molecule("O*.CC |m:1:2.3|"))
        let limited = CDKQueryMatcher.match(query, in: try molecule("CCO"), options: .init(maximumSearchSteps: 3))
        XCTAssertEqual(limited.verdict, .unresolved); XCTAssertFalse(limited.isExhaustive)
    }
}
