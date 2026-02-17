import XCTest
@testable import CDKPortExperimental

final class SmilesReactionParserPortTests: XCTestCase {

    private let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    // Mirrors CDK SmilesParserTest.testReaction.
    func testReactionCounts() throws {
        let reaction = try parser.parseReactionSmiles("O>>[H+].[OH-]")
        XCTAssertEqual(reaction.reactantCount, 1)
        XCTAssertEqual(reaction.productCount, 2)
    }

    // Mirrors CDK SmilesParserTest.noReactants/noProducts/noReaction/onlyAgents.
    func testReactionWithMissingSides() throws {
        let noReactants = try parser.parseReactionSmiles(">>C")
        XCTAssertEqual(noReactants.reactantCount, 0)
        XCTAssertEqual(noReactants.productCount, 1)

        let noProducts = try parser.parseReactionSmiles("C>>")
        XCTAssertEqual(noProducts.reactantCount, 1)
        XCTAssertEqual(noProducts.productCount, 0)

        let noReaction = try parser.parseReactionSmiles(">>")
        XCTAssertEqual(noReaction.reactantCount, 0)
        XCTAssertEqual(noReaction.productCount, 0)

        let onlyAgents = try parser.parseReactionSmiles(">C>")
        XCTAssertEqual(onlyAgents.reactantCount, 0)
        XCTAssertEqual(onlyAgents.agentCount, 1)
        XCTAssertEqual(onlyAgents.productCount, 0)
    }

    // Mirrors CDK SmilesParserTest.testReactionWithAgents.
    func testReactionWithAgents() throws {
        let reaction = try parser.parseReactionSmiles("CCO.CC(=O)O>[H+]>CC(=O)OCC.O")
        XCTAssertEqual(reaction.reactantCount, 2)
        XCTAssertEqual(reaction.agentCount, 1)
        XCTAssertEqual(reaction.productCount, 2)
        XCTAssertEqual(reaction.agents.first?.atomCount, 1)
    }

    // Mirrors CxSmilesParserTest.atomOrderingWithNonContiguousFragments.
    func testReactionFragmentGroupingAndAtomLabels() throws {
        let reaction = try parser.parseReactionSmiles("C.*.C>> |$;R1;$,f:0.2|")
        XCTAssertEqual(reaction.reactantCount, 2)
        XCTAssertEqual(reaction.productCount, 0)

        let atomCounts = reaction.reactants.map(\.atomCount).sorted()
        XCTAssertEqual(atomCounts, [1, 2])

        let pseudo = try XCTUnwrap(reaction.reactants.first(where: { $0.atomCount == 1 }))
        XCTAssertEqual(pseudo.atoms.first?.element, "R1")
    }

    // Mirrors CxSmilesParserTest.loadRacComponentsWithFragGrouping grouping aspect.
    func testReactionCxFragmentGroupingAgents() throws {
        let reaction = try parser.parseReactionSmiles("c1ccccc1[C@H](O)C>[Na+].[Cl-]>CC[C@H](O)[C@H](O)CCCCCC |f:1.2,r:3|")
        XCTAssertEqual(reaction.reactantCount, 1)
        XCTAssertEqual(reaction.agentCount, 1)
        XCTAssertEqual(reaction.productCount, 1)
        XCTAssertEqual(reaction.agents.first?.atomCount, 2)
    }

    // Mirrors CxSmilesParserTest.loadRacComponents shape with r:1 on products.
    func testReactionCxRacemicProductSelection() throws {
        let reaction = try parser.parseReactionSmiles("c1ccccc1[C@H](O)C>>CC[C@H](O)[C@H](O)CCCCCC |r:1|")
        XCTAssertEqual(reaction.reactantCount, 1)
        XCTAssertEqual(reaction.productCount, 1)
    }

    func testRejectsCxFragmentGroupingAcrossReactionSides() {
        XCTAssertThrowsError(try parser.parseReactionSmiles("C>>C |f:0.1|"))
    }

    func testRejectsMalformedReactionSyntax() {
        XCTAssertThrowsError(try parser.parseReactionSmiles("C>C"))
        XCTAssertThrowsError(try parser.parseReactionSmiles("C>>>O"))
    }
}
