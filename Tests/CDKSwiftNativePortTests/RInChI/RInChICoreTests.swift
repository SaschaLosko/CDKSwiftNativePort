import Foundation
import XCTest
@testable import CDKSwiftNativePort

final class RInChICoreTests: XCTestCase {
    func testInChIGeneratorExposesAuxInfo() throws {
        let molecule = try CDKSmilesParser().parseSmiles("CCO")
        let generator = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)

        XCTAssertEqual(generator.getStatus(), .success)
        XCTAssertTrue(try generator.getAuxInfo().hasPrefix("AuxInfo=1/"))
    }

    func testGeneratorReportsNullReactionInput() {
        let generator = CDKRInChIGenerator().generate(nil)
        XCTAssertEqual(generator.getStatus(), .error)
        XCTAssertEqual(generator.getMessages(), ["CDKReaction object provided as input is 'null'."])
    }

    func testDecompositionExtractsRolesAndDirection() {
        let decomposition = CDKRInChIDecomposition(
            rinchi: "RInChI=1.00.1S/C3H6O/c1-3(2)4/h1-2H3<>C4H10/c1-4(2)3/h4H,1-3H3<>H2O/h1H2/d="
        ).decompose()

        XCTAssertEqual(decomposition.getStatus(), .success)
        XCTAssertEqual(decomposition.getReactionDirection(), .bidirectional)
        XCTAssertEqual(
            decomposition.getComponents(),
            [
                CDKRInChIDecompositionComponent(
                    inchi: "InChI=1S/C3H6O/c1-3(2)4/h1-2H3",
                    auxInfo: "",
                    reactionRole: .reactant
                ),
                CDKRInChIDecompositionComponent(
                    inchi: "InChI=1S/C4H10/c1-4(2)3/h4H,1-3H3",
                    auxInfo: "",
                    reactionRole: .product
                ),
                CDKRInChIDecompositionComponent(
                    inchi: "InChI=1S/H2O/h1H2",
                    auxInfo: "",
                    reactionRole: .agent
                ),
            ]
        )
    }

    func testRInChIToReactionRoundTripsReferenceExample() throws {
        let input = "RInChI=1.00.1S/CH4/h1H4!ClH/h1H<>CH4O/c1-2/h2H,1H3!ClH/h1H<>H2O/h1H2/d+"
        let reaction = try XCTUnwrap(CDKRInChIToReaction(rinchi: input).getReaction())

        XCTAssertEqual(reaction.reactantCount, 2)
        XCTAssertEqual(reaction.productCount, 2)
        XCTAssertEqual(reaction.agentCount, 1)

        let regenerated = CDKRInChIGenerator().generate(reaction)
        XCTAssertEqual(regenerated.getStatus(), .success)
        XCTAssertEqual(regenerated.getRInChI(), input)
    }

    func testForceEquilibriumAndReverseChangeOnlyDirectionLayer() throws {
        let parser = CDKSmilesParser()
        let reaction = try parser.parseReactionSmiles("C=C.C=CC=C>>C1CCCCC1")
        let forward = CDKRInChIGenerator().generate(reaction)
        let equilibrium = CDKRInChIGenerator(options: CDKRInChIOptions(forceEquilibrium: true)).generate(reaction)
        let reversed = CDKRInChIGenerator().generate(CDKReactionManipulator.reverse(reaction))

        let forwardRInChI = try XCTUnwrap(forward.getRInChI())
        let equilibriumRInChI = try XCTUnwrap(equilibrium.getRInChI())
        let reversedRInChI = try XCTUnwrap(reversed.getRInChI())

        XCTAssertTrue(forwardRInChI.hasSuffix("/d+"))
        XCTAssertTrue(equilibriumRInChI.hasSuffix("/d="))
        XCTAssertTrue(reversedRInChI.hasSuffix("/d-"))

        XCTAssertEqual(String(forwardRInChI.dropLast()), String(equilibriumRInChI.dropLast()))
        XCTAssertEqual(String(forwardRInChI.dropLast()), String(reversedRInChI.dropLast()))
    }
}
