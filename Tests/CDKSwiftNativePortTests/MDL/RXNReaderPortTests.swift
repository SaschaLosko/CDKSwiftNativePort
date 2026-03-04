import XCTest
@testable import CDKSwiftNativePort

final class RXNReaderPortTests: XCTestCase {
    func testParsesBasicV2000ReactionBlock() throws {
        let text = """
        $RXN
        DemoReaction
          CDKSwiftNativePort

          1  1
        $MOL
        Reactant
          CDK

          1  0  0  0  0  0            999 V2000
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        M  END
        $MOL
        Product
          CDK

          1  0  0  0  0  0            999 V2000
            0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
        M  END
        """

        let molecules = try CDKRXNReader.read(text: text)
        XCTAssertEqual(molecules.count, 2)
        XCTAssertEqual(molecules[0].atoms.first?.element.uppercased(), "C")
        XCTAssertEqual(molecules[1].atoms.first?.element.uppercased(), "O")

        let reaction = try CDKRXNReader.readReaction(text: text)
        XCTAssertEqual(reaction.reactantCount, 1)
        XCTAssertEqual(reaction.productCount, 1)
        XCTAssertEqual(reaction.agentCount, 0)
        XCTAssertEqual(reaction.reactants.first?.atoms.first?.element.uppercased(), "C")
        XCTAssertEqual(reaction.products.first?.atoms.first?.element.uppercased(), "O")
    }

    func testRejectsTextWithoutReactionHeader() {
        XCTAssertThrowsError(try CDKRXNReader.read(text: "not an rxn"))
    }

    func testParsesAtomMapsInsideReactionParticipants() throws {
        let text = """
        $RXN
        MappedReaction
          CDKSwiftNativePort

          1  1
        $MOL
        Reactant
          CDK

          1  0  0  0  0  0            999 V2000
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  7  0  0
        M  END
        $MOL
        Product
          CDK

          1  0  0  0  0  0            999 V2000
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  7  0  0
        M  END
        """

        let reaction = try CDKRXNReader.readReaction(text: text)
        XCTAssertEqual(reaction.reactants.first?.atoms.first?.atomMapNumber, 7)
        XCTAssertEqual(reaction.products.first?.atoms.first?.atomMapNumber, 7)
    }

    func testParsesParticipantRolesAndOptionalStoichiometryFromMolHeaders() throws {
        let text = """
        $RXN
        StoichReaction
          CDKSwiftNativePort

          1  1  1
        $MOL 2
        Reactant
          CDK

          1  0  0  0  0  0            999 V2000
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        M  END
        $MOL
        Product
          CDK

          1  0  0  0  0  0            999 V2000
            1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
        M  END
        $MOL STOICH 0.5
        Agent
          CDK

          1  0  0  0  0  0            999 V2000
            2.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
        M  END
        """

        let reaction = try CDKRXNReader.readReaction(text: text)
        XCTAssertEqual(reaction.reactantParticipants.count, 1)
        XCTAssertEqual(reaction.productParticipants.count, 1)
        XCTAssertEqual(reaction.agentParticipants.count, 1)

        XCTAssertEqual(reaction.reactantParticipants.first?.role, .reactant)
        XCTAssertEqual(reaction.productParticipants.first?.role, .product)
        XCTAssertEqual(reaction.agentParticipants.first?.role, .agent)

        let reactantStoich = try XCTUnwrap(reaction.reactantParticipants.first?.stoichiometry)
        XCTAssertEqual(reactantStoich, 2.0, accuracy: 0.000_001)
        XCTAssertNil(reaction.productParticipants.first?.stoichiometry)
        let agentStoich = try XCTUnwrap(reaction.agentParticipants.first?.stoichiometry)
        XCTAssertEqual(agentStoich, 0.5, accuracy: 0.000_001)
    }
}
