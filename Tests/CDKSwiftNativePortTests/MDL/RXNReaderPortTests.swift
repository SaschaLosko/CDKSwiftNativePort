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
        XCTAssertEqual(reaction.name, "DemoReaction")
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

    func testParsesV3000ReactionBlockWithAgents() throws {
        let text = """
        $RXN V3000
        DemoV3000
          CDKSwiftNativePort

        M  V30 COUNTS 2 1 2
        M  V30 BEGIN REACTANT
        M  V30 BEGIN CTAB
        M  V30 COUNTS 2 1 0 0 0
        M  V30 BEGIN ATOM
        M  V30 1 C 0.0 0.0 0.0 1
        M  V30 2 O 1.0 0.0 0.0 2
        M  V30 END ATOM
        M  V30 BEGIN BOND
        M  V30 1 1 1 2
        M  V30 END BOND
        M  V30 END CTAB
        M  V30 BEGIN CTAB
        M  V30 COUNTS 1 0 0 0 0
        M  V30 BEGIN ATOM
        M  V30 1 Cl 0.0 0.0 0.0 0
        M  V30 END ATOM
        M  V30 END CTAB
        M  V30 END REACTANT
        M  V30 BEGIN PRODUCT
        M  V30 BEGIN CTAB
        M  V30 COUNTS 2 1 0 0 0
        M  V30 BEGIN ATOM
        M  V30 1 C 0.0 0.0 0.0 1
        M  V30 2 N 1.0 0.0 0.0 2
        M  V30 END ATOM
        M  V30 BEGIN BOND
        M  V30 1 2 1 2
        M  V30 END BOND
        M  V30 END CTAB
        M  V30 END PRODUCT
        M  V30 BEGIN AGENT
        M  V30 BEGIN CTAB
        M  V30 COUNTS 1 0 0 0 0
        M  V30 BEGIN ATOM
        M  V30 1 Al 0.0 0.0 0.0 0 CHG=3
        M  V30 END ATOM
        M  V30 END CTAB
        M  V30 BEGIN CTAB
        M  V30 COUNTS 1 0 0 0 0
        M  V30 BEGIN ATOM
        M  V30 1 Cl 0.0 0.0 0.0 0 CHG=-1
        M  V30 END ATOM
        M  V30 END CTAB
        M  V30 END AGENT
        M  END
        """

        let reaction = try CDKRXNReader.readReaction(text: text)
        XCTAssertEqual(reaction.name, "DemoV3000")
        XCTAssertEqual(reaction.reactantCount, 2)
        XCTAssertEqual(reaction.productCount, 1)
        XCTAssertEqual(reaction.agentCount, 2)
        XCTAssertEqual(reaction.reactants[0].atoms.compactMap(\.atomMapNumber), [1, 2])
        XCTAssertEqual(reaction.products[0].atoms.compactMap(\.atomMapNumber), [1, 2])
        XCTAssertEqual(reaction.agents[0].atoms.first?.charge, 3)
        XCTAssertEqual(reaction.agents[1].atoms.first?.charge, -1)
    }

    func testRelaxedModeInfersAgentsFromExtraEntitiesWhenCountsOmitThem() throws {
        let text = try loadUpstreamReactionFixture(named: "ethylesterification_countsLineHasNoAgents_countsLineMismatchMolfiles.rxn")

        let reaction = try CDKRXNReader.readReaction(text: text, options: .init(mode: .relaxed))

        XCTAssertEqual(reaction.reactantCount, 2)
        XCTAssertEqual(reaction.productCount, 2)
        XCTAssertEqual(reaction.agentCount, 1)
    }

    func testStrictModeRejectsAgentCountExtension() throws {
        let text = try loadUpstreamReactionFixture(named: "ethylesterification_countsLineHasAgents_countsLineMatchesMolfiles.rxn")

        XCTAssertThrowsError(try CDKRXNReader.readReaction(text: text, options: .init(mode: .strict))) { error in
            XCTAssertTrue(error.localizedDescription.contains(
                "RXN files uses agent count extension. This is not supported in mode STRICT"
            ))
        }
    }

    func testStrictModeRejectsExtraEntitiesWhenCountsOmitAgents() throws {
        let text = try loadUpstreamReactionFixture(named: "ethylesterification_countsLineHasNoAgents_countsLineMismatchMolfiles.rxn")

        XCTAssertThrowsError(try CDKRXNReader.readReaction(text: text, options: .init(mode: .strict))) { error in
            XCTAssertTrue(error.localizedDescription.contains(
                "Agents are not supported in mode STRICT. Found 5 molecular entities, but there are only 4 molecular entities declared on the counts line."
            ))
        }
    }

    func testParsesUpstreamV3000Resource() throws {
        let text = try loadUpstreamReactionFixture(named: "reaction_v3.rxn")

        let reaction = try CDKRXNReader.readReaction(text: text)

        XCTAssertEqual(reaction.reactantCount, 1)
        XCTAssertEqual(reaction.productCount, 1)
        XCTAssertEqual(reaction.reactants.first?.atomCount, 32)
        XCTAssertEqual(reaction.reactants.first?.bondCount, 29)
        XCTAssertEqual(reaction.products.first?.atomCount, 32)
        XCTAssertEqual(reaction.products.first?.bondCount, 29)
    }

    private func loadUpstreamReactionFixture(named name: String) throws -> String {
        let url = URL(fileURLWithPath: #filePath)
            .deletingLastPathComponent()
            .deletingLastPathComponent()
            .appendingPathComponent("Reaction/UpstreamReference/MDL")
            .appendingPathComponent(name)
        return try String(contentsOf: url, encoding: .utf8)
    }
}
