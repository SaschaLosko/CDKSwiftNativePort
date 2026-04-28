import XCTest
@testable import CDKSwiftNativePort

final class CMLReactionReaderPortTests: XCTestCase {
    func testParsesReactionFragmentWithAgent() throws {
        let reaction = try CDKCMLReactionReader.readReaction(text: CMLReactionFixtures.fragmentReaction)

        XCTAssertEqual(reaction.reactantCount, 1)
        XCTAssertEqual(reaction.productCount, 1)
        XCTAssertEqual(reaction.agentCount, 1)
        XCTAssertEqual(reaction.reactants[0].externalID, "react")
        XCTAssertEqual(reaction.products[0].externalID, "product")
        XCTAssertEqual(reaction.agents[0].externalID, "water")
    }

    func testParsesReactionPropertiesAndAtomIdentifiers() throws {
        let reaction = try CDKCMLReactionReader.readReaction(text: CMLReactionFixtures.reactionWithProperties)

        XCTAssertEqual(reaction.id, "reaction.2")
        XCTAssertEqual(reaction.properties["Ka"], "3")
        XCTAssertEqual(reaction.reactants[0].externalID, "react")
        XCTAssertEqual(reaction.products[0].externalID, "product")
        XCTAssertEqual(reaction.reactants[0].atoms[0].externalID, "a1")
        XCTAssertEqual(reaction.products[0].bonds[0].externalID, "b2")
        XCTAssertEqual(reaction.reactants[0].bondCount, 1)
        XCTAssertEqual(reaction.products[0].bonds[0].order, .double)
    }

    func testParsesReactionList() throws {
        let reactions = try CDKCMLReactionReader.readReactions(text: CMLReactionFixtures.reactionList)
        let list = try CDKCMLReactionReader.readReactionList(text: CMLReactionFixtures.reactionList)

        XCTAssertEqual(reactions.count, 2)
        XCTAssertEqual(reactions.map(\.id), ["r1", "r2"])
        XCTAssertEqual(reactions[0].reactants[0].externalID, "A")
        XCTAssertEqual(reactions[0].products[0].externalID, "B")
        XCTAssertEqual(reactions[1].reactants[0].externalID, "B")
        XCTAssertEqual(reactions[1].products[0].externalID, "C")
        XCTAssertTrue(reactions[0].reactants[0].atoms.isEmpty)
        XCTAssertFalse(list.isStepList)
        XCTAssertEqual(list.reactions.map(\.id), ["r1", "r2"])
    }

    func testParsesNestedReactionScheme() throws {
        let reactions = try CDKCMLReactionReader.readReactions(text: CMLReactionFixtures.reactionScheme)
        let scheme = try CDKCMLReactionReader.readReactionScheme(text: CMLReactionFixtures.reactionScheme)

        XCTAssertEqual(reactions.count, 4)
        XCTAssertEqual(reactions.map(\.id), ["r1", "r2", "r3", "r4"])
        XCTAssertEqual(reactions[0].reactants[0].externalID, "A")
        XCTAssertEqual(reactions[0].products[0].externalID, "B")
        XCTAssertEqual(reactions[3].reactants[0].externalID, "F")
        XCTAssertEqual(reactions[3].products[0].externalID, "G")
        XCTAssertEqual(scheme.entries.count, 2)
        guard case .scheme(let leftBranch) = scheme.entries[0] else {
            return XCTFail("Expected first nested entry to remain a reaction scheme.")
        }
        guard case .scheme(let rightBranch) = scheme.entries[1] else {
            return XCTFail("Expected second nested entry to remain a reaction scheme.")
        }
        XCTAssertEqual(leftBranch.flattenedReactions.map(\.id), ["r1", "r2"])
        XCTAssertEqual(rightBranch.flattenedReactions.map(\.id), ["r3", "r4"])
    }

    func testParsesReactionStepList() throws {
        let reactions = try CDKCMLReactionReader.readReactions(text: CMLReactionFixtures.reactionStepList)
        let list = try CDKCMLReactionReader.readReactionList(text: CMLReactionFixtures.reactionStepList)

        XCTAssertEqual(reactions.count, 3)
        XCTAssertEqual(reactions.map(\.id), ["r1", "r2", "r3"])
        XCTAssertEqual(reactions.map { $0.reactants[0].externalID }, ["A", "B", "C"])
        XCTAssertEqual(reactions.map { $0.products[0].externalID }, ["B", "C", "D"])
        XCTAssertTrue(list.isStepList)
        XCTAssertEqual(list.reactions.map(\.id), ["r1", "r2", "r3"])
    }

    func testParsesReactionSchemeStepList() throws {
        let reactions = try CDKCMLReactionReader.readReactions(text: CMLReactionFixtures.reactionSchemeStepList)
        let scheme = try CDKCMLReactionReader.readReactionScheme(text: CMLReactionFixtures.reactionSchemeStepList)

        XCTAssertEqual(reactions.count, 4)
        XCTAssertEqual(reactions.map(\.id), ["r1.1", "r1.2", "r2.1", "r2.2"])
        XCTAssertEqual(reactions.map { $0.reactants[0].externalID }, ["A", "B", "A", "D"])
        XCTAssertEqual(reactions.map { $0.products[0].externalID }, ["B", "C", "D", "E"])
        XCTAssertEqual(scheme.entries.count, 2)
        guard case .scheme(let nestedScheme) = scheme.entries[0] else {
            return XCTFail("Expected first scheme entry to remain nested.")
        }
        guard case .list(let stepList) = scheme.entries[1] else {
            return XCTFail("Expected second scheme entry to remain a step list.")
        }
        XCTAssertEqual(nestedScheme.flattenedReactions.map(\.id), ["r1.1", "r1.2"])
        XCTAssertTrue(stepList.isStepList)
        XCTAssertEqual(stepList.reactions.map(\.id), ["r2.1", "r2.2"])
    }

    func testParsesStepListContainingNestedSchemeWithoutFlatteningToRoot() throws {
        let list = try CDKCMLReactionReader.readReactionList(text: CMLReactionFixtures.reactionStepListWithNestedScheme)

        XCTAssertTrue(list.isStepList)
        XCTAssertEqual(list.flattenedReactions.map(\.id), ["r1", "r2", "r3"])
        XCTAssertEqual(list.entries.count, 2)
        guard case .scheme(let branchScheme) = list.entries[0] else {
            return XCTFail("Expected first step to preserve nested reaction scheme.")
        }
        guard case .reaction(let terminalReaction) = list.entries[1] else {
            return XCTFail("Expected second step to preserve direct reaction.")
        }
        XCTAssertEqual(branchScheme.flattenedReactions.map(\.id), ["r1", "r2"])
        XCTAssertEqual(terminalReaction.id, "r3")
    }

    func testParsesSharedMoleculeReferencesDefinedBeforeReaction() throws {
        let reaction = try CDKCMLReactionReader.readReaction(text: CMLReactionFixtures.sharedMoleculeListBeforeReaction)

        XCTAssertEqual(reaction.id, "shared-1")
        XCTAssertEqual(reaction.reactants[0].externalID, "A")
        XCTAssertEqual(reaction.reactants[0].atomCount, 2)
        XCTAssertEqual(reaction.reactants[0].atoms[0].externalID, "a1")
        XCTAssertEqual(reaction.products[0].externalID, "B")
        XCTAssertEqual(reaction.products[0].atomCount, 1)
        XCTAssertEqual(reaction.products[1].externalID, "C")
        XCTAssertEqual(reaction.products[1].atomCount, 1)
    }

    func testParsesSharedMoleculeReferencesDefinedAfterReaction() throws {
        let reaction = try CDKCMLReactionReader.readReaction(text: CMLReactionFixtures.sharedMoleculeListAfterReaction)

        XCTAssertEqual(reaction.id, "shared-1")
        XCTAssertEqual(reaction.reactants[0].externalID, "A")
        XCTAssertEqual(reaction.reactants[0].atomCount, 2)
        XCTAssertEqual(reaction.products[0].externalID, "B")
        XCTAssertEqual(reaction.products[0].atomCount, 1)
        XCTAssertEqual(reaction.products[1].externalID, "C")
        XCTAssertEqual(reaction.products[1].atomCount, 1)
    }

    func testParsesFormulaOnlyMoleculeSetReferences() throws {
        let reaction = try CDKCMLReactionReader.readReaction(text: CMLReactionFixtures.formulaMoleculeSetReaction)

        XCTAssertEqual(reaction.id, "react_1")
        XCTAssertEqual(reaction.reactants[0].externalID, "A")
        XCTAssertEqual(reaction.reactants[0].dataFieldValues(named: "Formula"), ["C 28 H 60 N 1"])
        XCTAssertEqual(reaction.products[0].externalID, "B")
        XCTAssertEqual(reaction.products[0].dataFieldValues(named: "Formula"), ["C 9 H 20 N 1"])
        XCTAssertEqual(reaction.products[1].externalID, "C")
        XCTAssertTrue(reaction.products[1].dataFields.isEmpty)
    }
}
