import XCTest
@testable import CDKSwiftNativePort

final class CMLReactionDisplayNameTests: XCTestCase {
    func testNamesEqualToGeneratedIDsSurviveRepeatedRoundTrips() throws {
        var reaction = try CDKFileImporter.readReaction(text: "CCO>O>CC=O", fileExtension: "rsmi")
        reaction.name = "Oxidation"
        reaction.reactantParticipants[0].molecule.name = "Ethanol"
        reaction.agentParticipants[0].molecule.name = "Water"
        reaction.productParticipants[0].molecule.name = "Acetaldehyde"
        for _ in 0..<2 {
            reaction = try CDKCMLReactionReader.readReaction(text: CDKCMLReactionWriter.write(reaction))
            XCTAssertEqual(reaction.name, "Oxidation")
            XCTAssertEqual(reaction.reactants.map(\.name), ["Ethanol"])
            XCTAssertEqual(reaction.agents.map(\.name), ["Water"])
            XCTAssertEqual(reaction.products.map(\.name), ["Acetaldehyde"])
            XCTAssertEqual(reaction.participants.map { $0.molecule.atomCount }, [3, 1, 3])
        }
    }

    func testDuplicateIDsDoNotMergeDistinctParticipantsAcrossReactions() throws {
        var first = try CDKFileImporter.readReaction(text: "CCO>O>CC=O", fileExtension: "rsmi")
        first.name = "First"
        first.reactantParticipants[0].molecule.externalID = "shared"
        first.agentParticipants[0].molecule.externalID = "shared"
        first.productParticipants[0].molecule.externalID = "shared_2"
        first.reactantParticipants[0].molecule.name = "Ethanol"
        first.agentParticipants[0].molecule.name = "Water"
        first.productParticipants[0].molecule.name = "Acetaldehyde"
        var second = first
        second.name = "Second"
        second.productParticipants[0].molecule.name = "Other product"
        let text = try CDKCMLReactionWriter.write([first, second])
        let reopened = try CDKCMLReactionReader.readReactions(text: text)
        XCTAssertEqual(reopened.map(\.name), ["First", "Second"])
        XCTAssertEqual(reopened[0].participants.map { $0.molecule.name }, ["Ethanol", "Water", "Acetaldehyde"])
        XCTAssertEqual(reopened[1].products[0].name, "Other product")
        let molecules = reopened.flatMap(\.participants).map(\.molecule)
        XCTAssertEqual(Set(molecules.compactMap(\.externalID)).count, 4)
        XCTAssertEqual(reopened[0].reactants[0].externalID, reopened[1].reactants[0].externalID)
        XCTAssertEqual(reopened[0].agents[0].externalID, reopened[1].agents[0].externalID)
        XCTAssertNotEqual(reopened[0].products[0].externalID, reopened[1].products[0].externalID)
        XCTAssertEqual(molecules.map(\.atomCount), [3, 1, 3, 3, 1, 3])
    }
}
