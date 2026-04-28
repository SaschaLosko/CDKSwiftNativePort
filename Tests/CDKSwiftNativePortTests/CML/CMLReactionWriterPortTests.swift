import Foundation
import XCTest
#if canImport(CoreGraphics)
import CoreGraphics
#endif
@testable import CDKSwiftNativePort

final class CMLReactionWriterPortTests: XCTestCase {
    func testWritesReactionCustomizationLikeCDK() throws {
        let reaction = CDKReaction(reactants: [Molecule(name: "react", externalID: "react")],
                                   agents: [Molecule(name: "agent", externalID: "agent")],
                                   products: [Molecule(name: "product", externalID: "product")],
                                   id: "reaction1")

        let text = try CDKCMLReactionWriter.write(reaction)

        XCTAssertTrue(text.contains("<reaction id=\"reaction1\""))
        XCTAssertTrue(text.contains("<molecule id=\"react\""))
        XCTAssertTrue(text.contains("<molecule id=\"product\""))
        XCTAssertTrue(text.contains("<molecule id=\"agent\""))
    }

    func testWritesReactionProperties() throws {
        let reaction = CDKReaction(reactants: [], agents: [], products: [],
                                   id: "r1",
                                   properties: ["blabla": "blabla2"])

        let text = try CDKCMLReactionWriter.write(reaction)
        XCTAssertTrue(text.contains("<scalar dictRef=\"cdk:reactionProperty\" title=\"blabla\""))
        XCTAssertTrue(text.contains(">blabla2</scalar>"))
    }

    func testRoundTripsReaction() throws {
        let reactant = Molecule(name: "Reactant",
                                externalID: "react",
                                atoms: [Atom(id: 1, externalID: "a1", element: "C", position: .zero)],
                                bonds: [])
        let product = Molecule(name: "Product",
                               externalID: "product",
                               atoms: [Atom(id: 1, externalID: "a2", element: "R", position: CGPoint(x: 1, y: 0), aliasLabel: "R")],
                               bonds: [])
        let agent = Molecule(name: "Water",
                             externalID: "water",
                             atoms: [Atom(id: 1, externalID: "a3", element: "H", position: CGPoint(x: 2, y: 0))],
                             bonds: [])
        let reaction = CDKReaction(reactants: [reactant],
                                   agents: [agent],
                                   products: [product],
                                   id: "reaction.1")

        let text = try CDKCMLReactionWriter.write(reaction)
        let roundTripped = try CDKCMLReactionReader.readReaction(text: text)

        XCTAssertEqual(roundTripped.id, "reaction.1")
        XCTAssertEqual(roundTripped.reactantCount, 1)
        XCTAssertEqual(roundTripped.agentCount, 1)
        XCTAssertEqual(roundTripped.productCount, 1)
        XCTAssertEqual(roundTripped.reactants[0].externalID, "react")
        XCTAssertEqual(roundTripped.products[0].externalID, "product")
        XCTAssertEqual(roundTripped.agents[0].externalID, "water")
        XCTAssertEqual(roundTripped.reactants[0].atoms[0].externalID, "a1")
    }

    func testWritesReactionListAndRoundTrips() throws {
        let first = CDKReaction(reactants: [Molecule(name: "A", externalID: "A")],
                                agents: [],
                                products: [Molecule(name: "B", externalID: "B")],
                                id: "r1")
        let second = CDKReaction(reactants: [Molecule(name: "B", externalID: "B")],
                                 agents: [],
                                 products: [Molecule(name: "C", externalID: "C")],
                                 id: "r2")

        let text = try CDKCMLReactionWriter.write([first, second])
        let roundTripped = try CDKCMLReactionReader.readReactions(text: text)

        XCTAssertTrue(text.contains("<reactionList "))
        XCTAssertEqual(roundTripped.count, 2)
        XCTAssertEqual(roundTripped.map(\.id), ["r1", "r2"])
    }

    func testWritesReactionSchemeAndRoundTripsHierarchy() throws {
        let first = CDKReaction(reactants: [Molecule(name: "A", externalID: "A")],
                                agents: [],
                                products: [Molecule(name: "B", externalID: "B")],
                                id: "r1")
        let second = CDKReaction(reactants: [Molecule(name: "B", externalID: "B")],
                                 agents: [],
                                 products: [Molecule(name: "C", externalID: "C")],
                                 id: "r2")
        let third = CDKReaction(reactants: [Molecule(name: "A", externalID: "A")],
                                agents: [],
                                products: [Molecule(name: "D", externalID: "D")],
                                id: "r3")

        let scheme = CDKReactionScheme(id: "rs0",
                                       entries: [
                                        .scheme(CDKReactionScheme(id: "rs1",
                                                                  entries: [.reaction(first), .reaction(second)])),
                                        .list(CDKReactionList(id: "rsl1",
                                                              reactions: [third],
                                                              isStepList: true))
                                       ])

        let text = try CDKCMLReactionWriter.write(scheme)
        let hierarchy = try CDKCMLReactionReader.readHierarchy(text: text)

        XCTAssertTrue(text.contains("<reactionScheme"))
        XCTAssertTrue(text.contains("<reactionStepList"))
        guard case .scheme(let parsedScheme) = hierarchy else {
            return XCTFail("Expected reaction scheme root after round-trip.")
        }
        XCTAssertEqual(parsedScheme.flattenedReactions.map(\.id), ["r1", "r2", "r3"])
        guard case .list(let parsedStepList) = parsedScheme.entries[1] else {
            return XCTFail("Expected nested step list to round-trip.")
        }
        XCTAssertTrue(parsedStepList.isStepList)
    }

    func testWritesStepListContainingNestedSchemeAndRoundTripsHierarchy() throws {
        let branchA = CDKReaction(reactants: [Molecule(name: "A", externalID: "A")],
                                  agents: [],
                                  products: [Molecule(name: "B", externalID: "B")],
                                  id: "r1")
        let branchB = CDKReaction(reactants: [Molecule(name: "B", externalID: "B")],
                                  agents: [],
                                  products: [Molecule(name: "C", externalID: "C")],
                                  id: "r2")
        let finalReaction = CDKReaction(reactants: [Molecule(name: "A", externalID: "A")],
                                        agents: [],
                                        products: [Molecule(name: "D", externalID: "D")],
                                        id: "r3")

        let stepList = CDKReactionList(id: "rsl-branch",
                                       entries: [
                                        .scheme(CDKReactionScheme(id: "rs-branch",
                                                                  entries: [.reaction(branchA), .reaction(branchB)])),
                                        .reaction(finalReaction)
                                       ],
                                       isStepList: true)

        let text = try CDKCMLReactionWriter.write(stepList)
        let parsedList = try CDKCMLReactionReader.readReactionList(text: text)

        XCTAssertTrue(text.contains("<reactionStepList"))
        XCTAssertTrue(text.contains("<reactionScheme"))
        XCTAssertEqual(parsedList.entries.count, 2)
        guard case .scheme(let parsedBranch) = parsedList.entries[0] else {
            return XCTFail("Expected nested scheme to survive step-list round-trip.")
        }
        guard case .reaction(let parsedFinalReaction) = parsedList.entries[1] else {
            return XCTFail("Expected terminal reaction to survive step-list round-trip.")
        }
        XCTAssertEqual(parsedBranch.flattenedReactions.map(\.id), ["r1", "r2"])
        XCTAssertEqual(parsedFinalReaction.id, "r3")
    }

    func testWritesFormulaOnlyParticipantMolecule() throws {
        var formulaOnly = Molecule(name: "A", externalID: "A")
        formulaOnly.appendDataFieldValue("C 28 H 60 N 1", named: "Formula")

        let reaction = CDKReaction(reactants: [formulaOnly], agents: [], products: [], id: "r-formula")
        let text = try CDKCMLReactionWriter.write(reaction)

        XCTAssertTrue(text.contains("<formula concise=\"C 28 H 60 N 1\""))

        let parsed = try CDKCMLReactionReader.readReaction(text: text)
        XCTAssertEqual(parsed.reactants[0].dataFieldValues(named: "Formula"), ["C 28 H 60 N 1"])
    }
}
