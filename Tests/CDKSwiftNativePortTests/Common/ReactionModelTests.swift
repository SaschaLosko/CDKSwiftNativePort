import XCTest
@testable import CDKSwiftNativePort

final class ReactionModelTests: XCTestCase {
    func testMoleculeArrayInitializerCreatesRoleParticipants() {
        let reactant = makeMolecule(name: "R", element: "C")
        let agent = makeMolecule(name: "A", element: "N")
        let product = makeMolecule(name: "P", element: "O")

        let reaction = CDKReaction(reactants: [reactant], agents: [agent], products: [product])

        XCTAssertEqual(reaction.reactantParticipants.count, 1)
        XCTAssertEqual(reaction.agentParticipants.count, 1)
        XCTAssertEqual(reaction.productParticipants.count, 1)

        XCTAssertEqual(reaction.reactantParticipants.first?.role, .reactant)
        XCTAssertEqual(reaction.agentParticipants.first?.role, .agent)
        XCTAssertEqual(reaction.productParticipants.first?.role, .product)
        XCTAssertTrue(reaction.participants.allSatisfy { $0.stoichiometry == nil })
    }

    func testParticipantInitializerNormalizesRolesBySectionAndRetainsStoichiometry() throws {
        let participant = CDKReactionParticipant(
            molecule: makeMolecule(name: "Input", element: "C"),
            role: .product,
            stoichiometry: 3.0
        )

        let reaction = CDKReaction(reactantParticipants: [participant],
                                   agentParticipants: [],
                                   productParticipants: [])

        XCTAssertEqual(reaction.reactantParticipants.count, 1)
        XCTAssertEqual(reaction.reactantParticipants.first?.role, .reactant)
        let reactantStoich = try XCTUnwrap(reaction.reactantParticipants.first?.stoichiometry)
        XCTAssertEqual(reactantStoich, 3.0, accuracy: 0.000_001)
    }

    func testFlatParticipantInitializerGroupsByRole() throws {
        let participants = [
            CDKReactionParticipant(molecule: makeMolecule(name: "R", element: "C"), role: .reactant, stoichiometry: 1.5),
            CDKReactionParticipant(molecule: makeMolecule(name: "A", element: "N"), role: .agent, stoichiometry: nil),
            CDKReactionParticipant(molecule: makeMolecule(name: "P", element: "O"), role: .product, stoichiometry: 2.0)
        ]

        let reaction = CDKReaction(participants: participants)

        XCTAssertEqual(reaction.reactantCount, 1)
        XCTAssertEqual(reaction.agentCount, 1)
        XCTAssertEqual(reaction.productCount, 1)
        let reactantStoich = try XCTUnwrap(reaction.reactantParticipants.first?.stoichiometry)
        XCTAssertEqual(reactantStoich, 1.5, accuracy: 0.000_001)
        XCTAssertNil(reaction.agentParticipants.first?.stoichiometry)
        let productStoich = try XCTUnwrap(reaction.productParticipants.first?.stoichiometry)
        XCTAssertEqual(productStoich, 2.0, accuracy: 0.000_001)
    }

    private func makeMolecule(name: String, element: String) -> Molecule {
        Molecule(name: name,
                 atoms: [Atom(id: 1, element: element, position: .zero)],
                 bonds: [])
    }
}
