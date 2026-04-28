import XCTest
@testable import CDKSwiftNativePort

final class ReactionSchemeManipulatorTests: XCTestCase {
    func testGetAllAtomContainersDeduplicatesAcrossNestedTopology() {
        let shared = makeMolecule(id: "B", element: "N")
        let start = makeMolecule(id: "A", element: "C")
        let end = makeMolecule(id: "C", element: "O")
        let extra = makeMolecule(id: "D", element: "S")

        let first = CDKReaction(reactants: [start], agents: [], products: [shared], id: "r1")
        let second = CDKReaction(reactants: [shared], agents: [], products: [end], id: "r2")
        let third = CDKReaction(reactants: [shared], agents: [], products: [extra], id: "r3")

        let scheme = CDKReactionScheme(id: "scheme-1",
                                       entries: [
                                        .reaction(first),
                                        .list(CDKReactionList(id: "steps",
                                                              entries: [.reaction(second)],
                                                              isStepList: true)),
                                        .scheme(CDKReactionScheme(id: "branch",
                                                                  entries: [.reaction(third)]))
                                       ])

        let molecules = CDKReactionSchemeManipulator.getAllAtomContainers(scheme)
        XCTAssertEqual(molecules.map(\.externalID), ["A", "B", "C", "D"])
    }

    func testGetAllIDsIncludesNestedSchemeAndListIDs() {
        let first = CDKReaction(reactants: [makeMolecule(id: "A", element: "C")],
                                agents: [],
                                products: [makeMolecule(id: "B", element: "N")],
                                id: "r1")
        let scheme = CDKReactionScheme(id: "scheme-1",
                                       entries: [
                                        .reaction(first),
                                        .list(CDKReactionList(id: "steps",
                                                              entries: [],
                                                              isStepList: true)),
                                        .scheme(CDKReactionScheme(id: "branch",
                                                                  entries: []))
                                       ])

        let ids = CDKReactionSchemeManipulator.getAllIDs(scheme)
        XCTAssertEqual(ids.prefix(6), ["scheme-1", "r1", "A", "B", "steps", "branch"])
    }

    func testCreateReactionSchemeBuildsUpstreamStyleBranching() throws {
        let molA = makeMolecule(id: "A", element: "C")
        let molB = makeMolecule(id: "B", element: "N")
        let molC = makeMolecule(id: "C", element: "O")
        let molD = makeMolecule(id: "D", element: "S")
        let molE = makeMolecule(id: "E", element: "P")

        let r1 = CDKReaction(reactants: [molA], agents: [], products: [molB], id: "r1")
        let r2 = CDKReaction(reactants: [molB], agents: [], products: [molC], id: "r2")
        let r3 = CDKReaction(reactants: [molB], agents: [], products: [molD], id: "r3")
        let r4 = CDKReaction(reactants: [molC], agents: [], products: [molE], id: "r4")

        let scheme = CDKReactionSchemeManipulator.createReactionScheme(
            CDKReactionSet(id: "set-1", members: [.reaction(r1), .reaction(r2), .reaction(r3), .reaction(r4)])
        )

        guard case .reaction(let rootReaction) = try XCTUnwrap(scheme.entries.first) else {
            return XCTFail("Expected the root scheme to start with a top reaction.")
        }
        XCTAssertEqual(rootReaction.id, "r1")

        guard case .scheme(let branch) = try XCTUnwrap(scheme.entries.dropFirst().first) else {
            return XCTFail("Expected nested branching scheme.")
        }

        let branchReactionIDs = branch.entries.compactMap { entry -> String? in
            guard case .reaction(let reaction) = entry else { return nil }
            return reaction.id
        }
        XCTAssertEqual(branchReactionIDs, ["r2", "r3"])

        let nestedSchemes = branch.entries.compactMap { entry -> CDKReactionScheme? in
            guard case .scheme(let scheme) = entry else { return nil }
            return scheme
        }
        XCTAssertEqual(nestedSchemes.count, 1)
        XCTAssertEqual(nestedSchemes[0].flattenedReactions.map(\.id), ["r4"])
    }

    func testExtractTopReactionsFindsOnlyRootSteps() {
        let molA = makeMolecule(id: "A", element: "C")
        let molB = makeMolecule(id: "B", element: "N")
        let molC = makeMolecule(id: "C", element: "O")

        let r1 = CDKReaction(reactants: [molA], agents: [], products: [molB], id: "r1")
        let r2 = CDKReaction(reactants: [molB], agents: [], products: [molC], id: "r2")
        let scheme = CDKReactionScheme(id: "scheme-1",
                                       entries: [.reaction(r1), .scheme(CDKReactionScheme(entries: [.reaction(r2)]))])

        let top = CDKReactionSchemeManipulator.extractTopReactions(scheme)
        XCTAssertEqual(top.flattenedReactions.map(\.id), ["r1"])
    }

    func testGetAtomContainerSetReturnsReactionPath() {
        let molA = makeMolecule(id: "A", element: "C")
        let molB = makeMolecule(id: "B", element: "N")
        let molC = makeMolecule(id: "C", element: "O")

        let r1 = CDKReaction(reactants: [molA], agents: [], products: [molB], id: "r1")
        let r2 = CDKReaction(reactants: [molB], agents: [], products: [molC], id: "r2")
        let scheme = CDKReactionSchemeManipulator.createReactionScheme(
            CDKReactionSet(members: [.reaction(r1), .reaction(r2)])
        )

        let paths = CDKReactionSchemeManipulator.getAtomContainerSet(originMol: molA,
                                                                     finalMol: molC,
                                                                     reactionScheme: scheme)
        XCTAssertEqual(paths.count, 1)
        XCTAssertEqual(paths[0].map(\.externalID), ["A", "B", "C"])
    }

    private func makeMolecule(id: String, element: String) -> Molecule {
        Molecule(name: id,
                 externalID: id,
                 atoms: [Atom(id: 1, element: element, position: .zero)],
                 bonds: [])
    }
}
