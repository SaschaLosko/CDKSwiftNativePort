import XCTest
@testable import CDKSwiftNativePort

final class ReactionManipulatorTests: XCTestCase {
    func testReverseSwapsSidesAndDirection() throws {
        let reactant = CDKReactionParticipant(molecule: makeDiatomic(name: "water",
                                                                     moleculeID: "water",
                                                                     leftElement: "H",
                                                                     rightElement: "O"),
                                              role: .reactant,
                                              stoichiometry: 3.0)
        let product = CDKReactionParticipant(molecule: makeDiatomic(name: "hydroxide",
                                                                    moleculeID: "hydroxide",
                                                                    leftElement: "O",
                                                                    rightElement: "H"),
                                             role: .product,
                                             stoichiometry: 1.0)
        let reaction = CDKReaction(reactantParticipants: [reactant],
                                   agentParticipants: [],
                                   productParticipants: [product],
                                   id: "rxn-1",
                                   direction: .backward,
                                   name: "Hydrolysis")

        let reversed = CDKReactionManipulator.reverse(reaction)

        XCTAssertEqual(reversed.direction, .forward)
        XCTAssertEqual(reversed.reactantCount, 1)
        XCTAssertEqual(reversed.productCount, 1)
        XCTAssertEqual(reversed.reactants.first?.externalID, "hydroxide")
        XCTAssertEqual(reversed.products.first?.externalID, "water")
        let productStoich = try XCTUnwrap(reversed.productParticipants.first?.stoichiometry)
        XCTAssertEqual(productStoich, 3.0, accuracy: 0.000_001)
    }

    func testCollectsIDsAndRelevantMoleculeLookups() {
        let reaction = CDKReaction(reactants: [makeDiatomic(name: "reactant",
                                                            moleculeID: "react",
                                                            leftElement: "C",
                                                            rightElement: "O")],
                                   agents: [makeDiatomic(name: "agent",
                                                         moleculeID: "agent",
                                                         leftElement: "N",
                                                         rightElement: "H")],
                                   products: [makeDiatomic(name: "product",
                                                           moleculeID: "prod",
                                                           leftElement: "C",
                                                           rightElement: "N")],
                                   id: "rxn-ids")

        let ids = CDKReactionManipulator.getAllIDs(reaction)
        XCTAssertEqual(ids.prefix(4), ["rxn-ids", "react", "react-a1", "react-a2"])
        XCTAssertTrue(ids.contains("react-b1"))
        XCTAssertTrue(ids.contains("agent"))
        XCTAssertTrue(ids.contains("prod"))

        let atom = try? XCTUnwrap(reaction.reactants[0].atoms.first)
        XCTAssertEqual(CDKReactionManipulator.getRelevantAtomContainer(reaction, atom: atom!), reaction.reactants[0])
        let bond = try? XCTUnwrap(reaction.products[0].bonds.first)
        XCTAssertEqual(CDKReactionManipulator.getRelevantAtomContainer(reaction, bond: bond!), reaction.products[0])
        XCTAssertEqual(CDKReactionManipulator.getAllAtomContainers(reaction).count, 3)
    }

    func testSetsAtomPropertiesAndCollectsChemObjects() {
        var reaction = CDKReaction(reactants: [makeDiatomic(name: "reactant",
                                                            moleculeID: "react",
                                                            leftElement: "C",
                                                            rightElement: "O")],
                                   agents: [],
                                   products: [makeDiatomic(name: "product",
                                                           moleculeID: "prod",
                                                           leftElement: "C",
                                                           rightElement: "N")],
                                   id: "rxn-props")

        CDKReactionManipulator.setAtomProperties(&reaction, key: "test", value: "ok")

        for molecule in CDKReactionManipulator.getAllAtomContainers(reaction) {
            for atom in molecule.atoms {
                XCTAssertEqual(atom.properties["test"], "ok")
            }
        }

        let objects = CDKReactionManipulator.getAllChemObjects(reaction)
        XCTAssertEqual(objects.count, 3)
    }

    func testFindsMappedAtomsAndBondsAcrossReactionSides() throws {
        let reactant = makeMappedEthene(name: "reactant", moleculeID: "reactant", startMap: 1)
        let product = makeMappedEthene(name: "product", moleculeID: "product", startMap: 1)
        let reaction = CDKReaction(reactants: [reactant], agents: [], products: [product])

        let mappedAtom = try XCTUnwrap(CDKReactionManipulator.getMappedChemObject(reaction,
                                                                                   reactant.atoms[0]))
        XCTAssertEqual(mappedAtom.atomMapNumber, 1)
        let mappedBond = try XCTUnwrap(CDKReactionManipulator.getMappedChemObject(reaction,
                                                                                   reactant.bonds[0]))
        XCTAssertEqual(Set([product.atom(id: mappedBond.a1)?.atomMapNumber,
                            product.atom(id: mappedBond.a2)?.atomMapNumber].compactMap { $0 }),
                       [1, 2])

        let mappedRefs = CDKReactionManipulator.findMappedBonds(reaction)
        XCTAssertEqual(mappedRefs.count, 2)
        XCTAssertEqual(Set(mappedRefs.map(\.role)), [.reactant, .product])
    }

    func testRemovesAtomsAndBondsFromParticipantMolecules() {
        var molecule = makeMappedEthene(name: "reactant", moleculeID: "reactant", startMap: 1)
        let removableAtom = molecule.atoms[0]
        let removableBond = molecule.bonds[0]
        var reaction = CDKReaction(reactants: [molecule], agents: [], products: [])

        XCTAssertTrue(CDKReactionManipulator.removeElectronContainer(&reaction, removableBond))
        XCTAssertEqual(reaction.reactants[0].bondCount, 0)

        molecule = makeMappedEthene(name: "reactant", moleculeID: "reactant", startMap: 1)
        reaction = CDKReaction(reactants: [molecule], agents: [], products: [])
        XCTAssertTrue(CDKReactionManipulator.removeAtomAndConnectedElectronContainers(&reaction, removableAtom))
        XCTAssertEqual(reaction.reactants[0].atomCount, 1)
        XCTAssertEqual(reaction.reactants[0].bondCount, 0)
    }

    func testInlinesReactionToMoleculeAndReconstructsIt() throws {
        let reaction = CDKReaction(
            reactantParticipants: [
                CDKReactionParticipant(molecule: makeMappedEthene(name: "reactant-a",
                                                                  moleculeID: "react-a",
                                                                  startMap: 1),
                                       role: .reactant,
                                       stoichiometry: 2.0),
                CDKReactionParticipant(molecule: makeDiatomic(name: "reactant-b",
                                                              moleculeID: "react-b",
                                                              leftElement: "N",
                                                              rightElement: "N"),
                                       role: .reactant)
            ],
            agentParticipants: [
                CDKReactionParticipant(molecule: makeDiatomic(name: "agent",
                                                              moleculeID: "agent",
                                                              leftElement: "Cl",
                                                              rightElement: "Cl"),
                                       role: .agent)
            ],
            productParticipants: [
                CDKReactionParticipant(molecule: makeMappedEthene(name: "product",
                                                                  moleculeID: "product",
                                                                  startMap: 1),
                                       role: .product)
            ],
            id: "inline-rxn",
            direction: .retroSynthetic,
            name: "Inline Reaction",
            properties: ["Ka": "3"]
        )

        let inlined = CDKReactionManipulator.toMolecule(reaction)
        XCTAssertEqual(Set(inlined.atoms.compactMap(\.reactionRole)), [.reactant, .agent, .product])
        XCTAssertEqual(Set(inlined.atoms.compactMap(\.componentGroupID)), [1, 2, 3, 4])
        XCTAssertEqual(inlined.externalID, "inline-rxn")
        XCTAssertEqual(inlined.dataFieldValues(named: "__cdkReactionDirection__"), ["retroSynthetic"])

        let reconstructed = try CDKReactionManipulator.toReaction(inlined)
        XCTAssertEqual(reconstructed.id, "inline-rxn")
        XCTAssertEqual(reconstructed.name, "Inline Reaction")
        XCTAssertEqual(reconstructed.direction, .retroSynthetic)
        XCTAssertEqual(reconstructed.properties["Ka"], "3")
        XCTAssertEqual(reconstructed.reactantCount, 2)
        XCTAssertEqual(reconstructed.agentCount, 1)
        XCTAssertEqual(reconstructed.productCount, 1)
    }

    func testPerceivesAndClearsAtomConfigurations() {
        var reaction = CDKReaction(reactants: [makeMappedEthene(name: "reactant",
                                                                moleculeID: "reactant",
                                                                startMap: 1)],
                                   agents: [],
                                   products: [])

        CDKReactionManipulator.perceiveAtomTypesAndConfigureAtoms(&reaction)
        let configuredAtom = reaction.reactants[0].atoms[0]
        XCTAssertEqual(configuredAtom.atomTypeName, "C.SP2")
        XCTAssertEqual(configuredAtom.maximumBondOrder, .double)
        XCTAssertEqual(configuredAtom.formalNeighbourCount, 1)
        XCTAssertEqual(configuredAtom.hybridization, .sp2)
        XCTAssertEqual(configuredAtom.valency, 4)
        XCTAssertEqual(configuredAtom.bondOrderSum ?? 0, 2.0, accuracy: 0.000_001)

        reaction.reactants[0].atoms[0].atomTypeName = "Preset"
        reaction.reactants[0].atoms[0].formalNeighbourCount = 99
        CDKReactionManipulator.perceiveAtomTypesAndConfigureUnsetProperties(&reaction)
        XCTAssertEqual(reaction.reactants[0].atoms[0].atomTypeName, "Preset")
        XCTAssertEqual(reaction.reactants[0].atoms[0].formalNeighbourCount, 99)
        XCTAssertEqual(reaction.reactants[0].atoms[0].maximumBondOrder, .double)

        CDKReactionManipulator.clearAtomConfigurations(&reaction)
        let clearedAtom = reaction.reactants[0].atoms[0]
        XCTAssertNil(clearedAtom.atomTypeName)
        XCTAssertNil(clearedAtom.maximumBondOrder)
        XCTAssertNil(clearedAtom.bondOrderSum)
        XCTAssertNil(clearedAtom.valency)
        XCTAssertNil(clearedAtom.formalNeighbourCount)
        XCTAssertNil(clearedAtom.hybridization)
    }

    func testReactionSetManipulatorFindsRelevantReactionsAndIDs() throws {
        let shared = makeDiatomic(name: "shared", moleculeID: "shared", leftElement: "C", rightElement: "O")
        let other = makeDiatomic(name: "other", moleculeID: "other", leftElement: "O", rightElement: "H")
        let product = makeDiatomic(name: "product", moleculeID: "product", leftElement: "C", rightElement: "N")

        let first = CDKReaction(reactants: [shared], agents: [], products: [product], id: "r1")
        let second = CDKReaction(reactants: [other], agents: [], products: [shared], id: "r2")
        let set = CDKReactionSet(id: "set-1",
                                 members: [.reaction(first), .reaction(second)])

        XCTAssertEqual(CDKReactionSetManipulator.getAtomCount(set), 8)
        XCTAssertEqual(CDKReactionSetManipulator.getBondCount(set), 4)
        XCTAssertEqual(CDKReactionSetManipulator.getAllMolecules(set).count, 3)
        XCTAssertEqual(CDKReactionSetManipulator.getReactionByReactionID(set, id: "r2")?.id, "r2")
        XCTAssertEqual(CDKReactionSetManipulator.getReactionByAtomContainerID(set, id: "shared")?.id, "r2")

        let relevant = CDKReactionSetManipulator.getRelevantReactions(set, molecule: shared)
        XCTAssertEqual(relevant.flattenedReactions.map(\.id), ["r1", "r2"])
        let asReactant = CDKReactionSetManipulator.getRelevantReactionsAsReactant(set, molecule: shared)
        XCTAssertEqual(asReactant.flattenedReactions.map(\.id), ["r1"])
        let asProduct = CDKReactionSetManipulator.getRelevantReactionsAsProduct(set, molecule: shared)
        XCTAssertEqual(asProduct.flattenedReactions.map(\.id), ["r2"])

        let atom = try XCTUnwrap(shared.atoms.first)
        XCTAssertEqual(CDKReactionSetManipulator.getRelevantReaction(set, atom: atom)?.id, "r1")
    }

    func testReactionSetManipulatorSetsAtomPropertiesAndCollectsChemObjects() {
        let shared = makeDiatomic(name: "shared", moleculeID: "shared", leftElement: "C", rightElement: "O")
        let product = makeDiatomic(name: "product", moleculeID: "product", leftElement: "C", rightElement: "N")

        let first = CDKReaction(reactants: [shared], agents: [], products: [product], id: "r1")
        var set = CDKReactionSet(id: "set-1", members: [.reaction(first)])

        CDKReactionSetManipulator.setAtomProperties(&set, key: "set-test", value: "ok")

        for molecule in CDKReactionSetManipulator.getAllAtomContainers(set) {
            for atom in molecule.atoms {
                XCTAssertEqual(atom.properties["set-test"], "ok")
            }
        }

        let objects = CDKReactionSetManipulator.getAllChemObjects(set)
        XCTAssertEqual(objects.count, 4)
    }

    private func makeMappedEthene(name: String, moleculeID: String, startMap: Int) -> Molecule {
        Molecule(name: name,
                 externalID: moleculeID,
                 atoms: [
                    Atom(id: 1,
                         externalID: "\(moleculeID)-a1",
                         element: "C",
                         position: .zero,
                         atomMapNumber: startMap),
                    Atom(id: 2,
                         externalID: "\(moleculeID)-a2",
                         element: "C",
                         position: CGPoint(x: 1, y: 0),
                         atomMapNumber: startMap + 1)
                 ],
                 bonds: [
                    Bond(id: 1,
                         externalID: "\(moleculeID)-b1",
                         a1: 1,
                         a2: 2,
                         order: .double)
                 ])
    }

    private func makeDiatomic(name: String,
                              moleculeID: String,
                              leftElement: String,
                              rightElement: String) -> Molecule {
        Molecule(name: name,
                 externalID: moleculeID,
                 atoms: [
                    Atom(id: 1,
                         externalID: "\(moleculeID)-a1",
                         element: leftElement,
                         position: .zero),
                    Atom(id: 2,
                         externalID: "\(moleculeID)-a2",
                         element: rightElement,
                         position: CGPoint(x: 1, y: 0))
                 ],
                 bonds: [
                    Bond(id: 1,
                         externalID: "\(moleculeID)-b1",
                         a1: 1,
                         a2: 2,
                         order: .single)
                 ])
    }
}
