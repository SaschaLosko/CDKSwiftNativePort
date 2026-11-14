import Foundation
import XCTest
@testable import CDKSwiftNativePort

final class CMLStereoRoundTripTests: XCTestCase {
    private func identity(_ molecule: Molecule) throws -> String {
        let result = try CDKInChIOfficialLibraryGenerator.generate(for: molecule)
        XCTAssertTrue(result.inchi.hasPrefix("InChI=1S/"))
        return result.inchi
    }

    func testTetrahedralAndAlkeneStereoSurviveMoleculeCML() throws {
        let inputs = [
            "N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O", "F[C@](Cl)(Br)I",
            "F/C=C/F", "F/C=C\\F", "C[C@H](O)[C@@H](N)C", "CC(O)C", "FC=CF",
        ]
        var identities: [String] = []
        for smiles in inputs {
            var source = try CDKSmilesParser().parseSmiles(smiles)
            // Stable stereo must not depend on stored array or drawing order.
            source.atoms.reverse(); source.bonds.reverse()
            let text = try CDKCMLWriter.write([source])
            let restored = try XCTUnwrap(CDKCMLReader.read(text: text).first)
            XCTAssertEqual(try identity(restored), try identity(source), smiles)
            if source.atoms.contains(where: { $0.chirality != .none }) { XCTAssertTrue(text.contains("<atomParity")) }
            if source.bonds.contains(where: { $0.doubleBondStereo != nil }) {
                XCTAssertTrue(text.contains("<bondStereo atomRefs4="))
            }
            identities.append(try identity(restored))
        }
        XCTAssertNotEqual(identities[0], identities[1])
        XCTAssertNotEqual(identities[3], identities[4])
        XCTAssertNotEqual(identities[3], identities[7])
    }

    func testReactionParticipantsPreserveStereoAndExternalReferences() throws {
        var first = try CDKSmilesParser().parseSmiles("N[C@@H](C)C(=O)O")
        for index in first.atoms.indices { first.atoms[index].externalID = "atom_\(index)" }
        let second = try CDKSmilesParser().parseSmiles("F/C=C\\F")
        let unspecified = try CDKSmilesParser().parseSmiles("FC=CF")
        let reaction = CDKReaction(reactants: [first], agents: [unspecified], products: [second])
        let text = try CDKCMLReactionWriter.write(reaction)
        let restored = try CDKCMLReactionReader.readReaction(text: text)
        XCTAssertEqual(try identity(restored.reactants[0]), try identity(first))
        XCTAssertEqual(try identity(restored.products[0]), try identity(second))
        XCTAssertEqual(try identity(restored.agents[0]), try identity(unspecified))
    }

    func testReversedWedgeEndpointIsRetained() throws {
        for stereo in [BondStereo.upReversed, .downReversed] {
            let source = Molecule(
                name: "Wedge",
                atoms: [Atom(id: 2, element: "C", position: .zero), Atom(id: 7, element: "F", position: .zero)],
                bonds: [Bond(id: 1, a1: 2, a2: 7, order: .single, stereo: stereo)])
            let text = try CDKCMLWriter.write([source])
            let restored = try XCTUnwrap(CDKCMLReader.read(text: text).first)
            XCTAssertEqual(restored.bonds.first?.stereo, stereo)
        }
    }

    func testMalformedStereoCannotSilentlyBecomeUnspecified() throws {
        var source = try CDKSmilesParser().parseSmiles("F/C=C/F")
        source.bonds[source.bonds.firstIndex { $0.order == .double }!].stereoReferenceAtomIDs = [1, 2, 3, 999]
        XCTAssertThrowsError(try CDKCMLWriter.write([source]))
        let valid = try CDKCMLWriter.write([CDKSmilesParser().parseSmiles("F[C@](Cl)(Br)I")])
        let broken = valid.replacingOccurrences(of: "<atomParity atomRefs4=", with: "<atomParity brokenRefs=")
        XCTAssertThrowsError(try CDKCMLReader.read(text: broken))
    }
}
