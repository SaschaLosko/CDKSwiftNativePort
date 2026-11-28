import XCTest

@testable import CDKSwiftNativePort

/// The writer derives `/` and `\` from the stored cis/trans configuration, so the
/// configuration survives edits, molfile input and generated wedges. The removal
/// case mirrors cdk/cdk #1337.
final class SmilesDirectionalStereoTests: XCTestCase {
    private let parser = CDKSmilesParser(flavor: .cdkDefault, generateCoordinates: false)

    /// Configurations as unordered pairs of reference neighbors and their relation.
    private func configurations(_ molecule: Molecule) -> Set<String> {
        let symbol = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0.element) })
        return Set(
            molecule.bonds.compactMap { bond -> String? in
                guard let configuration = bond.doubleBondStereo,
                    let references = bond.stereoReferenceAtomIDs, references.count == 4
                else { return nil }
                let ends = [symbol[references[0]]!, symbol[references[3]]!].sorted()
                return "\(ends[0])\(configuration == .cis ? "-cis-" : "-trans-")\(ends[1])"
            })
    }

    func testRoundTripsPreserveConfiguration() throws {
        for (input, expected) in [
            ("F/C=C/C", "F/C=C/C"),
            ("F\\C=C\\C", "F/C=C/C"),
            ("F/C=C\\C", "F/C=C\\C"),
            ("C(/F)=C/C", "C(/F)=C/C"),
            ("F/C(Cl)=C/C", "F/C(Cl)=C/C"),
            ("C/C=C/C=C/C", "C/C=C/C=C/C"),
            ("C/C=C\\C=C/C", "C/C=C\\C=C/C"),
            ("OC(=O)/C=C/C(=O)O", "OC(=O)/C=C/C(=O)O"),
            ("C(=C/C)\\Cl", "C(=C/C)\\Cl"),
            ("C1CCCCCC/C=C/1", "C1CCCCCC/C=C/1"),
            ("C/C=C/c1ccccc1", "C/C=C/c1ccccc1"),
        ] {
            let molecule = try parser.parseSmiles(input)
            let written = CDKSmilesGenerator().create(molecule)
            XCTAssertEqual(written, expected, input)
            XCTAssertEqual(configurations(try parser.parseSmiles(written)), configurations(molecule), input)
        }
    }

    func testConfigurationSurvivesCarrierRemoval() throws {
        for removesBondFirst in [false, true] {
            var reaction = try parser.parseReactionSmiles("F/C(Cl)=C/C>>C")
            let fluorine = reaction.reactantParticipants[0].molecule.atoms.first { $0.element == "F" }!
            if removesBondFirst {
                let bond = reaction.reactantParticipants[0].molecule.bonds.first {
                    $0.a1 == fluorine.id || $0.a2 == fluorine.id
                }!
                CDKReactionManipulator.removeElectronContainer(&reaction, bond)
            }
            CDKReactionManipulator.removeAtomAndConnectedElectronContainers(&reaction, fluorine)
            let written = CDKSmilesGenerator().create(reaction.reactantParticipants[0].molecule)
            XCTAssertEqual(written, "C(/Cl)=C/C")
            XCTAssertEqual(configurations(try parser.parseSmiles(written)), ["C-cis-Cl"])
        }
    }

    func testMolfileAlkeneIsWrittenWithConfiguration() throws {
        let molfile = """

              directional

              4  3  0  0  0  0  0  0  0  0999 V2000
                0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                1.0000    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                3.0000    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
              1  2  1  0
              2  3  2  0
              3  4  1  0
            M  END
            """
        XCTAssertEqual(CDKSmilesGenerator().create(try CDKMDLV2000Reader.read(text: molfile)), "C/C=C/C")
    }

    func testDepictionWedgesAreNotWrittenAsDirectionalBonds() throws {
        let molecule = try CDKSmilesParser().parseSmiles("C[C@H](O)CC")
        XCTAssertTrue(molecule.bonds.contains { $0.stereo != .none })
        XCTAssertEqual(CDKSmilesGenerator(flavor: [.isomeric]).create(molecule), "C[C@H](O)CC")
    }

    func testNonIsomericOutputHasNoDirectionalBonds() throws {
        XCTAssertEqual(CDKSmilesGenerator(flavor: []).create(try parser.parseSmiles("F/C=C/C")), "FC=CC")
    }
}
