import Foundation
import XCTest
#if canImport(CoreGraphics)
import CoreGraphics
#endif
@testable import CDKSwiftNativePort

final class MDLV2000WriterTests: XCTestCase {

    func testWritesEmptyV2000Molfile() throws {
        let text = try CDKMDLV2000Writer.write(Molecule(name: "Empty"))

        XCTAssertTrue(text.contains("Empty"))
        XCTAssertTrue(text.contains("V2000"))
        XCTAssertTrue(text.contains("  0  0  0  0  0  0"))
        XCTAssertTrue(text.contains("M  END"))
    }

    func testWritesAliasAndAtomValueLines() throws {
        let molecule = Molecule(
            name: "AliasAndValue",
            atoms: [
                Atom(id: 1, element: "C", position: .zero, aliasLabel: "Phenyl"),
                Atom(id: 2, element: "C", position: CGPoint(x: 1.5, y: 0), atomValue: "Carbon note")
            ],
            bonds: [Bond(id: 1, a1: 1, a2: 2, order: .single)]
        )

        let text = try CDKMDLV2000Writer.write(molecule)

        XCTAssertTrue(text.contains("A    1"))
        XCTAssertTrue(text.contains("Phenyl"))
        XCTAssertTrue(text.contains("V    2 Carbon note"))
    }

    func testWritesHydrogenIsotopesAndRGroupLines() throws {
        let molecule = Molecule(
            name: "Isotopes",
            atoms: [
                Atom(id: 1, element: "H", position: .zero, isotopeMassNumber: 2),
                Atom(id: 2, element: "H", position: CGPoint(x: 1.5, y: 0), isotopeMassNumber: 3),
                Atom(id: 3, element: "*", position: CGPoint(x: 3.0, y: 0), rGroupLabel: 12)
            ],
            bonds: []
        )

        let text = try CDKMDLV2000Writer.write(molecule)

        XCTAssertTrue(text.contains(" D "))
        XCTAssertTrue(text.contains(" T "))
        XCTAssertTrue(text.contains("R#"))
        XCTAssertTrue(text.contains("M  RGP  1   3  12"))
    }

    func testWritesQueryAtomListsAndQueryBondsWithTopology() throws {
        let molecule = Molecule(
            name: "Queries",
            atoms: [
                Atom(id: 1, element: "L", position: .zero, queryType: .anyAtom, atomList: ["C", "N"], atomListIsNegated: true),
                Atom(id: 2, element: "C", position: CGPoint(x: 1.5, y: 0)),
                Atom(id: 3, element: "C", position: CGPoint(x: 3.0, y: 0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single, queryType: .singleOrDouble, topology: .ring),
                Bond(id: 2, a1: 2, a2: 3, order: .single, queryType: .any, topology: .chain)
            ]
        )

        let text = try CDKMDLV2000Writer.write(molecule)

        XCTAssertTrue(text.contains("M  ALS   1  2 T C   N"))
        XCTAssertTrue(text.contains("  1 T    2   6   7"))
        XCTAssertTrue(text.contains("  1  2  5  0  0  1  0"))
        XCTAssertTrue(text.contains("  2  3  8  0  0  2  0"))
    }

    func testWritesParityChiralFlagValenceAndAtomMap() throws {
        let molecule = Molecule(
            name: "Parity",
            atoms: [
                Atom(id: 1,
                     element: "C",
                     position: .zero,
                     chirality: .clockwise,
                     explicitHydrogenCount: 0,
                     valenceOverride: 0,
                     cxStereoGroup: CDKCxSmilesParser.encodeStereoGroup(kind: "abs", number: 0),
                     atomMapNumber: 7)
            ],
            bonds: []
        )

        let text = try CDKMDLV2000Writer.write(molecule)

        XCTAssertTrue(text.contains("  1  0            999 V2000"))
        XCTAssertTrue(text.contains("C   0  0  1  0  0 15  0  0  0  7"))
    }

    func testPreservesExplicitDefaultValenceForChEBIIronSulfurCluster() throws {
        let source = """
        CHEBI:33740
          Marvin  03051214042D

          7  9  0  0  0  0            999 V2000
            0.4125    0.4125    0.0000 S   0  5  0  0  0  2  0  0  0  0  0  0
           -0.4125   -0.4125    0.0000 S   0  0  0  0  0  3  0  0  0  0  0  0
            0.4125   -0.4125    0.0000 Fe  0  0  0  0  0  0  0  0  0  0  0  0
           -0.4125    0.4125    0.0000 Fe  0  0  0  0  0  0  0  0  0  0  0  0
           -0.8690   -0.8690    0.0000 Fe  0  0  0  0  0  0  0  0  0  0  0  0
           -0.8690   -0.0440    0.0000 S   0  5  0  0  0  2  0  0  0  0  0  0
           -0.0440   -0.8690    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
          3  2  1  0  0  0  0
          3  1  1  0  0  0  0
          2  4  1  0  0  0  0
          1  4  1  0  0  0  0
          2  5  1  0  0  0  0
          3  7  1  0  0  0  0
          4  6  1  0  0  0  0
          6  5  1  0  0  0  0
          5  7  1  0  0  0  0
        M  CHG  2   1  -1   6  -1
        M  END
        """
        let molecule = try CDKMDLV2000Reader.read(text: source)

        XCTAssertEqual(molecule.atoms[0].valenceOverride, 2)
        XCTAssertEqual(molecule.atoms[5].valenceOverride, 2)

        let serialized = try CDKMDLV2000Writer.write(molecule)
        let reparsed = try CDKMDLV2000Reader.read(text: serialized)
        XCTAssertEqual(reparsed.atoms[0].valenceOverride, 2)
        XCTAssertEqual(reparsed.atoms[5].valenceOverride, 2)

        let identifiers = CDKMoleculeIdentifierService.compute(
            for: molecule,
            recalculateInChI: true)
        XCTAssertEqual(identifiers.inchi, "InChI=1S/3Fe.4S/q;;;;;2*-1")
        XCTAssertEqual(identifiers.inchiKey, "ZCPXJLJVCGFYTC-UHFFFAOYSA-N")
        XCTAssertEqual(identifiers.fixedHInchi, "InChI=1/3Fe.4S/q;;;;;2*-1")
        XCTAssertEqual(identifiers.fixedHInchiKey, "ZCPXJLJVCGFYTC-UHFFFAOYNA-N")
    }

    func testWritesSgroupsAndDataSgroups() throws {
        var molecule = Molecule(
            name: "Sgroups",
            atoms: [
                Atom(id: 1, element: "C", position: .zero),
                Atom(id: 2, element: "O", position: CGPoint(x: 1.5, y: 0)),
                Atom(id: 3, element: "C", position: CGPoint(x: 3.0, y: 0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single),
                Bond(id: 2, a1: 2, a2: 3, order: .single)
            ]
        )
        molecule.sgroups = [
            MoleculeSgroup(kind: .structureRepeatUnit,
                           keyword: "SRU",
                           atomIDs: [1, 2],
                           crossingBondIDs: [1],
                           subscriptText: "n",
                           roundBrackets: true,
                           connectivity: "HT",
                           brackets: [
                            MoleculeSgroupBracket(firstPoint: CGPoint(x: -1.0, y: -1.0),
                                                  secondPoint: CGPoint(x: -1.0, y: 1.0))
                           ]),
            MoleculeSgroup(kind: .data,
                           keyword: "DAT",
                           atomIDs: [2],
                           crossingBondIDs: [2],
                           dataFieldName: "FIELD",
                           dataValue: "Atom data",
                           dataOperator: ">=",
                           dataUnit: "u",
                           dataTag: "TG")
        ]

        let text = try CDKMDLV2000Writer.write(molecule)

        XCTAssertTrue(text.contains("M  STY  2   1 SRU   2 DAT"))
        XCTAssertTrue(text.contains("M  SAL   1  2   1   2"))
        XCTAssertTrue(text.contains("M  SBL   1  1   1"))
        XCTAssertTrue(text.contains("M  SMT   1 n"))
        XCTAssertTrue(text.contains("M  SCN  1   1 HT"))
        XCTAssertTrue(text.contains("M  SBT  1   1   1"))
        XCTAssertTrue(text.contains("M  SDT   2 FIELD"))
        XCTAssertTrue(text.contains("M  SED   2 Atom data"))
    }

    func testWritesDimensionFieldFor2DAnd3DScenes() throws {
        let planar = Molecule(
            name: "Planar",
            atoms: [Atom(id: 1, element: "C", position: CGPoint(x: 0.5, y: 0.5))],
            bonds: []
        )
        let planarText = try CDKMDLV2000Writer.write(planar)
        XCTAssertTrue(planarText.contains("2D"))

        let threeD = Molecule(
            name: "ThreeD",
            atoms: [Atom(id: 1, element: "C", position: CGPoint(x: 0.5, y: 0.5), zPosition: 0.1)],
            bonds: []
        )
        let threeDText = try CDKMDLV2000Writer.write(threeD)
        XCTAssertTrue(threeDText.contains("3D"))
    }

    func testWritesSignedZeroCoordinatesWithoutNegativeZeroArtifacts() throws {
        let molecule = Molecule(
            name: "SignedZero",
            atoms: [Atom(id: 1, element: "C", position: CGPoint(x: -0.0, y: 0.0), zPosition: -0.0)],
            bonds: []
        )

        let text = try CDKMDLV2000Writer.write(molecule)

        XCTAssertTrue(text.contains("    0.0000    0.0000    0.0000"))
        XCTAssertFalse(text.contains("-0.0000"))
    }

    func testRGroupPropertyWrapsAcrossMultipleLines() throws {
        let atoms = (1...19).map { index in
            Atom(id: index, element: "*", position: CGPoint(x: Double(index), y: 0), rGroupLabel: index)
        }
        let molecule = Molecule(name: "RGroups", atoms: atoms, bonds: [])

        let text = try CDKMDLV2000Writer.write(molecule)

        XCTAssertTrue(text.contains("M  RGP  8"))
        XCTAssertTrue(text.contains("M  RGP  3  17  17  18  18  19  19"))
    }

    func testRoundTripsAliasAtomValueSgroupsAndTopology() throws {
        var molecule = Molecule(
            name: "RoundTrip",
            atoms: [
                Atom(id: 1, element: "C", position: .zero, aliasLabel: "Phenyl"),
                Atom(id: 2, element: "L", position: CGPoint(x: 1.5, y: 0), queryType: .anyAtom, atomList: ["C", "N"], atomListIsNegated: false, atomValue: "Center"),
                Atom(id: 3, element: "*", position: CGPoint(x: 3.0, y: 0), rGroupLabel: 2)
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single, queryType: .singleOrAromatic, topology: .ring),
                Bond(id: 2, a1: 2, a2: 3, order: .single)
            ]
        )
        molecule.sgroups = [
            MoleculeSgroup(kind: .structureRepeatUnit,
                           keyword: "SRU",
                           atomIDs: [1, 2],
                           crossingBondIDs: [1],
                           subscriptText: "n")
        ]

        let text = try CDKMDLV2000Writer.write(molecule)
        let parsed = try CDKMDLV2000Reader.read(text: text)

        XCTAssertEqual(parsed.atoms[0].aliasLabel, "Phenyl")
        XCTAssertEqual(parsed.atoms[1].atomValue, "Center")
        XCTAssertEqual(parsed.atoms[1].atomList ?? [], ["C", "N"])
        XCTAssertEqual(parsed.bonds[0].queryType, .singleOrAromatic)
        XCTAssertEqual(parsed.bonds[0].topology, .ring)
        XCTAssertEqual(parsed.atoms[2].rGroupLabel, 2)
        XCTAssertEqual(parsed.sgroups.first?.keyword, "SRU")
    }
}
