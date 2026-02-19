import XCTest
@testable import CDKSwiftNativePort

final class MDLV2000ReaderPortTests: XCTestCase {

    func testParsesBasicV2000Molfile() throws {
        let molecule = try CDKMDLV2000Reader.read(text: ethanolMol)

        XCTAssertEqual(molecule.name, "Ethanol")
        XCTAssertEqual(molecule.atomCount, 3)
        XCTAssertEqual(molecule.bondCount, 2)
        XCTAssertEqual(molecule.atoms.map(\.element), ["C", "C", "O"])
    }

    func testParsesBondStereoCode() throws {
        let molecule = try CDKMDLV2000Reader.read(text: stereoMol)
        XCTAssertEqual(molecule.bondCount, 1)
        XCTAssertEqual(molecule.bonds[0].stereo, .up)
    }

    func testParsesAtomBlockChargeCode() throws {
        let molecule = try CDKMDLV2000Reader.read(text: chargedAtomMol)
        XCTAssertEqual(molecule.atomCount, 1)
        XCTAssertEqual(molecule.atoms[0].charge, -1)
    }

    func testParsesMChgPropertyLine() throws {
        let molecule = try CDKMDLV2000Reader.read(text: mChgMol)
        XCTAssertEqual(molecule.atomCount, 1)
        XCTAssertEqual(molecule.atoms[0].charge, 1)
    }

    func testParsesMIsotopePropertyLine() throws {
        let molecule = try CDKMDLV2000Reader.read(text: isotopeMol)
        XCTAssertEqual(molecule.atomCount, 1)
        XCTAssertEqual(molecule.atoms[0].isotopeMassNumber, 13)
    }

    func testParsesAtomListFromMALSRecord() throws {
        let molecule = try CDKMDLV2000Reader.read(text: atomListMol)
        let atom = try XCTUnwrap(molecule.atoms.first)

        XCTAssertEqual(atom.element, "L")
        XCTAssertEqual(atom.atomList ?? [], ["C", "N"])
        XCTAssertFalse(atom.atomListIsNegated)
        XCTAssertEqual(atom.queryType, .anyAtom)
    }

    func testParsesQueryBondTypes() throws {
        let molecule = try CDKMDLV2000Reader.read(text: queryBondMol)
        XCTAssertEqual(molecule.bondCount, 4)

        XCTAssertEqual(molecule.bonds[0].queryType, .singleOrDouble)
        XCTAssertEqual(molecule.bonds[1].queryType, .singleOrAromatic)
        XCTAssertEqual(molecule.bonds[2].queryType, .doubleOrAromatic)
        XCTAssertEqual(molecule.bonds[3].queryType, .any)
    }

    func testParsesAdditionalMAtomRecords() throws {
        let molecule = try CDKMDLV2000Reader.read(text: additionalMRecordsMol)
        let atom = try XCTUnwrap(molecule.atoms.first)

        XCTAssertEqual(atom.radical, 2)
        XCTAssertEqual(atom.rGroupLabel, 7)
        XCTAssertEqual(atom.attachmentPoint, 3)
        XCTAssertEqual(atom.substitutionCount, 2)
        XCTAssertEqual(atom.unsaturated, 1)
        XCTAssertEqual(atom.ringBondCount, 0)
    }

    func testParsesAtomQuerySymbolTypes() throws {
        let molecule = try CDKMDLV2000Reader.read(text: queryAtomSymbolMol)

        XCTAssertEqual(molecule.atomCount, 3)
        XCTAssertEqual(molecule.atoms[0].queryType, .anyNonHydrogen)
        XCTAssertEqual(molecule.atoms[1].queryType, .anyHetero)
        XCTAssertEqual(molecule.atoms[2].queryType, .anyAtom)
    }

    func testRejectsV3000InputWithoutCTAB() {
        XCTAssertThrowsError(try CDKMDLV2000Reader.read(text: v3000MolMissingCTAB))
    }

    func testRejectsTruncatedInput() {
        XCTAssertThrowsError(try CDKMDLV2000Reader.read(text: truncatedMol))
    }

    func testGeneratesLayoutWhenCoordinatesAreMissing() throws {
        let molecule = try CDKMDLV2000Reader.read(text: zeroCoordinateMol)
        let box = try XCTUnwrap(molecule.boundingBox())
        XCTAssertGreaterThan(box.width, 0.1)
    }

    private let ethanolMol = """
Ethanol
CDKSwiftNativePort

  3  2  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
"""

    private let stereoMol = """
Stereo
CDKSwiftNativePort

  2  1  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1  0  0  0
M  END
"""

    private let chargedAtomMol = """
ChargedAtom
CDKSwiftNativePort

  1  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
M  END
"""

    private let mChgMol = """
MChg
CDKSwiftNativePort

  1  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
M  CHG  1   1   1
M  END
"""

    private let isotopeMol = """
MIso
CDKSwiftNativePort

  1  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  ISO  1   1  13
M  END
"""

    private let atomListMol = """
AtomList
CDKSwiftNativePort

  1  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
M  ALS   1  2 F C N
M  END
"""

    private let queryBondMol = """
QueryBonds
CDKSwiftNativePort

  5  4  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  5  0  0  0  0
  2  3  6  0  0  0  0
  3  4  7  0  0  0  0
  4  5  8  0  0  0  0
M  END
"""

    private let additionalMRecordsMol = """
AdditionalMRecords
CDKSwiftNativePort

  1  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  RAD  1   1   2
M  RGP  1   1   7
M  APO  1   1   3
M  SUB  1   1   2
M  UNS  1   1   1
M  RBC  1   1   0
M  END
"""

    private let queryAtomSymbolMol = """
QueryAtomTypes
CDKSwiftNativePort

  3  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
    1.2000    0.0000    0.0000 Q   0  0  0  0  0  0  0  0  0  0  0  0
    2.4000    0.0000    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0
M  END
"""

    private let v3000MolMissingCTAB = """
V3000Example
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  END
"""

    private let truncatedMol = """
Truncated
CDKSwiftNativePort

  2  1  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
"""

    private let zeroCoordinateMol = """
ZeroCoordinates
CDKSwiftNativePort

  2  1  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
"""
}
