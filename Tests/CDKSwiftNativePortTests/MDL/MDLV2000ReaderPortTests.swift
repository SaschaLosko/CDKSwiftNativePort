import Foundation
import XCTest
#if canImport(CoreGraphics)
import CoreGraphics
#endif
@testable import CDKSwiftNativePort

final class MDLV2000ReaderPortTests: XCTestCase {

    func testParsesBasicV2000Molfile() throws {
        let molecule = try CDKMDLV2000Reader.read(text: ethanolMol)

        XCTAssertEqual(molecule.name, "Ethanol")
        XCTAssertEqual(molecule.atomCount, 3)
        XCTAssertEqual(molecule.bondCount, 2)
        XCTAssertEqual(molecule.atoms.map(\.element), ["C", "C", "O"])
    }

    func testParsesEmptyV2000Molfile() throws {
        let molecule = try CDKMDLV2000Reader.read(text: emptyMol)

        XCTAssertEqual(molecule.atomCount, 0)
        XCTAssertEqual(molecule.bondCount, 0)
        XCTAssertEqual(molecule.name, "Empty")
    }

    func testParsesAliasAndAtomValueLines() throws {
        let molecule = try CDKMDLV2000Reader.read(text: aliasAndValueMol)

        XCTAssertEqual(molecule.atomCount, 2)
        XCTAssertEqual(molecule.atoms[0].aliasLabel, "Phenyl")
        XCTAssertEqual(molecule.atoms[1].atomValue, "Carbon note")
    }

    func testParsesHydrogenIsotopesAndMassDiff() throws {
        let molecule = try CDKMDLV2000Reader.read(text: isotopeMol)

        XCTAssertEqual(molecule.atoms[0].element, "H")
        XCTAssertEqual(molecule.atoms[0].isotopeMassNumber, 2)
        XCTAssertEqual(molecule.atoms[1].element, "H")
        XCTAssertEqual(molecule.atoms[1].isotopeMassNumber, 3)
        XCTAssertEqual(molecule.atoms[2].isotopeMassNumber, 13)
    }

    func testParsesValenceParityAndChiralFlag() throws {
        let molecule = try CDKMDLV2000Reader.read(text: parityMol)
        let atom = try XCTUnwrap(molecule.atoms.first)

        XCTAssertEqual(atom.chirality, .clockwise)
        XCTAssertEqual(atom.valenceOverride, 0)
        XCTAssertEqual(atom.cxStereoGroup,
                       CDKCxSmilesParser.encodeStereoGroup(kind: "abs", number: 0))
    }

    func testParsesQueryBondsAndBondTopology() throws {
        let molecule = try CDKMDLV2000Reader.read(text: queryBondMol)

        XCTAssertEqual(molecule.bonds[0].queryType, .singleOrDouble)
        XCTAssertEqual(molecule.bonds[0].topology, .ring)
        XCTAssertEqual(molecule.bonds[1].queryType, .any)
        XCTAssertEqual(molecule.bonds[1].topology, .chain)
    }

    func testParsesLegacyAtomListAndMALSOverrides() throws {
        let molecule = try CDKMDLV2000Reader.read(text: atomListMol)
        let atom = try XCTUnwrap(molecule.atoms.first)

        XCTAssertEqual(atom.element, "L")
        XCTAssertEqual(atom.queryType, .anyAtom)
        XCTAssertEqual(atom.atomList ?? [], ["C", "N"])
        XCTAssertTrue(atom.atomListIsNegated)
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

    func testParsesSgroupsAndDataSgroups() throws {
        let molecule = try CDKMDLV2000Reader.read(text: sgroupMol)

        XCTAssertEqual(molecule.sgroups.count, 2)
        XCTAssertEqual(molecule.sgroups[0].keyword, "SRU")
        XCTAssertEqual(molecule.sgroups[0].atomIDs, [1, 2])
        XCTAssertEqual(molecule.sgroups[0].crossingBondIDs, [1])
        XCTAssertEqual(molecule.sgroups[0].subscriptText, "n")
        XCTAssertTrue(molecule.sgroups[0].roundBrackets)
        XCTAssertEqual(molecule.sgroups[0].connectivity, "HT")
        XCTAssertEqual(molecule.sgroups[1].keyword, "DAT")
        XCTAssertEqual(molecule.sgroups[1].dataFieldName, "FIELD")
        XCTAssertEqual(molecule.sgroups[1].dataValue, "Atom data")
        XCTAssertEqual(molecule.dataFieldValues(named: "FIELD"), ["Atom data"])
    }

    func testGeneratesLayoutWhenCoordinatesAreMissing() throws {
        let molecule = try CDKMDLV2000Reader.read(text: zeroCoordinateMol)
        let box = try XCTUnwrap(molecule.boundingBox())
        XCTAssertGreaterThan(box.width, 0.1)
    }

    private let ethanolMol = """
Ethanol
CDKSwiftNativePort

  3  2  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
"""

    private let emptyMol = """
Empty
CDKSwiftNativePort

  0  0  0  0  0  0            999 V2000
M  END
"""

    private let aliasAndValueMol = """
AliasAndValue
CDKSwiftNativePort

  2  1  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
A    1
Phenyl
V    2 Carbon note
M  END
"""

    private let isotopeMol = """
Isotopes
CDKSwiftNativePort

  3  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 D   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 T   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   1  0  0  0  0  0  0  0  0  0  0  0
M  END
"""

    private let parityMol = """
Parity
CDKSwiftNativePort

  1  0  0  0  1  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  1  0  0 15  0  0  0  0  0  0
M  END
"""

    private let queryBondMol = """
QueryBond
CDKSwiftNativePort

  3  2  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  5  0  0  1  0
  2  3  8  4  0  2  0
M  END
"""

    private let atomListMol = """
AtomList
CDKSwiftNativePort

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
  1 T    2   6   7
M  ALS   1  2 T C N
M  END
"""

    private let additionalMRecordsMol = """
AdditionalMRecords
CDKSwiftNativePort

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
M  RAD  1   1   2
M  RGP  1   1   7
M  APO  1   1   3
M  SUB  1   1   2
M  UNS  1   1   1
M  RBC  1   1   0
M  END
"""

    private let sgroupMol = """
Sgroup
CDKSwiftNativePort

  2  1  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  STY  2   1 SRU   2 DAT
M  SAL   1  2   1   2
M  SBL   1  1   1
M  SMT   1 n
M  SCN  1   1 HT
M  SDI   1  4   -1.0000   -1.0000   -1.0000    1.0000
M  SBT  1   1   1
M  SAL   2  1   2
M  SBL   2  1   1
M  SDT   2 FIELD
M  SED   2 Atom data
M  END
"""

    private let zeroCoordinateMol = """
ZeroCoordinates
CDKSwiftNativePort

  3  2  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
"""
}
