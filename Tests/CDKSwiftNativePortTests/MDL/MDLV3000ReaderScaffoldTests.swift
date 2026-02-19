import XCTest
@testable import CDKSwiftNativePort

final class MDLV3000ReaderScaffoldTests: XCTestCase {

    func testParsesV3000Scaffold() throws {
        let scaffold = try CDKMDLV3000Reader.parseScaffold(text: minimalV3000)

        XCTAssertEqual(scaffold.title, "V3000Example")
        XCTAssertEqual(scaffold.atomCount, 2)
        XCTAssertEqual(scaffold.bondCount, 1)
        XCTAssertEqual(scaffold.atomRecords.count, 2)
        XCTAssertEqual(scaffold.bondRecords.count, 1)
        XCTAssertEqual(scaffold.atomRecords[0].symbol, "C")
        XCTAssertEqual(scaffold.bondRecords[0].type, "1")
        XCTAssertEqual(scaffold.bondRecords[0].a1, 1)
        XCTAssertEqual(scaffold.bondRecords[0].a2, 2)
    }

    func testReadParsesBasicV3000IntoMolecule() throws {
        let molecule = try CDKMDLV3000Reader.read(text: minimalV3000)

        XCTAssertEqual(molecule.name, "V3000Example")
        XCTAssertEqual(molecule.atomCount, 2)
        XCTAssertEqual(molecule.bondCount, 1)
        XCTAssertEqual(molecule.atoms.map(\.element), ["C", "O"])
        XCTAssertEqual(molecule.bonds[0].order, .single)
        XCTAssertEqual(molecule.bonds[0].a1, 1)
        XCTAssertEqual(molecule.bonds[0].a2, 2)
    }

    func testMDLReaderDispatchesV3000ToMolecule() throws {
        let molecule = try CDKMDLReader.read(text: minimalV3000)
        XCTAssertEqual(molecule.atomCount, 2)
        XCTAssertEqual(molecule.bondCount, 1)
    }

    func testParsesV3000AtomAndBondAttributes() throws {
        let molecule = try CDKMDLV3000Reader.read(text: attributedV3000)
        XCTAssertEqual(molecule.atomCount, 4)
        XCTAssertEqual(molecule.bondCount, 3)

        let atom1 = try XCTUnwrap(molecule.atoms.first(where: { $0.id == 1 }))
        XCTAssertEqual(atom1.charge, 1)
        XCTAssertEqual(atom1.isotopeMassNumber, 13)
        XCTAssertEqual(atom1.explicitHydrogenCount, 0)

        let atom2 = try XCTUnwrap(molecule.atoms.first(where: { $0.id == 2 }))
        XCTAssertEqual(atom2.radical, 2)

        let atom3 = try XCTUnwrap(molecule.atoms.first(where: { $0.id == 3 }))
        XCTAssertEqual(atom3.element, "L")
        XCTAssertEqual(atom3.queryType, .anyAtom)
        XCTAssertEqual(atom3.atomList ?? [], ["N", "O"])

        let atom4 = try XCTUnwrap(molecule.atoms.first(where: { $0.id == 4 }))
        XCTAssertEqual(atom4.queryType, .anyNonHydrogen)

        let bond1 = try XCTUnwrap(molecule.bonds.first(where: { $0.id == 1 }))
        XCTAssertEqual(bond1.queryType, .singleOrAromatic)
        XCTAssertEqual(bond1.stereo, .up)

        let bond2 = try XCTUnwrap(molecule.bonds.first(where: { $0.id == 2 }))
        XCTAssertEqual(bond2.queryType, .doubleOrAromatic)

        let bond3 = try XCTUnwrap(molecule.bonds.first(where: { $0.id == 3 }))
        XCTAssertEqual(bond3.queryType, .any)
    }

    func testParsesContinuationLines() throws {
        let molecule = try CDKMDLV3000Reader.read(text: continuationV3000)
        let atom = try XCTUnwrap(molecule.atoms.first(where: { $0.id == 1 }))

        XCTAssertEqual(atom.element, "L")
        XCTAssertEqual(atom.atomList ?? [], ["C", "N", "O"])
        XCTAssertEqual(atom.queryType, .anyAtom)
    }

    func testGeneratesLayoutWhenCoordinatesAreMissing() throws {
        let molecule = try CDKMDLV3000Reader.read(text: zeroCoordinateV3000)
        let box = try XCTUnwrap(molecule.boundingBox())
        XCTAssertGreaterThan(box.width, 0.1)
    }

    func testParseScaffoldRejectsMissingCTAB() {
        XCTAssertThrowsError(try CDKMDLV3000Reader.parseScaffold(text: invalidV3000))
    }

    func testAppliesCollectionAndSGroupSemantics() throws {
        let scaffold = try CDKMDLV3000Reader.parseScaffold(text: collectionAndSGroupV3000)
        XCTAssertEqual(scaffold.collectionLines.count, 1)
        XCTAssertEqual(scaffold.sgroupLines.count, 2)

        let molecule = try CDKMDLV3000Reader.read(text: collectionAndSGroupV3000)
        let bond = try XCTUnwrap(molecule.bonds.first(where: { $0.id == 1 }))
        XCTAssertEqual(bond.order, .double)
        XCTAssertEqual(bond.stereo, .either)

        let atom = try XCTUnwrap(molecule.atoms.first(where: { $0.id == 1 }))
        XCTAssertEqual(atom.element, "Me")
        XCTAssertEqual(atom.attachmentPoint, 1)
        XCTAssertTrue(molecule.name.contains("AliasName"))
    }

    private let minimalV3000 = """
V3000Example
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0
M  V30 2 O 1.2000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
"""

    private let attributedV3000 = """
V3000Attributes
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0 CHG=1 MASS=13 HCOUNT=0
M  V30 2 O 1.2000 0.0000 0.0000 0 RAD=2
M  V30 3 [N,O] 2.4000 0.0000 0.0000 0
M  V30 4 A 3.6000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 6 1 2 CFG=1
M  V30 2 7 2 3
M  V30 3 8 3 4
M  V30 END BOND
M  V30 END CTAB
M  END
"""

    private let continuationV3000 = """
V3000Continuation
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 L 0.0000 0.0000 0.0000 0 ATOMLIST=(3 C N -
M  V30 O)
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
"""

    private let zeroCoordinateV3000 = """
V3000ZeroCoords
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0
M  V30 2 C 0.0000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
"""

    private let invalidV3000 = """
V3000Invalid
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  END
"""

    private let collectionAndSGroupV3000 = """
V3000CollectionAndSGroup
CDKSwiftNativePort

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0
M  V30 2 C 1.3000 0.0000 0.0000 0
M  V30 3 O 2.6000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS BONDS=(1 1)
M  V30 END COLLECTION
M  V30 BEGIN SGROUP
M  V30 1 SUP 0 ATOMS=(1 1) LABEL="Me" ATTCHPT=1
M  V30 2 DAT 0 ATOMS=(1 3) FIELDNAME="CDK" FIELDDATA="AliasName"
M  V30 END SGROUP
M  V30 END CTAB
M  END
"""
}
