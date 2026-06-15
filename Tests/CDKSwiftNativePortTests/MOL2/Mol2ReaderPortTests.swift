import XCTest
@testable import CDKSwiftNativePort

final class Mol2ReaderPortTests: XCTestCase {
    func testParsesAromaticMol2BondType() throws {
        let text = """
        @<TRIPOS>MOLECULE
        aromatic_pair
         2 1 0 0 0
        SMALL
        NO_CHARGES

        @<TRIPOS>ATOM
              1 C1          0.0000    0.0000    0.0000 C.ar       1 RING       0.0000
              2 C2          1.3900    0.0000    0.0000 C.ar       1 RING       0.0000
        @<TRIPOS>BOND
             1    1    2 ar
        """

        let molecules = try CDKMol2Reader.read(text: text)
        XCTAssertEqual(molecules.count, 1)

        let molecule = molecules[0]
        XCTAssertEqual(molecule.atomCount, 2)
        XCTAssertEqual(molecule.bondCount, 1)
        XCTAssertEqual(molecule.bonds[0].order, .aromatic)
        XCTAssertTrue(molecule.atoms.allSatisfy { $0.aromatic })
        XCTAssertEqual(try XCTUnwrap(molecule.atoms[0].zPosition), 0.0, accuracy: 0.00001)
    }

    func testReadsAndWritesZCoordinatesWithoutFlattening() throws {
        let text = """
        @<TRIPOS>MOLECULE
        zed
         2 1 0 0 0
        SMALL
        USER_CHARGES

        @<TRIPOS>ATOM
              1 C1         -1.2500    0.5000    2.7500 C.3        1 MOL        0.0000
              2 O2          0.1000   -0.2500   -1.5000 O.3        1 MOL        0.0000
        @<TRIPOS>BOND
             1    1    2 1
        """

        let molecule = try XCTUnwrap(CDKMol2Reader.read(text: text).first)

        XCTAssertEqual(try XCTUnwrap(molecule.atoms[0].zPosition), 2.75, accuracy: 0.00001)
        XCTAssertEqual(try XCTUnwrap(molecule.atoms[1].zPosition), -1.50, accuracy: 0.00001)

        let roundTripText = try CDKMol2Writer.write([molecule])
        let roundTripped = try XCTUnwrap(CDKMol2Reader.read(text: roundTripText).first)

        XCTAssertEqual(try XCTUnwrap(roundTripped.atoms[0].zPosition), 2.75, accuracy: 0.00001)
        XCTAssertEqual(try XCTUnwrap(roundTripped.atoms[1].zPosition), -1.50, accuracy: 0.00001)
    }

    func testReadsMultipleMoleculeBlocks() throws {
        let text = """
        @<TRIPOS>MOLECULE
        one
         1 0 0 0 0
        SMALL
        NO_CHARGES
        @<TRIPOS>ATOM
             1 He1 0.0 0.0 0.0 He 1 ONE 0.0

        @<TRIPOS>MOLECULE
        two
         1 0 0 0 0
        SMALL
        NO_CHARGES
        @<TRIPOS>ATOM
             1 Ne1 0.0 0.0 0.0 Ne 1 TWO 0.0
        """

        let molecules = try CDKMol2Reader.read(text: text)
        XCTAssertEqual(molecules.count, 2)
        XCTAssertEqual(molecules[0].name, "one")
        XCTAssertEqual(molecules[1].name, "two")
    }

    func testRejectsInputWithoutTriposMoleculeSection() {
        XCTAssertThrowsError(try CDKMol2Reader.read(text: "@<TRIPOS>ATOM\n"))
    }
}
