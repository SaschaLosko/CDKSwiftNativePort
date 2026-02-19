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
