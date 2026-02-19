import XCTest
@testable import CDKSwiftNativePort

final class PDBReaderPortTests: XCTestCase {
    func testReadsPDBWithConectBonds() throws {
        let text = """
        HEADER    TEST LIGAND
        ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C
        ATOM      2  O1  LIG A   1       1.220   0.000   0.000  1.00  0.00           O
        CONECT    1    2
        END
        """

        let molecules = try CDKPDBReader.read(text: text)
        XCTAssertEqual(molecules.count, 1)

        let molecule = molecules[0]
        XCTAssertEqual(molecule.atomCount, 2)
        XCTAssertEqual(molecule.bondCount, 1)
        XCTAssertEqual(molecule.atoms[1].element.uppercased(), "O")
    }

    func testInfersBondsWithoutConect() throws {
        let text = """
        ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C
        ATOM      2  C2  LIG A   1       1.450   0.000   0.000  1.00  0.00           C
        END
        """

        let molecules = try CDKPDBReader.read(text: text)
        XCTAssertEqual(molecules.count, 1)
        XCTAssertEqual(molecules[0].atomCount, 2)
        XCTAssertEqual(molecules[0].bondCount, 1)
    }

    func testThrowsOnMissingAtoms() {
        XCTAssertThrowsError(try CDKPDBReader.read(text: "HEADER only\nEND\n"))
    }
}
