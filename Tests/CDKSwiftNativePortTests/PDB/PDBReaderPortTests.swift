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

    func testReadsPDBZCoordinatesAndUsesThemFor3DScene() throws {
        let text = """
        HEADER    TEST 3D LIGAND
        MODEL        1
        ATOM      1  N   LIG A   1      -6.778  -1.424   4.200  1.00  0.00           N
        ATOM      2  CA  LIG A   1      -6.878  -0.708   2.896  1.00  0.00           C
        ATOM      3  C   LIG A   1      -5.557  -0.840   2.138  1.00  0.00           C
        ENDMDL
        MODEL        2
        ATOM      1  N   LIG A   1       9.000   9.000   9.000  1.00  0.00           N
        ENDMDL
        END
        """

        let molecule = try XCTUnwrap(CDKPDBReader.read(text: text).first)

        XCTAssertEqual(molecule.atomCount, 3)
        XCTAssertEqual(try XCTUnwrap(molecule.atoms[0].zPosition), 4.200, accuracy: 0.0001)
        XCTAssertEqual(try XCTUnwrap(molecule.atoms[1].zPosition), 2.896, accuracy: 0.0001)
        XCTAssertEqual(try XCTUnwrap(molecule.atoms[2].zPosition), 2.138, accuracy: 0.0001)

        let scene = CDKMetal3DSceneBuilder.build(molecule: molecule)
        XCTAssertTrue(scene.hasExplicit3DCoordinates)
        XCTAssertEqual(try XCTUnwrap(scene.atoms.first { $0.id == 1 }?.center.z), 4.200, accuracy: 0.0001)
        XCTAssertEqual(try XCTUnwrap(scene.atoms.first { $0.id == 2 }?.center.z), 2.896, accuracy: 0.0001)
    }

    func testWritesPDBZCoordinatesAndRoundTripsThem() throws {
        let molecule = Molecule(
            name: "3D fragment",
            atoms: [
                Atom(id: 1, element: "C", position: CGPoint(x: -1.25, y: 0.50), zPosition: 2.75),
                Atom(id: 2, element: "O", position: CGPoint(x: 0.10, y: -0.25), zPosition: -1.50),
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single),
            ])

        let text = try CDKPDBWriter.write([molecule])
        let parsed = try XCTUnwrap(CDKPDBReader.read(text: text).first)

        XCTAssertEqual(try XCTUnwrap(parsed.atoms[0].zPosition), 2.75, accuracy: 0.0001)
        XCTAssertEqual(try XCTUnwrap(parsed.atoms[1].zPosition), -1.50, accuracy: 0.0001)
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
