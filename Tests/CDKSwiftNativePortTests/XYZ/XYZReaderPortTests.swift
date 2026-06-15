import XCTest
@testable import CDKSwiftNativePort

final class XYZReaderPortTests: XCTestCase {
    func testReadsSingleXYZBlockAndInfersBonds() throws {
        let text = """
        3
        water
        O 0.0000 0.0000 0.0000
        H 0.9584 0.0000 0.0000
        H -0.2396 0.9271 0.0000
        """

        let molecules = try CDKXYZReader.read(text: text)
        XCTAssertEqual(molecules.count, 1)

        let molecule = molecules[0]
        XCTAssertEqual(molecule.name, "water")
        XCTAssertEqual(molecule.atomCount, 3)
        XCTAssertEqual(molecule.bondCount, 2)
        XCTAssertEqual(try XCTUnwrap(molecule.atoms[0].zPosition), 0.0, accuracy: 0.00001)
    }

    func testReadsAndWritesZCoordinatesWithoutFlattening() throws {
        let text = """
        2
        zed
        C -1.2500 0.5000 2.7500
        O 0.1000 -0.2500 -1.5000
        """

        let molecule = try XCTUnwrap(CDKXYZReader.read(text: text).first)

        XCTAssertEqual(try XCTUnwrap(molecule.atoms[0].zPosition), 2.75, accuracy: 0.00001)
        XCTAssertEqual(try XCTUnwrap(molecule.atoms[1].zPosition), -1.50, accuracy: 0.00001)

        let roundTripText = try CDKXYZWriter.write([molecule])
        let roundTripped = try XCTUnwrap(CDKXYZReader.read(text: roundTripText).first)

        XCTAssertEqual(try XCTUnwrap(roundTripped.atoms[0].zPosition), 2.75, accuracy: 0.00001)
        XCTAssertEqual(try XCTUnwrap(roundTripped.atoms[1].zPosition), -1.50, accuracy: 0.00001)
    }

    func testReadsNumericFirstXYZRowsAsXYZSymbol() throws {
        let text = """
        1
        numeric-first
        -1.2500 0.5000 2.7500 C
        """

        let molecule = try XCTUnwrap(CDKXYZReader.read(text: text).first)

        XCTAssertEqual(molecule.atoms[0].element, "C")
        XCTAssertEqual(molecule.atoms[0].position.x, -1.25, accuracy: 0.00001)
        XCTAssertEqual(molecule.atoms[0].position.y, 0.50, accuracy: 0.00001)
        XCTAssertEqual(try XCTUnwrap(molecule.atoms[0].zPosition), 2.75, accuracy: 0.00001)
    }

    func testBondPerceptionUsesZDistanceWhenCoordinatesAre3D() throws {
        let text = """
        2
        separated-on-z
        C 0.0 0.0 0.0
        O 0.0 0.0 4.0
        """

        let molecule = try XCTUnwrap(CDKXYZReader.read(text: text).first)

        XCTAssertEqual(molecule.bondCount, 0)
        XCTAssertEqual(try XCTUnwrap(molecule.atoms[1].zPosition), 4.0, accuracy: 0.00001)
    }

    func testReadsConcatenatedXYZBlocks() throws {
        let text = """
        1
        helium
        He 0.0 0.0 0.0
        1
        neon
        Ne 1.0 0.0 0.0
        """

        let molecules = try CDKXYZReader.read(text: text)
        XCTAssertEqual(molecules.count, 2)
        XCTAssertEqual(molecules[0].name, "helium")
        XCTAssertEqual(molecules[1].name, "neon")
    }

    func testThrowsOnTruncatedXYZBlock() {
        let text = """
        2
        bad
        C 0 0 0
        """

        XCTAssertThrowsError(try CDKXYZReader.read(text: text))
    }
}
