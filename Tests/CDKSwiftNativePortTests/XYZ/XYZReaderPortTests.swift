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
