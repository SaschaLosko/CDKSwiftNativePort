import XCTest
@testable import CDKSwiftNativePort

final class SmilesReaderPortTests: XCTestCase {
    func testReadsMultipleSmilesLinesWithNames() throws {
        let text = """
        CCO ethanol
        c1ccccc1 benzene
        """

        let molecules = try CDKSMILESReader.read(text: text)
        XCTAssertEqual(molecules.count, 2)
        XCTAssertEqual(molecules[0].name, "ethanol")
        XCTAssertEqual(molecules[1].name, "benzene")
        XCTAssertEqual(molecules[0].atomCount, 3)
        XCTAssertEqual(molecules[1].atomCount, 6)
    }

    func testIgnoresBlankAndCommentLines() throws {
        let text = """

        # comment
        // another comment
        C
        """

        let molecules = try CDKSMILESReader.read(text: text)
        XCTAssertEqual(molecules.count, 1)
        XCTAssertEqual(molecules[0].atomCount, 1)
    }

    func testThrowsOnNoMolecules() {
        let text = """

        # comment
        // comment
        """

        XCTAssertThrowsError(try CDKSMILESReader.read(text: text))
    }
}
