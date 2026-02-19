import XCTest
@testable import CDKSwiftNativePort

final class InChIReaderPortTests: XCTestCase {
    func testReadsMultipleInChILinesWithNames() throws {
        let text = """
        InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3 ethanol
        InChI=1S/CH4/h1H4 methane
        """

        let molecules = try CDKInChIReader.read(text: text)
        XCTAssertEqual(molecules.count, 2)
        XCTAssertEqual(molecules[0].name, "ethanol")
        XCTAssertEqual(molecules[1].name, "methane")
        XCTAssertGreaterThan(molecules[0].atomCount, 0)
        XCTAssertGreaterThan(molecules[1].atomCount, 0)
    }

    func testIgnoresBlankAndCommentLines() throws {
        let text = """
        # comment
        // comment
        InChI=1S/H2O/h1H2 water
        """

        let molecules = try CDKInChIReader.read(text: text)
        XCTAssertEqual(molecules.count, 1)
        XCTAssertEqual(molecules[0].name, "water")
    }

    func testThrowsOnNoMolecules() {
        XCTAssertThrowsError(try CDKInChIReader.read(text: "   \n# only comments\n"))
    }
}
