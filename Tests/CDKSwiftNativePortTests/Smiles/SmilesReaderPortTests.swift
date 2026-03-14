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

    func testReadsCxSmilesLineWithoutDroppingCxSuffix() throws {
        let text = "C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|"

        let molecules = try CDKSMILESReader.read(text: text)
        let molecule = try XCTUnwrap(molecules.first)

        XCTAssertEqual(molecules.count, 1)
        XCTAssertEqual(molecule.name, "")
        XCTAssertEqual(molecule.sgroups.count, 1)
        XCTAssertEqual(molecule.sgroups.first?.subscriptText, "1-3")
        XCTAssertEqual(molecule.atoms.filter { $0.rGroupMembership == nil && $0.symbolToDraw == "R1" }.count, 1)
    }

    func testThrowsOnNoMolecules() {
        let text = """

        # comment
        // comment
        """

        XCTAssertThrowsError(try CDKSMILESReader.read(text: text))
    }
}
