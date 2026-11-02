import XCTest
@testable import CDKSwiftNativePort

final class SmilesReaderPortTests: XCTestCase {
    func testGraphOnlyRecordsPreserveTitlesStereoAndFamilyDefinitions() throws {
        let text = """
            # retained record semantics without generated layout
            N[C@@H](C)C(=O)O alanine
            C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}| family
            """
        let graphs = try CDKSMILESReader.read(text: text, generateCoordinates: false)
        let drawn = try CDKSMILESReader.read(text: text)
        XCTAssertEqual(graphs.map(\.name), ["alanine", "family"])
        XCTAssertTrue(graphs[0].atoms.allSatisfy { $0.position == .zero })
        XCTAssertTrue(drawn[0].atoms.contains { $0.position != .zero })
        XCTAssertEqual(CDKMarkushMatcher.match(graphs[0], in: drawn[0]).verdict, .match)
        let opposite = try CDKSmilesParser(generateCoordinates: false).parseSmiles("N[C@H](C)C(=O)O")
        XCTAssertEqual(CDKMarkushMatcher.match(opposite, in: graphs[0]).verdict, .noMatch)
        XCTAssertEqual(CDKMarkushRelations.compare(graphs[1], drawn[1], relation: .equivalent).verdict, .holds)
        let coordinates = try CDKSMILESReader.read(text: "CO |(1,2,;3,4,)|", generateCoordinates: false)
        XCTAssertEqual(coordinates[0].atoms[0].position, CGPoint(x: 1, y: 2))
        XCTAssertEqual(coordinates[0].atoms[1].position, CGPoint(x: 3, y: 4))
    }

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
