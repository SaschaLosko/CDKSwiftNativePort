import XCTest
@testable import CDKSwiftNativePort

final class CxSmilesParserPortTests: XCTestCase {

    func testParsesCxAtomLabelsLayer() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C* |$;R1$|")

        XCTAssertEqual(molecule.atomCount, 2)
        XCTAssertEqual(molecule.atoms[0].element, "C")
        XCTAssertEqual(molecule.atoms[1].element, "R1")
    }

    func testParsesCxRacemicLayerWithoutFailure() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("CC[C@H](O)[C@H](O)CCCCCC |r|")

        XCTAssertEqual(molecule.atomCount, 12)
        XCTAssertTrue(molecule.atoms.contains(where: { $0.chirality != .none }))
    }

    func testParsesCxRacemicFragmentIndicesLayer() throws {
        let split = try CDKCxSmilesParser.split("CC[C@H](O)[C@H](O)CCCCCC |r:1|", enabled: true)
        XCTAssertEqual(split.state.racemicFragments, [1])
    }

    func testParsesCxFragmentLayerWithoutFailure() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C.*.C |f:0.2|")
        XCTAssertEqual(molecule.atomCount, 3)
    }

    func testParsesCxWithTitle() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("CC |$C;C$| ethane")
        XCTAssertEqual(molecule.name, "ethane")
    }

    func testRejectsMalformedCxLayer() {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        XCTAssertThrowsError(try parser.parseSmiles("CC |$;R1$"))
    }

    func testRejectsMalformedCxRacemicFragmentLayer() {
        XCTAssertThrowsError(try CDKCxSmilesParser.split("CC |r:a|", enabled: true))
    }
}
