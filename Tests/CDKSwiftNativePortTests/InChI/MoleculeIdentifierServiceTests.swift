import XCTest
@testable import CDKSwiftNativePort

final class MoleculeIdentifierServiceTests: XCTestCase {
    private let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    func testComputesSmilesAndIsoSmiles() throws {
        let molecule = try parser.parseSmiles("CCO")
        let identifiers = CDKMoleculeIdentifierService.compute(for: molecule)

        XCTAssertEqual(identifiers.smiles, "CCO")
        XCTAssertEqual(identifiers.isoSmiles, "CCO")
    }

    func testInchiFieldsReflectBackendAvailability() throws {
        let molecule = try parser.parseSmiles("CCO")
        let identifiers = CDKMoleculeIdentifierService.compute(for: molecule)

        XCTAssertTrue(identifiers.inchi.hasPrefix("InChI=1S/"))
        XCTAssertEqual(identifiers.inchiKey.count, 27)
        XCTAssertFalse(identifiers.inchi.hasPrefix("Unavailable"))
        XCTAssertFalse(identifiers.inchiKey.hasPrefix("Unavailable"))
    }
}
