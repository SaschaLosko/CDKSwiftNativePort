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
        XCTAssertTrue(identifiers.fixedHInchi.hasPrefix("InChI=1/"))
        XCTAssertEqual(identifiers.fixedHInchiKey.count, 27)
        XCTAssertFalse(identifiers.inchi.hasPrefix("Unavailable"))
        XCTAssertFalse(identifiers.inchiKey.hasPrefix("Unavailable"))
        XCTAssertFalse(identifiers.fixedHInchi.hasPrefix("Unavailable"))
        XCTAssertFalse(identifiers.fixedHInchiKey.hasPrefix("Unavailable"))
    }

    func testFixedHIdentifiersDistinguishNeutralAndZwitterionicForms() throws {
        let neutral = try parser.parseSmiles("NCC(=O)O")
        let zwitterion = try parser.parseSmiles("[NH3+]CC(=O)[O-]")

        let neutralIdentifiers = CDKMoleculeIdentifierService.compute(
            for: neutral,
            recalculateInChI: true)
        let zwitterionIdentifiers = CDKMoleculeIdentifierService.compute(
            for: zwitterion,
            recalculateInChI: true)

        XCTAssertEqual(neutralIdentifiers.inchi, zwitterionIdentifiers.inchi)
        XCTAssertEqual(neutralIdentifiers.inchiKey, zwitterionIdentifiers.inchiKey)
        XCTAssertNotEqual(neutralIdentifiers.fixedHInchi, zwitterionIdentifiers.fixedHInchi)
        XCTAssertNotEqual(neutralIdentifiers.fixedHInchiKey, zwitterionIdentifiers.fixedHInchiKey)
    }
}
