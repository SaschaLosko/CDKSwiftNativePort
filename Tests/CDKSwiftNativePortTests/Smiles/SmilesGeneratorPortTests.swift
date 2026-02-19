import XCTest
@testable import CDKSwiftNativePort

final class SmilesGeneratorPortTests: XCTestCase {
    private let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    func testGeneratesSmilesForEthanol() throws {
        let molecule = try parser.parseSmiles("CCO")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor: [.useAromaticSymbols, .strict])

        let smiles = generator.create(molecule)
        XCTAssertEqual(smiles, "CCO")
    }

    func testGeneratesIsomericSmilesWithChirality() throws {
        let molecule = try parser.parseSmiles("C[C@H](N)O")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
            flavor: [.useAromaticSymbols, .isomeric, .strict]
        )

        let isomeric = generator.create(molecule)
        XCTAssertTrue(isomeric.contains("@"))

        let reparsed = try parser.parseSmiles(isomeric)
        XCTAssertEqual(reparsed.atomCount, molecule.atomCount)
        XCTAssertEqual(reparsed.bondCount, molecule.bondCount)
    }

    func testGeneratesRingClosureNotation() throws {
        let molecule = try parser.parseSmiles("c1ccccc1")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
            flavor: [.useAromaticSymbols, .isomeric, .strict]
        )

        let smiles = generator.create(molecule)
        XCTAssertTrue(smiles.contains("1"))

        let reparsed = try parser.parseSmiles(smiles)
        XCTAssertEqual(reparsed.atomCount, 6)
        XCTAssertEqual(reparsed.bondCount, 6)
    }

    func testGeneratesDisconnectedComponents() throws {
        let molecule = try parser.parseSmiles("[Na+].[OH-]")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
            flavor: [.useAromaticSymbols, .isomeric, .strict]
        )

        let smiles = generator.create(molecule)
        XCTAssertTrue(smiles.contains("."))

        let reparsed = try parser.parseSmiles(smiles)
        XCTAssertEqual(reparsed.atomCount, 2)
        XCTAssertEqual(reparsed.bondCount, 0)
    }

    func testNonIsomericSmilesOmitsChiralityMarker() throws {
        let molecule = try parser.parseSmiles("C[C@H](N)O")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor: [.useAromaticSymbols, .strict])

        let smiles = generator.create(molecule)
        XCTAssertFalse(smiles.contains("@"))
    }
}
