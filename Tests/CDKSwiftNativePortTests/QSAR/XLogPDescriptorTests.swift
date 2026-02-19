import XCTest
@testable import CDKSwiftNativePort

final class XLogPDescriptorTests: XCTestCase {
    private let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    func testMannholdDescriptorMatchesReferenceFormula() throws {
        let benzene = try parser.parseSmiles("c1ccccc1")
        let value = CDKMannholdLogPDescriptor.calculate(for: benzene)
        XCTAssertEqual(value, 2.12, accuracy: 0.001)
    }

    func testXLogPReturnsValueForParsedMolecule() throws {
        let ethanol = try parser.parseSmiles("CCO")
        let value = CDKXLogPDescriptor.calculate(for: ethanol)
        XCTAssertNotNil(value)
    }

    func testHydrophobicChainRanksAboveAlcohol() throws {
        let octane = try parser.parseSmiles("CCCCCCCC")
        let ethanol = try parser.parseSmiles("CCO")

        let octaneValue = try XCTUnwrap(CDKXLogPDescriptor.calculate(for: octane))
        let ethanolValue = try XCTUnwrap(CDKXLogPDescriptor.calculate(for: ethanol))

        XCTAssertGreaterThan(octaneValue, ethanolValue)
    }

    func testChargedSpeciesHasLowerXLogPThanNeutralAnalog() throws {
        let neutral = try parser.parseSmiles("CN(C)C")
        let charged = try parser.parseSmiles("C[N+](C)(C)C")

        let neutralValue = try XCTUnwrap(CDKXLogPDescriptor.calculate(for: neutral))
        let chargedValue = try XCTUnwrap(CDKXLogPDescriptor.calculate(for: charged))

        XCTAssertLessThan(chargedValue, neutralValue)
    }

    func testXLogPHandlesChargedNitroMotif() throws {
        let nitro = try parser.parseSmiles("C[N+](=O)[O-]")
        let value = CDKXLogPDescriptor.calculate(for: nitro)
        XCTAssertNotNil(value)
    }
}
