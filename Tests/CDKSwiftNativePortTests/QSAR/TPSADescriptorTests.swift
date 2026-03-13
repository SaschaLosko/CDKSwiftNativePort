import XCTest
@testable import CDKSwiftNativePort

final class TPSADescriptorTests: XCTestCase {
    private let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    func testCDKReferenceCases() throws {
        try assertTPSA("O=C(O)CC", expected: 37.29)
        try assertTPSA("C=NC(CC#N)N(C)C", expected: 39.39)
        try assertTPSA("CCCN(=O)=O", expected: 45.82)
        try assertTPSA("C#N=CC(CNC)N1CC1", expected: 28.632)
        try assertTPSA("c1ccncc1", expected: 12.892)
        try assertTPSA("[H][N+]([H])(C)C", expected: 16.61, tolerance: 0.1)
        try assertTPSA("C(I)I", expected: 0.0, tolerance: 0.1)
        try assertTPSA("C(O)O", expected: 40.45, tolerance: 0.1)
    }

    func testDescriptorHandlesRingSystems() throws {
        let molecule = try parser.parseSmiles("C1CCCC1CCC2CCCNC2")
        let value = CDKTPSADescriptor.calculate(for: molecule, checkAromaticity: true)
        XCTAssertTrue(value.isFinite)
    }

    private func assertTPSA(_ smiles: String,
                            expected: Double,
                            tolerance: Double = 0.01,
                            file: StaticString = #filePath,
                            line: UInt = #line) throws {
        let molecule = try parser.parseSmiles(smiles)
        let observed = CDKTPSADescriptor.calculate(for: molecule, checkAromaticity: true)
        XCTAssertEqual(observed, expected, accuracy: tolerance,
                       "SMILES \(smiles) expected TPSA \(expected), got \(observed)",
                       file: file,
                       line: line)
    }
}
