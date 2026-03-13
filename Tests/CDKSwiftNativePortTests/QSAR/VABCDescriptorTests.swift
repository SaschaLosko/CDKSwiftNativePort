import XCTest
@testable import CDKSwiftNativePort

final class VABCDescriptorTests: XCTestCase {
    private let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    func testUnsupportedElementReturnsNil() throws {
        let ironChloride = try parser.parseSmiles("Cl[Fe]Cl")
        XCTAssertNil(CDKVABCDescriptor.calculate(for: ironChloride))
    }

    func testCDKReferenceCases() throws {
        try assertVABC("C", expected: 25.8524433266667)
        try assertVABC("[H]C([H])([H])[H]", expected: 25.8524433266667)
        try assertVABC("CC", expected: 43.1484279525333)
        try assertVABC("CCCC", expected: 77.7403972042667)
        try assertVABC("CC#N", expected: 48.8722707591)
        try assertVABC("CC(=O)O", expected: 58.0924226528555)
        try assertVABC("CC(F)(F)Cl", expected: 70.4946134235795)
        try assertVABC("S=C=S", expected: 57.5975740402667)
        try assertVABC("CCOP(=O)(OCC)OCC", expected: 167.320526666244)
        try assertVABC("c1ccccc1", expected: 81.1665316528)
        try assertVABC("c1cc2ccc3cccc4ccc(c1)c2c34", expected: 171.174708305067)
        try assertVABC("CN1CCCC1c2cccnc2", expected: 159.9875318718)
        try assertVABC("COc2ccc1[nH]c(nc1c2)S(=O)Cc3ncc(C)c(OC)c3C", expected: 292.23)
        try assertVABC("O=C1NS(=O)(=O)c2ccccc12", expected: 139.35)
        try assertVABC("Nc1ncnc2n(CCOCP(=O)(O)O)cnc12", expected: 199.84)
    }

    private func assertVABC(_ smiles: String,
                            expected: Double,
                            tolerance: Double = 0.01,
                            file: StaticString = #filePath,
                            line: UInt = #line) throws {
        let molecule = try parser.parseSmiles(smiles)
        let observed = try XCTUnwrap(CDKVABCDescriptor.calculate(for: molecule),
                                     "Expected VABC value for \(smiles)",
                                     file: file,
                                     line: line)
        XCTAssertEqual(observed, expected, accuracy: tolerance,
                       "SMILES \(smiles) expected VABC \(expected), got \(observed)",
                       file: file,
                       line: line)
    }
}
