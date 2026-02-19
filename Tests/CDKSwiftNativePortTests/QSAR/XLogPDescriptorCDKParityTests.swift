import XCTest
@testable import CDKSwiftNativePort

final class XLogPDescriptorCDKParityTests: XCTestCase {
    private let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    func testCDKReferenceCases() throws {
        try assertXLogP("Nc2ccc(S(=O)(=O)c1ccc(N)cc1)cc2", expected: 0.86, tolerance: 1.0) // no1596
        try assertXLogP("O=C(O)C(N)CCCN", expected: -3.30, tolerance: 0.1) // no367
        try assertXLogP("O=P(N1CC1)(N2CC2)N3CC3", expected: -1.19, tolerance: 0.1) // no1837
        try assertXLogP("c1cc2ccc3ccc4ccc5cccc6c(c1)c2c3c4c56", expected: 7.00, tolerance: 0.1) // no87
        try assertXLogP("S1C2N(C(=O)C2NC(=O)C(c2ccccc2)C(=O)O)C(C(=O)O)C1(C)C", expected: 1.84, tolerance: 0.1) // no1782
        try assertXLogP("C(#Cc1ccccc1)c1ccccc1", expected: 4.62, tolerance: 0.1) // no30
        try assertXLogP("FC(F)(F)c1ccc(cc1)C(=O)N", expected: 1.834, tolerance: 1.0) // no990
        try assertXLogP("Clc1cccc(c1)/C=C/[N+](=O)[O-]", expected: 2.809, tolerance: 1.0) // no1000
        try assertXLogP("CC(=O)OC1=CC=CC=C1C(=O)O", expected: 1.422, tolerance: 0.1) // aspirin
        try assertXLogP("O=C(OC)CNC(=O)c1ccc(N)cc1", expected: 0.31, tolerance: 1.0) // no1429
        try assertXLogP("O=[N+]([O-])c1ccc(cc1)CC(N)C(=O)O", expected: -1.487, tolerance: 1.0) // no1274
        try assertXLogP("O=C1NC(=O)C=CN1C1OC(CO)C(O)C1O", expected: -2.11, tolerance: 0.1) // no454
        try assertXLogP("O=C1N(C)C=CC(=O)N1C", expected: -0.59, tolerance: 0.1) // no498
        try assertXLogP("CCN(CC)CCCN(C2Cc1ccccc1C2)c3ccccc3", expected: 5.03, tolerance: 1.0) // Aprindine
        try assertXLogP("Brc1cc(Cl)c(O[P+]([S-])(OC)OC)cc1Cl", expected: 5.22, tolerance: 1.0) // no1844
        try assertXLogP("Clc1ccc2Sc3ccccc3N(CCCN3CCN(C)CC3)c2c1", expected: 4.56, tolerance: 1.0) // no1810
        try assertXLogP("[S+]([O-])(CCC1C(=O)N(N(c2ccccc2)C1=O)c1ccccc1)c1ccccc1", expected: 2.36, tolerance: 0.1) // no1822
    }

    func testCDKBenzeneAromaticityModes() throws {
        let benzene = try parser.parseSmiles("C1=CC=CC=C1")
        let aromatic = try XCTUnwrap(CDKXLogPDescriptor.calculate(for: benzene,
                                                                  checkAromaticity: true,
                                                                  salicylFlag: true))
        let nonAromatic = try XCTUnwrap(CDKXLogPDescriptor.calculate(for: benzene,
                                                                     checkAromaticity: false,
                                                                     salicylFlag: true))

        XCTAssertEqual(aromatic, 2.02, accuracy: 0.01)
        XCTAssertEqual(nonAromatic, 2.08, accuracy: 0.01)
    }

    private func assertXLogP(_ smiles: String,
                             expected: Double,
                             tolerance: Double,
                             file: StaticString = #filePath,
                             line: UInt = #line) throws {
        let molecule = try parser.parseSmiles(smiles)
        let observed = try XCTUnwrap(CDKXLogPDescriptor.calculate(for: molecule,
                                                                  checkAromaticity: true,
                                                                  salicylFlag: false),
                                     file: file,
                                     line: line)
        XCTAssertEqual(observed, expected, accuracy: tolerance,
                       "SMILES \(smiles) expected \(expected), got \(observed)",
                       file: file,
                       line: line)
    }
}
