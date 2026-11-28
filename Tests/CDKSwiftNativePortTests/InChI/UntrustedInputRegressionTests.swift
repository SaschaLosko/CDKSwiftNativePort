import Foundation
import XCTest

@testable import CDKSwiftNativePort

/// Damaged or hostile identifiers and molfiles must be rejected or read quickly,
/// never trap or allocate without bound. Found by seeded mutation fuzzing.
final class UntrustedInputRegressionTests: XCTestCase {
    func testInChIFormulaAboveTheInChIAtomLimitIsRejected() {
        for inchi in ["InChI=1S/C891604/c1-5-H-6(2)", "InChI=1S/C999999999H999999999", "InChI=1S/999999999CH4"] {
            XCTAssertThrowsError(try CDKFileImporter.readMolecules(text: inchi, fileExtension: "inchi"), inchi)
        }
    }

    func testOversizedInChINumbersNeitherTrapNorIterateWithoutBound() {
        let inputs = [
            "InChI=1S/C2H6/c1-99999999999999999999999",
            "InChI=1S/C2H6/c1-2/h1-999999999999H3",
            "InChI=1S/CH4/h1H4/i/hD9223372036854775807D9223372036854775807",
            "InChI=1S/CH4/h1H4/p999999999*+1",
            "InChI=1S/C2H6/c1-\u{0663}/h1-2H3",
        ]
        for inchi in inputs {
            _ = try? CDKFileImporter.readMolecules(text: inchi, fileExtension: "inchi")
        }
    }

    func testMolfileWithBondsButNoAtomsIsRejected() {
        let molfile = "\n  CDK\n\n  0  1  0  0  0  0  0  0  0  0999 V2000\n  1  2  1  0\nM  END\n"
        XCTAssertThrowsError(try CDKFileImporter.readMolecules(text: molfile, fileExtension: "mol"))
    }

    /// One explicit hydrogen must not switch off implicit hydrogens on the other
    /// atoms, and a ring double bond carries no cis/trans configuration.
    func testPartiallyExplicitHydrogensKeepImplicitHydrogensForInChI() throws {
        // Cyclohexene with one explicit hydrogen on C3.
        let molfile = """

              CDK

              7  7  0  0  0  0  0  0  0  0999 V2000
                0.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                0.0000   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
               -0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
               -0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                1.7320   -1.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
              1  2  2  0
              2  3  1  0
              3  4  1  0
              4  5  1  0
              5  6  1  0
              6  1  1  0
              3  7  1  0
            M  END

            """
        let molecule = try XCTUnwrap(CDKFileImporter.readMolecules(text: molfile, fileExtension: "mol").first)
        XCTAssertTrue(molecule.bonds.allSatisfy { $0.doubleBondStereo == nil })
        let identifiers = CDKMoleculeIdentifierService.compute(for: molecule)
        XCTAssertEqual(identifiers.inchi, "InChI=1S/C6H10/c1-2-4-6-5-3-1/h1-2H,3-6H2")
    }
}
