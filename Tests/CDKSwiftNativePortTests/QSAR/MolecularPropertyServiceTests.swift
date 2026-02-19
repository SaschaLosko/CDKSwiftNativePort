import XCTest
@testable import CDKSwiftNativePort

final class MolecularPropertyServiceTests: XCTestCase {
    private let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    func testComputesFormulaMassAndLipinskiCountsForEthanol() throws {
        let molecule = try parser.parseSmiles("CCO")
        let properties = CDKMoleculePropertyService.compute(for: molecule)

        XCTAssertEqual(properties.formula, "C2H6O")
        XCTAssertEqual(properties.heavyAtomCount, 3)
        XCTAssertEqual(properties.hBondDonorCount, 1)
        XCTAssertEqual(properties.hBondAcceptorCount, 1)
        XCTAssertEqual(properties.rotatableBondCount, 0)
        XCTAssertEqual(properties.ringCount, 0)

        XCTAssertEqual(properties.molarMass, properties.molecularWeight, accuracy: 1e-9)
        XCTAssertEqual(properties.molecularWeight, 46.068, accuracy: 0.03)
        XCTAssertEqual(properties.monoisotopicMass, 46.0419, accuracy: 0.05)

        XCTAssertNotNil(properties.ruleOfFive.xlogP)
        XCTAssertEqual(properties.ruleOfFive.isCompliant, true)
        XCTAssertEqual(properties.ruleOfFive.violations, 0)
        XCTAssertEqual(properties.ruleOfFive.statusText, "Pass (0 violations)")
    }

    func testAmideNitrogenIsNotCountedAsAcceptor() throws {
        let molecule = try parser.parseSmiles("CC(=O)N")
        let properties = CDKMoleculePropertyService.compute(for: molecule)

        XCTAssertEqual(properties.hBondDonorCount, 1)
        XCTAssertEqual(properties.hBondAcceptorCount, 1)
    }

    func testAromaticFormulaCountsImplicitHydrogenCorrectly() throws {
        let benzene = try parser.parseSmiles("c1ccccc1")
        let benzeneProperties = CDKMoleculePropertyService.compute(for: benzene)
        XCTAssertEqual(benzeneProperties.formula, "C6H6")

        let pyridine = try parser.parseSmiles("n1ccccc1")
        let pyridineProperties = CDKMoleculePropertyService.compute(for: pyridine)
        XCTAssertEqual(pyridineProperties.formula, "C5H5N")
    }

    func testRotatableBondAndRingCountDescriptors() throws {
        let butane = try parser.parseSmiles("CCCC")
        XCTAssertEqual(CDKRotatableBondsCountDescriptor.calculate(for: butane), 1)
        XCTAssertEqual(CDKRingCountDescriptor.calculate(for: butane), 0)

        let benzene = try parser.parseSmiles("c1ccccc1")
        XCTAssertEqual(CDKRotatableBondsCountDescriptor.calculate(for: benzene), 0)
        XCTAssertEqual(CDKRingCountDescriptor.calculate(for: benzene), 1)
    }

    func testRuleOfFiveDescriptorWithExplicitXlogP() {
        let pass = CDKRuleOfFiveDescriptor.evaluate(molecularWeight: 320.5,
                                                    hBondDonorCount: 2,
                                                    hBondAcceptorCount: 6,
                                                    xlogP: 3.8)
        XCTAssertEqual(pass.violations, 0)
        XCTAssertEqual(pass.evaluatedCriteriaCount, 4)
        XCTAssertEqual(pass.isCompliant, true)
        XCTAssertEqual(pass.statusText, "Pass (0 violations)")

        let fail = CDKRuleOfFiveDescriptor.evaluate(molecularWeight: 710.2,
                                                    hBondDonorCount: 7,
                                                    hBondAcceptorCount: 14,
                                                    xlogP: 7.1)
        XCTAssertEqual(fail.violations, 4)
        XCTAssertEqual(fail.evaluatedCriteriaCount, 4)
        XCTAssertEqual(fail.isCompliant, false)
        XCTAssertEqual(fail.statusText, "Fail (4 violations)")
    }

    func testExactMassRespectsIsotopicLabels() throws {
        let molecule = try parser.parseSmiles("[13CH4]")
        let exactMass = CDKExactMassDescriptor.calculate(for: molecule)
        XCTAssertGreaterThan(exactMass, 17.0)
        XCTAssertLessThan(exactMass, 18.0)
    }
}
