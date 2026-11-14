import Foundation
import XCTest
@testable import CDKSwiftNativePort

final class IsotopeTests: XCTestCase {
    func testPublishedCommonAndRareIsotopeValues() throws {
        // CDK 2.14.1 Isotopes static tables / BODR 10. These are independent
        // expected values, not values reconstructed from the Swift dataset.
        for (element, number, mass) in [
            ("H", 2, 2.014101778), ("H", 3, 3.016049278), ("C", 13, 13.00335484),
            ("N", 15, 15.0001089), ("O", 18, 17.999161), ("Cu", 64, 63.9297642),
            ("Tc", 99, 98.9062547), ("U", 238, 238.0507882),
        ] {
            XCTAssertEqual(try CDKIsotopes.mass(element: element, massNumber: number), mass, accuracy: 1e-9)
        }
        XCTAssertEqual(try CDKIsotopes.averageMass(atomicNumber: 6), 12.01074)
        XCTAssertEqual(try CDKIsotopes.majorIsotope(atomicNumber: 29)?.massNumber, 63)
        XCTAssertEqual(try CDKIsotopes.isotope(atomicNumber: 6, massNumber: 13)?.naturalAbundance, 1.07)
    }

    func testCompleteDatasetAndDeterministicElementLookup() throws {
        var total = 0
        for number in 1...118 {
            let isotopes = try CDKIsotopes.isotopes(atomicNumber: number)
            XCTAssertEqual(isotopes.map(\.massNumber), isotopes.map(\.massNumber).sorted())
            XCTAssertEqual(Set(isotopes.map(\.massNumber)).count, isotopes.count)
            for isotope in isotopes {
                XCTAssertEqual(isotope.atomicNumber, number)
                XCTAssertEqual(try CDKIsotopes.isotope(atomicNumber: number, massNumber: isotope.massNumber), isotope)
                XCTAssertLessThan(abs(isotope.exactMass - Double(isotope.massNumber)), 1)
            }
            total += isotopes.count
        }
        XCTAssertEqual(total, 3171)
        XCTAssertEqual(CDKIsotopes.atomicNumber(symbol: " Og "), 118)
        XCTAssertEqual(CDKIsotopes.atomicNumber(symbol: "cu"), 29)
        XCTAssertNil(CDKIsotopes.atomicNumber(symbol: "R1"))
        XCTAssertNil(CDKIsotopes.atomicNumber(symbol: ""))
    }

    func testMissingNaturalAndExplicitMassesRemainUnavailable() throws {
        XCTAssertNil(try CDKIsotopes.averageMass(atomicNumber: 43))
        XCTAssertNil(try CDKIsotopes.majorIsotope(atomicNumber: 43))
        XCTAssertNil(try CDKIsotopes.isotope(atomicNumber: 6, massNumber: 999))
        XCTAssertNil(try CDKIsotopes.isotope(atomicNumber: 0, massNumber: 0))
        XCTAssertNil(try CDKIsotopes.averageMass(atomicNumber: Int.max))
        for input in ["[999C]", "[Tc]", "*"] {
            let molecule = try CDKSmilesParser().parseSmiles(input)
            XCTAssertThrowsError(try CDKExactMassDescriptor.calculateValidated(for: molecule))
            XCTAssertThrowsError(try CDKMolecularWeightDescriptor.calculateValidated(for: molecule))
            XCTAssertTrue(CDKExactMassDescriptor.calculate(for: molecule).isNaN)
            XCTAssertNil(CDKRuleOfFiveDescriptor.evaluate(for: molecule, xlogP: 1).isCompliant)
        }
    }

    func testIsotopologuesAndImplicitHydrogens() throws {
        let labeled = try CDKSmilesParser().parseSmiles("[13CH4]")
        let ordinary = try CDKSmilesParser().parseSmiles("C")
        XCTAssertEqual(try CDKExactMassDescriptor.calculateValidated(for: labeled), 17.034654968, accuracy: 1e-9)
        XCTAssertEqual(
            try CDKExactMassDescriptor.calculateValidated(for: labeled, includeImplicitHydrogens: false),
            13.00335484, accuracy: 1e-9)
        XCTAssertEqual(
            CDKExactMassDescriptor.calculate(for: labeled) - CDKExactMassDescriptor.calculate(for: ordinary),
            1.00335484, accuracy: 1e-9)
        XCTAssertEqual(CDKMolecularWeightDescriptor.calculate(for: labeled), 13.00335484 + 4 * 1.00794, accuracy: 1e-9)
        XCTAssertEqual(CDKMolecularFormulaDescriptor.calculate(for: labeled), "CH4")
        XCTAssertEqual(labeled.atoms.first?.isotopeMassNumber, 13)
        let heavyWater = try CDKSmilesParser().parseSmiles("[2H]O[2H]")
        XCTAssertEqual(CDKExactMassDescriptor.calculate(for: heavyWater), 20.023118176, accuracy: 1e-9)
        XCTAssertEqual(CDKMolecularWeightDescriptor.calculate(for: heavyWater), 20.027603556, accuracy: 1e-9)
    }

    func testRareElementExactMassIsNotItsAverageMass() throws {
        let xenon = try CDKSmilesParser().parseSmiles("[Xe]")
        XCTAssertEqual(CDKExactMassDescriptor.calculate(for: xenon), 131.9041535, accuracy: 1e-7)
        XCTAssertEqual(CDKMolecularWeightDescriptor.calculate(for: xenon), 131.29248, accuracy: 1e-7)
        let technetium = try CDKSmilesParser().parseSmiles("[99Tc]")
        XCTAssertEqual(CDKExactMassDescriptor.calculate(for: technetium), 98.9062547, accuracy: 1e-7)
        XCTAssertEqual(CDKMolecularWeightDescriptor.calculate(for: technetium), 98.9062547, accuracy: 1e-7)
        let cml = try CDKFileExporter.write(molecule: technetium, as: .cml)
        XCTAssertTrue(cml.contains("elementType=\"Tc\""))
        XCTAssertFalse(cml.contains("elementType=\"R\""))
    }

    func testUnknownMassCannotProduceCompleteRuleOfFiveVerdict() {
        let result = CDKRuleOfFiveDescriptor.evaluate(
            molecularWeight: .nan, hBondDonorCount: 1,
            hBondAcceptorCount: 1, xlogP: 2)
        XCTAssertEqual(result.evaluatedCriteriaCount, 3)
        XCTAssertNil(result.isCompliant)
        XCTAssertTrue(result.statusText.contains("molecular weight unavailable"))
        let nonfinite = CDKRuleOfFiveDescriptor.evaluate(
            molecularWeight: 100, hBondDonorCount: 1,
            hBondAcceptorCount: 1, xlogP: .infinity)
        XCTAssertEqual(nonfinite.evaluatedCriteriaCount, 3)
        XCTAssertNil(nonfinite.isCompliant)
    }
}
