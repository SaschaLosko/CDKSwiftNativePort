import XCTest
@testable import CDKSwiftNativePort

final class InChIToStructureTests: XCTestCase {

    // Mirrors CDK InChI-to-structure coverage at a Swift/XCTest level.
    func testParsesEthanolConnectivityAndFormula() throws {
        let inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"

        let parser = CDKInChIGeneratorFactory.shared.getInChIToStructure(inchi)
        XCTAssertNotEqual(parser.getStatus(), .error)

        let molecule = try parser.getAtomContainer()
        XCTAssertEqual(molecule.atomCount, 3)
        XCTAssertEqual(molecule.bondCount, 2)

        let elements = molecule.atoms.map { $0.element }.sorted()
        XCTAssertEqual(elements, ["C", "C", "O"])
    }

    func testInfersConjugatedPatternForBenzeneInchi() throws {
        let inchi = "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"

        let parser = CDKInChIGeneratorFactory.shared.getInChIToStructure(inchi)
        let molecule = try parser.getAtomContainer()

        XCTAssertEqual(molecule.atomCount, 6)
        XCTAssertEqual(molecule.bondCount, 6)

        let doubleLike = molecule.bonds.filter { $0.order == .double || $0.order == .aromatic }.count
        XCTAssertGreaterThanOrEqual(doubleLike, 3)
    }

    func testAppliesTetrahedralLayerHints() throws {
        let inchi = "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1"

        let parser = CDKInChIGeneratorFactory.shared.getInChIToStructure(inchi)
        let molecule = try parser.getAtomContainer()

        let hasChiralCenter = molecule.atoms.contains { $0.chirality != .none }
        XCTAssertTrue(hasChiralCenter)
    }

    func testRejectsNonInchiInput() {
        let parser = CDKInChIGeneratorFactory.shared.getInChIToStructure("C1=CC=CC=C1")
        XCTAssertEqual(parser.getStatus(), .error)
        XCTAssertThrowsError(try parser.getAtomContainer())
    }

    func testParsesChargeLayerQ() throws {
        let inchi = "InChI=1S/NH3/h1H3/q+1"
        let parser = CDKInChIGeneratorFactory.shared.getInChIToStructure(inchi)
        let molecule = try parser.getAtomContainer()

        XCTAssertEqual(molecule.atomCount, 1)
        XCTAssertEqual(molecule.atoms.first?.charge, 1)
    }

    func testParsesIsotopeLayerI() throws {
        let inchi = "InChI=1S/CH4/h1H4/i1+1"
        let parser = CDKInChIGeneratorFactory.shared.getInChIToStructure(inchi)
        let molecule = try parser.getAtomContainer()

        XCTAssertEqual(molecule.atomCount, 1)
        XCTAssertEqual(molecule.atoms.first?.element.uppercased(), "C")
        XCTAssertEqual(molecule.atoms.first?.isotopeMassNumber, 13)
    }

    func testParsesMobileHydrogenLayer() throws {
        let inchi = "InChI=1S/C2H5NO2/c3-1-2(4)5/h1,3H2,(H,4,5)"
        let parser = CDKInChIGeneratorFactory.shared.getInChIToStructure(inchi)
        let molecule = try parser.getAtomContainer()

        let oxygenAtoms = molecule.atoms.filter { $0.element.uppercased() == "O" }
        XCTAssertEqual(oxygenAtoms.count, 2)

        let implicitHCounts = oxygenAtoms.map { molecule.implicitHydrogenCount(for: $0.id) }
        XCTAssertTrue(implicitHCounts.contains(where: { $0 > 0 }))
        XCTAssertTrue(implicitHCounts.contains(where: { $0 == 0 }))
    }

    func testParsesProtonLayerP() throws {
        let inchi = "InChI=1S/H2O/h1H2/p+1"
        let parser = CDKInChIGeneratorFactory.shared.getInChIToStructure(inchi)
        let molecule = try parser.getAtomContainer()

        XCTAssertEqual(molecule.atomCount, 1)
        XCTAssertEqual(molecule.atoms.first?.element.uppercased(), "O")
        XCTAssertEqual(molecule.atoms.first?.charge, 1)
    }

    func testParsesDoubleBondStereoLayerB() throws {
        let inchi = "InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+"
        let parser = CDKInChIGeneratorFactory.shared.getInChIToStructure(inchi)
        let molecule = try parser.getAtomContainer()

        let hasStereoDoubleBond = molecule.bonds.contains { $0.order == .double && $0.stereo == .either }
        XCTAssertTrue(hasStereoDoubleBond)
    }
}
