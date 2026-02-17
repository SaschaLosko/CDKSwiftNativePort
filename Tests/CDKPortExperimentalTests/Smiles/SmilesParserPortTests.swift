import XCTest
@testable import CDKPortExperimental

final class SmilesParserPortTests: XCTestCase {
    private let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    func testParsesSingleDoubleTripleBondOrders() throws {
        let single = try parse("C-C")
        XCTAssertEqual(single.atomCount, 2)
        XCTAssertEqual(single.bondCount, 1)
        XCTAssertEqual(single.bonds.first?.order, .single)

        let double = try parse("C=C")
        XCTAssertEqual(double.bondCount, 1)
        XCTAssertEqual(double.bonds.first?.order, .double)

        let triple = try parse("C#C")
        XCTAssertEqual(triple.bondCount, 1)
        XCTAssertEqual(triple.bonds.first?.order, .triple)
    }

    // Mirrors CDK SmilesParser behavior where plain "Co" tokenizes as C + aromatic o.
    func testParsesPlainCoAsTwoAtoms() throws {
        let molecule = try parse("Co")
        XCTAssertEqual(molecule.atomCount, 2)
        XCTAssertEqual(molecule.atoms[0].element, "C")
        XCTAssertEqual(molecule.atoms[1].element, "O")
    }

    func testParsesBracketIsotopeChargeAndExplicitHydrogen() throws {
        let isotope = try parse("[13C]")
        XCTAssertEqual(isotope.atomCount, 1)
        XCTAssertEqual(isotope.atoms.first?.element, "C")
        XCTAssertEqual(isotope.atoms.first?.isotopeMassNumber, 13)

        let hydroxide = try parse("[OH-]")
        XCTAssertEqual(hydroxide.atomCount, 1)
        XCTAssertEqual(hydroxide.bondCount, 0)
        XCTAssertEqual(hydroxide.atoms.first?.element, "O")
        XCTAssertEqual(hydroxide.atoms.first?.charge, -1)
        if let oxygen = hydroxide.atoms.first {
            XCTAssertEqual(hydroxide.implicitHydrogenCount(for: oxygen.id), 1)
        }
    }

    func testParsesDisconnectedIons() throws {
        let molecule = try parse("[Na+].[OH-]")

        XCTAssertEqual(molecule.atoms.filter { $0.element == "Na" }.count, 1)
        XCTAssertEqual(molecule.atoms.filter { $0.element == "O" }.count, 1)
        XCTAssertEqual(molecule.atoms.filter { $0.element == "H" }.count, 0)
        XCTAssertEqual(molecule.atomCount, 2)
        XCTAssertEqual(molecule.bondCount, 0)

        let naID = try XCTUnwrap(molecule.atoms.first(where: { $0.element == "Na" })?.id)
        let oxygenID = try XCTUnwrap(molecule.atoms.first(where: { $0.element == "O" })?.id)
        XCTAssertNil(molecule.bond(between: naID, and: oxygenID))
    }

    func testParsesWildcardAtoms() throws {
        let plain = try parse("*C")
        XCTAssertEqual(plain.atomCount, 2)
        XCTAssertEqual(plain.atoms.first?.element, "*")
        XCTAssertEqual(plain.atoms.last?.element, "C")

        let bracket = try parse("[*]C")
        XCTAssertEqual(bracket.atomCount, 2)
        XCTAssertEqual(bracket.atoms.first?.element, "*")
        XCTAssertEqual(bracket.atoms.last?.element, "C")
    }

    func testParsesPercentRingIndices() throws {
        // CDK regression area: ring closures above single digits.
        let molecule = try parse("C%10CCCCC%10")
        XCTAssertEqual(molecule.atomCount, 6)
        XCTAssertEqual(molecule.bondCount, 6)

        let leadingZero = try parse("C%02CCCCC%02")
        XCTAssertEqual(leadingZero.atomCount, 6)
        XCTAssertEqual(leadingZero.bondCount, 6)
    }

    func testParsesAromaticSeleniumBracketAtom() throws {
        let molecule = try parse("c1ncc[se]1")
        XCTAssertGreaterThanOrEqual(molecule.atomCount, 5)
        XCTAssertTrue(molecule.atoms.contains { $0.element == "Se" && $0.aromatic })
    }

    func testParsesChiralBracketAnnotationAndAssignsWedgeHash() throws {
        let molecule = try parse("C[C@H](N)O")
        XCTAssertTrue(molecule.atoms.contains { $0.chirality != .none })
        XCTAssertTrue(
            molecule.bonds.contains { bond in
                switch bond.stereo {
                case .up, .down, .upReversed, .downReversed:
                    return true
                case .none, .either:
                    return false
                }
            }
        )
    }

    func testDirectionalBondsMarkDoubleBondIsomericHint() throws {
        let molecule = try parse("F/C=C/F")
        let doubleBond = try XCTUnwrap(molecule.bonds.first(where: { $0.order == .double }))
        XCTAssertEqual(doubleBond.stereo, .either)
    }

    func testParsesHydrogenIsotopeAliasesDAndT() throws {
        let molecule = try parse("D.T")
        XCTAssertEqual(molecule.atomCount, 2)
        XCTAssertEqual(molecule.bondCount, 0)

        let masses = molecule.atoms.map(\.isotopeMassNumber).sorted { ($0 ?? 0) < ($1 ?? 0) }
        XCTAssertEqual(masses, [2, 3])
        XCTAssertTrue(molecule.atoms.allSatisfy { $0.element == "H" })
    }

    // CDK regression area: ring closure across disconnected notation.
    func testParsesRingClosureAcrossDotNotation() throws {
        let molecule = try parse("C1.O1")
        XCTAssertEqual(molecule.atomCount, 2)
        XCTAssertEqual(molecule.bondCount, 1)
    }

    // Mirrors CDK SmilesParserTest pyrrolylpyrrole_invalid.
    func testRejectsInvalidAromaticPyrrolylpyrrole() {
        XCTAssertThrowsError(try parse("c1cccn1c2cccn2"))
    }

    // Mirrors CDK SmilesParserTest pyrrolylpyrrole_valid.
    func testAcceptsValidAromaticPyrrolylpyrrole() throws {
        let molecule = try parse("c1cccn1c2ccc[nH]2")
        XCTAssertFalse(molecule.atoms.isEmpty)
    }

    // Mirrors CDK SmilesParserTest.testPyrrole3.
    func testRejectsInvalidPyrroleWithoutDonor() {
        XCTAssertThrowsError(try parse("n1cccc1"))
    }

    // Mirrors CDK SmilesParserTest.testPyrrole1.
    func testAcceptsPyrroleWithBracketHydrogen() throws {
        let molecule = try parse("[nH]1cccc1")
        XCTAssertEqual(molecule.atomCount, 5)
        XCTAssertEqual(molecule.bondCount, 5)
    }

    // Mirrors CDK SmilesParserTest.testPyrrole2.
    func testAcceptsPyrroleWithExplicitHydrogenBranch() throws {
        let molecule = try parse("n1([H])cccc1")
        XCTAssertEqual(molecule.atomCount, 6)
        XCTAssertEqual(molecule.atoms.filter { $0.element == "H" }.count, 1)
    }

    // Mirrors CDK SmilesParserTest.testImidazole3.
    func testRejectsInvalidImidazoleWithoutDonor() {
        XCTAssertThrowsError(try parse("n1cncc1"))
    }

    // Mirrors CDK SmilesParserTest.testImidazole4.
    func testAcceptsImidazoleWithBracketHydrogen() throws {
        let molecule = try parse("n1cc[nH]c1")
        XCTAssertEqual(molecule.atomCount, 5)
    }

    func testRejectsUnterminatedRingClosure() {
        XCTAssertThrowsError(try parse("C1CC"))
    }

    func testRejectsUnbalancedBranch() {
        XCTAssertThrowsError(try parse("C)"))
    }

    private func parse(_ smiles: String) throws -> Molecule {
        try parser.parseSmiles(smiles)
    }
}
