import XCTest
import CoreGraphics
@testable import CDKSwiftNativePort

final class CMLWriterPortTests: XCTestCase {
    func testWritesHydrogenCountFromImplicitAndExplicitHydrogens() throws {
        let molecule = Molecule(name: "Methane",
                                atoms: [
                                    Atom(id: 1, element: "C", position: .zero, explicitHydrogenCount: 3),
                                    Atom(id: 2, element: "H", position: CGPoint(x: 1, y: 0))
                                ],
                                bonds: [
                                    Bond(id: 1, a1: 1, a2: 2, order: .single)
                                ])

        let text = try CDKCMLWriter.write([molecule])
        XCTAssertTrue(text.contains("hydrogenCount=\"4\""))
    }

    func testOmitsZeroFormalChargeAndWritesIsotopeNumber() throws {
        let molecule = Molecule(name: "Carbon-13",
                                atoms: [
                                    Atom(id: 1, element: "C", position: .zero, charge: 0, isotopeMassNumber: 13)
                                ],
                                bonds: [])

        let text = try CDKCMLWriter.write([molecule])
        XCTAssertTrue(text.contains("isotopeNumber=\"13\""))
        XCTAssertFalse(text.contains("formalCharge=\"0\""))
    }

    func testWrites2DAnd3DCoordinatesWhenZPresent() throws {
        let molecule = Molecule(name: "3D",
                                atoms: [
                                    Atom(id: 1, element: "N", position: CGPoint(x: 1.3, y: 1.4), zPosition: 0.9)
                                ],
                                bonds: [])

        let text = try CDKCMLWriter.write([molecule])
        XCTAssertTrue(text.contains("x2=\"1.30000\""))
        XCTAssertTrue(text.contains("y2=\"1.40000\""))
        XCTAssertTrue(text.contains("x3=\"1.30000\""))
        XCTAssertTrue(text.contains("y3=\"1.40000\""))
        XCTAssertTrue(text.contains("z3=\"0.90000\""))
    }

    func testWritesBondStereoAndPseudoAtomRoundTrip() throws {
        let molecule = Molecule(name: "Pseudo",
                                atoms: [
                                    Atom(id: 1, element: "R", position: CGPoint(x: 0, y: 0), zPosition: 0.1, aliasLabel: "Glu55"),
                                    Atom(id: 2, element: "O", position: CGPoint(x: 1.2, y: 0.0), zPosition: 0.3)
                                ],
                                bonds: [
                                    Bond(id: 1, a1: 1, a2: 2, order: .single, stereo: .up)
                                ])

        let text = try CDKCMLWriter.write([molecule])
        XCTAssertTrue(text.contains("elementType=\"Du\""))
        XCTAssertTrue(text.contains("title=\"Glu55\""))
        XCTAssertTrue(text.contains("<bondStereo dictRef=\"cml:W\">W</bondStereo>"))

        let roundTripped = try XCTUnwrap(CDKCMLReader.read(text: text).first)
        XCTAssertEqual(roundTripped.atoms[0].aliasLabel, "Glu55")
        XCTAssertEqual(roundTripped.bonds[0].stereo, .up)
        XCTAssertEqual(try XCTUnwrap(roundTripped.atoms[0].zPosition), 0.1, accuracy: 0.00001)
        XCTAssertEqual(try XCTUnwrap(roundTripped.atoms[1].zPosition), 0.3, accuracy: 0.00001)
    }

    func testWritesAromaticBondDecorationAndRoundTrips() throws {
        let molecule = Molecule(name: "Aromatic",
                                atoms: [
                                    Atom(id: 1, element: "C", position: CGPoint(x: 0, y: 0), aromatic: true),
                                    Atom(id: 2, element: "C", position: CGPoint(x: 1.4, y: 0), aromatic: true)
                                ],
                                bonds: [
                                    Bond(id: 1, a1: 1, a2: 2, order: .aromatic)
                                ])

        let text = try CDKCMLWriter.write([molecule])
        XCTAssertTrue(text.contains("order=\"A\""))
        XCTAssertTrue(text.contains("<bondType dictRef=\"cdk:aromaticBond\" />"))
        XCTAssertTrue(text.contains("<scalar dictRef=\"cdk:aromaticAtom\" />"))

        let roundTripped = try XCTUnwrap(CDKCMLReader.read(text: text).first)
        XCTAssertEqual(roundTripped.bonds[0].order, .aromatic)
        XCTAssertTrue(roundTripped.atoms[0].aromatic)
        XCTAssertTrue(roundTripped.atoms[1].aromatic)
    }

    func testWritesMultipleMoleculesAndRoundTrips() throws {
        let first = Molecule(name: "First",
                             atoms: [Atom(id: 1, element: "C", position: .zero)],
                             bonds: [])
        let second = Molecule(name: "Second",
                              atoms: [Atom(id: 1, element: "O", position: CGPoint(x: 1, y: 0))],
                              bonds: [])

        let text = try CDKCMLWriter.write([first, second])
        XCTAssertEqual(text.components(separatedBy: "<molecule ").count - 1, 2)

        let roundTripped = try CDKCMLReader.read(text: text)
        XCTAssertEqual(roundTripped.count, 2)
        XCTAssertEqual(roundTripped[0].name, "First")
        XCTAssertEqual(roundTripped[1].name, "Second")
    }
}
