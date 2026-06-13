import XCTest
#if canImport(CoreGraphics)
import CoreGraphics
#endif
@testable import CDKSwiftNativePort

final class KekulizationTests: XCTestCase {
    func testKekulizesAromaticBenzeneRing() throws {
        let molecule = aromaticRing(elements: Array(repeating: "C", count: 6))

        let kekulized = try CDKKekulization.kekulized(molecule)

        XCTAssertEqual(kekulized.bonds.filter { $0.order == .double }.count, 3)
        XCTAssertEqual(kekulized.bonds.filter { $0.order == .single }.count, 3)
        XCTAssertFalse(kekulized.bonds.contains { $0.order == .aromatic })
        XCTAssertTrue(kekulized.atoms.allSatisfy(\.aromatic))
    }

    func testKekulizesPyrroleWithoutAssigningDoubleBondToNH() throws {
        var molecule = aromaticRing(elements: ["N", "C", "C", "C", "C"])
        molecule.atoms[0].explicitHydrogenCount = 1

        let kekulized = try CDKKekulization.kekulized(molecule)

        let nitrogenBondOrders = kekulized.bonds
            .filter { $0.a1 == 1 || $0.a2 == 1 }
            .map(\.order)
        XCTAssertEqual(nitrogenBondOrders, [.single, .single])
        XCTAssertEqual(kekulized.bonds.filter { $0.order == .double }.count, 2)
        XCTAssertEqual(kekulized.bonds.filter { $0.order == .single }.count, 3)
    }

    func testKekulizationRejectsOddUnavailableMatching() throws {
        let molecule = aromaticRing(elements: Array(repeating: "C", count: 5))

        XCTAssertThrowsError(try CDKKekulization.kekulized(molecule)) { error in
            guard case let ChemError.unsupported(message) = error else {
                return XCTFail("Expected ChemError.unsupported, got \(error).")
            }
            XCTAssertTrue(message.contains("Cannot assign Kekulé structure"))
        }
    }

    private func aromaticRing(elements: [String]) -> Molecule {
        let atoms = elements.enumerated().map { offset, element in
            Atom(
                id: offset + 1,
                element: element,
                position: CGPoint(x: Double(offset), y: 0),
                aromatic: true)
        }
        let bonds = elements.indices.map { index in
            Bond(
                id: index + 1,
                a1: index + 1,
                a2: index == elements.count - 1 ? 1 : index + 2,
                order: .aromatic)
        }
        return Molecule(atoms: atoms, bonds: bonds)
    }
}
