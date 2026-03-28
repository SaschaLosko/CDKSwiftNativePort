import Foundation
#if canImport(CoreGraphics)
import CoreGraphics
#endif
import XCTest
@testable import CDKSwiftNativePort

final class SelectionVisibilityTests: XCTestCase {
    func testHighlightSelectionDoesNotForceVisibilityWithoutSelectedAtoms() {
        let molecule = methaneLikeMolecule()

        XCTAssertFalse(
            CDKSelectionVisibility.shouldForceVisible(atom: molecule.atoms[0],
                                                      in: molecule,
                                                      highlightedAtomIDs: [],
                                                      highlightedBondIDs: [])
        )
    }

    func testHighlightSelectionRecognizesSelectedAtom() {
        let molecule = methaneLikeMolecule()
        let selectedAtoms: Set<Int> = [molecule.atoms[0].id]

        XCTAssertTrue(
            CDKSelectionVisibility.isSelected(atom: molecule.atoms[0],
                                              highlightedAtomIDs: selectedAtoms)
        )
    }

    func testSelectedAtomWithoutSelectedNeighborBondIsForcedVisible() {
        let molecule = methaneLikeMolecule()
        let selectedAtoms: Set<Int> = [molecule.atoms[0].id]

        XCTAssertTrue(
            CDKSelectionVisibility.shouldForceVisible(atom: molecule.atoms[0],
                                                      in: molecule,
                                                      highlightedAtomIDs: selectedAtoms,
                                                      highlightedBondIDs: [])
        )
    }

    func testSelectedAtomWithSelectedNeighborBondIsNotForcedVisible() {
        let molecule = methaneLikeMolecule()
        let selectedAtoms: Set<Int> = [molecule.atoms[0].id]
        let selectedBonds: Set<Int> = [molecule.bonds[0].id]

        XCTAssertFalse(
            CDKSelectionVisibility.shouldForceVisible(atom: molecule.atoms[0],
                                                      in: molecule,
                                                      highlightedAtomIDs: selectedAtoms,
                                                      highlightedBondIDs: selectedBonds)
        )
    }

    private func methaneLikeMolecule() -> Molecule {
        Molecule(
            name: "MethaneLike",
            atoms: [
                Atom(id: 1, element: "C", position: CGPoint(x: 0, y: 0)),
                Atom(id: 2, element: "H", position: CGPoint(x: 0, y: 1)),
                Atom(id: 3, element: "H", position: CGPoint(x: 0, y: -1)),
                Atom(id: 4, element: "H", position: CGPoint(x: 1, y: 0)),
                Atom(id: 5, element: "H", position: CGPoint(x: -1, y: 0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single),
                Bond(id: 2, a1: 1, a2: 3, order: .single),
                Bond(id: 3, a1: 1, a2: 4, order: .single),
                Bond(id: 4, a1: 1, a2: 5, order: .single)
            ]
        )
    }
}
