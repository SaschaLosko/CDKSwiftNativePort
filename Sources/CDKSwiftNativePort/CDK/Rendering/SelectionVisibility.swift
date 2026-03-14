import Foundation

// Port of CDK's SelectionVisibility.disconnected() behavior for depiction.
// Highlighted atoms are treated as selected and their labels are forced visible
// only when they are not adjacent to any highlighted bonds.
enum CDKSelectionVisibility {
    static func isSelected(atom: Atom,
                           highlightedAtomIDs: Set<Int>) -> Bool {
        highlightedAtomIDs.contains(atom.id)
    }

    static func isSelected(bond: Bond,
                           highlightedBondIDs: Set<Int>) -> Bool {
        highlightedBondIDs.contains(bond.id)
    }

    static func shouldForceVisible(atom: Atom,
                                   in molecule: Molecule,
                                   highlightedAtomIDs: Set<Int>,
                                   highlightedBondIDs: Set<Int>) -> Bool {
        guard isSelected(atom: atom, highlightedAtomIDs: highlightedAtomIDs) else {
            return false
        }
        return !molecule.bonds(forAtom: atom.id).contains { bond in
            isSelected(bond: bond, highlightedBondIDs: highlightedBondIDs)
        }
    }
}
