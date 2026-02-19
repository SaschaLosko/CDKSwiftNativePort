import Foundation

public enum CDKHBondAcceptorCountDescriptor {
    public static func calculate(for molecule: Molecule) -> Int {
        molecule.atoms.reduce(0) { partial, atom in
            let element = CDKDescriptorSupport.canonicalElementSymbol(atom.element).uppercased()
            let hydrogenCount = CDKDescriptorSupport.totalHydrogenCount(on: atom.id, in: molecule)

            let isAcceptor: Bool
            switch element {
            case "N":
                let aromaticDonorNitrogen = atom.aromatic && hydrogenCount > 0
                isAcceptor = atom.charge <= 0 && !aromaticDonorNitrogen && !isAmideLikeNitrogen(atomID: atom.id, in: molecule)
            case "O":
                isAcceptor = atom.charge <= 0
            case "S":
                isAcceptor = atom.charge <= 0
            case "P":
                isAcceptor = atom.charge <= 0 && hydrogenCount == 0
            default:
                isAcceptor = false
            }
            return partial + (isAcceptor ? 1 : 0)
        }
    }

    private static func isAmideLikeNitrogen(atomID: Int, in molecule: Molecule) -> Bool {
        for bond in molecule.bonds(forAtom: atomID) {
            if CDKDescriptorSupport.isAmideLikeCNBond(bond, in: molecule) {
                return true
            }
        }
        return false
    }
}
