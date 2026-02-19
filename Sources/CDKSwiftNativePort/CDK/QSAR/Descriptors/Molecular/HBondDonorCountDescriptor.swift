import Foundation

public enum CDKHBondDonorCountDescriptor {
    public static func calculate(for molecule: Molecule) -> Int {
        molecule.atoms.reduce(0) { partial, atom in
            let element = CDKDescriptorSupport.canonicalElementSymbol(atom.element).uppercased()
            guard element == "N" || element == "O" || element == "S" || element == "P" else {
                return partial
            }

            let hydrogenCount = CDKDescriptorSupport.totalHydrogenCount(on: atom.id, in: molecule)
            return partial + (hydrogenCount > 0 ? 1 : 0)
        }
    }
}
