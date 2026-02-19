import Foundation

public enum CDKHeavyAtomCountDescriptor {
    public static func calculate(for molecule: Molecule) -> Int {
        molecule.atoms.reduce(0) { partial, atom in
            let isHydrogen = CDKDescriptorSupport.canonicalElementSymbol(atom.element).uppercased() == "H"
            return partial + (isHydrogen ? 0 : 1)
        }
    }
}
