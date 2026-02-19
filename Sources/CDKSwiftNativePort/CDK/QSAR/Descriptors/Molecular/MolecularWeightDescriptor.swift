import Foundation

public enum CDKMolecularWeightDescriptor {
    public static func calculate(for molecule: Molecule,
                                 includeImplicitHydrogens: Bool = true) -> Double {
        var totalMass = 0.0
        for atom in molecule.atoms {
            totalMass += CDKDescriptorSupport.averageAtomicMass(for: atom)
            if includeImplicitHydrogens &&
                CDKDescriptorSupport.canonicalElementSymbol(atom.element).uppercased() != "H" {
                totalMass += Double(CDKDescriptorSupport.implicitHydrogenCount(on: atom.id, in: molecule)) *
                    CDKDescriptorSupport.averageAtomicMass(forElementSymbol: "H")
            }
        }
        return totalMass
    }
}
