import Foundation

public enum CDKExactMassDescriptor {
    public static func calculate(for molecule: Molecule,
                                 includeImplicitHydrogens: Bool = true) -> Double {
        var totalMass = 0.0
        let hydrogenExactMass = CDKDescriptorSupport.monoisotopicAtomicMass(forElementSymbol: "H")

        for atom in molecule.atoms {
            totalMass += CDKDescriptorSupport.monoisotopicAtomicMass(for: atom)
            if includeImplicitHydrogens &&
                CDKDescriptorSupport.canonicalElementSymbol(atom.element).uppercased() != "H" {
                totalMass += Double(CDKDescriptorSupport.implicitHydrogenCount(on: atom.id, in: molecule)) * hydrogenExactMass
            }
        }
        return totalMass
    }
}
