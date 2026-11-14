import Foundation

public enum CDKExactMassDescriptor {
    /// Returns NaN when physical mass data is unavailable. Use calculateValidated
    /// to receive an explicit unsupported-isotope or dataset error.
    public static func calculate(
        for molecule: Molecule, includeImplicitHydrogens: Bool = true
    ) -> Double {
        (try? calculateValidated(for: molecule, includeImplicitHydrogens: includeImplicitHydrogens)) ?? .nan
    }

    public static func calculateValidated(
        for molecule: Molecule, includeImplicitHydrogens: Bool = true
    ) throws -> Double {
        var totalMass = 0.0
        let hydrogenMass = try CDKIsotopes.mass(element: "H", average: false)
        for atom in molecule.atoms {
            totalMass += try CDKIsotopes.mass(
                element: atom.element, massNumber: atom.isotopeMassNumber, average: false)
            if includeImplicitHydrogens && CDKIsotopes.atomicNumber(symbol: atom.element) != 1 {
                totalMass +=
                    Double(CDKDescriptorSupport.implicitHydrogenCount(on: atom.id, in: molecule)) * hydrogenMass
            }
        }
        return totalMass
    }
}
