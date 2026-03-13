import Foundation

/// CDK port of `org.openscience.cdk.qsar.descriptors.molecular.VABCDescriptor`
/// backed by the original `org.openscience.cdk.geometry.volume.VABCVolume`
/// volume model and constants.
public enum CDKVABCDescriptor {
    public static func calculate(for molecule: Molecule) -> Double? {
        try? CDKVABCVolume.calculate(for: molecule)
    }
}

enum CDKVABCVolume {
    private enum VolumeError: Error {
        case unsupportedElement(String)
    }

    static func calculate(for molecule: Molecule) throws -> Double {
        var sum = 0.0
        var totalImplicitHydrogenCount = 0

        for atom in molecule.atoms {
            let symbol = CDKDescriptorSupport.canonicalElementSymbol(atom.element)
            guard let atomVolume = bondiiVolumes[symbol] else {
                throw VolumeError.unsupportedElement(symbol)
            }

            sum += atomVolume

            let implicitHydrogenCount = CDKDescriptorSupport.implicitHydrogenCount(on: atom.id, in: molecule)
            if implicitHydrogenCount > 0 {
                sum += Double(implicitHydrogenCount) * hydrogenVolume
                totalImplicitHydrogenCount += implicitHydrogenCount
            }
        }

        sum -= 5.92 * Double(molecule.bondCount + totalImplicitHydrogenCount)

        let ringBasis = CDKDescriptorRingSupport.smallestRingBasis(in: molecule)
        if !ringBasis.isEmpty {
            let aromaticAtomIDs = CDKDescriptorRingSupport.aromaticAtomIDs(in: molecule, perceiveAromaticity: true)
            let aromaticRingCount = ringBasis.reduce(0) { partial, cycle in
                partial + (cycle.atoms.allSatisfy { aromaticAtomIDs.contains($0) } ? 1 : 0)
            }
            let nonAromaticRingCount = ringBasis.count - aromaticRingCount

            sum -= 14.7 * Double(aromaticRingCount)
            sum -= 3.8 * Double(nonAromaticRingCount)
        }

        return sum
    }

    private static let hydrogenVolume = 7.2382293504

    // Bondi atom volumes used by CDK's VABCVolume implementation.
    private static let bondiiVolumes: [String: Double] = [
        "H": 7.2382293504,
        "C": 20.5795259250667,
        "N": 15.5985308577667,
        "O": 14.7102267005611,
        "Cl": 22.4492971208333,
        "Br": 26.5218483279667,
        "F": 13.3057882007064,
        "I": 32.5150310206656,
        "S": 24.4290240576,
        "P": 24.4290240576,
        "As": 26.5218483279667,
        "B": 40.48,
        "Se": 28.7309115245333,
        "Si": 38.7923854248
    ]
}
