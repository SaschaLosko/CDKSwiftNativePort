import Foundation

/// CDK port of `org.openscience.cdk.qsar.descriptors.molecular.MannholdLogPDescriptor`.
/// The descriptor estimates LogP from atom counts:
/// `logP = 1.46 + 0.11 * (#C) - 0.11 * (#hetero)`.
public enum CDKMannholdLogPDescriptor {
    public static func calculate(for molecule: Molecule) -> Double {
        var carbonCount = 0
        var heteroAtomCount = 0

        for atom in molecule.atoms {
            let symbol = CDKDescriptorSupport.canonicalElementSymbol(atom.element).uppercased()
            switch symbol {
            case "C":
                carbonCount += 1
            case "H":
                continue
            default:
                heteroAtomCount += 1
            }
        }

        return 1.46 + 0.11 * Double(carbonCount) - 0.11 * Double(heteroAtomCount)
    }
}
