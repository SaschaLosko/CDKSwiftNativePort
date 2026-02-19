import Foundation

public enum CDKRingCountDescriptor {
    public static func calculate(for molecule: Molecule) -> Int {
        let heavyAtomIDs = CDKDescriptorSupport.heavyAtomIDs(in: molecule)
        guard !heavyAtomIDs.isEmpty else { return 0 }

        let heavyBondCount = molecule.bonds.reduce(0) { partial, bond in
            let includesOnlyHeavyAtoms = heavyAtomIDs.contains(bond.a1) && heavyAtomIDs.contains(bond.a2)
            return partial + (includesOnlyHeavyAtoms ? 1 : 0)
        }
        let heavyAtomCount = heavyAtomIDs.count
        let componentCount = CDKDescriptorSupport.heavyAtomConnectedComponentCount(in: molecule)
        return max(0, heavyBondCount - heavyAtomCount + componentCount)
    }
}
