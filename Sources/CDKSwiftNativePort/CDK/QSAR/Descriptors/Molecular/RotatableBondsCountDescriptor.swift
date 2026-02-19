import Foundation

public enum CDKRotatableBondsCountDescriptor {
    public static func calculate(for molecule: Molecule) -> Int {
        molecule.bonds.reduce(0) { partial, bond in
            guard bond.order == .single else { return partial }
            guard !CDKDescriptorSupport.isBondInRing(bond, in: molecule) else { return partial }
            guard !CDKDescriptorSupport.isAmideLikeCNBond(bond, in: molecule) else { return partial }

            guard let atomA = molecule.atom(id: bond.a1), let atomB = molecule.atom(id: bond.a2) else {
                return partial
            }
            let symbolA = CDKDescriptorSupport.canonicalElementSymbol(atomA.element).uppercased()
            let symbolB = CDKDescriptorSupport.canonicalElementSymbol(atomB.element).uppercased()
            guard symbolA != "H", symbolB != "H" else { return partial }

            let heavyDegreeA = CDKDescriptorSupport.heavyDegree(of: bond.a1, in: molecule)
            let heavyDegreeB = CDKDescriptorSupport.heavyDegree(of: bond.a2, in: molecule)
            guard heavyDegreeA > 1, heavyDegreeB > 1 else { return partial }

            return partial + 1
        }
    }
}
