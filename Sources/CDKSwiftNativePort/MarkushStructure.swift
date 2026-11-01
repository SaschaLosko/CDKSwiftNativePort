import Foundation

extension Atom {
    /// A wildcard's display label is not an element assignment, even for labels such as B or Y.
    public var isUnresolvedPseudoAtom: Bool {
        get { properties["cdk:UnresolvedPseudoAtom"] == "true" }
        set { properties["cdk:UnresolvedPseudoAtom"] = newValue ? "true" : nil }
    }
}

public enum CDKMarkushStructure {
    static let elementSymbols: Set<String> = Set(
        "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og D T"
            .split(separator: " ").map(String.init))
    public static func requiresSpecifiedStructure(_ molecule: Molecule) -> Bool {
        molecule.name.trimmingCharacters(in: .whitespacesAndNewlines).hasPrefix("|")
            || molecule.sgroups.contains {
                [.structureRepeatUnit, .polymer, .extMulticenter].contains($0.kind)
            }
            || molecule.cxState?.linkNodes.isEmpty == false
            || molecule.cxState?.positionalVariations.isEmpty == false
            || !molecule.rGroupLogicDefinitions.isEmpty
            || molecule.atoms.contains { atom in
                !elementSymbols.contains(atom.element) || atom.isUnresolvedPseudoAtom
                    || atom.queryType != nil || atom.atomList?.isEmpty == false
                    || atom.atomListIsNegated || atom.rGroupLabel != nil
                    || atom.rGroupMembership != nil
                    || atom.attachmentPoint != nil || atom.substitutionCount != nil
                    || atom.unsaturated != nil
                    || atom.ringBondCount != nil
                    || ["*", "A", "Q", "R", "R#"].contains(atom.element.uppercased())
            }
            || molecule.bonds.contains { $0.queryType != nil || $0.topology != nil }
    }
}
