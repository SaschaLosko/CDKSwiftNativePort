import Foundation

enum CDKDepictionPreprocessor {
    static func prepareForRendering(molecule: Molecule, style: RenderStyle) -> Molecule {
        guard !style.showExplicitHydrogens else { return molecule }
        return suppressibleHydrogensCollapsed(in: molecule)
    }

    private static func suppressibleHydrogensCollapsed(in molecule: Molecule) -> Molecule {
        guard molecule.atoms.contains(where: isHydrogen) else { return molecule }

        var suppressedHydrogenIDs = Set<Int>()
        var removedHydrogenCountByNeighbor: [Int: Int] = [:]

        for atom in molecule.atoms where isHydrogen(atom) {
            guard isSuppressibleHydrogen(atomID: atom.id, in: molecule) else { continue }
            suppressedHydrogenIDs.insert(atom.id)
            if let bond = molecule.bonds(forAtom: atom.id).first {
                let neighborID = (bond.a1 == atom.id) ? bond.a2 : bond.a1
                removedHydrogenCountByNeighbor[neighborID, default: 0] += 1
            }
        }

        guard !suppressedHydrogenIDs.isEmpty else { return molecule }

        var prepared = molecule
        prepared.atoms = molecule.atoms.compactMap { atom in
            guard !suppressedHydrogenIDs.contains(atom.id) else { return nil }

            var updated = atom
            if let removed = removedHydrogenCountByNeighbor[atom.id], removed > 0 {
                updated.explicitHydrogenCount = max(0, updated.explicitHydrogenCount ?? 0) + removed
            }
            return updated
        }
        prepared.bonds = molecule.bonds.filter { bond in
            !suppressedHydrogenIDs.contains(bond.a1) && !suppressedHydrogenIDs.contains(bond.a2)
        }
        return prepared
    }

    private static func isSuppressibleHydrogen(atomID: Int, in molecule: Molecule) -> Bool {
        guard let atom = molecule.atom(id: atomID), isHydrogen(atom) else { return false }
        guard atom.charge == 0,
              atom.isotopeMassNumber == nil,
              atom.radical == nil,
              atom.queryType == nil,
              atom.atomList == nil else {
            return false
        }

        let attachedBonds = molecule.bonds(forAtom: atom.id)
        guard attachedBonds.count == 1, let bond = attachedBonds.first else { return false }
        guard bond.order == .single, bond.stereo == .none, bond.queryType == nil else { return false }

        let neighborID = (bond.a1 == atom.id) ? bond.a2 : bond.a1
        guard let neighbor = molecule.atom(id: neighborID), !isHydrogen(neighbor) else { return false }
        return true
    }

    private static func isHydrogen(_ atom: Atom) -> Bool {
        atom.element.trimmingCharacters(in: .whitespacesAndNewlines).uppercased() == "H"
    }
}
