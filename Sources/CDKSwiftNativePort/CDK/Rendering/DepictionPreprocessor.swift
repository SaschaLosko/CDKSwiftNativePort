#if canImport(CoreGraphics)
import CoreGraphics
#endif
import Foundation

enum CDKDepictionPreprocessor {
    static func prepareForRendering(molecule: Molecule, style: RenderStyle) -> Molecule {
        var prepared = displayShortcutsApplied(to: molecule)
        if !style.showExplicitHydrogens {
            prepared = suppressibleHydrogensCollapsed(in: prepared)
        }
        return prepared
    }

    private struct ShortcutState {
        var hiddenAtomIDs: Set<Int> = []
        var hiddenBondIDs: Set<Int> = []
        var forcedVisibleAtomIDs: Set<Int> = []
        var aliasLabels: [Int: String] = [:]
        var positionOverrides: [Int: CGPoint] = [:]
    }

    private static func displayShortcutsApplied(to molecule: Molecule) -> Molecule {
        guard !molecule.sgroups.isEmpty else { return molecule }

        var prepared = molecule
        var state = ShortcutState()

        for sgroup in molecule.sgroups {
            let keyword = sgroup.keyword?.trimmingCharacters(in: .whitespacesAndNewlines).uppercased() ?? ""
            switch keyword {
            case "SUP":
                guard !sgroup.expanded,
                      let label = sgroup.subscriptText?.trimmingCharacters(in: .whitespacesAndNewlines),
                      !label.isEmpty,
                      canContractAbbreviation(sgroup: sgroup, in: molecule) else {
                    continue
                }
                contractAbbreviation(sgroup: sgroup, label: label, molecule: molecule, state: &state)
            case "MUL":
                hideMultipleGroupParts(sgroup: sgroup, molecule: molecule, state: &state)
            default:
                if sgroup.kind == .extMulticenter {
                    hideExtMulticenterAtoms(sgroup: sgroup, molecule: molecule, state: &state)
                }
            }
        }

        guard !state.hiddenAtomIDs.isEmpty
                || !state.hiddenBondIDs.isEmpty
                || !state.aliasLabels.isEmpty
                || !state.positionOverrides.isEmpty else {
            return prepared
        }

        prepared.bonds = prepared.bonds.filter { !state.hiddenBondIDs.contains($0.id) }
        let bondedAtomIDs = Set(prepared.bonds.flatMap { [$0.a1, $0.a2] })
        let keptHiddenAtomIDs = state.hiddenAtomIDs.intersection(bondedAtomIDs)
        let keepAtomIDs = bondedAtomIDs
            .union(state.forcedVisibleAtomIDs)
            .union(prepared.highlightedAtomIDs)
            .union(keptHiddenAtomIDs)

        prepared.atoms = prepared.atoms.compactMap { atom in
            let isHidden = state.hiddenAtomIDs.contains(atom.id)
            if isHidden, !keepAtomIDs.contains(atom.id) {
                return nil
            }

            var updated = atom
            if isHidden {
                updated.element = ""
                updated.aliasLabel = nil
                updated.explicitHydrogenCount = 0
                updated.atomMapNumber = nil
                updated.rGroupLabel = nil
                updated.rGroupMembership = nil
            }
            if let alias = state.aliasLabels[atom.id] {
                updated.aliasLabel = alias
            }
            if let position = state.positionOverrides[atom.id] {
                updated.position = position
            }
            return updated
        }

        return prepared
    }

    private static func canContractAbbreviation(sgroup: MoleculeSgroup, in molecule: Molecule) -> Bool {
        let internalAtomIDs = Set(sgroup.atomIDs)
        guard !internalAtomIDs.isEmpty else { return false }

        let internalBondIDs = Set(molecule.bonds.compactMap { bond in
            internalAtomIDs.contains(bond.a1) && internalAtomIDs.contains(bond.a2) ? bond.id : nil
        })
        let highlightedAtomIDs = Set(molecule.highlightedAtomIDs)
        let highlightedBondIDs = Set(molecule.highlightedBondIDs)

        let highlightedInternalAtoms = internalAtomIDs.intersection(highlightedAtomIDs)
        let highlightedInternalBonds = internalBondIDs.intersection(highlightedBondIDs)
        if highlightedInternalAtoms.isEmpty && highlightedInternalBonds.isEmpty {
            return true
        }
        return highlightedInternalAtoms.count == internalAtomIDs.count
            && highlightedInternalBonds.count == internalBondIDs.count
    }

    private static func contractAbbreviation(sgroup: MoleculeSgroup,
                                             label: String,
                                             molecule: Molecule,
                                             state: inout ShortcutState) {
        let internalAtomIDs = Set(sgroup.atomIDs)
        guard !internalAtomIDs.isEmpty else { return }

        let crossingBonds = sgroup.crossingBondIDs.compactMap { crossingBondID in
            molecule.bonds.first(where: { $0.id == crossingBondID })
        }

        if crossingBonds.count > 1 {
            var sharedInternalAtomID: Int?
            for bond in crossingBonds {
                let internalAtomID: Int?
                if internalAtomIDs.contains(bond.a1) {
                    internalAtomID = bond.a1
                } else if internalAtomIDs.contains(bond.a2) {
                    internalAtomID = bond.a2
                } else {
                    internalAtomID = nil
                }
                if let internalAtomID {
                    if let sharedInternalAtomID, sharedInternalAtomID != internalAtomID {
                        return
                    }
                    sharedInternalAtomID = internalAtomID
                }
            }
        }

        for atomID in internalAtomIDs {
            state.hiddenAtomIDs.insert(atomID)
        }
        for bond in molecule.bonds where internalAtomIDs.contains(bond.a1) || internalAtomIDs.contains(bond.a2) {
            state.hiddenBondIDs.insert(bond.id)
        }

        if crossingBonds.isEmpty {
            guard let anchorAtomID = internalAtomIDs.sorted().first else { return }
            let points = internalAtomIDs.compactMap { molecule.atom(id: $0)?.position }
            let center = centroid(of: points) ?? molecule.atom(id: anchorAtomID)?.position ?? .zero
            state.hiddenAtomIDs.remove(anchorAtomID)
            state.forcedVisibleAtomIDs.insert(anchorAtomID)
            state.aliasLabels[anchorAtomID] = label
            state.positionOverrides[anchorAtomID] = center
            return
        }

        for bond in crossingBonds {
            state.hiddenBondIDs.remove(bond.id)
            if internalAtomIDs.contains(bond.a1) {
                state.hiddenAtomIDs.remove(bond.a1)
                state.forcedVisibleAtomIDs.insert(bond.a1)
                state.aliasLabels[bond.a1] = label
            }
            if internalAtomIDs.contains(bond.a2) {
                state.hiddenAtomIDs.remove(bond.a2)
                state.forcedVisibleAtomIDs.insert(bond.a2)
                state.aliasLabels[bond.a2] = label
            }
        }
    }

    private static func hideMultipleGroupParts(sgroup: MoleculeSgroup,
                                               molecule: Molecule,
                                               state: inout ShortcutState) {
        let atomIDs = Set(sgroup.atomIDs)
        let parentAtomIDs = Set(sgroup.parentAtomIDs)
        let crossingBondIDs = Set(sgroup.crossingBondIDs)
        guard !atomIDs.isEmpty else { return }

        for bond in molecule.bonds {
            if parentAtomIDs.contains(bond.a1), parentAtomIDs.contains(bond.a2) {
                continue
            }
            if atomIDs.contains(bond.a1) || atomIDs.contains(bond.a2) {
                state.hiddenBondIDs.insert(bond.id)
            }
        }
        for atomID in atomIDs where !parentAtomIDs.contains(atomID) {
            state.hiddenAtomIDs.insert(atomID)
        }
        for bondID in crossingBondIDs {
            state.hiddenBondIDs.remove(bondID)
        }
    }

    private static func hideExtMulticenterAtoms(sgroup: MoleculeSgroup,
                                                molecule: Molecule,
                                                state: inout ShortcutState) {
        let atomIDs = Set(sgroup.atomIDs)
        for crossingBondID in sgroup.crossingBondIDs {
            guard let bond = molecule.bonds.first(where: { $0.id == crossingBondID }) else { continue }
            if atomIDs.contains(bond.a1) {
                state.hiddenAtomIDs.insert(bond.a1)
            } else if atomIDs.contains(bond.a2) {
                state.hiddenAtomIDs.insert(bond.a2)
            }
        }
    }

    private static func centroid(of points: [CGPoint]) -> CGPoint? {
        guard !points.isEmpty else { return nil }
        let count = CGFloat(points.count)
        let total = points.reduce(CGPoint.zero) { partial, point in
            CGPoint(x: partial.x + point.x, y: partial.y + point.y)
        }
        return CGPoint(x: total.x / count, y: total.y / count)
    }

    private static func suppressibleHydrogensCollapsed(in molecule: Molecule) -> Molecule {
        guard molecule.atoms.contains(where: isHydrogen) else { return molecule }

        let highlightedAtomIDs = Set(molecule.highlightedAtomIDs)
        let highlightedBondIDs = Set(molecule.highlightedBondIDs)
        var suppressedHydrogenIDs = Set<Int>()
        var removedHydrogenCountByNeighbor: [Int: Int] = [:]

        for atom in molecule.atoms where isHydrogen(atom) {
            guard !highlightedAtomIDs.contains(atom.id),
                  isSuppressibleHydrogen(atomID: atom.id,
                                         in: molecule,
                                         highlightedBondIDs: highlightedBondIDs) else {
                continue
            }
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

    private static func isSuppressibleHydrogen(atomID: Int,
                                               in molecule: Molecule,
                                               highlightedBondIDs: Set<Int>) -> Bool {
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
        guard !highlightedBondIDs.contains(bond.id) else { return false }
        guard bond.order == .single, bond.stereo == .none, bond.queryType == nil else { return false }

        let neighborID = (bond.a1 == atom.id) ? bond.a2 : bond.a1
        guard let neighbor = molecule.atom(id: neighborID), !isHydrogen(neighbor) else { return false }
        return true
    }

    private static func isHydrogen(_ atom: Atom) -> Bool {
        atom.element.trimmingCharacters(in: .whitespacesAndNewlines).uppercased() == "H"
    }
}
