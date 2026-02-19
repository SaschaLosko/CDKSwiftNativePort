import Foundation

extension CDKSmilesParser {
    /// Swift counterpart of CDK's `parseReactionSmiles`.
    public func parseReactionSmiles(_ reactionSmiles: String) throws -> CDKReaction {
        let split = try CDKCxSmilesParser.split(reactionSmiles, enabled: flavor.contains(.cxsmiles))
        let parts = split.coreSmiles.split(separator: ">", omittingEmptySubsequences: false).map(String.init)
        guard parts.count == 3 else {
            throw ChemError.parseFailed("Reaction SMILES must contain exactly two '>' separators.")
        }

        let reactantComponents = try parseReactionSide(parts[0], side: .reactant)
        let agentComponents = try parseReactionSide(parts[1], side: .agent)
        let productComponents = try parseReactionSide(parts[2], side: .product)

        var components = reactantComponents + agentComponents + productComponents
        applyCxAtomLabels(to: &components, state: split.state)
        let grouped = try applyFragmentGrouping(to: components, state: split.state)

        let reactants = grouped.filter { $0.side == .reactant }.map(\.molecule)
        let agents = grouped.filter { $0.side == .agent }.map(\.molecule)
        let products = grouped.filter { $0.side == .product }.map(\.molecule)

        return CDKReaction(reactants: reactants, agents: agents, products: products, cxState: split.state)
    }
}

private enum CDKReactionSide {
    case reactant
    case agent
    case product
}

private struct CDKIndexedReactionComponent {
    var molecule: Molecule
    let side: CDKReactionSide
    let globalIndex: Int
}

private extension CDKSmilesParser {
    func parseReactionSide(_ raw: String, side: CDKReactionSide) throws -> [CDKIndexedReactionComponent] {
        let trimmed = raw.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { return [] }

        let fragments = trimmed.split(separator: ".", omittingEmptySubsequences: false).map(String.init)
        if fragments.contains(where: { $0.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty }) {
            throw ChemError.parseFailed("Invalid empty fragment in reaction side.")
        }

        return try fragments.map { fragment in
            CDKIndexedReactionComponent(
                molecule: try parseSmiles(fragment),
                side: side,
                globalIndex: -1
            )
        }
    }

    func applyCxAtomLabels(to components: inout [CDKIndexedReactionComponent], state: CDKCxSmilesState) {
        guard !state.atomLabels.isEmpty else { return }

        var componentRanges: [Range<Int>] = []
        var atomOffset = 0
        for component in components {
            let upper = atomOffset + component.molecule.atomCount
            componentRanges.append(atomOffset..<upper)
            atomOffset = upper
        }

        for (globalAtomIndex, label) in state.atomLabels {
            guard let compIndex = componentRanges.firstIndex(where: { $0.contains(globalAtomIndex) }) else { continue }
            let localAtomIndex = globalAtomIndex - componentRanges[compIndex].lowerBound
            guard localAtomIndex >= 0, localAtomIndex < components[compIndex].molecule.atoms.count else { continue }

            let atomID = components[compIndex].molecule.atoms[localAtomIndex].id
            guard let idx = components[compIndex].molecule.atoms.firstIndex(where: { $0.id == atomID }) else { continue }
            let old = components[compIndex].molecule.atoms[idx]
            components[compIndex].molecule.atoms[idx] = Atom(id: old.id,
                                                              element: label,
                                                              position: old.position,
                                                              charge: old.charge,
                                                              isotopeMassNumber: old.isotopeMassNumber,
                                                              aromatic: false,
                                                              chirality: old.chirality,
                                                              explicitHydrogenCount: old.explicitHydrogenCount,
                                                              queryType: old.queryType,
                                                              atomList: old.atomList,
                                                              atomListIsNegated: old.atomListIsNegated,
                                                              radical: old.radical,
                                                              rGroupLabel: old.rGroupLabel,
                                                              substitutionCount: old.substitutionCount,
                                                              unsaturated: old.unsaturated,
                                                              ringBondCount: old.ringBondCount,
                                                              attachmentPoint: old.attachmentPoint,
                                                              atomClass: old.atomClass,
                                                              atomMapNumber: old.atomMapNumber)
        }
    }

    func applyFragmentGrouping(to components: [CDKIndexedReactionComponent],
                               state: CDKCxSmilesState) throws -> [CDKIndexedReactionComponent] {
        guard !state.fragmentGroups.isEmpty else {
            return components.enumerated().map {
                CDKIndexedReactionComponent(molecule: $0.element.molecule,
                                            side: $0.element.side,
                                            globalIndex: $0.offset)
            }
        }

        let normalized = components.enumerated().map {
            CDKIndexedReactionComponent(molecule: $0.element.molecule,
                                        side: $0.element.side,
                                        globalIndex: $0.offset)
        }
        let byIndex = Dictionary(uniqueKeysWithValues: normalized.map { ($0.globalIndex, $0) })

        var groupsByIndex: [Int: [Int]] = [:]
        for rawGroup in state.fragmentGroups {
            let group = Array(Set(rawGroup)).sorted()
            guard !group.isEmpty else { continue }
            for idx in group where idx < 0 || idx >= normalized.count {
                throw ChemError.parseFailed("CXSMILES fragment-group index \(idx) is out of range.")
            }

            var perSide: [CDKReactionSide: [Int]] = [:]
            for idx in group {
                guard let side = byIndex[idx]?.side else { continue }
                perSide[side, default: []].append(idx)
            }

            for (_, sideGroupRaw) in perSide {
                let sideGroup = Array(Set(sideGroupRaw)).sorted()
                guard !sideGroup.isEmpty else { continue }
                for idx in sideGroup {
                    if groupsByIndex[idx] != nil {
                        throw ChemError.parseFailed("CXSMILES fragment index \(idx) appears in multiple groups.")
                    }
                    groupsByIndex[idx] = sideGroup
                }
            }
        }

        var emittedGroups: Set<[Int]> = []
        var out: [CDKIndexedReactionComponent] = []

        for component in normalized {
            if let group = groupsByIndex[component.globalIndex] {
                let first = group.min() ?? component.globalIndex
                guard component.globalIndex == first else { continue }
                if emittedGroups.contains(group) { continue }
                emittedGroups.insert(group)

                let molecules = group.compactMap { byIndex[$0]?.molecule }
                let side = byIndex[first]?.side ?? component.side
                let merged = mergeDisconnectedMolecules(molecules)
                out.append(CDKIndexedReactionComponent(molecule: merged, side: side, globalIndex: first))
            } else {
                out.append(component)
            }
        }

        return out
    }

    func mergeDisconnectedMolecules(_ molecules: [Molecule]) -> Molecule {
        var out = Molecule(name: "SMILES")
        var nextAtomID = 0
        var nextBondID = 0

        for molecule in molecules {
            var idMap: [Int: Int] = [:]

            for atom in molecule.atoms {
                nextAtomID += 1
                idMap[atom.id] = nextAtomID
                out.atoms.append(
                    Atom(id: nextAtomID,
                         element: atom.element,
                         position: atom.position,
                         charge: atom.charge,
                         isotopeMassNumber: atom.isotopeMassNumber,
                         aromatic: atom.aromatic,
                         chirality: atom.chirality,
                         explicitHydrogenCount: atom.explicitHydrogenCount,
                         queryType: atom.queryType,
                         atomList: atom.atomList,
                         atomListIsNegated: atom.atomListIsNegated,
                         radical: atom.radical,
                         rGroupLabel: atom.rGroupLabel,
                         substitutionCount: atom.substitutionCount,
                         unsaturated: atom.unsaturated,
                         ringBondCount: atom.ringBondCount,
                         attachmentPoint: atom.attachmentPoint,
                         atomClass: atom.atomClass,
                         atomMapNumber: atom.atomMapNumber)
                )
            }

            for bond in molecule.bonds {
                guard let a1 = idMap[bond.a1], let a2 = idMap[bond.a2] else { continue }
                nextBondID += 1
                out.bonds.append(Bond(id: nextBondID, a1: a1, a2: a2, order: bond.order, stereo: bond.stereo))
            }
        }

        return out
    }
}
