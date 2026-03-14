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

        let reactants = grouped.filter { $0.side == .reactant }.map {
            CDKReactionParticipant(molecule: $0.molecule, role: .reactant, stoichiometry: $0.stoichiometry)
        }
        let agents = grouped.filter { $0.side == .agent }.map {
            CDKReactionParticipant(molecule: $0.molecule, role: .agent, stoichiometry: $0.stoichiometry)
        }
        let products = grouped.filter { $0.side == .product }.map {
            CDKReactionParticipant(molecule: $0.molecule, role: .product, stoichiometry: $0.stoichiometry)
        }

        return CDKReaction(reactantParticipants: reactants,
                           agentParticipants: agents,
                           productParticipants: products,
                           name: split.title,
                           cxState: split.state)
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
    let stoichiometry: Double?
}

private extension CDKSmilesParser {
    func parseReactionSide(_ raw: String, side: CDKReactionSide) throws -> [CDKIndexedReactionComponent] {
        let trimmed = raw.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { return [] }

        var out: [CDKIndexedReactionComponent] = []
        var index = trimmed.startIndex

        while index < trimmed.endIndex {
            if trimmed[index] == "." {
                throw ChemError.parseFailed("Invalid empty fragment in reaction side.")
            }

            var fragmentEnd = index
            while fragmentEnd < trimmed.endIndex, trimmed[fragmentEnd] != "." {
                fragmentEnd = trimmed.index(after: fragmentEnd)
            }

            let fragment = String(trimmed[index..<fragmentEnd]).trimmingCharacters(in: .whitespacesAndNewlines)
            guard !fragment.isEmpty else {
                throw ChemError.parseFailed("Invalid empty fragment in reaction side.")
            }

            let molecule = try parseSmiles(fragment)
            out.append(
                CDKIndexedReactionComponent(
                    molecule: molecule,
                    side: side,
                    globalIndex: -1,
                    stoichiometry: nil
                )
            )

            guard fragmentEnd < trimmed.endIndex else { break }
            index = trimmed.index(after: fragmentEnd)
            if index == trimmed.endIndex {
                throw ChemError.parseFailed("Invalid empty fragment in reaction side.")
            }
        }

        return out
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
                                                              zPosition: old.zPosition,
                                                              charge: old.charge,
                                                              isotopeMassNumber: old.isotopeMassNumber,
                                                              aromatic: false,
                                                              chirality: old.chirality,
                                                              explicitHydrogenCount: old.explicitHydrogenCount,
                                                              queryType: old.queryType,
                                                              atomList: old.atomList,
                                                              atomListIsNegated: old.atomListIsNegated,
                                                              radical: old.radical,
                                                              radicalType: old.radicalType,
                                                              atomValue: old.atomValue,
                                                              rGroupLabel: old.rGroupLabel,
                                                              rGroupMembership: old.rGroupMembership,
                                                              componentGroupID: old.componentGroupID,
                                                              substitutionCount: old.substitutionCount,
                                                              unsaturated: old.unsaturated,
                                                              ringBondCount: old.ringBondCount,
                                                              attachmentPoint: old.attachmentPoint,
                                                              valenceOverride: old.valenceOverride,
                                                              cxStereoGroup: old.cxStereoGroup,
                                                              ligandOrderingAtomIDs: old.ligandOrderingAtomIDs,
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
                                            globalIndex: $0.offset,
                                            stoichiometry: $0.element.stoichiometry)
            }
        }

        let normalized = components.enumerated().map {
            CDKIndexedReactionComponent(molecule: $0.element.molecule,
                                        side: $0.element.side,
                                        globalIndex: $0.offset,
                                        stoichiometry: $0.element.stoichiometry)
        }
        let byIndex = Dictionary(uniqueKeysWithValues: normalized.map { ($0.globalIndex, $0) })

        var groupsByIndex: [Int: [Int]] = [:]
        for rawGroup in state.fragmentGroups {
            let group = Array(Set(rawGroup)).sorted()
            guard !group.isEmpty else { continue }
            for idx in group where idx < 0 || idx >= normalized.count {
                return normalized
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
                        return normalized
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
                let groupedStoichiometry = try mergedStoichiometry(for: group, byIndex: byIndex)
                let merged = mergeDisconnectedMolecules(molecules)
                out.append(CDKIndexedReactionComponent(molecule: merged,
                                                       side: side,
                                                       globalIndex: first,
                                                       stoichiometry: groupedStoichiometry))
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
            var bondIDMap: [Int: Int] = [:]

            for atom in molecule.atoms {
                nextAtomID += 1
                idMap[atom.id] = nextAtomID
                out.atoms.append(
                    Atom(id: nextAtomID,
                         element: atom.element,
                         position: atom.position,
                         zPosition: atom.zPosition,
                         charge: atom.charge,
                         isotopeMassNumber: atom.isotopeMassNumber,
                         aromatic: atom.aromatic,
                         chirality: atom.chirality,
                         explicitHydrogenCount: atom.explicitHydrogenCount,
                         queryType: atom.queryType,
                         atomList: atom.atomList,
                         atomListIsNegated: atom.atomListIsNegated,
                         radical: atom.radical,
                         radicalType: atom.radicalType,
                         atomValue: atom.atomValue,
                         rGroupLabel: atom.rGroupLabel,
                         rGroupMembership: atom.rGroupMembership,
                         componentGroupID: atom.componentGroupID,
                         substitutionCount: atom.substitutionCount,
                         unsaturated: atom.unsaturated,
                         ringBondCount: atom.ringBondCount,
                         attachmentPoint: atom.attachmentPoint,
                         valenceOverride: atom.valenceOverride,
                         cxStereoGroup: atom.cxStereoGroup,
                         ligandOrderingAtomIDs: atom.ligandOrderingAtomIDs,
                         atomClass: atom.atomClass,
                         atomMapNumber: atom.atomMapNumber)
                )
            }

            for bond in molecule.bonds {
                guard let a1 = idMap[bond.a1], let a2 = idMap[bond.a2] else { continue }
                nextBondID += 1
                bondIDMap[bond.id] = nextBondID
                out.bonds.append(Bond(id: nextBondID,
                                      a1: a1,
                                      a2: a2,
                                      order: bond.order,
                                      stereo: bond.stereo,
                                      queryType: bond.queryType,
                                      topology: bond.topology))
            }

            for sgroup in molecule.sgroups {
                let atomIDs = sgroup.atomIDs.compactMap { idMap[$0] }
                guard !atomIDs.isEmpty else { continue }
                let crossingBondIDs = sgroup.crossingBondIDs.compactMap { bondIDMap[$0] }
                out.sgroups.append(
                    MoleculeSgroup(kind: sgroup.kind,
                                   keyword: sgroup.keyword,
                                   atomIDs: atomIDs,
                                   crossingBondIDs: crossingBondIDs,
                                   subscriptText: sgroup.subscriptText,
                                   superscriptText: sgroup.superscriptText,
                                   roundBrackets: sgroup.roundBrackets,
                                   connectivity: sgroup.connectivity,
                                   dataFieldName: sgroup.dataFieldName,
                                   dataValue: sgroup.dataValue,
                                   dataOperator: sgroup.dataOperator,
                                   dataUnit: sgroup.dataUnit,
                                   dataTag: sgroup.dataTag,
                                   subtype: sgroup.subtype,
                                   parentAtomIDs: sgroup.parentAtomIDs.compactMap { idMap[$0] },
                                   componentNumber: sgroup.componentNumber,
                                   expanded: sgroup.expanded,
                                   brackets: sgroup.brackets,
                                   childGroupIndices: sgroup.childGroupIndices)
                )
            }
        }

        return out
    }

    func mergedStoichiometry(for group: [Int],
                             byIndex: [Int: CDKIndexedReactionComponent]) throws -> Double? {
        let specified = group.compactMap { byIndex[$0]?.stoichiometry }
        guard let first = specified.first else { return nil }
        if specified.contains(where: { $0 != first }) {
            throw ChemError.parseFailed("Conflicting stoichiometry values in CXSMILES fragment grouping.")
        }
        return first
    }
}
