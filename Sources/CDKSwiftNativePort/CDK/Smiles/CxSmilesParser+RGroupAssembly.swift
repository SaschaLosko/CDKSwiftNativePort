import Foundation

// Resolves definition molecules and their attachment points without changing
// lexical CX decoding or the public metadata-application order.
extension CDKCxSmilesParser {
    static func appendRGroupDefinitions(
        to molecule: inout Molecule,
        state: CDKCxSmilesState,
        parseDefinition: (String) throws -> Molecule
    ) throws {
        guard !state.rGroupDefinitions.isEmpty else { return }

        var seenNestedLabels = Set<String>()
        var nextComponentGroupID = (molecule.atoms.compactMap(\.componentGroupID).max() ?? 0) + 1

        for label in orderedRGroupLabels(from: state) {
            guard let definitions = state.rGroupDefinitions[label], !definitions.isEmpty else {
                continue
            }

            let rootAtoms = molecule.atoms.filter { rGroupDisplayLabel(for: $0) == label }
            let incomingOrders = try incomingBondOrders(for: rootAtoms, in: molecule, label: label)

            for definitionSmiles in definitions {
                var definition = try parseDefinition(definitionSmiles)
                ensureNonOverlappingDefinitions(in: &definition, seen: &seenNestedLabels)
                try ensureAttachmentPoints(in: &definition, matching: incomingOrders)
                nextComponentGroupID = assignComponentGroupIDs(
                    in: &definition,
                    nextComponentGroupID: nextComponentGroupID)

                for idx in definition.atoms.indices
                where definition.atoms[idx].rGroupMembership == nil {
                    definition.atoms[idx].rGroupMembership = label
                }

                merge(definition, into: &molecule)
            }
        }
    }

    private static func orderedRGroupLabels(from state: CDKCxSmilesState) -> [String] {
        var labels: [String] = []
        var seen = Set<String>()

        for label in state.rGroupOrder where seen.insert(label).inserted {
            labels.append(label)
        }

        for label in state.rGroupDefinitions.keys.sorted() where seen.insert(label).inserted {
            labels.append(label)
        }

        return labels
    }

    private static func incomingBondOrders(
        for atoms: [Atom],
        in molecule: Molecule,
        label: String
    ) throws -> [BondOrder]? {
        var reference: [BondOrder]? = nil

        for atom in atoms {
            let orders = molecule.bonds(forAtom: atom.id)
                .map(\.order)
                .sorted { $0.rawValue < $1.rawValue }

            if reference == nil {
                reference = orders
            } else if reference != orders {
                throw ChemError.parseFailed("R-Group \(label) has different incoming bond orders.")
            }
        }

        return reference
    }

    private static func assignComponentGroupIDs(
        in molecule: inout Molecule,
        nextComponentGroupID: Int
    ) -> Int {
        guard !molecule.atoms.isEmpty else { return nextComponentGroupID }

        let baseGroupID = max(1, nextComponentGroupID)
        var maxAssigned = baseGroupID

        for idx in molecule.atoms.indices {
            if let existing = molecule.atoms[idx].componentGroupID {
                molecule.atoms[idx].componentGroupID = baseGroupID + existing
                maxAssigned = max(maxAssigned, baseGroupID + existing)
            } else {
                molecule.atoms[idx].componentGroupID = baseGroupID
            }
        }

        return maxAssigned + 1
    }

    private static func ensureAttachmentPoints(
        in molecule: inout Molecule,
        matching bondOrders: [BondOrder]?
    ) throws {
        guard let bondOrders else { return }

        var attachmentPointIndices = molecule.atoms.indices.filter { idx in
            molecule.atoms[idx].attachmentPoint != nil
                && molecule.atoms[idx].rGroupMembership == nil
        }

        if attachmentPointIndices.isEmpty {
            guard !molecule.atoms.isEmpty else { return }

            let anchorAtomID = molecule.atoms[0].id
            var nextAtomID = molecule.atoms.map(\.id).max() ?? 0
            var nextBondID = molecule.bonds.map(\.id).max() ?? 0
            var hydrogenAdjustment = 0

            for (offset, order) in bondOrders.enumerated() {
                nextAtomID += 1
                nextBondID += 1
                hydrogenAdjustment += integerBondContribution(order)

                molecule.atoms.append(
                    Atom(
                        id: nextAtomID,
                        element: "*",
                        position: .zero,
                        explicitHydrogenCount: 0,
                        queryType: .anyAtom,
                        attachmentPoint: offset + 1)
                )
                molecule.bonds.append(
                    Bond(
                        id: nextBondID,
                        a1: anchorAtomID,
                        a2: nextAtomID,
                        order: order))
            }

            if let hCount = molecule.atoms[0].explicitHydrogenCount {
                molecule.atoms[0].explicitHydrogenCount = max(0, hCount - hydrogenAdjustment)
            }

            attachmentPointIndices = molecule.atoms.indices.filter { idx in
                molecule.atoms[idx].attachmentPoint != nil
                    && molecule.atoms[idx].rGroupMembership == nil
            }
        }

        if attachmentPointIndices.count != bondOrders.count {
            throw ChemError.parseFailed(
                "Number of R-group attachment points does not match incoming bond orders.")
        }
    }

    private static func integerBondContribution(_ order: BondOrder) -> Int {
        switch order {
        case .single, .aromatic:
            return 1
        case .double:
            return 2
        case .triple:
            return 3
        }
    }

    private static func ensureNonOverlappingDefinitions(
        in molecule: inout Molecule, seen: inout Set<String>
    ) {
        for idx in molecule.atoms.indices {
            guard let label = rGroupDisplayLabel(for: molecule.atoms[idx]),
                label.range(of: #"^R\d*$"#, options: .regularExpression) != nil
            else {
                continue
            }
            guard !seen.insert(label).inserted else { continue }

            let numericSeed = Int(label.dropFirst()) ?? 0
            var candidate = label
            var nextValue = max(0, numericSeed)

            while seen.contains(candidate) {
                nextValue += 1
                candidate = "R\(nextValue)"
            }

            seen.insert(candidate)
            renameRGroup(label, to: candidate, in: &molecule)
        }
    }

    private static func renameRGroup(
        _ oldLabel: String, to newLabel: String, in molecule: inout Molecule
    ) {
        for idx in molecule.atoms.indices {
            if rGroupDisplayLabel(for: molecule.atoms[idx]) == oldLabel {
                molecule.atoms[idx].element = newLabel
                molecule.atoms[idx].rGroupLabel = parseRGroupLabel(newLabel)
            }
            if molecule.atoms[idx].rGroupMembership == oldLabel {
                molecule.atoms[idx].rGroupMembership = newLabel
            }
        }
    }

    private static func merge(_ source: Molecule, into destination: inout Molecule) {
        var atomMap: [Int: Int] = [:]
        var bondMap: [Int: Int] = [:]
        var nextAtomID = destination.atoms.map(\.id).max() ?? 0
        var nextBondID = destination.bonds.map(\.id).max() ?? 0

        for atom in source.atoms {
            nextAtomID += 1
            atomMap[atom.id] = nextAtomID
        }
        for atom in source.atoms {
            var copy = RGroupEnumerationPlan.copyAtom(atom, id: atomMap[atom.id]!)
            copy.ligandOrderingAtomIDs = atom.ligandOrderingAtomIDs?.map { atomMap[$0] ?? $0 }
            destination.atoms.append(copy)
        }

        for bond in source.bonds {
            guard let a1 = atomMap[bond.a1], let a2 = atomMap[bond.a2] else { continue }
            nextBondID += 1
            bondMap[bond.id] = nextBondID
            destination.bonds.append(
                Bond(
                    id: nextBondID,
                    a1: a1,
                    a2: a2,
                    order: bond.order,
                    stereo: bond.stereo,
                    doubleBondStereo: bond.doubleBondStereo,
                    stereoReferenceAtomIDs: bond.stereoReferenceAtomIDs?.compactMap { atomMap[$0] },
                    queryType: bond.queryType,
                    topology: bond.topology,
                    coordinateBondReferenceAtomID: bond.coordinateBondReferenceAtomID.flatMap {
                        atomMap[$0]
                    },
                    reactingCenterStatus: bond.reactingCenterStatus))
        }

        for sgroup in source.sgroups {
            let atomIDs = sgroup.atomIDs.compactMap { atomMap[$0] }
            guard !atomIDs.isEmpty else { continue }
            let crossingBondIDs = sgroup.crossingBondIDs.compactMap { bondMap[$0] }
            destination.sgroups.append(
                MoleculeSgroup(
                    kind: sgroup.kind,
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
                    parentAtomIDs: sgroup.parentAtomIDs.compactMap { atomMap[$0] },
                    componentNumber: sgroup.componentNumber,
                    expanded: sgroup.expanded,
                    brackets: sgroup.brackets,
                    childGroupIndices: sgroup.childGroupIndices)
            )
        }
    }

    private static func rGroupDisplayLabel(for atom: Atom) -> String? {
        if atom.isUnresolvedPseudoAtom { return atom.aliasLabel ?? atom.element }
        let trimmed = atom.element.trimmingCharacters(in: .whitespacesAndNewlines)
        if trimmed.range(of: #"^R\d*$"#, options: .regularExpression) != nil {
            return trimmed
        }
        if let rGroupLabel = atom.rGroupLabel {
            return "R\(rGroupLabel)"
        }
        return nil
    }

    static func parseRGroupLabel(_ label: String) -> Int? {
        let trimmed = label.trimmingCharacters(in: .whitespacesAndNewlines)
        guard trimmed.uppercased().hasPrefix("R") else { return nil }
        let digits = trimmed.dropFirst()
        return digits.isEmpty ? 0 : Int(digits)
    }

}
