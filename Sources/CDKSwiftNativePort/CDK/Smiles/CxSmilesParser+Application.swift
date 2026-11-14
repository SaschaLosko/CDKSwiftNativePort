import Foundation

// Applies decoded metadata after the core SMILES graph has been parsed.
extension CDKCxSmilesParser {
    public static func apply(
        to molecule: inout Molecule,
        state: CDKCxSmilesState,
        parseDefinition: (String) throws -> Molecule
    ) throws {
        applyAtomLabels(to: &molecule, state: state)
        applyAtomValues(to: &molecule, state: state)
        applyCoordinates(to: &molecule, state: state)
        applyRadicals(to: &molecule, state: state)
        applyStereoMetadata(to: &molecule, state: state)
        applyLigandOrdering(to: &molecule, state: state)
        applyHighlights(to: &molecule, state: state)
        applyCoordinateBonds(to: &molecule, state: state)
        applySgroups(to: &molecule, state: state)
        applyLinkNodes(to: &molecule, state: state)
        try appendRGroupDefinitions(to: &molecule, state: state, parseDefinition: parseDefinition)
        molecule.rGroupLogicDefinitions = state.rGroupLogic ?? [:]
        molecule.cxState = state
    }

    public static func applyAtomLabels(to molecule: inout Molecule, state: CDKCxSmilesState) {
        guard !state.atomLabels.isEmpty else { return }
        guard !molecule.atoms.isEmpty else { return }

        for (atomIndex, rawLabel) in state.atomLabels {
            guard atomIndex >= 0, atomIndex < molecule.atoms.count else { continue }

            var existing = molecule.atoms[atomIndex]
            var label = rawLabel

            if label.hasSuffix("_p") {
                label.removeLast(2)
            }

            if label.hasPrefix("_AP") {
                existing.attachmentPoint = Int(label.dropFirst(3))
                existing.explicitHydrogenCount = 0
            } else {
                // Keep wildcard provenance even for aliases that resemble elements (B, Y).
                if existing.element == "*" { existing.isUnresolvedPseudoAtom = true }
                existing.element = label
                existing.aromatic = false
                if let rLabel = parseRGroupLabel(label) {
                    existing.rGroupLabel = rLabel
                }
            }

            molecule.atoms[atomIndex] = existing
        }
    }

    public static func applyAtomValues(to molecule: inout Molecule, state: CDKCxSmilesState) {
        guard !state.atomValues.isEmpty else { return }
        for (atomIndex, value) in state.atomValues where molecule.atoms.indices.contains(atomIndex) {
            molecule.atoms[atomIndex].atomValue = value
        }
    }

    public static func applyCoordinates(to molecule: inout Molecule, state: CDKCxSmilesState) {
        guard !state.atomCoordinates.isEmpty else { return }
        let limit = min(molecule.atoms.count, state.atomCoordinates.count)
        guard limit > 0 else { return }

        for index in 0..<limit {
            let coord = state.atomCoordinates[index]
            molecule.atoms[index].position = CGPoint(x: coord.x, y: coord.y)
            molecule.atoms[index].zPosition = state.has3DCoordinates ? coord.z : nil
        }
    }

    public static func applyRadicals(to molecule: inout Molecule, state: CDKCxSmilesState) {
        guard !state.atomRadicals.isEmpty else { return }
        for (atomIndex, radicalType) in state.atomRadicals
        where molecule.atoms.indices.contains(atomIndex) {
            molecule.atoms[atomIndex].radicalType = radicalType
            molecule.atoms[atomIndex].radical = radicalType.electronCount
        }
    }

    public static func applyStereoMetadata(to molecule: inout Molecule, state: CDKCxSmilesState) {
        guard state.racemic || !state.stereoGroups.isEmpty else { return }
        for idx in molecule.atoms.indices where molecule.atoms[idx].chirality != .none {
            if state.racemic {
                molecule.atoms[idx].cxStereoGroup = encodeStereoGroup(kind: "and", number: 1)
            }
        }
        for (atomIndex, group) in state.stereoGroups
        where molecule.atoms.indices.contains(atomIndex) {
            molecule.atoms[atomIndex].cxStereoGroup = group
        }
    }

    public static func applyLigandOrdering(to molecule: inout Molecule, state: CDKCxSmilesState) {
        guard !state.ligandOrdering.isEmpty else { return }
        for (atomIndex, neighborIndices) in state.ligandOrdering
        where molecule.atoms.indices.contains(atomIndex) {
            let neighborIDs = neighborIndices.compactMap { neighborIndex in
                molecule.atoms.indices.contains(neighborIndex)
                    ? molecule.atoms[neighborIndex].id : nil
            }
            molecule.atoms[atomIndex].ligandOrderingAtomIDs = neighborIDs
        }
    }

    public static func applyHighlights(to molecule: inout Molecule, state: CDKCxSmilesState) {
        if !state.atomHighlights.isEmpty {
            molecule.highlightedAtomIDs = state.atomHighlights.compactMap { atomIndex in
                molecule.atoms.indices.contains(atomIndex) ? molecule.atoms[atomIndex].id : nil
            }
        }
        if !state.bondHighlights.isEmpty {
            molecule.highlightedBondIDs = state.bondHighlights.compactMap { bondIndex in
                molecule.bonds.indices.contains(bondIndex) ? molecule.bonds[bondIndex].id : nil
            }
        }
    }

    public static func applyCoordinateBonds(to molecule: inout Molecule, state: CDKCxSmilesState) {
        guard let coordinateBonds = state.coordinateBonds, !coordinateBonds.isEmpty else { return }

        for (atomIndex, bondIndices) in coordinateBonds {
            guard molecule.atoms.indices.contains(atomIndex) else { continue }
            let atomID = molecule.atoms[atomIndex].id
            for bondIndex in bondIndices where molecule.bonds.indices.contains(bondIndex) {
                let bond = molecule.bonds[bondIndex]
                guard bond.a1 == atomID || bond.a2 == atomID else { continue }
                molecule.bonds[bondIndex].coordinateBondReferenceAtomID = atomID
            }
        }
    }

    public static func applySgroups(to molecule: inout Molecule, state: CDKCxSmilesState) {
        guard !state.sgroups.isEmpty || !state.positionalVariations.isEmpty else { return }

        var stateToMoleculeIndex: [Int: Int] = [:]
        var pendingChildIndices: [(stateIndex: Int, moleculeIndex: Int, childIndices: [Int])] = []

        for (index, definition) in state.sgroups.enumerated() {
            switch definition.kind {
            case .polymer:
                if let moleculeSgroup = buildPolymerSgroup(
                    from: definition, molecule: molecule, state: state)
                {
                    let moleculeIndex = molecule.sgroups.count
                    molecule.sgroups.append(moleculeSgroup)
                    stateToMoleculeIndex[index] = moleculeIndex
                    pendingChildIndices.append((index, moleculeIndex, definition.childIndices))
                }
            case .data:
                if definition.dataFieldName?.hasPrefix("cdk:") == true {
                    let fieldName = definition.dataFieldName ?? ""
                    let value = definition.dataValue ?? ""
                    molecule.dataFields[fieldName, default: []].append(value)
                    if !molecule.dataFieldOrder.contains(fieldName) {
                        molecule.dataFieldOrder.append(fieldName)
                    }
                }

                let atomIDs = definition.atomIndices.compactMap { atomIndex in
                    molecule.atoms.indices.contains(atomIndex) ? molecule.atoms[atomIndex].id : nil
                }
                let moleculeIndex = molecule.sgroups.count
                molecule.sgroups.append(
                    MoleculeSgroup(
                        kind: .data,
                        atomIDs: atomIDs,
                        dataFieldName: definition.dataFieldName,
                        dataValue: definition.dataValue,
                        dataOperator: definition.dataOperator,
                        dataUnit: definition.dataUnit,
                        dataTag: definition.dataTag)
                )
                stateToMoleculeIndex[index] = moleculeIndex
                pendingChildIndices.append((index, moleculeIndex, definition.childIndices))
            }
        }

        for pending in pendingChildIndices
        where molecule.sgroups.indices.contains(pending.moleculeIndex) {
            molecule.sgroups[pending.moleculeIndex].childGroupIndices = pending.childIndices
                .compactMap {
                    childIndex in
                    stateToMoleculeIndex[childIndex]
                }
        }

        for (atomIndex, endpoints) in state.positionalVariations
        where molecule.atoms.indices.contains(atomIndex) {
            let beginAtomID = molecule.atoms[atomIndex].id
            let endpointIDs = endpoints.compactMap { endpointIndex in
                molecule.atoms.indices.contains(endpointIndex)
                    ? molecule.atoms[endpointIndex].id : nil
            }
            let firstBondID = molecule.bonds(forAtom: beginAtomID).first?.id
            molecule.sgroups.append(
                MoleculeSgroup(
                    kind: .extMulticenter,
                    keyword: "m",
                    atomIDs: [beginAtomID] + endpointIDs,
                    crossingBondIDs: firstBondID.map { [$0] } ?? [])
            )
        }
    }

    private static func buildPolymerSgroup(
        from definition: CDKCxSmilesState.SgroupDefinition,
        molecule: Molecule,
        state: CDKCxSmilesState
    ) -> MoleculeSgroup? {
        var groupAtomIDs = Set(
            definition.atomIndices.compactMap { atomIndex in
                molecule.atoms.indices.contains(atomIndex) ? molecule.atoms[atomIndex].id : nil
            })
        guard !groupAtomIDs.isEmpty else { return nil }

        var crossingBondIDs = Set(
            definition.bondIndices.compactMap { bondIndex in
                molecule.bonds.indices.contains(bondIndex) ? molecule.bonds[bondIndex].id : nil
            })

        if groupAtomIDs.count == 1 && crossingBondIDs.isEmpty,
            let atomID = groupAtomIDs.first
        {
            let ringBondIDs = molecule.bonds(forAtom: atomID)
                .filter { isRingBond($0, in: molecule) }
                .map(\.id)
            if ringBondIDs.count == 2 {
                crossingBondIDs = Set(ringBondIDs)
            }
        }

        if crossingBondIDs.isEmpty {
            for atomID in Array(groupAtomIDs) {
                for bond in molecule.bonds(forAtom: atomID) {
                    let otherAtomID = bond.a1 == atomID ? bond.a2 : bond.a1
                    guard !groupAtomIDs.contains(otherAtomID) else { continue }

                    var isCrossing = true
                    if let otherIndex = molecule.indexOfAtom(id: otherAtomID),
                        let endpoints = state.positionalVariations[otherIndex]
                    {
                        let endpointIDs = Set(
                            endpoints.compactMap { endpointIndex in
                                molecule.atoms.indices.contains(endpointIndex)
                                    ? molecule.atoms[endpointIndex].id : nil
                            })
                        if !groupAtomIDs.isDisjoint(with: endpointIDs) {
                            isCrossing = false
                        }
                    }

                    if isCrossing {
                        crossingBondIDs.insert(bond.id)
                    }
                }
            }
        } else {
            var queue = Array(groupAtomIDs)
            var seen = groupAtomIDs
            while let atomID = queue.first {
                queue.removeFirst()
                for bond in molecule.bonds(forAtom: atomID) where !crossingBondIDs.contains(bond.id) {
                    let otherAtomID = bond.a1 == atomID ? bond.a2 : bond.a1
                    if seen.insert(otherAtomID).inserted {
                        queue.append(otherAtomID)
                        groupAtomIDs.insert(otherAtomID)
                    }
                }
            }
        }

        let kind = moleculeSgroupKind(for: definition.keyword)
        return MoleculeSgroup(
            kind: kind,
            keyword: definition.keyword,
            atomIDs: groupAtomIDs.sorted(),
            crossingBondIDs: crossingBondIDs.sorted(),
            subscriptText: definition.subscriptText,
            superscriptText: definition.superscriptText,
            roundBrackets: kind == .structureRepeatUnit,
            connectivity: definition.superscriptText)
    }

    private static func moleculeSgroupKind(for keyword: String?) -> MoleculeSgroup.Kind {
        switch keyword?.lowercased() {
        case "n":
            return .structureRepeatUnit
        case "mon", "mer", "co", "xl", "mod", "mix", "f", "any", "gen", "c", "grf", "alt", "ran",
            "blk":
            return .polymer
        default:
            return .generic
        }
    }

    public static func applyLinkNodes(to molecule: inout Molecule, state: CDKCxSmilesState) {
        guard !state.linkNodes.isEmpty else { return }

        for linkNode in state.linkNodes {
            guard molecule.atoms.indices.contains(linkNode.atomIndex) else { continue }

            var groupAtomIDs = Set([molecule.atoms[linkNode.atomIndex].id])
            var crossingBondIDs = Set(
                linkNode.bondIndices.compactMap { rawIndex in
                    molecule.bonds.indices.contains(rawIndex) ? molecule.bonds[rawIndex].id : nil
                })

            if let outer = linkNode.outerAtomIndices, !outer.isEmpty {
                let center = molecule.atoms[linkNode.atomIndex].id
                crossingBondIDs = Set(
                    outer.compactMap { index in
                        guard molecule.atoms.indices.contains(index) else { return nil }
                        let id = molecule.atoms[index].id
                        return molecule.bonds.first {
                            ($0.a1 == center && $0.a2 == id) || ($0.a2 == center && $0.a1 == id)
                        }?.id
                    })
            }

            if groupAtomIDs.count == 1 && crossingBondIDs.isEmpty {
                let atomID = molecule.atoms[linkNode.atomIndex].id
                let ringBondIDs = molecule.bonds(forAtom: atomID)
                    .filter { isRingBond($0, in: molecule) }
                    .map(\.id)
                if ringBondIDs.count == 2 {
                    crossingBondIDs = Set(ringBondIDs)
                }
            }

            if crossingBondIDs.isEmpty {
                for atomID in Array(groupAtomIDs) {
                    for bond in molecule.bonds(forAtom: atomID) {
                        let otherAtomID = bond.a1 == atomID ? bond.a2 : bond.a1
                        if !groupAtomIDs.contains(otherAtomID) {
                            crossingBondIDs.insert(bond.id)
                        }
                    }
                }
            } else {
                var queue = Array(groupAtomIDs)
                var seen = groupAtomIDs

                while let atomID = queue.first {
                    queue.removeFirst()
                    for bond in molecule.bonds(forAtom: atomID)
                    where !crossingBondIDs.contains(bond.id) {
                        let otherAtomID = bond.a1 == atomID ? bond.a2 : bond.a1
                        if seen.insert(otherAtomID).inserted {
                            queue.append(otherAtomID)
                            groupAtomIDs.insert(otherAtomID)
                        }
                    }
                }
            }

            molecule.sgroups.append(
                MoleculeSgroup(
                    kind: .structureRepeatUnit,
                    keyword: "n",
                    atomIDs: groupAtomIDs.sorted(),
                    crossingBondIDs: crossingBondIDs.sorted(),
                    subscriptText: "\(linkNode.lowerBound)-\(linkNode.upperBound)",
                    superscriptText: "ht",
                    roundBrackets: true,
                    connectivity: "ht")
            )
        }
    }

    private static func isRingBond(_ bond: Bond, in molecule: Molecule) -> Bool {
        struct EdgeKey: Hashable {
            let a: Int
            let b: Int

            init(_ u: Int, _ v: Int) {
                a = min(u, v)
                b = max(u, v)
            }
        }

        let ringEdges = Set(
            molecule.simpleCycles(maxSize: 12).flatMap { ring in
                guard let first = ring.first else { return [EdgeKey]() }
                let shifted = Array(ring.dropFirst()) + [first]
                return zip(ring, shifted).map { EdgeKey($0.0, $0.1) }
            })
        return ringEdges.contains(EdgeKey(bond.a1, bond.a2))
    }

}
