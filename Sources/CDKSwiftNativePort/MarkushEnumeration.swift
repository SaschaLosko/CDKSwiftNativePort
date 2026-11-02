import Foundation

public struct CDKMarkushConfiguration: Sendable {
    public let molecule: Molecule
    /// Counts in the source scaffold's repeat-Sgroup order.
    public let repeatCounts: [Int]
    public let positionVariationTargets: [Int: Int]
    public let stereoInversions: [Int]
    public let assignments: [CDKRGroupAssignment]
    public let sourceByMemberAtomID: [Int: CDKMarkushAtomSource]
}

/// Bounded expansion of finite, two-ended repeat units followed by R-group substitution.
/// This does not interpret textual definitions or perform family matching.
public enum CDKMarkushEnumerator {
    public static func configurations(
        of molecule: Molecule, maximumConfigurations: Int = 1_000, maximumAtoms: Int = 10_000
    ) throws -> [CDKMarkushConfiguration] {
        var result: [CDKMarkushConfiguration] = []
        try enumerate(
            molecule, maximumConfigurations: maximumConfigurations, maximumAtoms: maximumAtoms
        ) {
            result.append($0)
            return true
        }
        return result
    }

    @discardableResult
    public static func enumerate(
        _ molecule: Molecule, maximumConfigurations: Int = 1_000, maximumAtoms: Int = 10_000,
        checkpoint: () throws -> Void = {},
        _ visit: (CDKMarkushConfiguration) throws -> Bool
    ) throws -> Int {
        guard maximumConfigurations > 0, maximumAtoms > 0 else {
            throw CDKRGroupEnumerationError.invalidQuery("Limits must be positive.")
        }
        try Task.checkCancellation()
        try checkpoint()
        let atomIDs = Set(molecule.atoms.map(\.id))
        guard atomIDs.count == molecule.atoms.count,
            Set(molecule.bonds.map(\.id)).count == molecule.bonds.count,
            molecule.bonds.allSatisfy({ $0.a1 != $0.a2 && atomIDs.contains($0.a1) && atomIDs.contains($0.a2) })
        else {
            throw CDKRGroupEnumerationError.invalidQuery("Atom and bond IDs must describe a valid graph.")
        }
        let stereo = try MarkushStereoPlan(molecule)
        if stereo.hasMetadata {
            var total = 0
            var stopped = false
            try stereo.enumerate(checkpoint: checkpoint) { graph, inversions in
                try enumerate(
                    graph, maximumConfigurations: maximumConfigurations, maximumAtoms: maximumAtoms,
                    checkpoint: checkpoint
                ) { configuration in
                    guard total < maximumConfigurations else {
                        throw CDKRGroupEnumerationError.configurationLimitExceeded(maximumConfigurations)
                    }
                    total += 1
                    stopped = try !visit(
                        .init(
                            molecule: configuration.molecule,
                            repeatCounts: configuration.repeatCounts,
                            positionVariationTargets: configuration.positionVariationTargets,
                            stereoInversions: inversions,
                            assignments: configuration.assignments,
                            sourceByMemberAtomID: configuration.sourceByMemberAtomID))
                    return !stopped
                }
                return !stopped
            }
            return total
        }
        let positions = try MarkushPositionPlan(molecule)
        if !positions.sites.isEmpty {
            var total = 0
            var stopped = false
            try positions.enumerate(checkpoint: checkpoint) { graph, targets in
                try enumerate(
                    graph, maximumConfigurations: maximumConfigurations, maximumAtoms: maximumAtoms,
                    checkpoint: checkpoint
                ) { configuration in
                    guard total < maximumConfigurations else {
                        throw CDKRGroupEnumerationError.configurationLimitExceeded(maximumConfigurations)
                    }
                    total += 1
                    stopped = try !visit(
                        .init(
                            molecule: configuration.molecule,
                            repeatCounts: configuration.repeatCounts, positionVariationTargets: targets,
                            stereoInversions: configuration.stereoInversions,
                            assignments: configuration.assignments,
                            sourceByMemberAtomID: configuration.sourceByMemberAtomID))
                    return !stopped
                }
                return !stopped
            }
            return total
        }
        let molecule = try CDKMarkushVariableLabels.resolve(molecule).molecule
        let repeats = try RepeatExpansionPlan(molecule)
        var count = 0
        var stopped = false
        func deliver(_ member: Molecule, counts: [Int], assignments: [CDKRGroupAssignment], origins: [Int: Int]) throws
            -> Bool
        {
            try Task.checkCancellation()
            try checkpoint()
            guard !CDKMarkushStructure.requiresSpecifiedStructure(member) else {
                throw CDKRGroupEnumerationError.unsupportedFeature(
                    "A generated configuration still contains unresolved variables or queries.")
            }
            guard member.atoms.count <= maximumAtoms else {
                throw CDKRGroupEnumerationError.atomLimitExceeded(maximumAtoms)
            }
            guard count < maximumConfigurations else {
                throw CDKRGroupEnumerationError.configurationLimitExceeded(maximumConfigurations)
            }
            count += 1
            var sources: [Int: CDKMarkushAtomSource] = [:]
            for atom in member.atoms {
                if let original = origins[atom.id] {
                    sources[atom.id] = .init(
                        atomID: original, siteAtomID: nil, rGroupNumber: nil, alternativeIndex: nil)
                }
            }
            for assignment in assignments {
                for (source, target) in assignment.memberAtomIDByFragmentAtomID {
                    sources[target] = .init(
                        atomID: origins[source] ?? source,
                        siteAtomID: origins[assignment.siteAtomID] ?? assignment.siteAtomID,
                        rGroupNumber: assignment.rGroupNumber, alternativeIndex: assignment.alternativeIndex)
                }
            }
            stopped = try !visit(
                CDKMarkushConfiguration(
                    molecule: member, repeatCounts: counts, positionVariationTargets: [:], stereoInversions: [],
                    assignments: assignments, sourceByMemberAtomID: sources))
            return !stopped
        }
        func expand(_ index: Int, _ graph: Molecule, _ counts: [Int], _ origins: [Int: Int]) throws {
            try Task.checkCancellation()
            try checkpoint()
            guard !stopped else { return }
            guard index < repeats.units.count else {
                var expanded = graph
                expanded.cxState = nil
                expanded.sgroups.removeAll { $0.kind == .structureRepeatUnit }
                if expanded.atoms.contains(where: {
                    $0.rGroupLabel != nil || $0.rGroupMembership != nil
                })
                    || !expanded.rGroupLogicDefinitions.isEmpty
                {
                    let query = try CDKRGroupQueryManipulator.toRGroupQuery(expanded)
                    try query.enumerateConfigurationsWithAssignments(
                        maximumConfigurations: maximumConfigurations, maximumAtoms: maximumAtoms, checkpoint: checkpoint
                    ) {
                        try deliver($0.molecule, counts: counts, assignments: $0.assignments, origins: origins)
                    }
                } else {
                    _ = try deliver(expanded, counts: counts, assignments: [], origins: origins)
                }
                return
            }
            let unit = repeats.units[index]
            for repetitions in unit.range {
                var expandedOrigins = origins
                let expanded = try unit.expand(
                    graph, count: repetitions, maximumAtoms: maximumAtoms, origins: &expandedOrigins)
                try expand(index + 1, expanded, counts + [repetitions], expandedOrigins)
                if stopped { return }
            }
        }
        try expand(0, molecule, [], Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0.id) }))
        return count
    }
}

private struct RepeatExpansionPlan {
    struct Unit {
        let atomIDs: Set<Int>
        let crossingBondIDs: [Int]
        let range: ClosedRange<Int>

        func expand(_ molecule: Molecule, count: Int, maximumAtoms: Int, origins: inout [Int: Int]) throws -> Molecule {
            if count == 0 {
                let boundaries = crossingBondIDs.compactMap { id in molecule.bonds.first { $0.id == id } }
                guard boundaries.count == 2 else {
                    throw CDKRGroupEnumerationError.invalidQuery("Missing repeat boundary.")
                }
                let endpoints = boundaries.map { atomIDs.contains($0.a1) ? $0.a2 : $0.a1 }
                guard endpoints[0] != endpoints[1],
                    !molecule.bonds.contains(where: {
                        Set([$0.a1, $0.a2]) == Set(endpoints)
                    })
                else {
                    throw CDKRGroupEnumerationError.unsupportedFeature(
                        "Omitting this repeat unit would create a self bond or duplicate bond.")
                }
                var output = molecule
                output.atoms.removeAll { atomIDs.contains($0.id) }
                output.bonds.removeAll { atomIDs.contains($0.a1) || atomIDs.contains($0.a2) }
                output.bonds.append(
                    RGroupEnumerationPlan.copyBond(
                        boundaries[0], id: boundaries[0].id, a1: endpoints[0], a2: endpoints[1],
                        mapReference: { $0 }))
                for id in atomIDs { origins.removeValue(forKey: id) }
                return output
            }
            guard molecule.atoms.count <= maximumAtoms,
                atomIDs.count <= (maximumAtoms - molecule.atoms.count) / max(1, count - 1) || count == 1
            else {
                throw CDKRGroupEnumerationError.atomLimitExceeded(maximumAtoms)
            }
            guard count > 1 else { return molecule }
            let boundaries = crossingBondIDs.compactMap { id in molecule.bonds.first { $0.id == id }
            }
            guard boundaries.count == 2 else {
                throw CDKRGroupEnumerationError.invalidQuery("Missing repeat boundary.")
            }
            func inner(_ bond: Bond) -> Int { atomIDs.contains(bond.a1) ? bond.a1 : bond.a2 }
            func outer(_ bond: Bond) -> Int { atomIDs.contains(bond.a1) ? bond.a2 : bond.a1 }
            let first = boundaries[0]
            let last = boundaries[1]
            let sourceAtoms = molecule.atoms.filter { atomIDs.contains($0.id) }
            let sourceBonds = molecule.bonds.filter {
                atomIDs.contains($0.a1) && atomIDs.contains($0.a2)
            }
            var output = molecule
            var nextAtom = molecule.atoms.map(\.id).max() ?? 0
            var nextBond = molecule.bonds.map(\.id).max() ?? 0
            var maps = [Dictionary(uniqueKeysWithValues: atomIDs.map { ($0, $0) })]
            for _ in 1..<count {
                var map: [Int: Int] = [:]
                for atom in sourceAtoms {
                    guard nextAtom < Int.max else {
                        throw CDKRGroupEnumerationError.invalidQuery("Atom IDs overflow during repetition.")
                    }
                    nextAtom += 1
                    map[atom.id] = nextAtom
                    origins[nextAtom] = origins[atom.id] ?? atom.id
                }
                maps.append(map)
            }
            for copy in 1..<count {
                try Task.checkCancellation()
                let map = maps[copy]
                for atom in sourceAtoms {
                    var cloned = RGroupEnumerationPlan.copyAtom(atom, id: map[atom.id]!)
                    cloned.ligandOrderingAtomIDs = atom.ligandOrderingAtomIDs?.map { id in
                        if id == outer(first) { return maps[copy - 1][inner(last)]! }
                        if id == outer(last) {
                            return copy + 1 < count ? maps[copy + 1][inner(first)]! : id
                        }
                        return map[id] ?? id
                    }
                    output.atoms.append(cloned)
                }
                for bond in sourceBonds {
                    guard nextBond < Int.max else {
                        throw CDKRGroupEnumerationError.invalidQuery("Bond IDs overflow during repetition.")
                    }
                    nextBond += 1
                    output.bonds.append(
                        RGroupEnumerationPlan.copyBond(
                            bond, id: nextBond,
                            a1: map[bond.a1]!, a2: map[bond.a2]!, mapReference: { map[$0] ?? $0 }))
                }
                guard nextBond < Int.max else {
                    throw CDKRGroupEnumerationError.invalidQuery("Bond IDs overflow during repetition.")
                }
                nextBond += 1
                output.bonds.append(
                    Bond(
                        id: nextBond, a1: maps[copy - 1][inner(last)]!, a2: map[inner(first)]!,
                        order: .single))
            }
            // Keep the right boundary ID stable for adjacent, independently expanded units.
            let rightIndex = output.bonds.firstIndex { $0.id == last.id }!
            let lastInner = maps[count - 1][inner(last)]!
            output.bonds[rightIndex] = RGroupEnumerationPlan.copyBond(
                last, id: last.id,
                a1: atomIDs.contains(last.a1) ? lastInner : last.a1,
                a2: atomIDs.contains(last.a2) ? lastInner : last.a2, mapReference: { $0 })
            for index in output.atoms.indices {
                if output.atoms[index].id == outer(last) {
                    output.atoms[index].ligandOrderingAtomIDs = output.atoms[index]
                        .ligandOrderingAtomIDs?.map { $0 == inner(last) ? lastInner : $0 }
                } else if output.atoms[index].id == inner(last) {
                    output.atoms[index].ligandOrderingAtomIDs = output.atoms[index]
                        .ligandOrderingAtomIDs?.map {
                            $0 == outer(last) ? maps[1][inner(first)]! : $0
                        }
                }
            }
            return output
        }
    }
    let units: [Unit]

    init(_ molecule: Molecule) throws {
        guard !molecule.name.trimmingCharacters(in: .whitespacesAndNewlines).hasPrefix("|"),
            !molecule.sgroups.contains(where: { [.polymer, .extMulticenter].contains($0.kind) }),
            molecule.cxState?.positionalVariations.isEmpty != false
        else {
            throw CDKRGroupEnumerationError.unsupportedFeature(
                "Uninterpreted CX extensions, polymers and multicenter attachments are not supported."
            )
        }
        let allIDs = Set(molecule.atoms.map(\.id))
        guard allIDs.count == molecule.atoms.count,
            Set(molecule.bonds.map(\.id)).count == molecule.bonds.count,
            molecule.bonds.allSatisfy({
                allIDs.contains($0.a1) && allIDs.contains($0.a2) && $0.a1 != $0.a2
            })
        else {
            throw CDKRGroupEnumerationError.invalidQuery(
                "Duplicate identifiers or invalid endpoints.")
        }
        for node in molecule.cxState?.linkNodes ?? [] {
            guard molecule.atoms.indices.contains(node.atomIndex) else {
                throw CDKRGroupEnumerationError.invalidQuery("Link node references a missing atom.")
            }
            if let outer = node.outerAtomIndices, !outer.isEmpty {
                let center = molecule.atoms[node.atomIndex].id
                guard outer.count == 2, Set(outer).count == 2,
                    outer.allSatisfy({ index in
                        molecule.atoms.indices.contains(index)
                            && molecule.bonds.contains {
                                ($0.a1 == center && $0.a2 == molecule.atoms[index].id)
                                    || ($0.a2 == center && $0.a1 == molecule.atoms[index].id)
                            }
                    })
                else {
                    throw CDKRGroupEnumerationError.invalidQuery("Invalid link-node outer atoms.")
                }
            }
        }
        var seen: Set<Int> = []
        var units: [Unit] = []
        for group in molecule.sgroups where group.kind == .structureRepeatUnit {
            let ids = Set(group.atomIDs)
            let bounds = (group.subscriptText ?? "").split(
                separator: "-", omittingEmptySubsequences: false)
            guard (1...2).contains(bounds.count), let lower = bounds.first.flatMap({ Int($0) }),
                let upper = bounds.last.flatMap({ Int($0) }), lower >= 0, upper >= lower,
                upper <= 1_000,
                !ids.isEmpty, ids.isSubset(of: allIDs), seen.isDisjoint(with: ids),
                group.childGroupIndices.isEmpty,
                group.connectivity == nil || group.connectivity == "ht"
            else {
                throw CDKRGroupEnumerationError.unsupportedFeature(
                    "Repeat units require disjoint atoms, head-to-tail connectivity and explicit counts from 0 through 1000."
                )
            }
            let crossing = molecule.bonds.filter { ids.contains($0.a1) != ids.contains($0.a2) }
            guard crossing.count == 2, Set(crossing.map(\.id)) == Set(group.crossingBondIDs),
                crossing.allSatisfy({
                    $0.order == .single && $0.queryType == nil && $0.stereo == .none
                }),
                molecule.atoms.filter({ ids.contains($0.id) }).allSatisfy({
                    $0.rGroupMembership == nil && !$0.aromatic
                })
            else {
                throw CDKRGroupEnumerationError.unsupportedFeature(
                    "Repeat expansion requires two explicit single-bond boundaries on the scaffold."
                )
            }
            let boundaryAtoms = Set(crossing.flatMap { [$0.a1, $0.a2] }).union(ids)
            guard
                molecule.atoms.filter({ boundaryAtoms.contains($0.id) }).allSatisfy({
                    $0.chirality == .none && $0.cxStereoGroup == nil
                }),
                molecule.bonds.allSatisfy({ bond in
                    !Set(bond.stereoReferenceAtomIDs ?? []).isDisjoint(with: ids)
                        ? false
                        : (!ids.contains(bond.a1) && !ids.contains(bond.a2))
                            || (bond.doubleBondStereo == nil && bond.stereo == .none)
                })
            else {
                throw CDKRGroupEnumerationError.unsupportedFeature(
                    "Stereogenic repeat units or boundaries require separate stereo expansion.")
            }
            var reached: Set<Int> = [ids.first!]
            var pending = Array(reached)
            while let id = pending.popLast() {
                for bond in molecule.bonds(forAtom: id) {
                    let neighbor = bond.a1 == id ? bond.a2 : bond.a1
                    if ids.contains(neighbor), reached.insert(neighbor).inserted {
                        pending.append(neighbor)
                    }
                }
            }
            guard reached == ids else {
                throw CDKRGroupEnumerationError.invalidQuery("Disconnected repeat unit.")
            }
            seen.formUnion(ids)
            units.append(
                Unit(atomIDs: ids, crossingBondIDs: group.crossingBondIDs, range: lower...upper))
        }
        guard units.count <= 16 else {
            throw CDKRGroupEnumerationError.unsupportedFeature("More than 16 repeat units.")
        }
        self.units = units
    }
}
