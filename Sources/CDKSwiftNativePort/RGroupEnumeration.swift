// Derived from CDK 2.13 RGroupQuery by Mark Rijnbeek (Copyright 2010)
// and the explicit attachment-point support by John Mayfield (Copyright 2026).
// SPDX-License-Identifier: LGPL-2.1-or-later
import Foundation

public enum CDKRGroupEnumerationError: Error, Equatable, Sendable {
    case invalidQuery(String)
    case unsupportedFeature(String)
    case configurationLimitExceeded(Int)
}

extension Atom {
    /// CDK's REST_H query annotation on a neighbor of an omitted R-group.
    /// This is a matching constraint, not a hydrogen count.
    public var rGroupRestH: Bool? {
        get { properties["cdk:RestH"].flatMap(Bool.init) }
        set { properties["cdk:RestH"] = newValue.map(String.init) }
    }
}

extension CDKRGroupQuery {
    /// Enumerates positional assignments, including symmetry-equivalent duplicates, like CDK 2.13.
    /// The limit is explicit: exceeding it throws rather than returning a truncated family.
    public func getAllConfigurations(maximumConfigurations: Int = 100_000) throws -> [Molecule] {
        var result: [Molecule] = []
        try enumerateConfigurations(maximumConfigurations: maximumConfigurations) {
            result.append($0)
            return true
        }
        return result
    }

    /// Streams configurations in R-number, occurrence, site and alternative order.
    /// Return false to stop early. Returns the number delivered, including the final callback.
    /// Checks Swift task cancellation between search steps. No configurations are deduplicated.
    @discardableResult
    public func enumerateConfigurations(
        maximumConfigurations: Int = 100_000,
        checkpoint: () throws -> Void = {},
        _ visit: (Molecule) throws -> Bool
    ) throws -> Int {
        guard maximumConfigurations > 0 else {
            throw CDKRGroupEnumerationError.invalidQuery(
                "The configuration limit must be positive.")
        }
        try Task.checkCancellation()
        try checkpoint()
        let plan = try RGroupEnumerationPlan(query: self)
        guard !plan.groups.isEmpty else { return 0 }
        guard plan.groups.allSatisfy({ !$0.occurrences.isEmpty }) else { return 0 }
        var assignments: [Int: Int] = [:]  // Root site ID -> alternative index, or -1 for absence.
        var presence: [Int: Bool] = [:]
        var count = 0
        var stopped = false

        func conditionsMet() -> Bool {
            for group in plan.groups where presence[group.number] == true {
                let required = group.definition.requiredRGroupNumber
                if required != 0 {
                    if presence[required] == false { return false }
                    if plan.groups.first(where: { $0.number == required })?.occurrences.allSatisfy({
                        $0 == 0
                    }
                    ) == true {
                        return false
                    }
                }
            }
            return true
        }

        func searchGroup(_ level: Int) throws {
            try Task.checkCancellation()
            try checkpoint()
            guard !stopped else { return }
            if level == plan.groups.count {
                guard conditionsMet() else { return }
                guard count < maximumConfigurations else {
                    throw CDKRGroupEnumerationError.configurationLimitExceeded(
                        maximumConfigurations)
                }
                let molecule = try plan.configuration(assignments)
                count += 1
                stopped = try !visit(molecule)
                return
            }
            let group = plan.groups[level]
            func searchSite(_ index: Int, _ remaining: Int) throws {
                try Task.checkCancellation()
                try checkpoint()
                guard !stopped else { return }
                let available = group.sites.count - index
                guard remaining >= 0, remaining <= available else { return }
                if index == group.sites.count {
                    try searchGroup(level + 1)
                    return
                }
                let id = group.sites[index].id
                if remaining < available {
                    assignments[id] = -1
                    try searchSite(index + 1, remaining)
                }
                if remaining > 0 {
                    for alternative in group.definition.rGroups.indices {
                        assignments[id] = alternative
                        try searchSite(index + 1, remaining - 1)
                        if stopped { return }
                    }
                }
            }
            for occurrence in group.occurrences {
                presence[group.number] = occurrence > 0
                if conditionsMet() { try searchSite(0, occurrence) }
                if stopped { return }
            }
            presence[group.number] = nil
        }
        try searchGroup(0)
        return count
    }
}

struct RGroupEnumerationPlan {
    struct Group {
        var number: Int
        var definition: CDKRGroupList
        var sites: [Atom]
        var occurrences: [Int]
    }
    struct Fragment {
        var source: Molecule
        var atomMap: [Int: Int]
        var dummyIDs: Set<Int>
        var dummyRootBonds: [Int: Bond]
        var offset: CGPoint
    }
    var query: CDKRGroupQuery
    var groups: [Group] = []
    var attachments: [Int: [Bond]] = [:]
    var removedRootIDs: Set<Int>

    init(query: CDKRGroupQuery) throws {
        self.query = query
        removedRootIDs = Set(
            query.rootStructure.atoms.filter {
                $0.rGroupLabel != nil || $0.element == "R" || $0.element == "R#"
            }.map(\.id))
        try Self.validateGraph(query.rootStructure)
        guard query.areSubstituentsDefined(), query.areRootAtomsDefined() else {
            throw CDKRGroupEnumerationError.invalidQuery(
                "Missing R-group definitions or scaffold sites.")
        }
        guard query.rGroupQueryAtoms().count <= 256 else {
            throw CDKRGroupEnumerationError.unsupportedFeature("More than 256 R-group sites.")
        }
        let positions = Dictionary(
            uniqueKeysWithValues: query.rootStructure.atoms.enumerated().map { ($1.id, $0) })
        for number in query.rGroupDefinitions.keys.sorted() {
            let definition = query.rGroupDefinitions[number]!
            guard number >= 0, definition.rGroupNumber == number,
                CDKRGroupList.isValidOccurrenceSyntax(definition.occurrence),
                definition.requiredRGroupNumber == 0
                    || query.rGroupDefinitions[definition.requiredRGroupNumber] != nil
            else {
                throw CDKRGroupEnumerationError.invalidQuery(
                    "Invalid definition or dependency for R\(number).")
            }
            let sites = query.rGroupQueryAtoms(number)
            for site in sites {
                let connected = query.rootStructure.bonds(forAtom: site.id)
                if let order = site.ligandOrderingAtomIDs {
                    guard order.count == connected.count, Set(order).count == order.count,
                        Set(order) == Set(connected.map { $0.a1 == site.id ? $0.a2 : $0.a1 })
                    else {
                        throw CDKRGroupEnumerationError.invalidQuery(
                            "Invalid attachment order at atom \(site.id).")
                    }
                    attachments[site.id] = order.map { partner in
                        connected.first { $0.a1 == partner || $0.a2 == partner }!
                    }
                } else {
                    attachments[site.id] = connected.sorted {
                        positions[$0.a1 == site.id ? $0.a2 : $0.a1]! < positions[
                            $1.a1 == site.id ? $1.a2 : $1.a1]!
                    }
                }
            }
            for fragment in definition.rGroups {
                try Self.validateGraph(fragment.group)
                guard fragment.group.atoms.allSatisfy({ $0.rGroupLabel == nil }) else {
                    throw CDKRGroupEnumerationError.unsupportedFeature(
                        "Nested R-group definitions.")
                }
                // Validate all declared attachment points, even if no configuration uses them.
                for ordinal in 1...2 { _ = try Self.attachment(fragment, ordinal: ordinal) }
            }
            groups.append(
                Group(
                    number: number, definition: definition, sites: sites,
                    occurrences: definition.matchOccurrence(maxAttachments: sites.count)))
        }
    }

    static func validateGraph(_ molecule: Molecule) throws {
        let ids = Set(molecule.atoms.map(\.id))
        guard ids.count == molecule.atoms.count,
            Set(molecule.bonds.map(\.id)).count == molecule.bonds.count,
            molecule.atoms.allSatisfy({ Set($0.ligandOrderingAtomIDs ?? []).isSubset(of: ids) }),
            molecule.bonds.allSatisfy({ Set($0.stereoReferenceAtomIDs ?? []).isSubset(of: ids) }),
            molecule.bonds.allSatisfy({
                ids.contains($0.a1) && ids.contains($0.a2) && $0.a1 != $0.a2
            })
        else {
            throw CDKRGroupEnumerationError.invalidQuery("Duplicate IDs or invalid bond endpoints.")
        }
        guard
            !molecule.sgroups.contains(where: {
                [.structureRepeatUnit, .polymer, .extMulticenter].contains($0.kind)
            }),
            molecule.cxState?.linkNodes.isEmpty != false,
            !molecule.name.trimmingCharacters(in: .whitespacesAndNewlines).hasPrefix("|")
        else {
            throw CDKRGroupEnumerationError.unsupportedFeature(
                "Repeat units, multicenter attachments or uninterpreted CX extensions require separate expansion."
            )
        }
    }

    static func isDummy(_ atom: Atom) -> Bool {
        atom.element == "*" && (atom.attachmentPoint ?? 0) > 0
    }

    /// Returns the real atom plus an optional explicit AP dummy atom.
    static func attachment(_ fragment: CDKRGroup, ordinal: Int) throws -> (atom: Int, dummy: Int?)? {
        let explicit =
            ordinal == 1
            ? fragment.firstAttachmentPointAtomID : fragment.secondAttachmentPointAtomID
        let candidates = fragment.group.atoms.filter {
            $0.attachmentPoint == ordinal || $0.attachmentPoint == 3
        }
        guard explicit != nil || candidates.count <= 1 else {
            throw CDKRGroupEnumerationError.invalidQuery("Ambiguous attachment point \(ordinal).")
        }
        guard let id = explicit ?? candidates.first?.id else { return nil }
        guard let atom = fragment.group.atom(id: id) else {
            throw CDKRGroupEnumerationError.invalidQuery("Attachment references a missing atom.")
        }
        if !isDummy(atom) { return (id, nil) }
        let bonds = fragment.group.bonds(forAtom: id)
        guard bonds.count == 1,
            let neighbor = fragment.group.atom(id: bonds[0].a1 == id ? bonds[0].a2 : bonds[0].a1),
            !isDummy(neighbor)
        else {
            throw CDKRGroupEnumerationError.invalidQuery(
                "An explicit attachment dummy must have one real neighbor.")
        }
        return (neighbor.id, id)
    }

    func configuration(_ assignments: [Int: Int]) throws -> Molecule {
        var result = query.rootStructure
        var nextAtomID = result.atoms.map(\.id).max() ?? 0
        var nextBondID = result.bonds.map(\.id).max() ?? 0
        var replacements: [Int: [Int: Int]] = [:]  // Bond -> root site -> replacement atom.
        var fragments: [Fragment] = []
        func nextID(_ value: inout Int) throws -> Int {
            let increment = value.addingReportingOverflow(1)
            guard !increment.overflow else {
                throw CDKRGroupEnumerationError.invalidQuery("Atom or bond ID overflow.")
            }
            value = increment.partialValue
            return value
        }
        for group in groups {
            for site in group.sites {
                guard let index = assignments[site.id], index >= 0 else { continue }
                let chosen = group.definition.rGroups[index]
                let dummies = Set(chosen.group.atoms.filter(Self.isDummy).map(\.id))
                var map: [Int: Int] = [:]
                for atom in chosen.group.atoms where !dummies.contains(atom.id) {
                    map[atom.id] = try nextID(&nextAtomID)
                }
                var dummyBonds: [Int: Bond] = [:]
                for (ordinal, bond) in attachments[site.id, default: []].enumerated() {
                    // CDK's lenient RGfile behavior retains an unattached fragment when APO is absent.
                    guard let attachment = try Self.attachment(chosen, ordinal: min(ordinal + 1, 2))
                    else {
                        continue
                    }
                    replacements[bond.id, default: [:]][site.id] = map[attachment.atom]!
                    if let dummy = attachment.dummy { dummyBonds[dummy] = bond }
                }
                let first = try Self.attachment(chosen, ordinal: 1).flatMap {
                    chosen.group.atom(id: $0.atom)
                }
                let offset =
                    first.map {
                        CGPoint(
                            x: site.position.x - $0.position.x, y: site.position.y - $0.position.y)
                    } ?? .zero
                fragments.append(
                    Fragment(
                        source: chosen.group, atomMap: map, dummyIDs: dummies,
                        dummyRootBonds: dummyBonds, offset: offset))
            }
        }
        func rootEndpoint(_ id: Int, bond: Bond) -> Int? {
            removedRootIDs.contains(id) ? replacements[bond.id]?[id] : id
        }
        func rootReference(_ id: Int, center: Int) -> Int? {
            guard removedRootIDs.contains(id) else { return id }
            guard let bond = query.rootStructure.bond(between: id, and: center) else { return nil }
            return replacements[bond.id]?[id]
        }
        result.atoms.removeAll { removedRootIDs.contains($0.id) }
        for index in result.atoms.indices {
            let id = result.atoms[index].id
            result.atoms[index].ligandOrderingAtomIDs = result.atoms[index].ligandOrderingAtomIDs?
                .map {
                    rootReference($0, center: id) ?? id  // implicit H occupies the omitted ligand's position
                }
        }
        result.bonds = result.bonds.compactMap { bond in
            guard let a1 = rootEndpoint(bond.a1, bond: bond),
                let a2 = rootEndpoint(bond.a2, bond: bond)
            else { return nil }
            return Self.copyBond(
                bond, id: bond.id, a1: a1, a2: a2,
                mapReference: { ref in
                    if ref == bond.a1 { return a1 }
                    if ref == bond.a2 { return a2 }
                    return rootReference(ref, center: bond.a1)
                        ?? rootReference(ref, center: bond.a2)
                })
        }
        for fragment in fragments {
            func reference(_ id: Int) -> Int? {
                if let mapped = fragment.atomMap[id] { return mapped }
                guard let bond = fragment.dummyRootBonds[id] else { return nil }
                // The endpoint attached to this fragment is the one in its atom map.
                let ownIDs = Set(fragment.atomMap.values)
                let a1 = rootEndpoint(bond.a1, bond: bond)
                let a2 = rootEndpoint(bond.a2, bond: bond)
                return a1.map { ownIDs.contains($0) } == true ? a2 : a1
            }
            for atom in fragment.source.atoms where !fragment.dummyIDs.contains(atom.id) {
                var copied = Self.copyAtom(atom, id: fragment.atomMap[atom.id]!)
                copied.position.x += fragment.offset.x
                copied.position.y += fragment.offset.y
                copied.rGroupMembership = nil
                copied.componentGroupID = nil
                copied.attachmentPoint = nil
                copied.ligandOrderingAtomIDs = atom.ligandOrderingAtomIDs?.map {
                    reference($0) ?? copied.id
                }
                result.atoms.append(copied)
            }
            for bond in fragment.source.bonds {
                guard let a1 = fragment.atomMap[bond.a1], let a2 = fragment.atomMap[bond.a2] else {
                    continue
                }
                result.bonds.append(
                    Self.copyBond(
                        bond, id: try nextID(&nextBondID), a1: a1, a2: a2, mapReference: reference))
            }
        }
        for group in groups {
            for site in group.sites where assignments[site.id] == -1 {
                for bond in attachments[site.id, default: []] {
                    let other = bond.a1 == site.id ? bond.a2 : bond.a1
                    guard let id = rootEndpoint(other, bond: bond),
                        let index = result.indexOfAtom(id: id)
                    else { continue }
                    result.atoms[index].rGroupRestH = group.definition.restH
                    if let h = result.atoms[index].explicitHydrogenCount {
                        let updated = h.addingReportingOverflow(bond.order.rawValue)
                        guard h >= 0, !updated.overflow else {
                            throw CDKRGroupEnumerationError.invalidQuery("Invalid hydrogen count.")
                        }
                        result.atoms[index].explicitHydrogenCount = updated.partialValue
                    }
                }
            }
        }
        result.rGroupLogicDefinitions = [:]
        result.sgroups = []  // Like CDK; variable Sgroups were rejected before enumeration.
        result.cxState = nil  // Its source-indexed tables no longer refer to the assembled graph.
        result.highlightedAtomIDs.removeAll { removedRootIDs.contains($0) }
        let bondIDs = Set(result.bonds.map(\.id))
        result.highlightedBondIDs.removeAll { !bondIDs.contains($0) }
        return result
    }

    static func copyBond(_ bond: Bond, id: Int, a1: Int, a2: Int, mapReference: (Int) -> Int?)
        -> Bond
    {
        let references = bond.stereoReferenceAtomIDs?.compactMap(mapReference)
        let validStereo = references?.count == bond.stereoReferenceAtomIDs?.count
        return Bond(
            id: id, externalID: bond.externalID, a1: a1, a2: a2, order: bond.order,
            stereo: bond.stereo,
            doubleBondStereo: validStereo ? bond.doubleBondStereo : nil,
            stereoReferenceAtomIDs: validStereo ? references : nil, queryType: bond.queryType,
            topology: bond.topology,
            coordinateBondReferenceAtomID: bond.coordinateBondReferenceAtomID.flatMap(mapReference),
            reactingCenterStatus: bond.reactingCenterStatus)
    }

    static func copyAtom(_ atom: Atom, id: Int) -> Atom {
        Atom(
            id: id, externalID: atom.externalID, element: atom.element, position: atom.position,
            zPosition: atom.zPosition, charge: atom.charge,
            isotopeMassNumber: atom.isotopeMassNumber,
            aromatic: atom.aromatic, chirality: atom.chirality,
            explicitHydrogenCount: atom.explicitHydrogenCount,
            queryType: atom.queryType, atomList: atom.atomList,
            atomListIsNegated: atom.atomListIsNegated,
            radical: atom.radical, radicalType: atom.radicalType, atomValue: atom.atomValue,
            rGroupLabel: atom.rGroupLabel, rGroupMembership: atom.rGroupMembership,
            componentGroupID: atom.componentGroupID,
            reactionRole: atom.reactionRole, substitutionCount: atom.substitutionCount,
            unsaturated: atom.unsaturated,
            ringBondCount: atom.ringBondCount, attachmentPoint: atom.attachmentPoint,
            valenceOverride: atom.valenceOverride,
            cxStereoGroup: atom.cxStereoGroup, ligandOrderingAtomIDs: atom.ligandOrderingAtomIDs,
            atomClass: atom.atomClass, atomMapNumber: atom.atomMapNumber,
            aliasLabel: atom.aliasLabel,
            properties: atom.properties, atomTypeName: atom.atomTypeName,
            maximumBondOrder: atom.maximumBondOrder,
            bondOrderSum: atom.bondOrderSum, valency: atom.valency,
            formalNeighbourCount: atom.formalNeighbourCount,
            hybridization: atom.hybridization)
    }
}
