import Foundation

public struct CDKReactionAtomReference: Equatable, Hashable, Sendable {
    public let role: CDKReactionRole
    public let participantIndex: Int
    public let atom: Atom

    public init(role: CDKReactionRole, participantIndex: Int, atom: Atom) {
        self.role = role
        self.participantIndex = participantIndex
        self.atom = atom
    }
}

public struct CDKReactionBondReference: Equatable, Hashable, Sendable {
    public let role: CDKReactionRole
    public let participantIndex: Int
    public let bond: Bond

    public init(role: CDKReactionRole, participantIndex: Int, bond: Bond) {
        self.role = role
        self.participantIndex = participantIndex
        self.bond = bond
    }
}

public enum CDKReactionChemObject: Equatable {
    case reaction(CDKReaction)
    case reactionSet(CDKReactionSet)
    case reactionList(CDKReactionList)
    case reactionScheme(CDKReactionScheme)
    case molecule(Molecule)
}

public enum CDKReactionManipulator {
    private static let propertyPrefix = "__cdkReactionProperty__:"
    private static let directionFieldName = "__cdkReactionDirection__"

    public static func getAtomCount(_ reaction: CDKReaction) -> Int {
        getAllMolecules(reaction).reduce(0) { $0 + $1.atomCount }
    }

    public static func getBondCount(_ reaction: CDKReaction) -> Int {
        getAllMolecules(reaction).reduce(0) { $0 + $1.bondCount }
    }

    public static func getAllMolecules(_ reaction: CDKReaction) -> [Molecule] {
        getAllReactants(reaction) + getAllAgents(reaction) + getAllProducts(reaction)
    }

    public static func getAllReactants(_ reaction: CDKReaction) -> [Molecule] {
        reaction.reactants
    }

    public static func getAllAgents(_ reaction: CDKReaction) -> [Molecule] {
        reaction.agents
    }

    public static func getAllProducts(_ reaction: CDKReaction) -> [Molecule] {
        reaction.products
    }

    public static func getAllAtomContainers(_ reaction: CDKReaction) -> [Molecule] {
        getAllMolecules(reaction)
    }

    public static func reverse(_ reaction: CDKReaction) -> CDKReaction {
        CDKReaction(reactantParticipants: reaction.productParticipants,
                    agentParticipants: reaction.agentParticipants,
                    productParticipants: reaction.reactantParticipants,
                    id: reaction.id,
                    direction: reaction.direction.reversed,
                    name: reaction.name,
                    properties: reaction.properties,
                    cxState: reaction.cxState)
    }

    public static func getAllIDs(_ reaction: CDKReaction) -> [String] {
        var ids: [String] = []
        if let id = reaction.id {
            ids.append(id)
        }
        for molecule in getAllMolecules(reaction) {
            appendIDs(from: molecule, into: &ids)
        }
        return ids
    }

    public static func getRelevantAtomContainer(_ reaction: CDKReaction, atom: Atom) -> Molecule? {
        getAllMolecules(reaction).first { molecule in
            molecule.atoms.contains(atom)
        }
    }

    public static func getRelevantAtomContainer(_ reaction: CDKReaction, bond: Bond) -> Molecule? {
        getAllMolecules(reaction).first { molecule in
            molecule.bonds.contains(bond)
        }
    }

    public static func setAtomProperties(_ reaction: inout CDKReaction,
                                         key: String,
                                         value: String) {
        mutateAllParticipants(in: &reaction) { molecule in
            var updated = molecule
            updated.atoms = updated.atoms.map { atom in
                var atom = atom
                atom.properties[key] = value
                return atom
            }
            return updated
        }
    }

    public static func getAllChemObjects(_ reaction: CDKReaction) -> [CDKReactionChemObject] {
        [.reaction(reaction)] + getAllAtomContainers(reaction).map(CDKReactionChemObject.molecule)
    }

    @discardableResult
    public static func removeAtomAndConnectedElectronContainers(_ reaction: inout CDKReaction,
                                                                _ atom: Atom) -> Bool {
        removeMatchingObject(from: &reaction) { molecule in
            molecule.atoms.contains(atom)
        } mutation: { molecule in
            guard molecule.atoms.contains(atom) else { return molecule }
            var updated = molecule
            updated.atoms.removeAll { $0 == atom }
            updated.bonds.removeAll { $0.a1 == atom.id || $0.a2 == atom.id }
            updated.sgroups = updated.sgroups.compactMap { sgroup in
                let atomIDs = sgroup.atomIDs.filter { $0 != atom.id }
                let crossingBondIDs = sgroup.crossingBondIDs.filter { bondID in
                    updated.bonds.contains(where: { $0.id == bondID })
                }
                guard !atomIDs.isEmpty else { return nil }
                return MoleculeSgroup(kind: sgroup.kind,
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
                                      parentAtomIDs: sgroup.parentAtomIDs.filter { $0 != atom.id },
                                      componentNumber: sgroup.componentNumber,
                                      expanded: sgroup.expanded,
                                      brackets: sgroup.brackets,
                                      childGroupIndices: sgroup.childGroupIndices)
            }
            updated.highlightedAtomIDs.removeAll { $0 == atom.id }
            updated.highlightedBondIDs.removeAll { bondID in
                !updated.bonds.contains(where: { $0.id == bondID })
            }
            return updated
        }
    }

    @discardableResult
    public static func removeElectronContainer(_ reaction: inout CDKReaction,
                                               _ bond: Bond) -> Bool {
        removeMatchingObject(from: &reaction) { molecule in
            molecule.bonds.contains(bond)
        } mutation: { molecule in
            guard molecule.bonds.contains(bond) else { return molecule }
            var updated = molecule
            updated.bonds.removeAll { $0 == bond }
            updated.sgroups = updated.sgroups.map { sgroup in
                MoleculeSgroup(kind: sgroup.kind,
                               keyword: sgroup.keyword,
                               atomIDs: sgroup.atomIDs,
                               crossingBondIDs: sgroup.crossingBondIDs.filter { $0 != bond.id },
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
                               parentAtomIDs: sgroup.parentAtomIDs,
                               componentNumber: sgroup.componentNumber,
                               expanded: sgroup.expanded,
                               brackets: sgroup.brackets,
                               childGroupIndices: sgroup.childGroupIndices)
            }
            updated.highlightedBondIDs.removeAll { $0 == bond.id }
            return updated
        }
    }

    public static func getMappedChemObject(_ reaction: CDKReaction, _ atom: Atom) -> Atom? {
        guard let mapNumber = atom.atomMapNumber, mapNumber > 0 else { return nil }
        let sourceRole = participantReference(in: reaction, atom: atom)?.role
        return participantGroups(forSearchFrom: sourceRole, in: reaction)
            .flatMap(\.participants)
            .flatMap(\.molecule.atoms)
            .first { $0.atomMapNumber == mapNumber }
    }

    public static func getMappedChemObject(_ reaction: CDKReaction, _ bond: Bond) -> Bond? {
        guard let source = participantReference(in: reaction, bond: bond),
              let sourceMap = mappedBondKey(for: source.bond, in: source.molecule) else {
            return nil
        }

        for group in participantGroups(forSearchFrom: source.role, in: reaction) {
            for participant in group.participants {
                if let matched = participant.molecule.bonds.first(where: {
                    mappedBondKey(for: $0, in: participant.molecule) == sourceMap
                }) {
                    return matched
                }
            }
        }
        return nil
    }

    public static func findMappedBonds(_ reaction: CDKReaction) -> Set<CDKReactionBondReference> {
        let leftParticipants = indexedParticipants(for: reaction.reactantParticipants, role: .reactant)
            + indexedParticipants(for: reaction.agentParticipants, role: .agent)
        let rightParticipants = indexedParticipants(for: reaction.productParticipants, role: .product)

        let leftKeys = Set(leftParticipants.flatMap { participant in
            participant.molecule.bonds.compactMap { mappedBondKey(for: $0, in: participant.molecule) }
        })
        let rightKeys = Set(rightParticipants.flatMap { participant in
            participant.molecule.bonds.compactMap { mappedBondKey(for: $0, in: participant.molecule) }
        })
        let shared = leftKeys.intersection(rightKeys)
        guard !shared.isEmpty else { return [] }

        var references = Set<CDKReactionBondReference>()
        for participant in leftParticipants + rightParticipants {
            for bond in participant.molecule.bonds {
                guard let key = mappedBondKey(for: bond, in: participant.molecule),
                      shared.contains(key) else {
                    continue
                }
                references.insert(CDKReactionBondReference(role: participant.role,
                                                           participantIndex: participant.index,
                                                           bond: bond))
            }
        }
        return references
    }

    public static func toMolecule(_ reaction: CDKReaction) -> Molecule {
        var merged = Molecule(name: normalizedReactionName(reaction.name, fallback: "Reaction"),
                              externalID: reaction.id)
        merged.cxState = reaction.cxState

        if !reaction.properties.isEmpty {
            for key in reaction.properties.keys.sorted() {
                guard let value = reaction.properties[key] else { continue }
                merged.appendDataFieldValue(value, named: propertyPrefix + key)
            }
        }
        merged.appendDataFieldValue(reaction.direction.rawValue, named: directionFieldName)

        var nextAtomID = 0
        var nextBondID = 0
        var nextComponentGroupID = 1
        let orderedParticipants = indexedParticipants(for: reaction.reactantParticipants, role: .reactant)
            + indexedParticipants(for: reaction.agentParticipants, role: .agent)
            + indexedParticipants(for: reaction.productParticipants, role: .product)

        for participant in orderedParticipants {
            appendParticipant(participant.participant,
                              to: &merged,
                              role: participant.role,
                              componentGroupID: nextComponentGroupID,
                              nextAtomID: &nextAtomID,
                              nextBondID: &nextBondID)
            nextComponentGroupID += 1
        }

        return merged
    }

    public static func toReaction(_ molecule: Molecule) throws -> CDKReaction {
        let groupedAtoms = Dictionary(grouping: molecule.atoms) { $0.componentGroupID ?? -1 }
        var groups: [(groupID: Int, role: CDKReactionRole, atoms: [Atom], bonds: [Bond], sgroups: [MoleculeSgroup])] = []

        for groupID in groupedAtoms.keys.sorted() {
            guard groupID > 0, let groupAtoms = groupedAtoms[groupID], !groupAtoms.isEmpty else {
                throw ChemError.parseFailed("Inlined reaction molecule contains atoms without a reaction group.")
            }
            let roles = Set(groupAtoms.compactMap(\.reactionRole))
            guard roles.count == 1, let role = roles.first else {
                throw ChemError.parseFailed("Inlined reaction molecule contains a group without a unique reaction role.")
            }

            let atomIDSet = Set(groupAtoms.map(\.id))
            let groupBonds = molecule.bonds.filter { atomIDSet.contains($0.a1) || atomIDSet.contains($0.a2) }
            if groupBonds.contains(where: { !atomIDSet.contains($0.a1) || !atomIDSet.contains($0.a2) }) {
                throw ChemError.parseFailed("Inlined reaction molecule contains a bond crossing reaction groups.")
            }
            let groupBondIDs = Set(groupBonds.map(\.id))
            let groupSgroups = molecule.sgroups.compactMap { sgroup -> MoleculeSgroup? in
                let sgroupAtomIDs = Set(sgroup.atomIDs)
                guard !sgroupAtomIDs.isDisjoint(with: atomIDSet) else { return nil }
                guard sgroupAtomIDs.isSubset(of: atomIDSet) else {
                    return nil
                }
                return MoleculeSgroup(kind: sgroup.kind,
                                      keyword: sgroup.keyword,
                                      atomIDs: sgroup.atomIDs,
                                      crossingBondIDs: sgroup.crossingBondIDs.filter { groupBondIDs.contains($0) },
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
                                      parentAtomIDs: sgroup.parentAtomIDs.filter { atomIDSet.contains($0) },
                                      componentNumber: sgroup.componentNumber,
                                      expanded: sgroup.expanded,
                                      brackets: sgroup.brackets,
                                      childGroupIndices: sgroup.childGroupIndices)
            }

            groups.append((groupID, role, groupAtoms, groupBonds, groupSgroups))
        }

        let sortedGroups = groups.sorted { $0.groupID < $1.groupID }
        var reactants: [CDKReactionParticipant] = []
        var agents: [CDKReactionParticipant] = []
        var products: [CDKReactionParticipant] = []

        for group in sortedGroups {
            let participantMolecule = moleculeFromInlineGroup(group)
            let participant = CDKReactionParticipant(molecule: participantMolecule, role: group.role)
            switch group.role {
            case .reactant:
                reactants.append(participant)
            case .agent:
                agents.append(participant)
            case .product:
                products.append(participant)
            }
        }

        let propertyPairs: [(String, String)] = molecule.dataFields.compactMap { key, values in
            guard key.hasPrefix(propertyPrefix) else { return nil }
            return (String(key.dropFirst(propertyPrefix.count)), values.joined(separator: "\n"))
        }
        let properties = Dictionary(uniqueKeysWithValues: propertyPairs)
        let direction = molecule.dataFieldValues(named: directionFieldName).first
            .flatMap(CDKReactionDirection.init(rawValue:)) ?? .forward
        let reactionName = molecule.name.trimmingCharacters(in: .whitespacesAndNewlines)

        return CDKReaction(reactantParticipants: reactants,
                           agentParticipants: agents,
                           productParticipants: products,
                           id: molecule.externalID,
                           direction: direction,
                           name: reactionName.isEmpty ? nil : reactionName,
                           properties: properties,
                           cxState: molecule.cxState)
    }

    public static func perceiveAtomTypesAndConfigureAtoms(_ reaction: inout CDKReaction) {
        mutateAllParticipants(in: &reaction) { molecule in
            configured(molecule: molecule, onlyUnset: false)
        }
    }

    public static func perceiveAtomTypesAndConfigureUnsetProperties(_ reaction: inout CDKReaction) {
        mutateAllParticipants(in: &reaction) { molecule in
            configured(molecule: molecule, onlyUnset: true)
        }
    }

    public static func clearAtomConfigurations(_ reaction: inout CDKReaction) {
        mutateAllParticipants(in: &reaction) { molecule in
            var updated = molecule
            updated.atoms = updated.atoms.map { atom in
                var atom = atom
                atom.atomTypeName = nil
                atom.maximumBondOrder = nil
                atom.bondOrderSum = nil
                atom.valency = nil
                atom.formalNeighbourCount = nil
                atom.hybridization = nil
                return atom
            }
            return updated
        }
    }

    private static func normalizedReactionName(_ raw: String?, fallback: String) -> String {
        let trimmed = (raw ?? "").trimmingCharacters(in: .whitespacesAndNewlines)
        return trimmed.isEmpty ? fallback : trimmed
    }

    private static func appendIDs(from molecule: Molecule, into ids: inout [String]) {
        if let externalID = molecule.externalID {
            ids.append(externalID)
        }
        ids.append(contentsOf: molecule.atoms.compactMap(\.externalID))
        ids.append(contentsOf: molecule.bonds.compactMap(\.externalID))
    }

    private static func indexedParticipants(for participants: [CDKReactionParticipant],
                                            role: CDKReactionRole)
        -> [(index: Int, role: CDKReactionRole, participant: CDKReactionParticipant, molecule: Molecule)]
    {
        participants.enumerated().map { index, participant in
            (index, role, participant, participant.molecule)
        }
    }

    private static func participantReference(in reaction: CDKReaction,
                                             atom: Atom)
        -> (role: CDKReactionRole, index: Int, molecule: Molecule)?
    {
        for participant in indexedParticipants(for: reaction.reactantParticipants, role: .reactant) {
            if participant.molecule.atoms.contains(atom) {
                return (participant.role, participant.index, participant.molecule)
            }
        }
        for participant in indexedParticipants(for: reaction.agentParticipants, role: .agent) {
            if participant.molecule.atoms.contains(atom) {
                return (participant.role, participant.index, participant.molecule)
            }
        }
        for participant in indexedParticipants(for: reaction.productParticipants, role: .product) {
            if participant.molecule.atoms.contains(atom) {
                return (participant.role, participant.index, participant.molecule)
            }
        }
        return nil
    }

    private static func participantReference(in reaction: CDKReaction,
                                             bond: Bond)
        -> (role: CDKReactionRole, index: Int, molecule: Molecule, bond: Bond)?
    {
        for participant in indexedParticipants(for: reaction.reactantParticipants, role: .reactant) {
            if participant.molecule.bonds.contains(bond) {
                return (participant.role, participant.index, participant.molecule, bond)
            }
        }
        for participant in indexedParticipants(for: reaction.agentParticipants, role: .agent) {
            if participant.molecule.bonds.contains(bond) {
                return (participant.role, participant.index, participant.molecule, bond)
            }
        }
        for participant in indexedParticipants(for: reaction.productParticipants, role: .product) {
            if participant.molecule.bonds.contains(bond) {
                return (participant.role, participant.index, participant.molecule, bond)
            }
        }
        return nil
    }

    private static func participantGroups(forSearchFrom sourceRole: CDKReactionRole?,
                                          in reaction: CDKReaction)
        -> [(role: CDKReactionRole, participants: [(index: Int, role: CDKReactionRole, participant: CDKReactionParticipant, molecule: Molecule)])]
    {
        switch sourceRole {
        case .product:
            return [
                (.reactant, indexedParticipants(for: reaction.reactantParticipants, role: .reactant)),
                (.agent, indexedParticipants(for: reaction.agentParticipants, role: .agent))
            ]
        default:
            return [
                (.product, indexedParticipants(for: reaction.productParticipants, role: .product))
            ]
        }
    }

    private static func mappedBondKey(for bond: Bond, in molecule: Molecule) -> [Int]? {
        guard let a1 = molecule.atom(id: bond.a1)?.atomMapNumber,
              let a2 = molecule.atom(id: bond.a2)?.atomMapNumber,
              a1 > 0, a2 > 0 else {
            return nil
        }
        return [a1, a2].sorted()
    }

    private static func removeMatchingObject(from reaction: inout CDKReaction,
                                             contains: (Molecule) -> Bool,
                                             mutation: (Molecule) -> Molecule) -> Bool {
        let reactantRemoved = mutateParticipants(&reaction.reactantParticipants, contains: contains, mutation: mutation)
        let agentRemoved = mutateParticipants(&reaction.agentParticipants, contains: contains, mutation: mutation)
        let productRemoved = mutateParticipants(&reaction.productParticipants, contains: contains, mutation: mutation)
        return reactantRemoved || agentRemoved || productRemoved
    }

    private static func mutateParticipants(_ participants: inout [CDKReactionParticipant],
                                           contains: (Molecule) -> Bool,
                                           mutation: (Molecule) -> Molecule) -> Bool {
        var changed = false
        for index in participants.indices where contains(participants[index].molecule) {
            participants[index].molecule = mutation(participants[index].molecule)
            changed = true
        }
        return changed
    }

    private static func appendParticipant(_ participant: CDKReactionParticipant,
                                          to molecule: inout Molecule,
                                          role: CDKReactionRole,
                                          componentGroupID: Int,
                                          nextAtomID: inout Int,
                                          nextBondID: inout Int) {
        var atomIDMap: [Int: Int] = [:]
        var bondIDMap: [Int: Int] = [:]

        for atom in participant.molecule.atoms {
            nextAtomID += 1
            atomIDMap[atom.id] = nextAtomID
            molecule.atoms.append(
                Atom(id: nextAtomID,
                     externalID: atom.externalID,
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
                     componentGroupID: componentGroupID,
                     reactionRole: role,
                     substitutionCount: atom.substitutionCount,
                     unsaturated: atom.unsaturated,
                     ringBondCount: atom.ringBondCount,
                     attachmentPoint: atom.attachmentPoint,
                     valenceOverride: atom.valenceOverride,
                     cxStereoGroup: atom.cxStereoGroup,
                     ligandOrderingAtomIDs: atom.ligandOrderingAtomIDs,
                     atomClass: atom.atomClass,
                     atomMapNumber: atom.atomMapNumber,
                     aliasLabel: atom.aliasLabel,
                     properties: atom.properties,
                     atomTypeName: atom.atomTypeName,
                     maximumBondOrder: atom.maximumBondOrder,
                     bondOrderSum: atom.bondOrderSum,
                     valency: atom.valency,
                     formalNeighbourCount: atom.formalNeighbourCount,
                     hybridization: atom.hybridization)
            )
        }

        for bond in participant.molecule.bonds {
            guard let a1 = atomIDMap[bond.a1], let a2 = atomIDMap[bond.a2] else { continue }
            nextBondID += 1
            bondIDMap[bond.id] = nextBondID
            molecule.bonds.append(Bond(id: nextBondID,
                                       externalID: bond.externalID,
                                       a1: a1,
                                       a2: a2,
                                       order: bond.order,
                                       stereo: bond.stereo,
                                       queryType: bond.queryType,
                                       topology: bond.topology))
        }

        for sgroup in participant.molecule.sgroups {
            let atomIDs = sgroup.atomIDs.compactMap { atomIDMap[$0] }
            guard !atomIDs.isEmpty else { continue }
            let crossingBondIDs = sgroup.crossingBondIDs.compactMap { bondIDMap[$0] }
            molecule.sgroups.append(
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
                               parentAtomIDs: sgroup.parentAtomIDs.compactMap { atomIDMap[$0] },
                               componentNumber: sgroup.componentNumber,
                               expanded: sgroup.expanded,
                               brackets: sgroup.brackets,
                               childGroupIndices: sgroup.childGroupIndices)
            )
        }
    }

    private static func moleculeFromInlineGroup(_ group: (groupID: Int,
                                                         role: CDKReactionRole,
                                                         atoms: [Atom],
                                                         bonds: [Bond],
                                                         sgroups: [MoleculeSgroup])) -> Molecule {
        var participant = Molecule(name: "Participant")
        participant.atoms = group.atoms.map { atom in
            Atom(id: atom.id,
                 externalID: atom.externalID,
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
                 componentGroupID: nil,
                 reactionRole: nil,
                 substitutionCount: atom.substitutionCount,
                 unsaturated: atom.unsaturated,
                 ringBondCount: atom.ringBondCount,
                 attachmentPoint: atom.attachmentPoint,
                 valenceOverride: atom.valenceOverride,
                 cxStereoGroup: atom.cxStereoGroup,
                 ligandOrderingAtomIDs: atom.ligandOrderingAtomIDs,
                 atomClass: atom.atomClass,
                 atomMapNumber: atom.atomMapNumber,
                 aliasLabel: atom.aliasLabel,
                 properties: atom.properties,
                 atomTypeName: atom.atomTypeName,
                 maximumBondOrder: atom.maximumBondOrder,
                 bondOrderSum: atom.bondOrderSum,
                 valency: atom.valency,
                 formalNeighbourCount: atom.formalNeighbourCount,
                 hybridization: atom.hybridization)
        }
        participant.bonds = group.bonds
        participant.sgroups = group.sgroups
        return participant
    }

    private static func mutateAllParticipants(in reaction: inout CDKReaction,
                                              transform: (Molecule) -> Molecule) {
        reaction.reactantParticipants = reaction.reactantParticipants.map { participant in
            var participant = participant
            participant.molecule = transform(participant.molecule)
            return participant
        }
        reaction.agentParticipants = reaction.agentParticipants.map { participant in
            var participant = participant
            participant.molecule = transform(participant.molecule)
            return participant
        }
        reaction.productParticipants = reaction.productParticipants.map { participant in
            var participant = participant
            participant.molecule = transform(participant.molecule)
            return participant
        }
    }

    private static func configured(molecule: Molecule, onlyUnset: Bool) -> Molecule {
        var updated = molecule
        updated.atoms = molecule.atoms.map { atom in
            configure(atom: atom, in: molecule, onlyUnset: onlyUnset)
        }
        return updated
    }

    private static func configure(atom: Atom,
                                  in molecule: Molecule,
                                  onlyUnset: Bool) -> Atom {
        var atom = atom
        let bonds = molecule.bonds(forAtom: atom.id)
        let neighbourCount = molecule.neighbors(of: atom.id).count
        let maxOrder = bonds.max(by: { $0.order.rawValue < $1.order.rawValue })?.order
        let orderSum = bonds.reduce(0.0) { $0 + $1.order.valenceContribution }
        let hybridization = inferredHybridization(for: atom,
                                                  neighbourCount: neighbourCount,
                                                  maxBondOrder: maxOrder,
                                                  bondOrderSum: orderSum)
        let typeName = inferredAtomTypeName(for: atom, hybridization: hybridization)
        let valency = inferredValency(for: atom, bondOrderSum: orderSum)

        if !onlyUnset || atom.formalNeighbourCount == nil {
            atom.formalNeighbourCount = neighbourCount
        }
        if !onlyUnset || atom.maximumBondOrder == nil {
            atom.maximumBondOrder = maxOrder
        }
        if !onlyUnset || atom.bondOrderSum == nil {
            atom.bondOrderSum = orderSum
        }
        if !onlyUnset || atom.hybridization == nil {
            atom.hybridization = hybridization
        }
        if !onlyUnset || atom.atomTypeName == nil {
            atom.atomTypeName = typeName
        }
        if !onlyUnset || atom.valency == nil {
            atom.valency = valency
        }
        return atom
    }

    private static func inferredHybridization(for atom: Atom,
                                              neighbourCount: Int,
                                              maxBondOrder: BondOrder?,
                                              bondOrderSum: Double) -> AtomHybridization {
        if atom.element.uppercased() == "H" {
            return .s
        }
        if maxBondOrder == .triple || bondOrderSum >= 2.9 && neighbourCount <= 2 {
            return .sp1
        }
        if atom.aromatic || maxBondOrder == .double || bondOrderSum > Double(neighbourCount) {
            return neighbourCount == 3 ? .planar3 : .sp2
        }
        return .sp3
    }

    private static func inferredAtomTypeName(for atom: Atom,
                                             hybridization: AtomHybridization) -> String {
        "\(atom.element.uppercased()).\(hybridization.rawValue.uppercased())"
    }

    private static func inferredValency(for atom: Atom,
                                        bondOrderSum: Double) -> Int {
        if let override = atom.valenceOverride {
            return override
        }

        let target: Double
        switch atom.element.uppercased() {
        case "C":
            target = atom.aromatic ? 3.0 : 4.0
        case "N":
            target = atom.aromatic ? 3.0 : (atom.charge > 0 ? 4.0 : 3.0)
        case "O":
            target = atom.charge > 0 ? 3.0 : (atom.charge < 0 ? 1.0 : 2.0)
        case "S":
            target = atom.charge > 0 ? 3.0 : 2.0
        case "P":
            target = atom.charge > 0 ? 4.0 : 3.0
        case "B":
            target = atom.charge < 0 ? 4.0 : 3.0
        case "F", "CL", "BR", "I":
            target = 1.0
        case "H":
            target = 1.0
        default:
            target = max(0.0, ceil(bondOrderSum))
        }

        return Int(max(0.0, round(target)))
    }
}

public enum CDKReactionSetManipulator {
    public static func getAtomCount(_ set: CDKReactionSet) -> Int {
        set.flattenedReactions.reduce(0) { $0 + CDKReactionManipulator.getAtomCount($1) }
    }

    public static func getBondCount(_ set: CDKReactionSet) -> Int {
        set.flattenedReactions.reduce(0) { $0 + CDKReactionManipulator.getBondCount($1) }
    }

    public static func getAllMolecules(_ set: CDKReactionSet) -> [Molecule] {
        var unique: [Molecule] = []
        for molecule in set.flattenedReactions.flatMap(CDKReactionManipulator.getAllMolecules) {
            if !unique.contains(molecule) {
                unique.append(molecule)
            }
        }
        return unique
    }

    public static func getAllAtomContainers(_ set: CDKReactionSet) -> [Molecule] {
        getAllMolecules(set)
    }

    public static func getAllIDs(_ set: CDKReactionSet) -> [String] {
        var ids: [String] = []
        if let id = set.id {
            ids.append(id)
        }
        for reaction in set.flattenedReactions {
            ids.append(contentsOf: CDKReactionManipulator.getAllIDs(reaction))
        }
        return ids
    }

    public static func setAtomProperties(_ set: inout CDKReactionSet,
                                         key: String,
                                         value: String) {
        _ = mutateReactions(in: &set) { reaction in
            var updated = reaction
            CDKReactionManipulator.setAtomProperties(&updated, key: key, value: value)
            return (updated, true)
        }
    }

    public static func getAllChemObjects(_ set: CDKReactionSet) -> [CDKReactionChemObject] {
        var objects: [CDKReactionChemObject] = [.reactionSet(set)]
        for reaction in set.flattenedReactions {
            objects.append(contentsOf: CDKReactionManipulator.getAllChemObjects(reaction))
        }
        return objects
    }

    public static func getRelevantReaction(_ set: CDKReactionSet, atom: Atom) -> CDKReaction? {
        set.flattenedReactions.first { reaction in
            CDKReactionManipulator.getRelevantAtomContainer(reaction, atom: atom) != nil
        }
    }

    public static func getRelevantReaction(_ set: CDKReactionSet, bond: Bond) -> CDKReaction? {
        set.flattenedReactions.first { reaction in
            CDKReactionManipulator.getRelevantAtomContainer(reaction, bond: bond) != nil
        }
    }

    public static func getRelevantAtomContainer(_ set: CDKReactionSet, atom: Atom) -> Molecule? {
        set.flattenedReactions.lazy.compactMap { reaction in
            CDKReactionManipulator.getRelevantAtomContainer(reaction, atom: atom)
        }.first
    }

    public static func getRelevantAtomContainer(_ set: CDKReactionSet, bond: Bond) -> Molecule? {
        set.flattenedReactions.lazy.compactMap { reaction in
            CDKReactionManipulator.getRelevantAtomContainer(reaction, bond: bond)
        }.first
    }

    @discardableResult
    public static func removeAtomAndConnectedElectronContainers(_ set: inout CDKReactionSet,
                                                                _ atom: Atom) -> Bool {
        mutateReactions(in: &set) { reaction in
            var updated = reaction
            let changed = CDKReactionManipulator.removeAtomAndConnectedElectronContainers(&updated, atom)
            return (updated, changed)
        }
    }

    @discardableResult
    public static func removeElectronContainer(_ set: inout CDKReactionSet,
                                               _ bond: Bond) -> Bool {
        mutateReactions(in: &set) { reaction in
            var updated = reaction
            let changed = CDKReactionManipulator.removeElectronContainer(&updated, bond)
            return (updated, changed)
        }
    }

    public static func getRelevantReactions(_ set: CDKReactionSet, molecule: Molecule) -> CDKReactionSet {
        let matches = set.flattenedReactions.filter { reaction in
            reaction.reactants.contains(molecule) || reaction.products.contains(molecule)
        }
        return CDKReactionSet(id: set.id, name: set.name, members: matches.map(CDKReactionSetMember.reaction), properties: set.properties)
    }

    public static func getRelevantReactionsAsReactant(_ set: CDKReactionSet, molecule: Molecule) -> CDKReactionSet {
        let matches = set.flattenedReactions.filter { $0.reactants.contains(molecule) }
        return CDKReactionSet(id: set.id, name: set.name, members: matches.map(CDKReactionSetMember.reaction), properties: set.properties)
    }

    public static func getRelevantReactionsAsProduct(_ set: CDKReactionSet, molecule: Molecule) -> CDKReactionSet {
        let matches = set.flattenedReactions.filter { $0.products.contains(molecule) }
        return CDKReactionSet(id: set.id, name: set.name, members: matches.map(CDKReactionSetMember.reaction), properties: set.properties)
    }

    public static func getReactionByAtomContainerID(_ set: CDKReactionSet, id: String) -> CDKReaction? {
        if let productMatch = set.flattenedReactions.first(where: { reaction in
            reaction.products.contains { $0.externalID == id }
        }) {
            return productMatch
        }
        return set.flattenedReactions.first { reaction in
            reaction.reactants.contains { $0.externalID == id }
        }
    }

    public static func getReactionByReactionID(_ set: CDKReactionSet, id: String) -> CDKReaction? {
        set.flattenedReactions.first { $0.id == id }
    }

    @discardableResult
    private static func mutateReactions(in set: inout CDKReactionSet,
                                        update: (CDKReaction) -> (reaction: CDKReaction, changed: Bool)) -> Bool {
        var changed = false
        set.members = set.members.map { member in
            mutate(member: member, update: update, changed: &changed)
        }
        return changed
    }

    private static func mutate(member: CDKReactionSetMember,
                               update: (CDKReaction) -> (reaction: CDKReaction, changed: Bool),
                               changed: inout Bool) -> CDKReactionSetMember {
        switch member {
        case .reaction(let reaction):
            let result = update(reaction)
            changed = changed || result.changed
            return .reaction(result.reaction)
        case .list(let list):
            return .list(mutate(list: list, update: update, changed: &changed))
        case .scheme(let scheme):
            return .scheme(mutate(scheme: scheme, update: update, changed: &changed))
        }
    }

    private static func mutate(list: CDKReactionList,
                               update: (CDKReaction) -> (reaction: CDKReaction, changed: Bool),
                               changed: inout Bool) -> CDKReactionList {
        var updated = list
        updated.entries = list.entries.map { entry in
            switch entry {
            case .reaction(let reaction):
                let result = update(reaction)
                changed = changed || result.changed
                return .reaction(result.reaction)
            case .list(let nested):
                return .list(mutate(list: nested, update: update, changed: &changed))
            case .scheme(let nested):
                return .scheme(mutate(scheme: nested, update: update, changed: &changed))
            }
        }
        return updated
    }

    private static func mutate(scheme: CDKReactionScheme,
                               update: (CDKReaction) -> (reaction: CDKReaction, changed: Bool),
                               changed: inout Bool) -> CDKReactionScheme {
        var updated = scheme
        updated.entries = scheme.entries.map { entry in
            switch entry {
            case .reaction(let reaction):
                let result = update(reaction)
                changed = changed || result.changed
                return .reaction(result.reaction)
            case .list(let nested):
                return .list(mutate(list: nested, update: update, changed: &changed))
            case .scheme(let nested):
                return .scheme(mutate(scheme: nested, update: update, changed: &changed))
            }
        }
        return updated
    }
}
