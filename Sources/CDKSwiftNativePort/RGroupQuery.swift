import Foundation

public struct MoleculeRGroupLogic: Hashable, Codable, Sendable {
    public var occurrence: String
    public var requiredRGroupNumber: Int
    public var restH: Bool

    public init(occurrence: String = CDKRGroupList.defaultOccurrence,
                requiredRGroupNumber: Int = 0,
                restH: Bool = false) {
        self.occurrence = CDKRGroupList.normalizedOccurrence(occurrence)
        self.requiredRGroupNumber = requiredRGroupNumber
        self.restH = restH
    }
}

public struct CDKRGroup: Hashable, Codable, Sendable {
    public var group: Molecule
    public var firstAttachmentPointAtomID: Int?
    public var secondAttachmentPointAtomID: Int?

    public init(group: Molecule,
                firstAttachmentPointAtomID: Int? = nil,
                secondAttachmentPointAtomID: Int? = nil) {
        self.group = group
        self.firstAttachmentPointAtomID = firstAttachmentPointAtomID
        self.secondAttachmentPointAtomID = secondAttachmentPointAtomID
    }
}

public struct CDKRGroupList: Hashable, Codable, Sendable {
    public static let defaultOccurrence = ">0"

    public var rGroupNumber: Int
    public var restH: Bool
    public var occurrence: String
    public var requiredRGroupNumber: Int
    public var rGroups: [CDKRGroup]

    public init(rGroupNumber: Int,
                restH: Bool = false,
                occurrence: String = defaultOccurrence,
                requiredRGroupNumber: Int = 0,
                rGroups: [CDKRGroup] = []) throws {
        guard (1...32).contains(rGroupNumber) else {
            throw ChemError.parseFailed("R-group number must be between 1 and 32.")
        }
        guard CDKRGroupList.isValidOccurrenceSyntax(occurrence) else {
            throw ChemError.parseFailed("Invalid occurrence line: \(occurrence)")
        }
        self.rGroupNumber = rGroupNumber
        self.restH = restH
        self.occurrence = CDKRGroupList.normalizedOccurrence(occurrence)
        self.requiredRGroupNumber = requiredRGroupNumber
        self.rGroups = rGroups
    }

    public static func normalizedOccurrence(_ occurrence: String) -> String {
        let trimmed = occurrence.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { return defaultOccurrence }
        return trimmed.replacingOccurrences(of: " ", with: "")
    }

    public static func isValidOccurrenceSyntax(_ occurrence: String) -> Bool {
        let normalized = normalizedOccurrence(occurrence)
        let conditions = normalized.split(separator: ",").map(String.init)
        guard !conditions.isEmpty else { return false }

        for condition in conditions {
            if condition.allSatisfy(\.isNumber) {
                guard let value = Int(condition), value >= 0 else { return false }
                continue
            }

            if condition.hasPrefix("<") {
                let suffix = String(condition.dropFirst())
                guard let value = Int(suffix), value > 0 else { return false }
                continue
            }

            if condition.hasPrefix(">") {
                let suffix = String(condition.dropFirst())
                guard Int(suffix) != nil else { return false }
                continue
            }

            let parts = condition.split(separator: "-", maxSplits: 1).map(String.init)
            if parts.count == 2,
               let lower = Int(parts[0]),
               let upper = Int(parts[1]),
               lower >= 0,
               upper >= lower {
                continue
            }

            return false
        }

        return true
    }

    public func matchOccurrence(maxAttachments: Int) -> [Int] {
        guard maxAttachments >= 0 else { return [] }

        var matched = Set<Int>()
        for condition in CDKRGroupList.normalizedOccurrence(occurrence).split(separator: ",").map(String.init) {
            if condition.allSatisfy(\.isNumber), let exact = Int(condition) {
                if exact <= maxAttachments {
                    matched.insert(exact)
                }
                continue
            }

            if condition.hasPrefix("<"), let upperBound = Int(condition.dropFirst()) {
                if upperBound > 0 {
                    for value in 0..<min(upperBound, maxAttachments + 1) {
                        matched.insert(value)
                    }
                }
                continue
            }

            if condition.hasPrefix(">"), let lowerBound = Int(condition.dropFirst()) {
                let start = max(0, lowerBound + 1)
                if start <= maxAttachments {
                    for value in start...maxAttachments {
                        matched.insert(value)
                    }
                }
                continue
            }

            let parts = condition.split(separator: "-", maxSplits: 1).map(String.init)
            if parts.count == 2,
               let lower = Int(parts[0]),
               let upper = Int(parts[1]),
               lower <= upper {
                for value in lower...min(upper, maxAttachments) {
                    matched.insert(value)
                }
            }
        }

        return matched.sorted()
    }
}

public struct CDKRGroupQuery: Hashable, Codable, Sendable {
    public var rootStructure: Molecule
    public var rGroupDefinitions: [Int: CDKRGroupList]

    public init(rootStructure: Molecule = Molecule(name: "Root structure"),
                rGroupDefinitions: [Int: CDKRGroupList] = [:]) {
        self.rootStructure = rootStructure
        self.rGroupDefinitions = rGroupDefinitions
    }

    public static func isValidRGroupQueryLabel(_ label: String) -> Bool {
        let trimmed = label.trimmingCharacters(in: .whitespacesAndNewlines).uppercased()
        guard trimmed.hasPrefix("R"), let value = Int(trimmed.dropFirst()) else { return false }
        return (1...32).contains(value)
    }

    public func rGroupQueryAtoms(_ rGroupNumber: Int? = nil) -> [Atom] {
        rootStructure.atoms.filter { atom in
            guard let label = atom.rGroupLabel else { return false }
            if let rGroupNumber {
                return label == rGroupNumber
            }
            return true
        }
    }

    public func areSubstituentsDefined() -> Bool {
        for atom in rGroupQueryAtoms() {
            guard let rLabel = atom.rGroupLabel,
                  let definition = rGroupDefinitions[rLabel],
                  !definition.rGroups.isEmpty else {
                return false
            }
        }
        return true
    }

    public func areRootAtomsDefined() -> Bool {
        for number in rGroupDefinitions.keys {
            guard rGroupQueryAtoms(number).isEmpty == false else { return false }
        }
        return true
    }
}

public enum CDKRGroupQueryManipulator {
    public static func canConvertToQuery(_ molecule: Molecule) -> Bool {
        guard let query = try? toRGroupQuery(molecule) else { return false }
        return !query.rootStructure.atoms.isEmpty && query.areRootAtomsDefined() && query.areSubstituentsDefined()
    }

    public static func setRGroupLabel(_ label: String, on molecule: inout Molecule) {
        for index in molecule.atoms.indices {
            molecule.atoms[index].rGroupMembership = label
        }
    }

    public static func toFlatMolecule(_ query: CDKRGroupQuery) -> Molecule {
        var molecule = query.rootStructure
        molecule.rGroupLogicDefinitions = Dictionary(uniqueKeysWithValues: query.rGroupDefinitions.map { key, value in
            (key, MoleculeRGroupLogic(occurrence: value.occurrence,
                                      requiredRGroupNumber: value.requiredRGroupNumber,
                                      restH: value.restH))
        })

        var nextAtomID = molecule.atoms.map(\.id).max() ?? 0
        var nextBondID = molecule.bonds.map(\.id).max() ?? 0
        var nextComponentGroupID = (molecule.atoms.compactMap(\.componentGroupID).max() ?? 0) + 1

        for rGroupNumber in query.rGroupDefinitions.keys.sorted() {
            guard let definition = query.rGroupDefinitions[rGroupNumber] else { continue }
            let label = "R\(rGroupNumber)"

            for rGroup in definition.rGroups {
                let componentGroupID = nextComponentGroupID
                nextComponentGroupID += 1

                var atomIDMap: [Int: Int] = [:]
                for atom in rGroup.group.atoms {
                    nextAtomID += 1
                    atomIDMap[atom.id] = nextAtomID
                    var copied = atom
                    copied = Atom(id: nextAtomID,
                                  externalID: copied.externalID,
                                  element: copied.element,
                                  position: copied.position,
                                  zPosition: copied.zPosition,
                                  charge: copied.charge,
                                  isotopeMassNumber: copied.isotopeMassNumber,
                                  aromatic: copied.aromatic,
                                  chirality: copied.chirality,
                                  explicitHydrogenCount: copied.explicitHydrogenCount,
                                  queryType: copied.queryType,
                                  atomList: copied.atomList,
                                  atomListIsNegated: copied.atomListIsNegated,
                                  radical: copied.radical,
                                  radicalType: copied.radicalType,
                                  atomValue: copied.atomValue,
                                  rGroupLabel: copied.rGroupLabel,
                                  rGroupMembership: copied.rGroupMembership ?? label,
                                  componentGroupID: componentGroupID,
                                  substitutionCount: copied.substitutionCount,
                                  unsaturated: copied.unsaturated,
                                  ringBondCount: copied.ringBondCount,
                                  attachmentPoint: copied.attachmentPoint,
                                  valenceOverride: copied.valenceOverride,
                                  cxStereoGroup: copied.cxStereoGroup,
                                  ligandOrderingAtomIDs: copied.ligandOrderingAtomIDs,
                                  atomClass: copied.atomClass,
                                  atomMapNumber: copied.atomMapNumber,
                                  aliasLabel: copied.aliasLabel)
                    molecule.atoms.append(copied)
                }

                for bond in rGroup.group.bonds {
                    guard let a1 = atomIDMap[bond.a1], let a2 = atomIDMap[bond.a2] else { continue }
                    nextBondID += 1
                    molecule.bonds.append(Bond(id: nextBondID,
                                               externalID: bond.externalID,
                                               a1: a1,
                                               a2: a2,
                                               order: bond.order,
                                               stereo: bond.stereo,
                                               queryType: bond.queryType,
                                               topology: bond.topology))
                }

                let sourceToFlatIndex = Dictionary(uniqueKeysWithValues: molecule.atoms.enumerated().map { ($1.id, $0) })
                for atom in rGroup.group.atoms {
                    guard let mappedID = atomIDMap[atom.id],
                          let atomIndex = sourceToFlatIndex[mappedID] else { continue }
                    if atom.id == rGroup.firstAttachmentPointAtomID {
                        molecule.atoms[atomIndex].attachmentPoint = (atom.id == rGroup.secondAttachmentPointAtomID) ? 3 : 1
                    } else if atom.id == rGroup.secondAttachmentPointAtomID {
                        molecule.atoms[atomIndex].attachmentPoint = 2
                    }
                }
            }
        }

        return molecule
    }

    public static func toRGroupQuery(_ molecule: Molecule) throws -> CDKRGroupQuery {
        let rootAtomIDs = Set(molecule.atoms.compactMap { atom in
            atom.rGroupMembership == nil ? atom.id : nil
        })
        var rootStructure = extractSubmolecule(from: molecule,
                                               atomIDs: rootAtomIDs,
                                               name: molecule.name,
                                               preserveDataFields: true)
        rootStructure.rGroupLogicDefinitions = [:]

        var definitions: [Int: CDKRGroupList] = [:]
        let groupedByLabel = Dictionary(grouping: molecule.atoms.compactMap { atom -> (String, Atom)? in
            guard let membership = atom.rGroupMembership else { return nil }
            return (membership, atom)
        }, by: { $0.0 })

        for label in groupedByLabel.keys.sorted(by: compareRGroupLabels) {
            guard let rGroupNumber = parseRGroupNumber(label) else {
                throw ChemError.parseFailed("Unsupported R-group label '\(label)' for RGfile export.")
            }
            let atomsForLabel = groupedByLabel[label, default: []].map(\.1)
            let groupedByComponent = Dictionary(grouping: atomsForLabel, by: \.componentGroupID)
            var rGroups: [CDKRGroup] = []

            for componentGroupID in groupedByComponent.keys.compactMap({ $0 }).sorted() {
                guard let atoms = groupedByComponent[componentGroupID], !atoms.isEmpty else { continue }
                let atomIDs = Set(atoms.map(\.id))
                let group = extractSubmolecule(from: molecule,
                                               atomIDs: atomIDs,
                                               name: label,
                                               preserveDataFields: false)
                rGroups.append(CDKRGroup(group: group,
                                         firstAttachmentPointAtomID: attachmentPointAtomID(in: group, matching: [1, 3]),
                                         secondAttachmentPointAtomID: attachmentPointAtomID(in: group, matching: [2, 3])))
            }

            if let ungroupedAtoms = groupedByComponent[nil], !ungroupedAtoms.isEmpty {
                let components = connectedComponents(atomIDs: Set(ungroupedAtoms.map(\.id)), molecule: molecule)
                for component in components {
                    let group = extractSubmolecule(from: molecule,
                                                   atomIDs: component,
                                                   name: label,
                                                   preserveDataFields: false)
                    rGroups.append(CDKRGroup(group: group,
                                             firstAttachmentPointAtomID: attachmentPointAtomID(in: group, matching: [1, 3]),
                                             secondAttachmentPointAtomID: attachmentPointAtomID(in: group, matching: [2, 3])))
                }
            }

            let logic = molecule.rGroupLogicDefinitions[rGroupNumber] ?? MoleculeRGroupLogic()
            definitions[rGroupNumber] = try CDKRGroupList(rGroupNumber: rGroupNumber,
                                                          restH: logic.restH,
                                                          occurrence: logic.occurrence,
                                                          requiredRGroupNumber: logic.requiredRGroupNumber,
                                                          rGroups: rGroups)
        }

        return CDKRGroupQuery(rootStructure: rootStructure, rGroupDefinitions: definitions)
    }

    private static func parseRGroupNumber(_ label: String) -> Int? {
        let trimmed = label.trimmingCharacters(in: .whitespacesAndNewlines).uppercased()
        guard trimmed.hasPrefix("R"), trimmed.count > 1 else { return nil }
        return Int(trimmed.dropFirst())
    }

    private static func compareRGroupLabels(_ lhs: String, _ rhs: String) -> Bool {
        let lhsNumber = parseRGroupNumber(lhs) ?? .max
        let rhsNumber = parseRGroupNumber(rhs) ?? .max
        if lhsNumber != rhsNumber { return lhsNumber < rhsNumber }
        return lhs < rhs
    }

    private static func attachmentPointAtomID(in molecule: Molecule, matching values: Set<Int>) -> Int? {
        molecule.atoms.first(where: { atom in
            guard let attachmentPoint = atom.attachmentPoint else { return false }
            return values.contains(attachmentPoint)
        })?.id
    }

    private static func connectedComponents(atomIDs: Set<Int>, molecule: Molecule) -> [Set<Int>] {
        guard !atomIDs.isEmpty else { return [] }

        var adjacency: [Int: Set<Int>] = [:]
        for atomID in atomIDs {
            adjacency[atomID] = []
        }
        for bond in molecule.bonds where atomIDs.contains(bond.a1) && atomIDs.contains(bond.a2) {
            adjacency[bond.a1, default: []].insert(bond.a2)
            adjacency[bond.a2, default: []].insert(bond.a1)
        }

        var visited = Set<Int>()
        var components: [Set<Int>] = []

        for atomID in atomIDs.sorted() where !visited.contains(atomID) {
            var stack = [atomID]
            var component = Set<Int>()
            visited.insert(atomID)

            while let current = stack.popLast() {
                component.insert(current)
                for neighbor in adjacency[current, default: []].sorted() where !visited.contains(neighbor) {
                    visited.insert(neighbor)
                    stack.append(neighbor)
                }
            }

            components.append(component)
        }

        return components
    }

    private static func extractSubmolecule(from source: Molecule,
                                           atomIDs: Set<Int>,
                                           name: String,
                                           preserveDataFields: Bool) -> Molecule {
        let atoms = source.atoms.filter { atomIDs.contains($0.id) }
        let bondIDs = Set(source.bonds.compactMap { bond in
            atomIDs.contains(bond.a1) && atomIDs.contains(bond.a2) ? bond.id : nil
        })
        let bonds = source.bonds.filter { bondIDs.contains($0.id) }

        var result = Molecule(name: name,
                              externalID: source.externalID,
                              atoms: atoms,
                              bonds: bonds)
        result.sgroups = filteredSgroups(from: source.sgroups,
                                         atomIDs: atomIDs,
                                         bondIDs: bondIDs)
        result.highlightedAtomIDs = source.highlightedAtomIDs.filter { atomIDs.contains($0) }
        result.highlightedBondIDs = source.highlightedBondIDs.filter { bondIDs.contains($0) }
        result.cxState = nil
        if preserveDataFields {
            result.dataFields = source.dataFields
            result.dataFieldOrder = source.dataFieldOrder
            result.rGroupLogicDefinitions = source.rGroupLogicDefinitions
        }
        return result
    }

    private static func filteredSgroups(from sgroups: [MoleculeSgroup],
                                        atomIDs: Set<Int>,
                                        bondIDs: Set<Int>) -> [MoleculeSgroup] {
        let includedIndices = sgroups.indices.filter { index in
            let sgroup = sgroups[index]
            return Set(sgroup.atomIDs).isSubset(of: atomIDs)
                && Set(sgroup.crossingBondIDs).isSubset(of: bondIDs)
                && Set(sgroup.parentAtomIDs).isSubset(of: atomIDs)
        }

        let includedIndexSet = Set(includedIndices)
        let indexMap = Dictionary(uniqueKeysWithValues: includedIndices.enumerated().map { ($1, $0) })

        return includedIndices.map { index in
            var sgroup = sgroups[index]
            sgroup.childGroupIndices = sgroup.childGroupIndices
                .filter { includedIndexSet.contains($0) }
                .compactMap { indexMap[$0] }
            return sgroup
        }
    }
}
