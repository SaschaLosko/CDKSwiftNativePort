import Foundation

public struct CDKCxSmilesState: Equatable, Hashable, Codable, Sendable {
    public enum WedgeDisplay: String, Equatable, Hashable, Codable, Sendable {
        case wedgeBegin
        case wedgedHashBegin
    }

    public struct BondDisplayEntry: Equatable, Hashable, Codable, Sendable {
        public var atomIndex: Int
        public var bondIndex: Int
        public var display: WedgeDisplay

        public init(atomIndex: Int, bondIndex: Int, display: WedgeDisplay) {
            self.atomIndex = atomIndex
            self.bondIndex = bondIndex
            self.display = display
        }
    }

    public struct LinkNode: Equatable, Hashable, Codable, Sendable {
        public var atomIndex: Int
        public var lowerBound: Int
        public var upperBound: Int
        public var bondIndices: [Int]

        public init(atomIndex: Int,
                    lowerBound: Int,
                    upperBound: Int,
                    bondIndices: [Int] = []) {
            self.atomIndex = atomIndex
            self.lowerBound = lowerBound
            self.upperBound = upperBound
            self.bondIndices = bondIndices
        }
    }

    public struct SgroupDefinition: Equatable, Hashable, Codable, Sendable {
        public enum Kind: String, Codable, Hashable, Sendable {
            case polymer
            case data
        }

        public var kind: Kind
        public var keyword: String?
        public var atomIndices: [Int]
        public var bondIndices: [Int]
        public var subscriptText: String?
        public var superscriptText: String?
        public var dataFieldName: String?
        public var dataValue: String?
        public var dataOperator: String?
        public var dataUnit: String?
        public var dataTag: String?
        public var childIndices: [Int]

        public init(kind: Kind,
                    keyword: String? = nil,
                    atomIndices: [Int] = [],
                    bondIndices: [Int] = [],
                    subscriptText: String? = nil,
                    superscriptText: String? = nil,
                    dataFieldName: String? = nil,
                    dataValue: String? = nil,
                    dataOperator: String? = nil,
                    dataUnit: String? = nil,
                    dataTag: String? = nil,
                    childIndices: [Int] = []) {
            self.kind = kind
            self.keyword = keyword
            self.atomIndices = atomIndices
            self.bondIndices = bondIndices
            self.subscriptText = subscriptText
            self.superscriptText = superscriptText
            self.dataFieldName = dataFieldName
            self.dataValue = dataValue
            self.dataOperator = dataOperator
            self.dataUnit = dataUnit
            self.dataTag = dataTag
            self.childIndices = childIndices
        }
    }

    public var atomLabels: [Int: String] = [:]
    public var atomValues: [Int: String] = [:]
    public var atomCoordinates: [CxCoordinate] = []
    public var has3DCoordinates: Bool = false
    public var fragmentGroups: [[Int]] = []
    public var racemic: Bool = false
    public var racemicFragments: [Int] = []
    public var positionalVariations: [Int: [Int]] = [:]
    public var ligandOrdering: [Int: [Int]] = [:]
    public var stereoGroups: [Int: Int] = [:]
    public var atomRadicals: [Int: CxRadicalType] = [:]
    public var sgroups: [SgroupDefinition] = []
    public var atomHighlights: [Int] = []
    public var bondHighlights: [Int] = []
    public var bondDisplays: [BondDisplayEntry] = []
    public var linkNodes: [LinkNode] = []
    public var rGroupDefinitions: [String: [String]] = [:]
    public var rGroupOrder: [String] = []

    public init(atomLabels: [Int: String] = [:],
                atomValues: [Int: String] = [:],
                atomCoordinates: [CxCoordinate] = [],
                has3DCoordinates: Bool = false,
                fragmentGroups: [[Int]] = [],
                racemic: Bool = false,
                racemicFragments: [Int] = [],
                positionalVariations: [Int: [Int]] = [:],
                ligandOrdering: [Int: [Int]] = [:],
                stereoGroups: [Int: Int] = [:],
                atomRadicals: [Int: CxRadicalType] = [:],
                sgroups: [SgroupDefinition] = [],
                atomHighlights: [Int] = [],
                bondHighlights: [Int] = [],
                bondDisplays: [BondDisplayEntry] = [],
                linkNodes: [LinkNode] = [],
                rGroupDefinitions: [String: [String]] = [:],
                rGroupOrder: [String] = []) {
        self.atomLabels = atomLabels
        self.atomValues = atomValues
        self.atomCoordinates = atomCoordinates
        self.has3DCoordinates = has3DCoordinates
        self.fragmentGroups = fragmentGroups
        self.racemic = racemic
        self.racemicFragments = racemicFragments
        self.positionalVariations = positionalVariations
        self.ligandOrdering = ligandOrdering
        self.stereoGroups = stereoGroups
        self.atomRadicals = atomRadicals
        self.sgroups = sgroups
        self.atomHighlights = atomHighlights
        self.bondHighlights = bondHighlights
        self.bondDisplays = bondDisplays
        self.linkNodes = linkNodes
        self.rGroupDefinitions = rGroupDefinitions
        self.rGroupOrder = rGroupOrder
    }
}

public enum CDKCxSmilesParser {
    public struct SplitResult: Equatable {
        public let coreSmiles: String
        public let title: String?
        public let state: CDKCxSmilesState

        public init(coreSmiles: String, title: String?, state: CDKCxSmilesState) {
            self.coreSmiles = coreSmiles
            self.title = title
            self.state = state
        }
    }

    public static func split(_ input: String, enabled: Bool) throws -> SplitResult {
        let trimmed = input.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { throw ChemError.emptyInput }

        guard enabled else {
            return splitTitleField(from: trimmed)
        }

        let basic = splitTitleField(from: trimmed)
        guard let title = basic.title, title.hasPrefix("|") else {
            return basic
        }

        guard let secondPipe = findClosingCxPipe(in: title, startingAt: title.index(after: title.startIndex)) else {
            return basic
        }

        let cxBody = String(title[title.index(after: title.startIndex)..<secondPipe])
        guard let state = try? parseLayers(cxBody) else {
            return basic
        }

        let trailingStart = title.index(after: secondPipe)
        let trailing = String(title[trailingStart...]).trimmingCharacters(in: .whitespacesAndNewlines)
        return SplitResult(coreSmiles: basic.coreSmiles,
                           title: trailing,
                           state: state)
    }

    private static func splitTitleField(from input: String) -> SplitResult {
        let separators = CharacterSet.whitespacesAndNewlines
        guard let firstSeparatorRange = input.rangeOfCharacter(from: separators) else {
            return SplitResult(coreSmiles: input, title: nil, state: CDKCxSmilesState())
        }

        let core = String(input[..<firstSeparatorRange.lowerBound]).trimmingCharacters(in: separators)
        let trailing = String(input[firstSeparatorRange.upperBound...]).trimmingCharacters(in: separators)
        return SplitResult(coreSmiles: core,
                           title: trailing.isEmpty ? nil : trailing,
                           state: CDKCxSmilesState())
    }

    public static func apply(to molecule: inout Molecule,
                             state: CDKCxSmilesState,
                             parseDefinition: (String) throws -> Molecule) throws {
        applyAtomLabels(to: &molecule, state: state)
        applyAtomValues(to: &molecule, state: state)
        applyCoordinates(to: &molecule, state: state)
        applyRadicals(to: &molecule, state: state)
        applyStereoMetadata(to: &molecule, state: state)
        applyLigandOrdering(to: &molecule, state: state)
        applyHighlights(to: &molecule, state: state)
        applySgroups(to: &molecule, state: state)
        applyLinkNodes(to: &molecule, state: state)
        try appendRGroupDefinitions(to: &molecule, state: state, parseDefinition: parseDefinition)
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
        for (atomIndex, radicalType) in state.atomRadicals where molecule.atoms.indices.contains(atomIndex) {
            molecule.atoms[atomIndex].radicalType = radicalType
            molecule.atoms[atomIndex].radical = radicalType.electronCount
        }
    }

    public static func applyStereoMetadata(to molecule: inout Molecule, state: CDKCxSmilesState) {
        guard state.racemic || !state.stereoGroups.isEmpty else { return }
        for idx in molecule.atoms.indices where molecule.atoms[idx].chirality != .none {
            if state.racemic {
                molecule.atoms[idx].cxStereoGroup = 1
            }
        }
        for (atomIndex, group) in state.stereoGroups where molecule.atoms.indices.contains(atomIndex) {
            molecule.atoms[atomIndex].cxStereoGroup = group
        }
    }

    public static func applyLigandOrdering(to molecule: inout Molecule, state: CDKCxSmilesState) {
        guard !state.ligandOrdering.isEmpty else { return }
        for (atomIndex, neighborIndices) in state.ligandOrdering where molecule.atoms.indices.contains(atomIndex) {
            let neighborIDs = neighborIndices.compactMap { neighborIndex in
                molecule.atoms.indices.contains(neighborIndex) ? molecule.atoms[neighborIndex].id : nil
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

    public static func applySgroups(to molecule: inout Molecule, state: CDKCxSmilesState) {
        guard !state.sgroups.isEmpty || !state.positionalVariations.isEmpty else { return }

        var stateToMoleculeIndex: [Int: Int] = [:]
        var pendingChildIndices: [(stateIndex: Int, moleculeIndex: Int, childIndices: [Int])] = []

        for (index, definition) in state.sgroups.enumerated() {
            switch definition.kind {
            case .polymer:
                if let moleculeSgroup = buildPolymerSgroup(from: definition, molecule: molecule, state: state) {
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
                    MoleculeSgroup(kind: .data,
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

        for pending in pendingChildIndices where molecule.sgroups.indices.contains(pending.moleculeIndex) {
            molecule.sgroups[pending.moleculeIndex].childGroupIndices = pending.childIndices.compactMap { childIndex in
                stateToMoleculeIndex[childIndex]
            }
        }

        for (atomIndex, endpoints) in state.positionalVariations where molecule.atoms.indices.contains(atomIndex) {
            let beginAtomID = molecule.atoms[atomIndex].id
            let endpointIDs = endpoints.compactMap { endpointIndex in
                molecule.atoms.indices.contains(endpointIndex) ? molecule.atoms[endpointIndex].id : nil
            }
            let firstBondID = molecule.bonds(forAtom: beginAtomID).first?.id
            molecule.sgroups.append(
                MoleculeSgroup(kind: .extMulticenter,
                               keyword: "m",
                               atomIDs: [beginAtomID] + endpointIDs,
                               crossingBondIDs: firstBondID.map { [$0] } ?? [])
            )
        }
    }

    private static func appendRGroupDefinitions(to molecule: inout Molecule,
                                                state: CDKCxSmilesState,
                                                parseDefinition: (String) throws -> Molecule) throws {
        guard !state.rGroupDefinitions.isEmpty else { return }

        var seenNestedLabels = Set<String>()
        var nextComponentGroupID = (molecule.atoms.compactMap(\.componentGroupID).max() ?? 0) + 1

        for label in orderedRGroupLabels(from: state) {
            guard let definitions = state.rGroupDefinitions[label], !definitions.isEmpty else { continue }

            let rootAtoms = molecule.atoms.filter { rGroupDisplayLabel(for: $0) == label }
            let incomingOrders = try incomingBondOrders(for: rootAtoms, in: molecule, label: label)

            for definitionSmiles in definitions {
                var definition = try parseDefinition(definitionSmiles)
                ensureNonOverlappingDefinitions(in: &definition, seen: &seenNestedLabels)
                try ensureAttachmentPoints(in: &definition, matching: incomingOrders)
                nextComponentGroupID = assignComponentGroupIDs(in: &definition,
                                                               nextComponentGroupID: nextComponentGroupID)

                for idx in definition.atoms.indices where definition.atoms[idx].rGroupMembership == nil {
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

    private static func incomingBondOrders(for atoms: [Atom],
                                           in molecule: Molecule,
                                           label: String) throws -> [BondOrder]? {
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

    private static func assignComponentGroupIDs(in molecule: inout Molecule,
                                                nextComponentGroupID: Int) -> Int {
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

    private static func ensureAttachmentPoints(in molecule: inout Molecule,
                                               matching bondOrders: [BondOrder]?) throws {
        guard let bondOrders else { return }

        var attachmentPointIndices = molecule.atoms.indices.filter { idx in
            molecule.atoms[idx].attachmentPoint != nil && molecule.atoms[idx].rGroupMembership == nil
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
                    Atom(id: nextAtomID,
                         element: "*",
                         position: .zero,
                         explicitHydrogenCount: 0,
                         queryType: .anyAtom,
                         attachmentPoint: offset + 1)
                )
                molecule.bonds.append(Bond(id: nextBondID,
                                           a1: anchorAtomID,
                                           a2: nextAtomID,
                                           order: order))
            }

            if let hCount = molecule.atoms[0].explicitHydrogenCount {
                molecule.atoms[0].explicitHydrogenCount = max(0, hCount - hydrogenAdjustment)
            }

            attachmentPointIndices = molecule.atoms.indices.filter { idx in
                molecule.atoms[idx].attachmentPoint != nil && molecule.atoms[idx].rGroupMembership == nil
            }
        }

        if attachmentPointIndices.count != bondOrders.count {
            throw ChemError.parseFailed("Number of R-group attachment points does not match incoming bond orders.")
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

    private static func ensureNonOverlappingDefinitions(in molecule: inout Molecule, seen: inout Set<String>) {
        for idx in molecule.atoms.indices {
            guard let label = rGroupDisplayLabel(for: molecule.atoms[idx]),
                  label.range(of: #"^R\d*$"#, options: .regularExpression) != nil else {
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

    private static func renameRGroup(_ oldLabel: String, to newLabel: String, in molecule: inout Molecule) {
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
            destination.atoms.append(
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

        for bond in source.bonds {
            guard let a1 = atomMap[bond.a1], let a2 = atomMap[bond.a2] else { continue }
            nextBondID += 1
            bondMap[bond.id] = nextBondID
            destination.bonds.append(Bond(id: nextBondID,
                                          a1: a1,
                                          a2: a2,
                                          order: bond.order,
                                          stereo: bond.stereo,
                                          queryType: bond.queryType,
                                          topology: bond.topology))
        }

        for sgroup in source.sgroups {
            let atomIDs = sgroup.atomIDs.compactMap { atomMap[$0] }
            guard !atomIDs.isEmpty else { continue }
            let crossingBondIDs = sgroup.crossingBondIDs.compactMap { bondMap[$0] }
            destination.sgroups.append(
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
                               parentAtomIDs: sgroup.parentAtomIDs.compactMap { atomMap[$0] },
                               componentNumber: sgroup.componentNumber,
                               expanded: sgroup.expanded,
                               brackets: sgroup.brackets,
                               childGroupIndices: sgroup.childGroupIndices)
            )
        }
    }

    private static func rGroupDisplayLabel(for atom: Atom) -> String? {
        let trimmed = atom.element.trimmingCharacters(in: .whitespacesAndNewlines)
        if trimmed.range(of: #"^R\d*$"#, options: .regularExpression) != nil {
            return trimmed
        }
        if let rGroupLabel = atom.rGroupLabel {
            return "R\(rGroupLabel)"
        }
        return nil
    }

    private static func parseRGroupLabel(_ label: String) -> Int? {
        let trimmed = label.trimmingCharacters(in: .whitespacesAndNewlines)
        guard trimmed.uppercased().hasPrefix("R") else { return nil }
        let digits = trimmed.dropFirst()
        return digits.isEmpty ? nil : Int(digits)
    }

    private static func findClosingCxPipe(in input: String, startingAt start: String.Index) -> String.Index? {
        var cursor = start
        var braceDepth = 0

        while cursor < input.endIndex {
            let ch = input[cursor]
            if ch == "{" {
                braceDepth += 1
            } else if ch == "}" {
                braceDepth = max(0, braceDepth - 1)
            } else if ch == "|" && braceDepth == 0 {
                return cursor
            }
            cursor = input.index(after: cursor)
        }

        return nil
    }

    private static func parseLayers(_ raw: String) throws -> CDKCxSmilesState {
        var state = CDKCxSmilesState()
        guard !raw.isEmpty else { return state }

        let layers = splitTopLevelLayers(raw)
        for layer in layers {
            let token = layer.trimmingCharacters(in: .whitespacesAndNewlines)
            guard !token.isEmpty else { continue }

            if token.hasPrefix("$") {
                guard token.hasSuffix("$"), token.count >= 2 else {
                    throw ChemError.parseFailed("Malformed CXSMILES atom metadata layer.")
                }
                if token.hasPrefix("$_AV:") {
                    parseAtomMetadataLayer(String(token.dropFirst(5).dropLast()),
                                           into: &state.atomValues,
                                           stripLeadingRGroupUnderscore: false)
                } else {
                    parseAtomMetadataLayer(String(token.dropFirst().dropLast()),
                                           into: &state.atomLabels,
                                           stripLeadingRGroupUnderscore: true)
                }
                continue
            }

            if token.hasPrefix("(") {
                try parseCoordinates(token, into: &state)
                continue
            }

            if token == "r" {
                state.racemic = true
                continue
            }

            if token.hasPrefix("r:") {
                state.racemicFragments = try parseIntegerList(String(token.dropFirst(2)), separator: ",")
                continue
            }

            if token.hasPrefix("f:") {
                try parseFragmentGroups(String(token.dropFirst(2)), into: &state)
                continue
            }

            if token.hasPrefix("m:") {
                try parseIndexedIntegerListMap(String(token.dropFirst(2)), separator: ".", into: &state.positionalVariations)
                continue
            }

            if token.hasPrefix("LO:") {
                try parseIndexedIntegerListMap(String(token.dropFirst(3)), separator: ".", into: &state.ligandOrdering)
                continue
            }

            if token.hasPrefix("SgD:") {
                try parseDataSgroup(token, into: &state)
                continue
            }

            if token.hasPrefix("Sg:") {
                try parsePolymerSgroup(token, into: &state)
                continue
            }

            if token.hasPrefix("SgH:") {
                try parseSgroupHierarchy(token, into: &state)
                continue
            }

            if token.hasPrefix("^") {
                try parseRadicals(token, into: &state)
                continue
            }

            if token.hasPrefix("&") {
                try parseStereoGroups(token, into: &state, kind: "and")
                continue
            }

            if token.hasPrefix("o") {
                try parseStereoGroups(token, into: &state, kind: "or")
                continue
            }

            if token.hasPrefix("a:") {
                try parseStereoGroups(token, into: &state, kind: "abs")
                continue
            }

            if token.hasPrefix("ha:") {
                state.atomHighlights = try parseIntegerList(String(token.dropFirst(3)), separator: ",")
                continue
            }

            if token.hasPrefix("hb:") {
                state.bondHighlights = try parseIntegerList(String(token.dropFirst(3)), separator: ",")
                continue
            }

            if token.hasPrefix("wD:") {
                try parseBondDisplays(String(token.dropFirst(3)),
                                      display: .wedgeBegin,
                                      into: &state)
                continue
            }

            if token.hasPrefix("wU:") {
                try parseBondDisplays(String(token.dropFirst(3)),
                                      display: .wedgedHashBegin,
                                      into: &state)
                continue
            }

            if token.hasPrefix("RG:") {
                try parseRGroups(token, into: &state)
                continue
            }

            if token.hasPrefix("LN:") {
                try parseLinkNodes(token, into: &state)
                continue
            }

            if token.hasPrefix("c:") ||
                token.hasPrefix("t:") ||
                token.hasPrefix("ctu:") ||
                token.hasPrefix("lp:") ||
                token.hasPrefix("C:") ||
                token.hasPrefix("H:") {
                continue
            }

            throw ChemError.parseFailed("Unsupported or malformed CXSMILES layer: \(token)")
        }

        return state
    }

    private static func parseAtomMetadataLayer(_ content: String,
                                               into destination: inout [Int: String],
                                               stripLeadingRGroupUnderscore: Bool) {
        let entries = splitAtomMetadataEntries(content)
        for (idx, rawEntry) in entries.enumerated() {
            guard !rawEntry.isEmpty else { continue }
            var value = unescape(rawEntry)
            if stripLeadingRGroupUnderscore,
               value.hasPrefix("_R"),
               value.count > 2 {
                value.removeFirst()
            }
            destination[idx] = value
        }
    }

    private static func splitAtomMetadataEntries(_ content: String) -> [String] {
        let characters = Array(content)
        var entries: [String] = []
        var current = ""
        var index = 0

        while index < characters.count {
            let ch = characters[index]
            if ch == "&",
               index + 1 < characters.count,
               characters[index + 1] == "#" {
                current.append(ch)
                index += 1
                current.append(characters[index])
                index += 1
                while index < characters.count {
                    let entityChar = characters[index]
                    current.append(entityChar)
                    index += 1
                    if entityChar == ";" {
                        break
                    }
                }
                continue
            }

            if ch == ";" {
                entries.append(current)
                current.removeAll(keepingCapacity: true)
            } else {
                current.append(ch)
            }
            index += 1
        }

        entries.append(current)
        return entries
    }

    private static func parseFragmentGroups(_ body: String,
                                            into state: inout CDKCxSmilesState) throws {
        guard !body.isEmpty else { return }
        for group in body.split(separator: ",", omittingEmptySubsequences: false) where !group.isEmpty {
            let ids = try parseIntegerList(String(group), separator: ".")
            guard !ids.isEmpty else {
                throw ChemError.parseFailed("Malformed CXSMILES fragment-group layer.")
            }
            state.fragmentGroups.append(ids)
        }
    }

    private static func parseCoordinates(_ token: String,
                                         into state: inout CDKCxSmilesState) throws {
        guard token.hasPrefix("("), token.hasSuffix(")") else {
            throw ChemError.parseFailed("Malformed CXSMILES coordinate layer.")
        }
        let body = String(token.dropFirst().dropLast())
        guard !body.isEmpty else { return }

        for rawEntry in body.split(separator: ";", omittingEmptySubsequences: false) where !rawEntry.isEmpty {
            let parts = String(rawEntry).split(separator: ",", omittingEmptySubsequences: false)
            guard parts.count >= 3 else {
                throw ChemError.parseFailed("Malformed CXSMILES coordinate tuple.")
            }
            let x = try parseCxDouble(String(parts[0]))
            let y = try parseCxDouble(String(parts[1]))
            let z = try parseCxDouble(String(parts[2]))
            state.has3DCoordinates = state.has3DCoordinates || abs(z) > 0.000_000_1
            state.atomCoordinates.append(CxCoordinate(x: x, y: y, z: z))
        }
    }

    private static func parsePolymerSgroup(_ token: String,
                                           into state: inout CDKCxSmilesState) throws {
        let body = String(token.dropFirst(3))
        let parts = body.split(separator: ":", omittingEmptySubsequences: false).map(String.init)
        guard parts.count >= 2 else {
            throw ChemError.parseFailed("Malformed CXSMILES polymer Sgroup layer.")
        }

        let keyword = parts[0]
        let atomIndices = try parseIntegerList(parts[1], separator: ",")
        let subscriptText: String
        if parts.count >= 3 {
            let candidate = unescape(parts[2])
            subscriptText = candidate.isEmpty ? keyword : candidate
        } else {
            subscriptText = keyword
        }
        let superscriptText = parts.count >= 4 ? unescape(parts[3]) : ""

        state.sgroups.append(
            CDKCxSmilesState.SgroupDefinition(
                kind: .polymer,
                keyword: keyword,
                atomIndices: atomIndices,
                subscriptText: subscriptText,
                superscriptText: superscriptText
            )
        )
    }

    private static func parseDataSgroup(_ token: String,
                                        into state: inout CDKCxSmilesState) throws {
        let body = String(token.dropFirst(4))
        let parts = body.split(separator: ":", omittingEmptySubsequences: false).map(String.init)
        guard parts.count >= 3 else {
            throw ChemError.parseFailed("Malformed CXSMILES data Sgroup layer.")
        }

        let atomIndices = try parseIntegerList(parts[0], separator: ",", allowEmpty: true)
        var fieldName = unescape(parts[1])
        if fieldName.hasPrefix("cdk.") {
            fieldName = fieldName.replacingOccurrences(of: "cdk.", with: "cdk:")
        }

        state.sgroups.append(
            CDKCxSmilesState.SgroupDefinition(
                kind: .data,
                atomIndices: atomIndices,
                dataFieldName: fieldName,
                dataValue: unescape(parts[safe: 2] ?? ""),
                dataOperator: unescape(parts[safe: 3] ?? ""),
                dataUnit: unescape(parts[safe: 4] ?? ""),
                dataTag: unescape(parts[safe: 5] ?? "")
            )
        )
    }

    private static func parseSgroupHierarchy(_ token: String,
                                             into state: inout CDKCxSmilesState) throws {
        let body = String(token.dropFirst(4))
        guard !body.isEmpty else { return }

        for rawEntry in body.split(separator: ",", omittingEmptySubsequences: false) where !rawEntry.isEmpty {
            let parts = rawEntry.split(separator: ":", omittingEmptySubsequences: false)
            guard parts.count == 2,
                  let parentIndex = Int(parts[0]),
                  state.sgroups.indices.contains(parentIndex) else {
                throw ChemError.parseFailed("Malformed CXSMILES Sgroup hierarchy layer.")
            }
            let childIndices = try parseIntegerList(String(parts[1]), separator: ".")
            state.sgroups[parentIndex].childIndices.append(contentsOf: childIndices)
        }
    }

    private static func parseIndexedIntegerListMap(_ body: String,
                                                   separator: Character,
                                                   into destination: inout [Int: [Int]]) throws {
        guard !body.isEmpty else { return }

        for rawEntry in body.split(separator: ",", omittingEmptySubsequences: false) where !rawEntry.isEmpty {
            let parts = rawEntry.split(separator: ":", omittingEmptySubsequences: false)
            guard parts.count == 2,
                  let key = Int(parts[0]) else {
                throw ChemError.parseFailed("Malformed CXSMILES indexed list layer.")
            }
            destination[key] = try parseIntegerList(String(parts[1]), separator: separator)
        }
    }

    private static func parseStereoGroups(_ token: String,
                                          into state: inout CDKCxSmilesState,
                                          kind: String) throws {
        let number: Int
        let atomListBody: String

        switch kind {
        case "and", "or":
            guard let colon = token.firstIndex(of: ":"),
                  colon > token.startIndex else {
                throw ChemError.parseFailed("Malformed CXSMILES stereo-group layer.")
            }
            let numberText = String(token[token.index(after: token.startIndex)..<colon])
            guard let parsed = Int(numberText) else {
                throw ChemError.parseFailed("Malformed CXSMILES stereo-group identifier.")
            }
            number = parsed
            atomListBody = String(token[token.index(after: colon)...])
        case "abs":
            guard token.hasPrefix("a:") else {
                throw ChemError.parseFailed("Malformed CXSMILES absolute stereo-group layer.")
            }
            number = 0
            atomListBody = String(token.dropFirst(2))
        default:
            throw ChemError.parseFailed("Unsupported stereo-group kind.")
        }

        let encodedGroup = encodeStereoGroup(kind: kind, number: number)
        for atomIndex in try parseIntegerList(atomListBody, separator: ",") {
            state.stereoGroups[atomIndex] = encodedGroup
        }
    }

    private static func parseRadicals(_ token: String,
                                      into state: inout CDKCxSmilesState) throws {
        guard let colon = token.firstIndex(of: ":"),
              colon > token.startIndex else {
            throw ChemError.parseFailed("Malformed CXSMILES radical layer.")
        }
        let kindText = String(token[token.index(after: token.startIndex)..<colon])
        guard let rawValue = Int(kindText),
              let radicalType = CxRadicalType(rawValue: rawValue) else {
            throw ChemError.parseFailed("Malformed CXSMILES radical type.")
        }
        let atomIndices = try parseIntegerList(String(token[token.index(after: colon)...]), separator: ",")
        for atomIndex in atomIndices {
            state.atomRadicals[atomIndex] = radicalType
        }
    }

    private static func parseBondDisplays(_ body: String,
                                          display: CDKCxSmilesState.WedgeDisplay,
                                          into state: inout CDKCxSmilesState) throws {
        guard !body.isEmpty else { return }
        for rawEntry in body.split(separator: ",", omittingEmptySubsequences: false) where !rawEntry.isEmpty {
            let parts = rawEntry.split(separator: ".", omittingEmptySubsequences: false)
            guard parts.count == 2,
                  let atomIndex = Int(parts[0]),
                  let bondIndex = Int(parts[1]) else {
                throw ChemError.parseFailed("Malformed CXSMILES wedge layer.")
            }
            state.bondDisplays.append(.init(atomIndex: atomIndex, bondIndex: bondIndex, display: display))
        }
    }

    private static func parseRGroups(_ token: String, into state: inout CDKCxSmilesState) throws {
        let body = Array(token.dropFirst(3))
        var index = 0
        var currentLabel: String?

        while index < body.count {
            let ch = body[index]

            if ch == "_" {
                index += 1
                let start = index
                while index < body.count && (body[index].isLetter || body[index].isNumber) {
                    index += 1
                }
                guard start < index else {
                    throw ChemError.parseFailed("Malformed CXSMILES R-group label.")
                }
                guard index < body.count, body[index] == "=" else {
                    throw ChemError.parseFailed("Malformed CXSMILES R-group assignment.")
                }

                currentLabel = String(body[start..<index])
                if let currentLabel, !state.rGroupOrder.contains(currentLabel) {
                    state.rGroupOrder.append(currentLabel)
                }
                index += 1
                continue
            }

            if ch == "{" {
                guard let currentLabel else {
                    throw ChemError.parseFailed("Malformed CXSMILES R-group definition.")
                }

                index += 1
                let definitionStart = index
                var depth = 1

                while index < body.count, depth > 0 {
                    if body[index] == "{" {
                        depth += 1
                    } else if body[index] == "}" {
                        depth -= 1
                    }

                    if depth == 0 { break }
                    index += 1
                }

                guard depth == 0 else {
                    throw ChemError.parseFailed("Unterminated CXSMILES R-group definition.")
                }

                let definition = String(body[definitionStart..<index]).trimmingCharacters(in: .whitespacesAndNewlines)
                state.rGroupDefinitions[currentLabel, default: []].append(definition)
                index += 1
                if index < body.count, body[index] == "," {
                    index += 1
                }
                continue
            }

            if ch == "," || ch.isWhitespace {
                index += 1
                continue
            }

            throw ChemError.parseFailed("Malformed CXSMILES R-group layer.")
        }
    }

    private static func parseLinkNodes(_ token: String, into state: inout CDKCxSmilesState) throws {
        let body = String(token.dropFirst(3))
        guard !body.isEmpty else { return }

        for entry in body.split(separator: ",") {
            let parts = entry.split(separator: ":")
            guard parts.count == 2 || parts.count == 3,
                  let atomIndex = Int(parts[0]) else {
                throw ChemError.parseFailed("Malformed CXSMILES link-node layer.")
            }

            let values = parts[1].split(separator: ".").compactMap { Int($0) }
            guard values.count >= 2 else {
                throw ChemError.parseFailed("Malformed CXSMILES link-node bounds.")
            }

            let bondIndices = parts.count == 3
                ? parts[2].split(separator: ".").compactMap { Int($0) }
                : Array(values.dropFirst(2))
            state.linkNodes.append(
                CDKCxSmilesState.LinkNode(atomIndex: atomIndex,
                                          lowerBound: values[0],
                                          upperBound: values[1],
                                          bondIndices: bondIndices)
            )
        }
    }

    private static func splitTopLevelLayers(_ input: String) -> [String] {
        let characters = Array(input)
        var out: [String] = []
        var current = ""
        var inAtomLabel = false
        var parenDepth = 0
        var braceDepth = 0

        for (index, ch) in characters.enumerated() {
            if ch == "$" {
                inAtomLabel.toggle()
                current.append(ch)
                continue
            }

            if !inAtomLabel {
                let trimmed = current.trimmingCharacters(in: .whitespacesAndNewlines)
                let allowsInternalCommas = allowsInternalCommas(in: trimmed)
                if ch == "(" {
                    parenDepth += 1
                } else if ch == ")" {
                    parenDepth = max(0, parenDepth - 1)
                } else if ch == "{" {
                    braceDepth += 1
                } else if ch == "}" {
                    braceDepth = max(0, braceDepth - 1)
                } else if ch == "," &&
                            parenDepth == 0 &&
                            braceDepth == 0 &&
                            (!allowsInternalCommas || startsNewLayer(
                                (index + 1) < characters.count
                                    ? String(characters[(index + 1)...]).trimmingCharacters(in: .whitespacesAndNewlines)
                                    : ""
                            )) {
                    out.append(current)
                    current.removeAll(keepingCapacity: true)
                    continue
                }
            }

            current.append(ch)
        }

        out.append(current)
        return out
    }

    private static func startsNewLayer(_ text: String) -> Bool {
        guard !text.isEmpty else { return false }
        if text.hasPrefix("$") { return true }
        for prefix in [
            "RG:", "LN:", "LO:", "f:", "r:", "ha:", "hb:",
            "m:", "Sg:", "SgD:", "SgH:", "c:", "t:", "ctu:",
            "lp:", "C:", "H:", "wD:", "wU:", "&", "o", "a:", "^"
        ] {
            if text.hasPrefix(prefix) {
                return true
            }
        }
        return text == "r" || text.hasPrefix("r,") || text.hasPrefix("r ")
    }

    private static func allowsInternalCommas(in token: String) -> Bool {
        guard !token.isEmpty else { return false }
        return [
            "RG:", "LN:", "LO:", "f:", "r:", "ha:", "hb:",
            "m:", "Sg:", "SgD:", "SgH:", "c:", "t:", "ctu:",
            "lp:", "C:", "H:", "wD:", "wU:", "&", "o", "a:", "^"
        ].contains { token.hasPrefix($0) } || token == "r"
    }

    private static func buildPolymerSgroup(from definition: CDKCxSmilesState.SgroupDefinition,
                                           molecule: Molecule,
                                           state: CDKCxSmilesState) -> MoleculeSgroup? {
        var groupAtomIDs = Set(definition.atomIndices.compactMap { atomIndex in
            molecule.atoms.indices.contains(atomIndex) ? molecule.atoms[atomIndex].id : nil
        })
        guard !groupAtomIDs.isEmpty else { return nil }

        var crossingBondIDs = Set(definition.bondIndices.compactMap { bondIndex in
            molecule.bonds.indices.contains(bondIndex) ? molecule.bonds[bondIndex].id : nil
        })

        if groupAtomIDs.count == 1 && crossingBondIDs.isEmpty,
           let atomID = groupAtomIDs.first {
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
                       let endpoints = state.positionalVariations[otherIndex] {
                        let endpointIDs = Set(endpoints.compactMap { endpointIndex in
                            molecule.atoms.indices.contains(endpointIndex) ? molecule.atoms[endpointIndex].id : nil
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
        return MoleculeSgroup(kind: kind,
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
        case "mon", "mer", "co", "xl", "mod", "mix", "f", "any", "gen", "c", "grf", "alt", "ran", "blk":
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
            var crossingBondIDs = Set(linkNode.bondIndices.compactMap { rawIndex in
                molecule.bonds.indices.contains(rawIndex) ? molecule.bonds[rawIndex].id : nil
            })

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
                    for bond in molecule.bonds(forAtom: atomID) where !crossingBondIDs.contains(bond.id) {
                        let otherAtomID = bond.a1 == atomID ? bond.a2 : bond.a1
                        if seen.insert(otherAtomID).inserted {
                            queue.append(otherAtomID)
                            groupAtomIDs.insert(otherAtomID)
                        }
                    }
                }
            }

            molecule.sgroups.append(
                MoleculeSgroup(kind: .structureRepeatUnit,
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

        let ringEdges = Set(molecule.simpleCycles(maxSize: 12).flatMap { ring in
            guard let first = ring.first else { return [EdgeKey]() }
            let shifted = Array(ring.dropFirst()) + [first]
            return zip(ring, shifted).map { EdgeKey($0.0, $0.1) }
        })
        return ringEdges.contains(EdgeKey(bond.a1, bond.a2))
    }

    static func encodeStereoGroup(kind: String, number: Int) -> Int {
        let kindCode: Int
        switch kind {
        case "abs":
            kindCode = 1
        case "or":
            kindCode = 2
        case "and":
            kindCode = 3
        default:
            kindCode = 0
        }
        return (kindCode << 16) | max(0, number)
    }

    static func decodeStereoGroup(_ encoded: Int) -> (kind: String, number: Int)? {
        let kindCode = encoded >> 16
        let number = encoded & 0xFFFF
        switch kindCode {
        case 1:
            return ("abs", 0)
        case 2:
            return ("or", number)
        case 3:
            return ("and", number)
        default:
            return nil
        }
    }

    private static func parseIntegerList(_ body: String,
                                         separator: Character,
                                         allowEmpty: Bool = false) throws -> [Int] {
        if body.isEmpty {
            return allowEmpty ? [] : []
        }
        let parts = body.split(separator: separator, omittingEmptySubsequences: allowEmpty == false)
        if !allowEmpty && parts.isEmpty {
            return []
        }

        var values: [Int] = []
        values.reserveCapacity(parts.count)
        for part in parts where !part.isEmpty {
            guard let value = Int(part) else {
                throw ChemError.parseFailed("Malformed CXSMILES integer list.")
            }
            values.append(value)
        }

        if !allowEmpty && values.isEmpty && !body.isEmpty {
            throw ChemError.parseFailed("Malformed CXSMILES integer list.")
        }
        return values
    }

    private static func parseCxDouble(_ text: String) throws -> Double {
        let trimmed = text.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { return 0 }
        if trimmed == "." || trimmed == "-" || trimmed == "+" {
            throw ChemError.parseFailed("Malformed CXSMILES coordinate value.")
        }
        guard let value = Double(trimmed) else {
            throw ChemError.parseFailed("Malformed CXSMILES coordinate value.")
        }
        return value
    }

    // Minimal numeric HTML entity decoding used by CXSMILES escaped labels.
    private static func unescape(_ input: String) -> String {
        var out = input
        while let start = out.range(of: "&#"), let semi = out[start.upperBound...].firstIndex(of: ";") {
            let numeric = out[start.upperBound..<semi]
            if let value = Int(numeric), let scalar = UnicodeScalar(value) {
                out.replaceSubrange(start.lowerBound...semi, with: String(Character(scalar)))
            } else {
                break
            }
        }
        return out
    }
}

private extension Array {
    subscript(safe index: Int) -> Element? {
        indices.contains(index) ? self[index] : nil
    }
}
