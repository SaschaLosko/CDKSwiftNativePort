import Foundation

public struct CDKCxSmilesState: Equatable {
    public struct LinkNode: Equatable {
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

    public var atomLabels: [Int: String] = [:]
    public var fragmentGroups: [[Int]] = []
    public var racemic: Bool = false
    public var racemicFragments: [Int] = []
    public var linkNodes: [LinkNode] = []
    public var rGroupDefinitions: [String: [String]] = [:]
    public var rGroupOrder: [String] = []

    public init(atomLabels: [Int: String] = [:],
                fragmentGroups: [[Int]] = [],
                racemic: Bool = false,
                racemicFragments: [Int] = [],
                linkNodes: [LinkNode] = [],
                rGroupDefinitions: [String: [String]] = [:],
                rGroupOrder: [String] = []) {
        self.atomLabels = atomLabels
        self.fragmentGroups = fragmentGroups
        self.racemic = racemic
        self.racemicFragments = racemicFragments
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

        guard enabled, let firstPipe = trimmed.firstIndex(of: "|") else {
            return splitTitleField(from: trimmed)
        }

        let afterFirst = trimmed.index(after: firstPipe)
        guard let secondPipe = findClosingCxPipe(in: trimmed, startingAt: afterFirst) else {
            throw ChemError.parseFailed("Unterminated CXSMILES layer (missing closing '|').")
        }

        let core = String(trimmed[..<firstPipe]).trimmingCharacters(in: .whitespacesAndNewlines)
        let cxBody = String(trimmed[afterFirst..<secondPipe])
        let trailing = String(trimmed[trimmed.index(after: secondPipe)...]).trimmingCharacters(in: .whitespacesAndNewlines)

        guard !core.isEmpty else {
            throw ChemError.parseFailed("Missing core SMILES before CXSMILES layer.")
        }
        if trailing.contains("|") {
            throw ChemError.parseFailed("Malformed CXSMILES tail.")
        }

        let state = try parseLayers(cxBody)
        let title = trailing.isEmpty ? nil : trailing
        return SplitResult(coreSmiles: core, title: title, state: state)
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
        applyLinkNodes(to: &molecule, state: state)
        try appendRGroupDefinitions(to: &molecule, state: state, parseDefinition: parseDefinition)
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
                     rGroupMembership: atom.rGroupMembership,
                     componentGroupID: atom.componentGroupID,
                     substitutionCount: atom.substitutionCount,
                     unsaturated: atom.unsaturated,
                     ringBondCount: atom.ringBondCount,
                     attachmentPoint: atom.attachmentPoint,
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
                                          queryType: bond.queryType))
        }

        for sgroup in source.sgroups {
            let atomIDs = sgroup.atomIDs.compactMap { atomMap[$0] }
            guard !atomIDs.isEmpty else { continue }
            let crossingBondIDs = sgroup.crossingBondIDs.compactMap { bondMap[$0] }
            destination.sgroups.append(
                MoleculeSgroup(kind: sgroup.kind,
                               atomIDs: atomIDs,
                               crossingBondIDs: crossingBondIDs,
                               subscriptText: sgroup.subscriptText,
                               superscriptText: sgroup.superscriptText,
                               roundBrackets: sgroup.roundBrackets)
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
                    throw ChemError.parseFailed("Malformed CXSMILES atom-label layer.")
                }
                parseAtomLabels(token, into: &state)
                continue
            }

            if token == "r" {
                state.racemic = true
                continue
            }

            if token.hasPrefix("r:") {
                let body = String(token.dropFirst(2))
                let values = body.split(separator: ",").compactMap { Int($0) }
                if values.isEmpty && !body.isEmpty {
                    throw ChemError.parseFailed("Malformed CXSMILES racemic-fragment layer.")
                }
                state.racemicFragments = values
                continue
            }

            if token.hasPrefix("f:") {
                let body = String(token.dropFirst(2))
                if body.isEmpty { continue }
                let groups = body.split(separator: ",")
                for group in groups {
                    let ids = group.split(separator: ".").compactMap { Int($0) }
                    if ids.isEmpty {
                        throw ChemError.parseFailed("Malformed CXSMILES fragment-group layer.")
                    }
                    state.fragmentGroups.append(ids)
                }
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
        }

        return state
    }

    private static func parseAtomLabels(_ token: String, into state: inout CDKCxSmilesState) {
        let content = String(token.dropFirst().dropLast())
        let entries = content.split(separator: ";", omittingEmptySubsequences: false)
        for (idx, rawEntry) in entries.enumerated() {
            guard !rawEntry.isEmpty else { continue }
            state.atomLabels[idx] = unescape(String(rawEntry))
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
                let allowsInternalCommas = trimmed.hasPrefix("RG:")
                    || trimmed.hasPrefix("LN:")
                    || trimmed.hasPrefix("LO:")
                    || trimmed.hasPrefix("f:")
                    || trimmed.hasPrefix("r:")
                    || trimmed.hasPrefix("ha:")
                    || trimmed.hasPrefix("hb:")
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
        for prefix in ["RG:", "LN:", "LO:", "f:", "r:", "ha:", "hb:"] {
            if text.hasPrefix(prefix) {
                return true
            }
        }
        return text == "r" || text.hasPrefix("r,") || text.hasPrefix("r ")
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
                               atomIDs: groupAtomIDs.sorted(),
                               crossingBondIDs: crossingBondIDs.sorted(),
                               subscriptText: "\(linkNode.lowerBound)-\(linkNode.upperBound)",
                               superscriptText: "ht",
                               roundBrackets: true)
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
