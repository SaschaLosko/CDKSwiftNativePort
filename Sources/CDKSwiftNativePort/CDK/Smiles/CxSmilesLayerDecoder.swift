import Foundation

/// Decodes CX metadata without mutating a molecule. Unsupported layers still
/// fall back to the original title at the public parser boundary.
enum CDKCxSmilesLayerDecoder {
    static func findClosingCxPipe(in input: String, startingAt start: String.Index)
        -> String
        .Index?
    {
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

    static func decode(_ raw: String) throws -> CDKCxSmilesState {
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
                    parseAtomMetadataLayer(
                        String(token.dropFirst(5).dropLast()),
                        into: &state.atomValues,
                        stripLeadingRGroupUnderscore: false)
                } else {
                    parseAtomMetadataLayer(
                        String(token.dropFirst().dropLast()),
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
                state.racemicFragments = try parseIntegerList(
                    String(token.dropFirst(2)), separator: ",")
                continue
            }

            if token.hasPrefix("f:") {
                try parseFragmentGroups(String(token.dropFirst(2)), into: &state)
                continue
            }

            if token.hasPrefix("m:") {
                try parseIndexedIntegerListMap(
                    String(token.dropFirst(2)), separator: ".", into: &state.positionalVariations)
                continue
            }

            if token.hasPrefix("LO:") {
                try parseIndexedIntegerListMap(
                    String(token.dropFirst(3)), separator: ".", into: &state.ligandOrdering)
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
                state.atomHighlights = try parseIntegerList(
                    String(token.dropFirst(3)), separator: ",")
                continue
            }

            if token.hasPrefix("hb:") {
                state.bondHighlights = try parseIntegerList(
                    String(token.dropFirst(3)), separator: ",")
                continue
            }

            if token.hasPrefix("wD:") {
                try parseBondDisplays(
                    String(token.dropFirst(3)),
                    display: .wedgeBegin,
                    into: &state)
                continue
            }

            if token.hasPrefix("wU:") {
                try parseBondDisplays(
                    String(token.dropFirst(3)),
                    display: .wedgedHashBegin,
                    into: &state)
                continue
            }

            if token.hasPrefix("RG:") {
                try parseRGroups(token, into: &state)
                continue
            }

            if token.hasPrefix("LOG={") {
                try parseRGroupLogic(token, into: &state)
                continue
            }

            if token.hasPrefix("LN:") {
                try parseLinkNodes(token, into: &state)
                continue
            }

            if token.hasPrefix("C:") {
                try parseCoordinateBonds(token, into: &state)
                continue
            }

            if token.hasPrefix("c:") || token.hasPrefix("t:") || token.hasPrefix("ctu:")
                || token.hasPrefix("lp:") || token.hasPrefix("H:")
            {
                continue
            }

            throw ChemError.parseFailed("Unsupported or malformed CXSMILES layer: \(token)")
        }

        return state
    }

    private static func parseAtomMetadataLayer(
        _ content: String,
        into destination: inout [Int: String],
        stripLeadingRGroupUnderscore: Bool
    ) {
        let entries = splitAtomMetadataEntries(content)
        for (idx, rawEntry) in entries.enumerated() {
            guard !rawEntry.isEmpty else { continue }
            var value = unescape(rawEntry)
            if stripLeadingRGroupUnderscore,
                value.hasPrefix("_R"),
                value.count > 2
            {
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
                characters[index + 1] == "#"
            {
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

    private static func parseFragmentGroups(
        _ body: String,
        into state: inout CDKCxSmilesState
    ) throws {
        guard !body.isEmpty else { return }
        for group in body.split(separator: ",", omittingEmptySubsequences: false)
        where !group.isEmpty {
            let ids = try parseIntegerList(String(group), separator: ".")
            guard !ids.isEmpty else {
                throw ChemError.parseFailed("Malformed CXSMILES fragment-group layer.")
            }
            state.fragmentGroups.append(ids)
        }
    }

    private static func parseCoordinates(
        _ token: String,
        into state: inout CDKCxSmilesState
    ) throws {
        guard token.hasPrefix("("), token.hasSuffix(")") else {
            throw ChemError.parseFailed("Malformed CXSMILES coordinate layer.")
        }
        let body = String(token.dropFirst().dropLast())
        guard !body.isEmpty else { return }

        for rawEntry in body.split(separator: ";", omittingEmptySubsequences: false)
        where !rawEntry.isEmpty {
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

    private static func parsePolymerSgroup(
        _ token: String,
        into state: inout CDKCxSmilesState
    ) throws {
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

    private static func parseDataSgroup(
        _ token: String,
        into state: inout CDKCxSmilesState
    ) throws {
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

    private static func parseSgroupHierarchy(
        _ token: String,
        into state: inout CDKCxSmilesState
    ) throws {
        let body = String(token.dropFirst(4))
        guard !body.isEmpty else { return }

        for rawEntry in body.split(separator: ",", omittingEmptySubsequences: false)
        where !rawEntry.isEmpty {
            let parts = rawEntry.split(separator: ":", omittingEmptySubsequences: false)
            guard parts.count == 2,
                let parentIndex = Int(parts[0]),
                state.sgroups.indices.contains(parentIndex)
            else {
                throw ChemError.parseFailed("Malformed CXSMILES Sgroup hierarchy layer.")
            }
            let childIndices = try parseIntegerList(String(parts[1]), separator: ".")
            state.sgroups[parentIndex].childIndices.append(contentsOf: childIndices)
        }
    }

    private static func parseIndexedIntegerListMap(
        _ body: String,
        separator: Character,
        into destination: inout [Int: [Int]]
    ) throws {
        guard !body.isEmpty else { return }

        for rawEntry in body.split(separator: ",", omittingEmptySubsequences: false)
        where !rawEntry.isEmpty {
            let parts = rawEntry.split(separator: ":", omittingEmptySubsequences: false)
            guard parts.count == 2,
                let key = Int(parts[0])
            else {
                throw ChemError.parseFailed("Malformed CXSMILES indexed list layer.")
            }
            destination[key] = try parseIntegerList(String(parts[1]), separator: separator)
        }
    }

    private static func parseStereoGroups(
        _ token: String,
        into state: inout CDKCxSmilesState,
        kind: String
    ) throws {
        let number: Int
        let atomListBody: String

        switch kind {
        case "and", "or":
            guard let colon = token.firstIndex(of: ":"),
                colon > token.startIndex
            else {
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

        let encodedGroup = CDKCxSmilesParser.encodeStereoGroup(kind: kind, number: number)
        for atomIndex in try parseIntegerList(atomListBody, separator: ",") {
            state.stereoGroups[atomIndex] = encodedGroup
        }
    }

    private static func parseRadicals(
        _ token: String,
        into state: inout CDKCxSmilesState
    ) throws {
        guard let colon = token.firstIndex(of: ":"),
            colon > token.startIndex
        else {
            throw ChemError.parseFailed("Malformed CXSMILES radical layer.")
        }
        let kindText = String(token[token.index(after: token.startIndex)..<colon])
        guard let rawValue = Int(kindText),
            let radicalType = CxRadicalType(rawValue: rawValue)
        else {
            throw ChemError.parseFailed("Malformed CXSMILES radical type.")
        }
        let atomIndices = try parseIntegerList(
            String(token[token.index(after: colon)...]), separator: ",")
        for atomIndex in atomIndices {
            state.atomRadicals[atomIndex] = radicalType
        }
    }

    private static func parseBondDisplays(
        _ body: String,
        display: CDKCxSmilesState.WedgeDisplay,
        into state: inout CDKCxSmilesState
    ) throws {
        guard !body.isEmpty else { return }
        for rawEntry in body.split(separator: ",", omittingEmptySubsequences: false)
        where !rawEntry.isEmpty {
            let parts = rawEntry.split(separator: ".", omittingEmptySubsequences: false)
            guard parts.count == 2,
                let atomIndex = Int(parts[0]),
                let bondIndex = Int(parts[1])
            else {
                throw ChemError.parseFailed("Malformed CXSMILES wedge layer.")
            }
            state.bondDisplays.append(
                .init(atomIndex: atomIndex, bondIndex: bondIndex, display: display))
        }
    }

    private static func parseCoordinateBonds(_ token: String, into state: inout CDKCxSmilesState)
        throws
    {
        let body = String(token.dropFirst(2))
        guard !body.isEmpty else {
            throw ChemError.parseFailed("Malformed CXSMILES coordinate-bond layer.")
        }

        var parsed: [Int: [Int]] = [:]
        for rawEntry in body.split(separator: ",", omittingEmptySubsequences: false) {
            let fields = rawEntry.split(separator: ".", omittingEmptySubsequences: false)
            guard fields.count == 2,
                let atomIndex = Int(fields[0]), atomIndex >= 0,
                let bondIndex = Int(fields[1]), bondIndex >= 0
            else {
                throw ChemError.parseFailed("Malformed CXSMILES coordinate-bond entry.")
            }
            parsed[atomIndex, default: []].append(bondIndex)
        }
        state.coordinateBonds = parsed
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

                let definition = String(body[definitionStart..<index]).trimmingCharacters(
                    in: .whitespacesAndNewlines)
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

    private static func parseRGroupLogic(_ token: String, into state: inout CDKCxSmilesState) throws {
        guard token.hasSuffix("}") else {
            throw ChemError.parseFailed("Malformed CXSMILES LOG layer.")
        }
        var logic: [Int: MoleculeRGroupLogic] = [:]
        for entry in token.dropFirst(5).dropLast().split(
            separator: ".", omittingEmptySubsequences: false)
        {
            let pair = entry.split(separator: ":", maxSplits: 1, omittingEmptySubsequences: false)
            guard pair.count == 2, pair[0].hasPrefix("_R"), let number = Int(pair[0].dropFirst(2)),
                (0...32).contains(number), logic[number] == nil
            else {
                throw ChemError.parseFailed("Invalid CXSMILES R-group logic label.")
            }
            let fields = pair[1].split(separator: ";", omittingEmptySubsequences: false)
            guard fields.count == 3, fields[1].isEmpty || fields[1] == "H",
                CDKRGroupList.isValidOccurrenceSyntax(String(fields[2]))
            else {
                throw ChemError.parseFailed("Invalid CXSMILES R-group condition.")
            }
            let required: Int
            if fields[0].isEmpty {
                required = 0
            } else if fields[0].hasPrefix("_R"), let value = Int(fields[0].dropFirst(2)),
                (1...32).contains(value)
            {
                required = value
            } else {
                throw ChemError.parseFailed("Invalid CXSMILES R-group dependency.")
            }
            logic[number] = MoleculeRGroupLogic(
                occurrence: String(fields[2]), requiredRGroupNumber: required,
                restH: fields[1] == "H",
                source: .cxsmiles)
        }
        state.rGroupLogic = logic
    }

    private static func parseLinkNodes(_ token: String, into state: inout CDKCxSmilesState) throws {
        for entry in token.dropFirst(3).split(separator: ",", omittingEmptySubsequences: false) {
            let parts = entry.split(separator: ":", omittingEmptySubsequences: false)
            guard parts.count == 2 || parts.count == 3, let atomIndex = Int(parts[0]),
                atomIndex >= 0
            else {
                throw ChemError.parseFailed("Malformed CXSMILES link-node layer.")
            }
            let fields = parts[1].split(separator: ".", omittingEmptySubsequences: false)
            let values = fields.compactMap { Int($0) }
            guard values.count == fields.count, values.count == 2 || values.count == 4,
                values[0] >= 1, values[1] >= values[0]
            else {
                throw ChemError.parseFailed("Malformed CXSMILES link-node bounds.")
            }
            var node = CDKCxSmilesState.LinkNode(
                atomIndex: atomIndex, lowerBound: values[0], upperBound: values[1])
            if parts.count == 3 {
                // Preserve the legacy Swift port's explicit bond-index variant for existing documents.
                let raw = parts[2].split(separator: ".", omittingEmptySubsequences: false)
                node.bondIndices = raw.compactMap { Int($0) }
                guard node.bondIndices.count == raw.count, node.bondIndices.allSatisfy({ $0 >= 0 })
                else {
                    throw ChemError.parseFailed("Invalid link-node bond indices.")
                }
            } else {
                node.outerAtomIndices = Array(values.dropFirst(2))
                guard node.outerAtomIndices!.allSatisfy({ $0 >= 0 }) else {
                    throw ChemError.parseFailed("Invalid link-node outer atoms.")
                }
            }
            state.linkNodes.append(node)
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
                } else if ch == "," && parenDepth == 0 && braceDepth == 0
                    && (!allowsInternalCommas
                        || startsNewLayer(
                            (index + 1) < characters.count
                                ? String(characters[(index + 1)...]).trimmingCharacters(
                                    in: .whitespacesAndNewlines)
                                : ""
                        ))
                {
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
            "RG:", "LOG={", "LN:", "LO:", "f:", "r:", "ha:", "hb:",
            "m:", "Sg:", "SgD:", "SgH:", "c:", "t:", "ctu:",
            "lp:", "C:", "H:", "wD:", "wU:", "&", "o", "a:", "^",
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
            "RG:", "LOG={", "LN:", "LO:", "f:", "r:", "ha:", "hb:",
            "m:", "Sg:", "SgD:", "SgH:", "c:", "t:", "ctu:",
            "lp:", "C:", "H:", "wD:", "wU:", "&", "o", "a:", "^",
        ].contains { token.hasPrefix($0) } || token == "r"
    }

    private static func parseIntegerList(
        _ body: String,
        separator: Character,
        allowEmpty: Bool = false
    ) throws -> [Int] {
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
        while let start = out.range(of: "&#"),
            let semi = out[start.upperBound...].firstIndex(of: ";")
        {
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

extension Array {
    fileprivate subscript(safe index: Int) -> Element? {
        indices.contains(index) ? self[index] : nil
    }
}
