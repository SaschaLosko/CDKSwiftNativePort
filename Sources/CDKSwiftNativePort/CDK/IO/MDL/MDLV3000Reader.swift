import Foundation
#if canImport(CoreGraphics)
import CoreGraphics
#endif

public struct CDKMDLV3000Scaffold: Hashable, Codable {
    public struct AtomRecord: Hashable, Codable {
        public let index: Int
        public let symbol: String
        public let x: Double?
        public let y: Double?
        public let z: Double?
        public let raw: String
    }

    public struct BondRecord: Hashable, Codable {
        public let index: Int
        public let type: String
        public let a1: Int
        public let a2: Int
        public let raw: String
    }

    public let title: String
    public let countsLine: String
    public let atomCount: Int
    public let bondCount: Int
    public let sgroupCount: Int
    public let chiralFlag: Int
    public let atomRecords: [AtomRecord]
    public let bondRecords: [BondRecord]
    public let collectionLines: [String]
    public let sgroupLines: [String]
}

public enum CDKMDLV3000Reader {
    private static let markushDefinitionSubtype = "MARKUSHRGROUP"

    private struct BondParseResult {
        let bond: Bond
        let positionalVariationSgroup: MoleculeSgroup?
    }

    private struct RGroupBlock {
        let number: Int
        let memberCTABs: [[String]]
    }

    public static func read(text: String) throws -> Molecule {
        let normalized = text
            .replacingOccurrences(of: "\r\n", with: "\n")
            .replacingOccurrences(of: "\r", with: "\n")
        return try read(lines: normalized.components(separatedBy: "\n"))
    }

    public static func read(lines: [String]) throws -> Molecule {
        let trimmed = dropTrailingEmptyLines(lines)
        let scaffold = try parseScaffold(lines: trimmed)
        let rgroupBlocks = try parseRGroupBlocks(lines: trimmed)

        var atoms: [Atom] = []
        atoms.reserveCapacity(scaffold.atomRecords.count)
        var seenAtomIDs = Set<Int>()

        for record in scaffold.atomRecords {
            let atom = try parseAtom(rawLine: record.raw)
            guard seenAtomIDs.insert(atom.id).inserted else {
                throw ChemError.parseFailed("Duplicate V3000 atom index \(atom.id).")
            }
            atoms.append(atom)
        }

        if scaffold.atomCount > 0, atoms.count != scaffold.atomCount {
            throw ChemError.parseFailed(
                "V3000 atom count mismatch (counts=\(scaffold.atomCount), parsed=\(atoms.count))."
            )
        }

        let atomIDs = Set(atoms.map(\.id))

        var bonds: [Bond] = []
        bonds.reserveCapacity(scaffold.bondRecords.count)
        var seenBondIDs = Set<Int>()
        var positionalVariationGroups: [MoleculeSgroup] = []

        for record in scaffold.bondRecords {
            let parsed = try parseBond(rawLine: record.raw, validAtomIDs: atomIDs)
            let bond = parsed.bond
            guard seenBondIDs.insert(bond.id).inserted else {
                throw ChemError.parseFailed("Duplicate V3000 bond index \(bond.id).")
            }
            bonds.append(bond)
            if let sgroup = parsed.positionalVariationSgroup {
                positionalVariationGroups.append(sgroup)
            }
        }

        if scaffold.bondCount > 0, bonds.count != scaffold.bondCount {
            throw ChemError.parseFailed(
                "V3000 bond count mismatch (counts=\(scaffold.bondCount), parsed=\(bonds.count))."
            )
        }

        let title = scaffold.title.trimmingCharacters(in: .whitespacesAndNewlines)
        var molecule = Molecule(
            name: title.isEmpty ? "Molecule" : title,
            atoms: atoms.sorted(by: { $0.id < $1.id }),
            bonds: bonds.sorted(by: { $0.id < $1.id })
        )

        molecule.sgroups.append(contentsOf: positionalVariationGroups)
        applyCollectionSemantics(scaffold.collectionLines, chiralFlag: scaffold.chiralFlag, to: &molecule)
        applySGroupSemantics(scaffold.sgroupLines, to: &molecule)
        try applyRGroupBlocks(rgroupBlocks, to: &molecule)
        CDKSDFDataFieldParser.applyParsedFields(from: trimmed, to: &molecule)

        if let box = molecule.boundingBox(), box.width <= 0.0001 && box.height <= 0.0001 {
            molecule = Depiction2DGenerator.generate(for: molecule)
        }

        return molecule
    }

    public static func parseScaffold(text: String) throws -> CDKMDLV3000Scaffold {
        let normalized = text
            .replacingOccurrences(of: "\r\n", with: "\n")
            .replacingOccurrences(of: "\r", with: "\n")
        return try parseScaffold(lines: normalized.components(separatedBy: "\n"))
    }

    public static func parseScaffold(lines: [String]) throws -> CDKMDLV3000Scaffold {
        let trimmed = dropTrailingEmptyLines(lines)
        if trimmed.isEmpty {
            throw ChemError.parseFailed("Expected a header line, but found nothing.")
        }
        guard trimmed.count >= 4 else {
            throw ChemError.parseFailed("Molfile too short.")
        }

        let countsLine = trimmed[3]
        guard countsLine.uppercased().contains("V3000") else {
            throw ChemError.parseFailed("Not an MDL V3000 molfile counts line.")
        }

        let unfolded = unfoldV30Continuations(trimmed)

        guard let beginCTAB = unfolded.firstIndex(where: { $0.trimmingCharacters(in: .whitespaces) == "M  V30 BEGIN CTAB" }),
              let endCTAB = unfolded.firstIndex(where: { $0.trimmingCharacters(in: .whitespaces) == "M  V30 END CTAB" }),
              beginCTAB < endCTAB else {
            throw ChemError.parseFailed("V3000 CTAB block missing.")
        }

        enum Section {
            case none
            case atom
            case bond
            case collection
            case sgroup
        }

        var section: Section = .none
        var atomCount = -1
        var bondCount = -1
        var sgroupCount = -1
        var chiralFlag = 0
        var atomRecords: [CDKMDLV3000Scaffold.AtomRecord] = []
        var bondRecords: [CDKMDLV3000Scaffold.BondRecord] = []
        var collectionLines: [String] = []
        var sgroupLines: [String] = []

        for index in (beginCTAB + 1)..<endCTAB {
            let rawLine = unfolded[index]
            let line = rawLine.trimmingCharacters(in: .whitespaces)

            if line == "M  V30 BEGIN ATOM" {
                section = .atom
                continue
            }
            if line == "M  V30 END ATOM" {
                section = .none
                continue
            }
            if line == "M  V30 BEGIN BOND" {
                section = .bond
                continue
            }
            if line == "M  V30 END BOND" {
                section = .none
                continue
            }
            if line == "M  V30 BEGIN COLLECTION" {
                section = .collection
                continue
            }
            if line == "M  V30 END COLLECTION" {
                section = .none
                continue
            }
            if line == "M  V30 BEGIN SGROUP" {
                section = .sgroup
                continue
            }
            if line == "M  V30 END SGROUP" {
                section = .none
                continue
            }

            if line.hasPrefix("M  V30 COUNTS") {
                let fields = tokenizeV30(line)
                if fields.count >= 3 {
                    atomCount = Int(fields[1]) ?? atomCount
                    bondCount = Int(fields[2]) ?? bondCount
                    if fields.count >= 4 {
                        sgroupCount = Int(fields[3]) ?? sgroupCount
                    }
                    if fields.count >= 6 {
                        chiralFlag = Int(fields[5]) ?? chiralFlag
                    }
                }
                continue
            }

            guard line.hasPrefix("M  V30 ") else { continue }

            switch section {
            case .atom:
                let parsed = try parseAtomRecord(rawLine)
                atomRecords.append(parsed)
            case .bond:
                let parsed = try parseBondRecord(rawLine)
                bondRecords.append(parsed)
            case .collection:
                collectionLines.append(rawLine)
            case .sgroup:
                sgroupLines.append(rawLine)
            case .none:
                break
            }
        }

        if atomCount < 0 { atomCount = atomRecords.count }
        if bondCount < 0 { bondCount = bondRecords.count }
        if sgroupCount < 0 { sgroupCount = sgroupLines.count }

        let title = unfolded[0].trimmingCharacters(in: .whitespacesAndNewlines)

        return CDKMDLV3000Scaffold(
            title: title,
            countsLine: countsLine,
            atomCount: atomCount,
            bondCount: bondCount,
            sgroupCount: sgroupCount,
            chiralFlag: chiralFlag,
            atomRecords: atomRecords,
            bondRecords: bondRecords,
            collectionLines: collectionLines,
            sgroupLines: sgroupLines
        )
    }

    private static func parseRGroupBlocks(lines: [String]) throws -> [RGroupBlock] {
        let unfolded = unfoldV30Continuations(lines)
        guard let beginCTAB = unfolded.firstIndex(where: { $0.trimmingCharacters(in: .whitespaces) == "M  V30 BEGIN CTAB" }),
              let endCTAB = unfolded.firstIndex(where: { $0.trimmingCharacters(in: .whitespaces) == "M  V30 END CTAB" }),
              beginCTAB < endCTAB else {
            return []
        }

        var blocks: [RGroupBlock] = []
        var index = endCTAB + 1

        while index < unfolded.count {
            let line = unfolded[index].trimmingCharacters(in: .whitespaces)
            guard line.hasPrefix("M  V30 BEGIN RGROUP") else {
                index += 1
                continue
            }

            let fields = tokenizeV30(line)
            guard fields.count >= 3, let number = Int(fields[2]) else {
                throw ChemError.parseFailed("Invalid V3000 RGROUP header: \(unfolded[index])")
            }

            index += 1
            var memberCTABs: [[String]] = []

            while index < unfolded.count {
                let current = unfolded[index].trimmingCharacters(in: .whitespaces)
                if current == "M  V30 END RGROUP" {
                    break
                }

                if current == "M  V30 BEGIN CTAB" {
                    var ctabLines: [String] = [unfolded[index]]
                    index += 1

                    while index < unfolded.count {
                        ctabLines.append(unfolded[index])
                        if unfolded[index].trimmingCharacters(in: .whitespaces) == "M  V30 END CTAB" {
                            break
                        }
                        index += 1
                    }

                    guard ctabLines.last?.trimmingCharacters(in: .whitespaces) == "M  V30 END CTAB" else {
                        throw ChemError.parseFailed("Unterminated V3000 RGROUP member CTAB.")
                    }

                    memberCTABs.append(ctabLines)
                }

                index += 1
            }

            guard index < unfolded.count,
                  unfolded[index].trimmingCharacters(in: .whitespaces) == "M  V30 END RGROUP" else {
                throw ChemError.parseFailed("Unterminated V3000 RGROUP block.")
            }

            blocks.append(RGroupBlock(number: number, memberCTABs: memberCTABs))
            index += 1
        }

        return blocks
    }

    private static func parseAtom(rawLine: String) throws -> Atom {
        let fields = tokenizeV30(rawLine)
        guard fields.count >= 6,
              let id = Int(fields[0]),
              let x = Double(fields[2]),
              let y = Double(fields[3]) else {
            throw ChemError.parseFailed("Invalid V3000 atom line: \(rawLine)")
        }

        let symbolToken = fields[1]
        let (attributes, bareAttributes) = parseAttributes(fields.dropFirst(6))

        let atomAtomMapping = (Int(fields[5]) ?? 0) > 0 ? Int(fields[5]) : nil
        var element = normalizeElementToken(symbolToken)
        var aromatic = isAromaticSymbol(symbolToken)
        var queryType = queryType(for: symbolToken)
        var atomList: [String]? = nil
        var atomListNegated = false
        var isotopeMassNumber = parseIntAttribute("MASS", in: attributes)

        if symbolToken.uppercased() == "D" {
            element = "H"
            isotopeMassNumber = isotopeMassNumber ?? 2
            aromatic = false
        } else if symbolToken.uppercased() == "T" {
            element = "H"
            isotopeMassNumber = isotopeMassNumber ?? 3
            aromatic = false
        }

        if let symbolList = parseSymbolAtomList(symbolToken) {
            atomList = symbolList.entries
            atomListNegated = symbolList.negated
            queryType = .anyAtom
            element = "L"
            aromatic = false
        }

        if let attrListValue = attributes["ATOMLIST"],
           let parsed = parseAtomListValue(attrListValue) {
            atomList = parsed.entries
            atomListNegated = parsed.negated || parseBoolAttribute("NOT", attributes: attributes, bare: bareAttributes)
            queryType = .anyAtom
            element = "L"
            aromatic = false
        }

        if symbolToken.uppercased() == "ANY" {
            element = "*"
            queryType = .anyAtom
            aromatic = false
        }

        var atom = Atom(
            id: id,
            element: element,
            position: CGPoint(x: x, y: y),
            zPosition: Double(fields[4]),
            charge: parseIntAttribute("CHG", in: attributes) ?? 0,
            isotopeMassNumber: isotopeMassNumber,
            aromatic: aromatic,
            chirality: .none,
            explicitHydrogenCount: parseHydrogenCount(attributes),
            queryType: queryType,
            atomList: atomList,
            atomListIsNegated: atomListNegated,
            radical: parseIntAttribute("RAD", in: attributes),
            rGroupLabel: parseRGroup(attributes),
            substitutionCount: parseIntAttribute("SUBST", in: attributes) ?? parseIntAttribute("SUB", in: attributes),
            unsaturated: parseIntAttribute("UNSAT", in: attributes) ?? parseIntAttribute("UNS", in: attributes),
            ringBondCount: parseIntAttribute("RBCNT", in: attributes) ?? parseIntAttribute("RBC", in: attributes),
            attachmentPoint: parseIntAttribute("ATTCHPT", in: attributes) ?? parseIntAttribute("APO", in: attributes),
            valenceOverride: parseIntAttribute("VAL", in: attributes),
            atomClass: atomAtomMapping,
            atomMapNumber: atomAtomMapping
        )

        if let cfg = parseIntAttribute("CFG", in: attributes) {
            if cfg == 1 {
                atom.chirality = .clockwise
            } else if cfg == 2 {
                atom.chirality = .anticlockwise
            }
        }

        return atom
    }

    private static func parseBond(rawLine: String, validAtomIDs: Set<Int>) throws -> BondParseResult {
        let fields = tokenizeV30(rawLine)
        guard fields.count >= 4,
              let id = Int(fields[0]),
              let a1 = Int(fields[2]),
              let a2 = Int(fields[3]) else {
            throw ChemError.parseFailed("Invalid V3000 bond line: \(rawLine)")
        }

        guard validAtomIDs.contains(a1), validAtomIDs.contains(a2) else {
            throw ChemError.parseFailed("V3000 bond references unknown atom index.")
        }

        let (order, queryType) = try bondFromTypeToken(fields[1], rawLine: rawLine)
        let (attributes, _) = parseAttributes(fields.dropFirst(4))

        let stereo = stereoFromCFG(attributes["CFG"])
        let topology = topologyFromToken(attributes["TOPO"])
        let positionalVariation = positionalVariationSgroup(bondID: id, beginAtomID: a1, attributes: attributes)

        return BondParseResult(
            bond: Bond(
                id: id,
                a1: a1,
                a2: a2,
                order: order,
                stereo: stereo,
                queryType: queryType,
                topology: topology
            ),
            positionalVariationSgroup: positionalVariation
        )
    }

    private static func parseAtomRecord(_ rawLine: String) throws -> CDKMDLV3000Scaffold.AtomRecord {
        let fields = tokenizeV30(rawLine)
        guard fields.count >= 5,
              let index = Int(fields[0]) else {
            throw ChemError.parseFailed("Invalid V3000 atom record: \(rawLine)")
        }

        let symbol = fields[1]
        let x = Double(fields[2])
        let y = Double(fields[3])
        let z = Double(fields[4])

        return CDKMDLV3000Scaffold.AtomRecord(
            index: index,
            symbol: symbol,
            x: x,
            y: y,
            z: z,
            raw: rawLine
        )
    }

    private static func parseBondRecord(_ rawLine: String) throws -> CDKMDLV3000Scaffold.BondRecord {
        let fields = tokenizeV30(rawLine)
        guard fields.count >= 4,
              let index = Int(fields[0]),
              let a1 = Int(fields[2]),
              let a2 = Int(fields[3]) else {
            throw ChemError.parseFailed("Invalid V3000 bond record: \(rawLine)")
        }

        return CDKMDLV3000Scaffold.BondRecord(
            index: index,
            type: fields[1],
            a1: a1,
            a2: a2,
            raw: rawLine
        )
    }

    private static func parseAttributes(_ tokens: ArraySlice<String>) -> ([String: String], Set<String>) {
        var attributes: [String: String] = [:]
        var bare: Set<String> = []

        for token in tokens {
            if let eq = token.firstIndex(of: "=") {
                let key = String(token[..<eq]).uppercased()
                let value = String(token[token.index(after: eq)...])
                attributes[key] = value
            } else {
                bare.insert(token.uppercased())
            }
        }

        return (attributes, bare)
    }

    private static func parseIntAttribute(_ key: String, in attributes: [String: String]) -> Int? {
        guard let raw = attributes[key.uppercased()]?.trimmingCharacters(in: .whitespacesAndNewlines), !raw.isEmpty else {
            return nil
        }
        return Int(raw)
    }

    private static func parseBoolAttribute(_ key: String,
                                           attributes: [String: String],
                                           bare: Set<String>) -> Bool {
        let upper = key.uppercased()
        if bare.contains(upper) {
            return true
        }

        guard let raw = attributes[upper]?.trimmingCharacters(in: .whitespacesAndNewlines).uppercased() else {
            return false
        }
        return raw == "1" || raw == "T" || raw == "TRUE" || raw == "YES"
    }

    private static func parseHydrogenCount(_ attributes: [String: String]) -> Int? {
        guard let hCount = parseIntAttribute("HCOUNT", in: attributes) else {
            return nil
        }
        return hCount >= 0 ? hCount : nil
    }

    private static func parseRGroup(_ attributes: [String: String]) -> Int? {
        if let single = parseIntAttribute("RGROUP", in: attributes) {
            return single
        }

        guard let listValue = attributes["RGROUPS"] else {
            return nil
        }

        let values = parseIntegerListValue(listValue)
        return values.first
    }

    private static func parseIntegerListValue(_ value: String) -> [Int] {
        var payload = value.trimmingCharacters(in: .whitespacesAndNewlines)
        if payload.hasPrefix("(") && payload.hasSuffix(")") {
            payload.removeFirst()
            payload.removeLast()
        }
        payload = payload.replacingOccurrences(of: ",", with: " ")

        let fields = payload.split(whereSeparator: \.isWhitespace).map(String.init)
        guard !fields.isEmpty else { return [] }

        var start = 0
        var expected: Int? = nil
        if let n = Int(fields[0]), n >= 0 {
            expected = n
            start = 1
        }

        var values: [Int] = fields.dropFirst(start).compactMap(Int.init)
        if let expected, values.count > expected {
            values = Array(values.prefix(expected))
        }
        return values
    }

    private static func parseAtomListValue(_ value: String) -> (entries: [String], negated: Bool)? {
        var payload = value.trimmingCharacters(in: .whitespacesAndNewlines)
        var negated = false

        if payload.uppercased().hasPrefix("NOT") {
            negated = true
            payload = String(payload.dropFirst(3)).trimmingCharacters(in: .whitespacesAndNewlines)
        }

        if payload.hasPrefix("(") && payload.hasSuffix(")") {
            payload.removeFirst()
            payload.removeLast()
        }
        if payload.hasPrefix("[") && payload.hasSuffix("]") {
            payload.removeFirst()
            payload.removeLast()
        }

        payload = payload.replacingOccurrences(of: ",", with: " ")

        let fields = payload.split(whereSeparator: \.isWhitespace).map(String.init)
        guard !fields.isEmpty else { return nil }

        var start = 0
        var expectedCount: Int? = nil
        if let n = Int(fields[0]), n >= 0 {
            expectedCount = n
            start = 1
        }

        var entries = fields.dropFirst(start).map(normalizeElementToken)
        if let expectedCount, entries.count > expectedCount {
            entries = Array(entries.prefix(expectedCount))
        }
        guard !entries.isEmpty else { return nil }

        return (entries: entries, negated: negated)
    }

    private static func parseSymbolAtomList(_ symbolToken: String) -> (entries: [String], negated: Bool)? {
        var token = symbolToken.trimmingCharacters(in: .whitespacesAndNewlines)
        var negated = false

        if token.uppercased().hasPrefix("NOT[") {
            negated = true
            token = String(token.dropFirst(3))
        }

        guard token.hasPrefix("[") && token.hasSuffix("]") else {
            return nil
        }

        token.removeFirst()
        token.removeLast()
        token = token.replacingOccurrences(of: ",", with: " ")

        let entries = token
            .split(whereSeparator: \.isWhitespace)
            .map(String.init)
            .map(normalizeElementToken)

        guard !entries.isEmpty else { return nil }
        return (entries: entries, negated: negated)
    }

    private static func queryType(for symbolToken: String) -> AtomQueryType? {
        switch symbolToken.uppercased() {
        case "*", "ANY":
            return .anyAtom
        case "A":
            return .anyNonHydrogen
        case "Q":
            return .anyHetero
        default:
            return nil
        }
    }

    private static func isAromaticSymbol(_ symbolToken: String) -> Bool {
        let token = symbolToken.trimmingCharacters(in: .whitespacesAndNewlines)
        guard token.allSatisfy({ $0.isLetter }) else { return false }
        let lower = token.lowercased()
        return token == lower && aromaticElementTokens.contains(lower)
    }

    private static func normalizeElementToken(_ token: String) -> String {
        let trimmed = token.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { return "C" }

        if trimmed == "*" {
            return "*"
        }

        if trimmed.uppercased() == "ANY" {
            return "*"
        }

        if trimmed.allSatisfy({ $0.isLetter }) {
            guard let first = trimmed.first else { return trimmed }
            let head = String(first).uppercased()
            let tail = String(trimmed.dropFirst()).lowercased()
            return head + tail
        }

        return trimmed
    }

    private static func bondFromTypeToken(_ token: String, rawLine: String) throws -> (BondOrder, BondQueryType?) {
        switch token.uppercased() {
        case "1", "S":
            return (.single, nil)
        case "2", "D":
            return (.double, nil)
        case "3", "T":
            return (.triple, nil)
        case "4", "AR":
            return (.aromatic, nil)
        case "5", "SD":
            return (.single, .singleOrDouble)
        case "6", "SA":
            return (.single, .singleOrAromatic)
        case "7", "DA":
            return (.double, .doubleOrAromatic)
        case "8", "ANY", "Q":
            return (.single, .any)
        case "9", "10":
            throw ChemError.parseFailed("Error while parsing bond type: Unsupported bond type: \(token), line='\(rawLine)'")
        default:
            if let numeric = Int(token) {
                switch numeric {
                case 1: return (.single, nil)
                case 2: return (.double, nil)
                case 3: return (.triple, nil)
                case 4: return (.aromatic, nil)
                case 5: return (.single, .singleOrDouble)
                case 6: return (.single, .singleOrAromatic)
                case 7: return (.double, .doubleOrAromatic)
                case 8: return (.single, .any)
                case 9, 10:
                    throw ChemError.parseFailed("Error while parsing bond type: Unsupported bond type: \(token), line='\(rawLine)'")
                default:
                    throw ChemError.parseFailed("Error while parsing bond type: Invalid bond type: \(token), line='\(rawLine)'")
                }
            }
            throw ChemError.parseFailed("Error while parsing bond type: Invalid bond type: \(token), line='\(rawLine)'")
        }
    }

    private static func topologyFromToken(_ raw: String?) -> BondTopology? {
        guard let raw = raw?.trimmingCharacters(in: .whitespacesAndNewlines),
              let value = Int(raw) else {
            return nil
        }
        switch value {
        case 1:
            return .ring
        case 2:
            return .chain
        default:
            return nil
        }
    }

    private static func positionalVariationSgroup(bondID: Int,
                                                  beginAtomID: Int,
                                                  attributes: [String: String]) -> MoleculeSgroup? {
        guard let attach = attributes["ATTACH"]?.uppercased(),
              attach == "ANY" else {
            return nil
        }
        let endpointIDs = parseIntegerListValue(attributes["ENDPTS"] ?? "")
        guard !endpointIDs.isEmpty else { return nil }
        return MoleculeSgroup(kind: .extMulticenter,
                              keyword: "m",
                              atomIDs: [beginAtomID] + endpointIDs,
                              crossingBondIDs: [bondID])
    }

    private static func stereoFromCFG(_ rawCFG: String?) -> BondStereo {
        guard let raw = rawCFG?.trimmingCharacters(in: .whitespacesAndNewlines),
              let cfg = Int(raw) else {
            return .none
        }

        switch cfg {
        case 1: return .up
        case 2: return .either
        case 3: return .down
        default: return .none
        }
    }

    private static func applyCollectionSemantics(_ lines: [String],
                                                 chiralFlag: Int,
                                                 to molecule: inout Molecule) {
        if chiralFlag == 1 {
            for idx in molecule.atoms.indices where molecule.atoms[idx].chirality != .none {
                molecule.atoms[idx].cxStereoGroup = CDKCxSmilesParser.encodeStereoGroup(kind: "abs", number: 0)
            }
        } else {
            for idx in molecule.atoms.indices where molecule.atoms[idx].chirality != .none {
                molecule.atoms[idx].cxStereoGroup = CDKCxSmilesParser.encodeStereoGroup(kind: "or", number: 1)
            }
        }

        guard !lines.isEmpty else { return }

        var atomIndexByID: [Int: Int] = [:]
        for idx in molecule.atoms.indices {
            atomIndexByID[molecule.atoms[idx].id] = idx
        }
        var bondIndexByID: [Int: Int] = [:]
        for idx in molecule.bonds.indices {
            bondIndexByID[molecule.bonds[idx].id] = idx
        }

        for rawLine in lines {
            let fields = tokenizeV30(rawLine)
            guard !fields.isEmpty else { continue }

            let collectionType = fields[0].uppercased()
            let (attributes, _) = parseAttributes(fields.dropFirst())
            let atomIDs = parseIntegerListValue(attributes["ATOMS"] ?? "")
            let bondIDs = parseIntegerListValue(attributes["BONDS"] ?? "")

            if collectionType.contains("HILITE") {
                molecule.highlightedAtomIDs.append(contentsOf: atomIDs.filter { atomIndexByID[$0] != nil })
                molecule.highlightedBondIDs.append(contentsOf: bondIDs.filter { bondIndexByID[$0] != nil })
                continue
            }

            if collectionType.contains("STE") {
                for bondID in bondIDs {
                    guard let idx = bondIndexByID[bondID], molecule.bonds[idx].order == .double else { continue }
                    molecule.bonds[idx].stereo = .either
                }

                let stereoGroup: Int?
                if collectionType.contains("ABS") {
                    stereoGroup = CDKCxSmilesParser.encodeStereoGroup(kind: "abs", number: 0)
                } else if collectionType.contains("RAC"), let number = trailingCollectionNumber(in: collectionType) {
                    stereoGroup = CDKCxSmilesParser.encodeStereoGroup(kind: "or", number: number)
                } else if collectionType.contains("REL"), let number = trailingCollectionNumber(in: collectionType) {
                    stereoGroup = CDKCxSmilesParser.encodeStereoGroup(kind: "and", number: number)
                } else {
                    stereoGroup = nil
                }

                for atomID in atomIDs {
                    guard let idx = atomIndexByID[atomID] else { continue }
                    if molecule.atoms[idx].chirality == .none && molecule.neighbors(of: atomID).count >= 3 {
                        molecule.atoms[idx].chirality = .clockwise
                    }
                    if let stereoGroup {
                        molecule.atoms[idx].cxStereoGroup = stereoGroup
                    }
                }
            }
        }
    }

    private static func applyRGroupBlocks(_ blocks: [RGroupBlock], to molecule: inout Molecule) throws {
        guard !blocks.isEmpty else { return }

        var nextComponentGroupID = (molecule.atoms.compactMap(\.componentGroupID).max() ?? 0) + 1

        for block in blocks.sorted(by: { $0.number < $1.number }) {
            let label = "R\(block.number)"

            for memberCTAB in block.memberCTABs {
                var member = try read(lines: wrappedMolfileLines(for: memberCTAB, title: label))
                let componentGroupID = nextComponentGroupID
                nextComponentGroupID += 1

                for idx in member.atoms.indices {
                    if member.atoms[idx].rGroupMembership == nil {
                        member.atoms[idx].rGroupMembership = label
                    }
                    member.atoms[idx].componentGroupID = componentGroupID
                }

                merge(member, into: &molecule)
            }
        }
    }

    private static func applySGroupSemantics(_ lines: [String], to molecule: inout Molecule) {
        guard !lines.isEmpty else { return }

        var atomIndexByID: [Int: Int] = [:]
        for idx in molecule.atoms.indices {
            atomIndexByID[molecule.atoms[idx].id] = idx
        }
        var bondIndexByID: [Int: Int] = [:]
        for idx in molecule.bonds.indices {
            bondIndexByID[molecule.bonds[idx].id] = idx
        }
        var nextMarkushComponentGroupID = (molecule.atoms.compactMap(\.componentGroupID).max() ?? 0) + 1

        var parsedSgroups: [(sourceIndex: Int, parentIndex: Int?, sgroup: MoleculeSgroup)] = []

        for rawLine in lines {
            let fields = tokenizeV30(rawLine)
            guard fields.count >= 3, let sourceIndex = Int(fields[0]) else { continue }

            let sgroupType = fields[1].uppercased()
            let (attributeArrays, _) = parseRepeatedAttributes(fields.dropFirst(3))
            let attributes = firstAttributes(from: attributeArrays)
            let atomIDs = parseIntegerListValue(attributes["ATOMS"] ?? "")
            let parentAtomIDs = parseIntegerListValue(attributes["PATOMS"] ?? "")
            let crossingBondKey = sgroupType == "DAT" ? "CBONDS" : "XBONDS"
            let crossingBondIDs = parseIntegerListValue(attributes[crossingBondKey] ?? "")
            let parentIndex = parseIntAttribute("PARENT", in: attributes)
            let subtype = parseQuotedOrBareString(attributes["SUBTYPE"] ?? "")?.uppercased()

            if let attchpt = parseIntAttribute("ATTCHPT", in: attributes) {
                for atomID in atomIDs.isEmpty ? parentAtomIDs : atomIDs {
                    guard let idx = atomIndexByID[atomID] else { continue }
                    molecule.atoms[idx].attachmentPoint = attchpt
                }
            }

            if sgroupType == "SUP",
               subtype == markushDefinitionSubtype,
               let rawLabel = attributes["LABEL"],
               let label = parseQuotedOrBareString(rawLabel), !label.isEmpty {
                let componentGroupID = parseIntAttribute("COMPNO", in: attributes) ?? nextMarkushComponentGroupID
                if parseIntAttribute("COMPNO", in: attributes) == nil {
                    nextMarkushComponentGroupID += 1
                }
                for atomID in atomIDs {
                    guard let idx = atomIndexByID[atomID] else { continue }
                    molecule.atoms[idx].rGroupMembership = label
                    molecule.atoms[idx].componentGroupID = componentGroupID
                }
                continue
            }

            if sgroupType == "SUP",
               atomIDs.count == 1,
               let rawLabel = attributes["LABEL"],
               let label = parseQuotedOrBareString(rawLabel), !label.isEmpty,
               let anchorAtomID = atomIDs.first, let anchorIndex = atomIndexByID[anchorAtomID] {
                molecule.atoms[anchorIndex].element = label
                molecule.atoms[anchorIndex].aromatic = false
                if label.uppercased().hasPrefix("R"),
                   let rValue = Int(label.dropFirst()) {
                    molecule.atoms[anchorIndex].rGroupLabel = rValue
                }
            }

            if sgroupType == "DAT" {
                let fieldName = parseQuotedOrBareString(attributes["FIELDNAME"] ?? "") ?? "DATA"
                let fieldData = parseQuotedOrBareString(attributes["FIELDDATA"] ?? attributes["DATA"] ?? "")
                molecule.ensureDataField(named: fieldName)
                if let fieldData, !fieldData.isEmpty {
                    molecule.appendDataFieldValue(fieldData, named: fieldName)
                }
            }

            let subscriptText: String?
            if sgroupType == "MUL" {
                subscriptText = parseQuotedOrBareString(attributes["MULT"] ?? "")
            } else {
                subscriptText = parseQuotedOrBareString(attributes["LABEL"] ?? "")
            }

            parsedSgroups.append((
                sourceIndex: sourceIndex,
                parentIndex: parentIndex,
                sgroup: MoleculeSgroup(
                    kind: sgroupKind(for: sgroupType),
                    keyword: sgroupType,
                    atomIDs: atomIDs,
                    crossingBondIDs: crossingBondIDs.filter { bondIndexByID[$0] != nil },
                    subscriptText: subscriptText,
                    superscriptText: nil,
                    roundBrackets: attributes["BRKTYP"]?.uppercased() == "PAREN",
                    connectivity: parseQuotedOrBareString(attributes["CONNECT"] ?? ""),
                    dataFieldName: parseQuotedOrBareString(attributes["FIELDNAME"] ?? ""),
                    dataValue: parseQuotedOrBareString(attributes["FIELDDATA"] ?? attributes["DATA"] ?? ""),
                    dataOperator: nil,
                    dataUnit: nil,
                    dataTag: nil,
                    subtype: parseQuotedOrBareString(attributes["SUBTYPE"] ?? ""),
                    parentAtomIDs: parentAtomIDs,
                    componentNumber: parseIntAttribute("COMPNO", in: attributes),
                    expanded: attributes["ESTATE"]?.uppercased() == "E",
                    brackets: parseBrackets(attributeArrays["BRKXYZ"] ?? []),
                    childGroupIndices: []
                )
            ))
        }

        var sourceToMoleculeIndex: [Int: Int] = [:]
        for parsed in parsedSgroups {
            sourceToMoleculeIndex[parsed.sourceIndex] = molecule.sgroups.count
            molecule.sgroups.append(parsed.sgroup)
        }

        for parsed in parsedSgroups {
            guard let parentSourceIndex = parsed.parentIndex,
                  let parentMoleculeIndex = sourceToMoleculeIndex[parentSourceIndex],
                  let childMoleculeIndex = sourceToMoleculeIndex[parsed.sourceIndex],
                  molecule.sgroups.indices.contains(parentMoleculeIndex) else {
                continue
            }
            molecule.sgroups[parentMoleculeIndex].childGroupIndices.append(childMoleculeIndex)
        }
    }

    private static func sgroupKind(for sgroupType: String) -> MoleculeSgroup.Kind {
        switch sgroupType {
        case "SRU":
            return .structureRepeatUnit
        case "DAT":
            return .data
        case "COP", "CRO", "GRA", "MOD", "MON", "MER":
            return .polymer
        default:
            return .generic
        }
    }

    private static func parseBrackets(_ values: [String]) -> [MoleculeSgroupBracket] {
        values.compactMap { value in
            let numbers = parseDoubleListValue(value)
            guard numbers.count >= 5 else { return nil }
            return MoleculeSgroupBracket(
                firstPoint: CGPoint(x: numbers[0], y: numbers[1]),
                secondPoint: CGPoint(x: numbers[3], y: numbers[4])
            )
        }
    }

    private static func parseDoubleListValue(_ value: String) -> [Double] {
        var payload = value.trimmingCharacters(in: .whitespacesAndNewlines)
        if payload.hasPrefix("(") && payload.hasSuffix(")") {
            payload.removeFirst()
            payload.removeLast()
        }
        payload = payload.replacingOccurrences(of: ",", with: " ")
        let fields = payload.split(whereSeparator: \.isWhitespace)
        guard !fields.isEmpty else { return [] }
        let start = Int(fields[0]) != nil ? 1 : 0
        return fields.dropFirst(start).compactMap { Double(String($0)) }
    }

    private static func trailingCollectionNumber(in collectionType: String) -> Int? {
        let digits = collectionType.reversed().prefix { $0.isNumber }.reversed()
        guard !digits.isEmpty else { return nil }
        return Int(String(digits))
    }

    private static func parseRepeatedAttributes(_ tokens: ArraySlice<String>) -> ([String: [String]], Set<String>) {
        var attributes: [String: [String]] = [:]
        var bare: Set<String> = []

        for token in tokens {
            if let eq = token.firstIndex(of: "=") {
                let key = String(token[..<eq]).uppercased()
                let value = String(token[token.index(after: eq)...])
                attributes[key, default: []].append(value)
            } else {
                bare.insert(token.uppercased())
            }
        }

        return (attributes, bare)
    }

    private static func firstAttributes(from attributes: [String: [String]]) -> [String: String] {
        attributes.compactMapValues { $0.first }
    }

    private static func parseQuotedOrBareString(_ raw: String) -> String? {
        let value = raw.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !value.isEmpty else { return nil }
        if value.hasPrefix("\""), value.hasSuffix("\""), value.count >= 2 {
            return String(value.dropFirst().dropLast())
        }
        return value
    }

    private static func wrappedMolfileLines(for ctabLines: [String], title: String) -> [String] {
        [title, "", "", "  0  0  0     0  0            999 V3000"] + ctabLines + ["M  END"]
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
                     atomMapNumber: atom.atomMapNumber,
                     aliasLabel: atom.aliasLabel)
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

    private static func tokenizeV30(_ line: String) -> [String] {
        let trimmed = line.trimmingCharacters(in: .whitespaces)
        guard trimmed.hasPrefix("M  V30 ") else { return [] }
        let payload = String(trimmed.dropFirst(7))
        return tokenizeV30Payload(payload)
    }

    private static func tokenizeV30Payload(_ payload: String) -> [String] {
        var tokens: [String] = []
        var current = ""
        var parenDepth = 0
        var inQuotes = false

        for scalar in payload.unicodeScalars {
            let character = Character(scalar)

            if character == "\"" {
                inQuotes.toggle()
                current.append(character)
                continue
            }

            if !inQuotes {
                if character == "(" {
                    parenDepth += 1
                    current.append(character)
                    continue
                }
                if character == ")" {
                    parenDepth = max(0, parenDepth - 1)
                    current.append(character)
                    continue
                }

                if character.isWhitespace && parenDepth == 0 {
                    if !current.isEmpty {
                        tokens.append(current)
                        current.removeAll(keepingCapacity: true)
                    }
                    continue
                }
            }

            current.append(character)
        }

        if !current.isEmpty {
            tokens.append(current)
        }

        return tokens
    }

    private static func unfoldV30Continuations(_ lines: [String]) -> [String] {
        var result: [String] = []
        result.reserveCapacity(lines.count)

        var pendingPayload: String? = nil

        for rawLine in lines {
            let trimmed = rawLine.trimmingCharacters(in: .whitespaces)
            guard trimmed.hasPrefix("M  V30 ") else {
                if let pending = pendingPayload {
                    result.append("M  V30 \(pending)")
                    pendingPayload = nil
                }
                result.append(rawLine)
                continue
            }

            var payload = String(trimmed.dropFirst(7))
            if let pending = pendingPayload {
                payload = pending + payload
            }

            if payload.trimmingCharacters(in: .whitespacesAndNewlines).hasSuffix("-") {
                payload = removeContinuationMarker(payload)
                pendingPayload = payload
            } else {
                result.append("M  V30 \(payload)")
                pendingPayload = nil
            }
        }

        if let pending = pendingPayload {
            result.append("M  V30 \(pending)")
        }

        return result
    }

    private static func removeContinuationMarker(_ payload: String) -> String {
        var trimmed = payload.trimmingCharacters(in: .whitespacesAndNewlines)
        if trimmed.hasSuffix("-") {
            trimmed.removeLast()
        }
        return trimmed + " "
    }

    private static func dropTrailingEmptyLines(_ lines: [String]) -> [String] {
        var out = lines
        while let last = out.last, last.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty {
            out.removeLast()
        }
        return out
    }

    private static let aromaticElementTokens: Set<String> = ["b", "c", "n", "o", "p", "s", "se", "as"]
}
