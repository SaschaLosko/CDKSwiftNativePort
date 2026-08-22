import Foundation

#if canImport(CoreGraphics)
    import CoreGraphics
#endif

public enum CDKMDLV2000Reader {
    private struct SgroupDraft {
        var sourceIndex: Int
        var keyword: String?
        var atomIDs: [Int] = []
        var crossingBondIDs: [Int] = []
        var subscriptText: String?
        var roundBrackets = false
        var connectivity: String?
        var dataFieldName: String?
        var dataValueSegments: [String] = []
        var dataOperator: String?
        var dataUnit: String?
        var dataTag: String?
        var subtype: String?
        var parentAtomIDs: [Int] = []
        var componentNumber: Int?
        var expanded = false
        var brackets: [MoleculeSgroupBracket] = []
        var parentSourceIndices: [Int] = []
    }

    public static func read(text: String) throws -> Molecule {
        let normalized =
            text
            .replacingOccurrences(of: "\r\n", with: "\n")
            .replacingOccurrences(of: "\r", with: "\n")
        return try read(lines: normalized.components(separatedBy: "\n"))
    }

    public static func read(lines: [String]) throws -> Molecule {
        let trimmed = dropTrailingEmptyLines(lines)
        guard trimmed.count >= 4 else {
            throw ChemError.parseFailed("Molfile too short.")
        }

        let countsLine = trimmed[3]
        if countsLine.uppercased().contains("V3000") {
            return try CDKMDLV3000Reader.read(text: trimmed.joined(separator: "\n"))
        }

        let atomCount = parseCountsValue(in: countsLine, at: 0)
        let bondCount = parseCountsValue(in: countsLine, at: 1)
        let chiralFlag = parseCountsValue(in: countsLine, at: 4)

        guard atomCount >= 0 else {
            throw ChemError.parseFailed("Atom count missing/invalid.")
        }
        guard bondCount >= 0 else {
            throw ChemError.parseFailed("Bond count missing/invalid.")
        }

        let firstAtomLine = 4
        let firstBondLine = firstAtomLine + atomCount
        let firstPropertiesLine = firstBondLine + bondCount

        guard trimmed.count >= firstPropertiesLine else {
            throw ChemError.parseFailed("Molfile truncated.")
        }

        var atoms: [Atom] = []
        atoms.reserveCapacity(atomCount)

        for atomIndex in 0..<atomCount {
            let line = trimmed[firstAtomLine + atomIndex]
            atoms.append(try parseAtomLine(line, atomID: atomIndex + 1))
        }

        var bonds: [Bond] = []
        bonds.reserveCapacity(bondCount)

        for bondIndex in 0..<bondCount {
            let line = trimmed[firstBondLine + bondIndex]
            bonds.append(try parseBondLine(line, bondID: bondIndex + 1, atomCount: atomCount))
        }

        let title = trimmed[0].trimmingCharacters(in: .whitespacesAndNewlines)
        var molecule = Molecule(name: title.isEmpty ? "Molecule" : title, atoms: atoms, bonds: bonds)
        var sgroupDrafts: [Int: SgroupDraft] = [:]

        var foundMEnd = false
        var lineIndex = firstPropertiesLine
        while lineIndex < trimmed.count {
            let line = trimmed[lineIndex]
            if line.hasPrefix("M  END") {
                foundMEnd = true
                break
            }

            lineIndex += consumePropertyLine(
                trimmed,
                index: lineIndex,
                molecule: &molecule,
                sgroupDrafts: &sgroupDrafts)
        }

        guard foundMEnd else {
            throw ChemError.parseFailed("Molfile missing M  END record.")
        }

        applySgroupDrafts(sgroupDrafts, to: &molecule)
        applyChiralFlag(chiralFlag, to: &molecule)
        CDKSDFDataFieldParser.applyParsedFields(from: trimmed, to: &molecule)

        if molecule.atomCount > 0,
            let box = molecule.boundingBox(),
            box.width <= 0.0001 && box.height <= 0.0001
        {
            molecule = Depiction2DGenerator.generate(for: molecule)
        }

        return molecule
    }

    private static func parseAtomLine(_ line: String, atomID: Int) throws -> Atom {
        let parts = line.split(whereSeparator: \.isWhitespace).map(String.init)
        guard parts.count >= 4 else {
            throw ChemError.parseFailed("Invalid atom line at atom index \(atomID).")
        }

        let x = parseFixedDouble(line, start: 0, length: 10) ?? Double(parts[0]) ?? 0
        let y = parseFixedDouble(line, start: 10, length: 10) ?? Double(parts[1]) ?? 0
        let z = parseFixedDouble(line, start: 20, length: 10) ?? Double(parts[2]) ?? 0

        let rawSymbol = parseAtomSymbol(line: line, parts: parts)
        var element = normalizeElementToken(rawSymbol)
        var aromatic = isAromaticSymbolToken(rawSymbol)
        var query = queryType(for: rawSymbol)
        var isotopeMassNumber: Int? = nil

        if rawSymbol.uppercased() == "D" {
            element = "H"
            aromatic = false
            isotopeMassNumber = 2
        } else if rawSymbol.uppercased() == "T" {
            element = "H"
            aromatic = false
            isotopeMassNumber = 3
        }

        let massDiff = parseSignedFixedInt(line, start: 34, length: 2) ?? parseInt(parts, index: 4) ?? 0
        if massDiff != 0, isotopeMassNumber == nil {
            let baseMass = majorMassNumber(for: element)
            if baseMass > 0 {
                isotopeMassNumber = baseMass + massDiff
            }
        }

        let chargeCode = parseFixedInt(line, start: 36, length: 3) ?? parseInt(parts, index: 5) ?? 0
        let parity = parseFixedInt(line, start: 39, length: 3) ?? parseInt(parts, index: 6) ?? 0
        let hydrogenCountCode =
            parseFixedInt(line, start: 42, length: 3) ?? parseInt(parts, index: 7) ?? 0
        let valenceCode = parseFixedInt(line, start: 48, length: 3) ?? parseInt(parts, index: 9) ?? 0
        let atomMapFixed = parseFixedInt(line, start: 60, length: 3)
        let atomMapToken = parseInt(parts, index: 12)
        let atomMapValue = [atomMapFixed, atomMapToken].compactMap { $0 }.first(where: { $0 > 0 })

        let (charge, radical) = decodeChargeAndRadical(chargeCode)

        let chirality: AtomChirality
        switch parity {
        case 1:
            chirality = .clockwise
        case 2:
            chirality = .anticlockwise
        default:
            chirality = .none
        }

        let explicitHydrogenCount: Int?
        if query != nil && hydrogenCountCode > 0 {
            explicitHydrogenCount = max(0, hydrogenCountCode - 1)
        } else {
            explicitHydrogenCount = nil
        }

        let valenceOverride: Int?
        if valenceCode == 15 {
            valenceOverride = 0
        } else if valenceCode > 0 {
            valenceOverride = valenceCode
        } else {
            valenceOverride = nil
        }

        if rawSymbol.uppercased() == "ANY" {
            element = "*"
            query = .anyAtom
            aromatic = false
        }

        return Atom(
            id: atomID,
            element: element,
            position: CGPoint(x: x, y: y),
            zPosition: z,
            charge: charge,
            isotopeMassNumber: isotopeMassNumber,
            aromatic: aromatic,
            chirality: chirality,
            explicitHydrogenCount: explicitHydrogenCount,
            queryType: query,
            radical: radical,
            valenceOverride: valenceOverride,
            atomMapNumber: atomMapValue)
    }

    private static func parseBondLine(
        _ line: String,
        bondID: Int,
        atomCount: Int
    ) throws -> Bond {
        let parts = line.split(whereSeparator: \.isWhitespace).map(String.init)
        guard parts.count >= 3 else {
            throw ChemError.parseFailed("Invalid bond line at bond index \(bondID).")
        }

        let a1 = parseFixedInt(line, start: 0, length: 3) ?? parseInt(parts, index: 0) ?? 0
        let a2 = parseFixedInt(line, start: 3, length: 3) ?? parseInt(parts, index: 1) ?? 0
        let typeCode = parseFixedInt(line, start: 6, length: 3) ?? parseInt(parts, index: 2) ?? 1
        let stereoCode = parseFixedInt(line, start: 9, length: 3) ?? parseInt(parts, index: 3) ?? 0
        let topologyCode = parseFixedInt(line, start: 15, length: 3) ?? parseInt(parts, index: 5) ?? 0
        let reactingCenterCode =
            parseFixedInt(line, start: 18, length: 3) ?? parseInt(parts, index: 6) ?? 0

        guard a1 == 0 || (1...atomCount).contains(a1),
            a2 == 0 || (1...atomCount).contains(a2)
        else {
            throw ChemError.parseFailed("Bond references unknown atom index.")
        }
        guard a1 > 0 && a2 > 0 else {
            throw ChemError.parseFailed("Bond references unknown atom index.")
        }

        let (order, queryType) = try bondFromTypeCode(typeCode, rawLine: line)
        return Bond(
            id: bondID,
            a1: a1,
            a2: a2,
            order: order,
            stereo: stereoFromMolfile(stereoCode, bondTypeCode: typeCode),
            queryType: queryType,
            topology: topologyFromMolfile(topologyCode),
            reactingCenterStatus: reactingCenterCode == 0
                ? nil
                : CDKMDLReactingCenterStatus(rawValue: reactingCenterCode))
    }

    private static func consumePropertyLine(
        _ lines: [String],
        index: Int,
        molecule: inout Molecule,
        sgroupDrafts: inout [Int: SgroupDraft]
    ) -> Int {
        let line = lines[index]
        let trimmed = line.trimmingCharacters(in: .whitespaces)

        if isAliasLine(trimmed) {
            guard index + 1 < lines.count,
                let atomID = parseTrailingInt(in: trimmed)
            else {
                return 1
            }
            applyAlias(lines[index + 1], atomID: atomID, to: &molecule)
            return 2
        }

        if trimmed.hasPrefix("V  ") {
            applyAtomValueLine(trimmed, to: &molecule)
            return 1
        }

        if trimmed.hasPrefix("G  ") {
            guard index + 1 < lines.count,
                let atomID = parseFixedInt(trimmed, start: 3, length: 3)
            else {
                return 1
            }
            applyAlias(lines[index + 1], atomID: atomID, to: &molecule)
            return 2
        }

        if isLegacyAtomListLine(trimmed) {
            applyLegacyAtomListLine(trimmed, to: &molecule)
            return 1
        }

        guard trimmed.hasPrefix("M  ") else {
            return 1
        }

        let fields = trimmed.split(whereSeparator: \.isWhitespace).map(String.init)
        guard fields.count >= 2 else { return 1 }

        switch fields[1] {
        case "CHG":
            applyAtomIntPairs(fields, to: &molecule) { atom, value in
                atom.charge = value
            }
        case "ISO":
            applyAtomIntPairs(fields, to: &molecule) { atom, value in
                atom.isotopeMassNumber = value > 0 ? value : nil
            }
        case "MAP":
            applyAtomIntPairs(fields, to: &molecule) { atom, value in
                atom.atomMapNumber = value > 0 ? value : nil
            }
        case "RAD":
            applyAtomIntPairs(fields, to: &molecule) { atom, value in
                atom.radical = value > 0 ? value : nil
            }
        case "RGP":
            applyAtomIntPairs(fields, to: &molecule) { atom, value in
                atom.rGroupLabel = (0...32).contains(value) ? value : nil
                if atom.rGroupLabel != nil {
                    atom.aliasLabel = nil
                }
            }
        case "ALS":
            applyMALSLine(fields, to: &molecule)
        case "APO":
            applyAtomIntPairs(fields, to: &molecule) { atom, value in
                atom.attachmentPoint = value
            }
        case "SUB":
            applyAtomIntPairs(fields, to: &molecule) { atom, value in
                atom.substitutionCount = value
            }
        case "UNS":
            applyAtomIntPairs(fields, to: &molecule) { atom, value in
                atom.unsaturated = value
            }
        case "RBC":
            applyAtomIntPairs(fields, to: &molecule) { atom, value in
                atom.ringBondCount = value
            }
        case "STY":
            applySTypeLine(fields, drafts: &sgroupDrafts)
        case "SAL":
            applySgroupIndexList(fields, drafts: &sgroupDrafts, kind: .atoms)
        case "SBL":
            applySgroupIndexList(fields, drafts: &sgroupDrafts, kind: .bonds)
        case "SPL":
            applySPLLine(fields, drafts: &sgroupDrafts)
        case "SCN":
            applySgroupRepeatedTextLine(fields, drafts: &sgroupDrafts) { draft, value in
                draft.connectivity = value
            }
        case "SDI":
            applySDILine(fields, drafts: &sgroupDrafts)
        case "SMT":
            applySingleSgroupTextLine(trimmed, drafts: &sgroupDrafts) { draft, value in
                draft.subscriptText = value
            }
        case "SBT":
            applySgroupRepeatedIntLine(fields, drafts: &sgroupDrafts) { draft, value in
                draft.roundBrackets = value == 1
            }
        case "SST":
            applySgroupRepeatedTextLine(fields, drafts: &sgroupDrafts) { draft, value in
                draft.subtype = value
            }
        case "SDS":
            applySDSLine(fields, drafts: &sgroupDrafts)
        case "SPA":
            applySgroupIndexList(fields, drafts: &sgroupDrafts, kind: .parentAtoms)
        case "SNC":
            applySgroupRepeatedIntLine(fields, drafts: &sgroupDrafts) { draft, value in
                draft.componentNumber = value
            }
        case "SDT":
            applySDTLine(trimmed, drafts: &sgroupDrafts)
        case "SCD", "SED":
            applySgroupDataLine(trimmed, drafts: &sgroupDrafts)
        default:
            break
        }

        return 1
    }

    private enum SgroupIndexListKind {
        case atoms
        case bonds
        case parentAtoms
    }

    private static func applySTypeLine(_ fields: [String], drafts: inout [Int: SgroupDraft]) {
        guard let count = parseInt(fields, index: 2), count > 0 else { return }
        for pair in 0..<count {
            let sourceIndex = 3 + pair * 2
            let keywordIndex = sourceIndex + 1
            guard keywordIndex < fields.count,
                let source = Int(fields[sourceIndex])
            else { continue }
            var draft = drafts[source] ?? SgroupDraft(sourceIndex: source)
            draft.keyword = fields[keywordIndex].trimmingCharacters(in: .whitespacesAndNewlines)
                .uppercased()
            drafts[source] = draft
        }
    }

    private static func applySgroupIndexList(
        _ fields: [String],
        drafts: inout [Int: SgroupDraft],
        kind: SgroupIndexListKind
    ) {
        guard let source = parseInt(fields, index: 2),
            let count = parseInt(fields, index: 3),
            count > 0
        else { return }

        let values = fields.dropFirst(4).prefix(count).compactMap(Int.init)
        guard !values.isEmpty else { return }

        var draft = drafts[source] ?? SgroupDraft(sourceIndex: source)
        switch kind {
        case .atoms:
            draft.atomIDs.append(contentsOf: values)
            draft.atomIDs = uniqueInts(draft.atomIDs)
        case .bonds:
            draft.crossingBondIDs.append(contentsOf: values)
            draft.crossingBondIDs = uniqueInts(draft.crossingBondIDs)
        case .parentAtoms:
            draft.parentAtomIDs.append(contentsOf: values)
            draft.parentAtomIDs = uniqueInts(draft.parentAtomIDs)
        }
        drafts[source] = draft
    }

    private static func applySPLLine(_ fields: [String], drafts: inout [Int: SgroupDraft]) {
        guard let count = parseInt(fields, index: 2), count > 0 else { return }
        for pair in 0..<count {
            let childIndex = 3 + pair * 2
            let parentIndex = childIndex + 1
            guard parentIndex < fields.count,
                let child = Int(fields[childIndex]),
                let parent = Int(fields[parentIndex])
            else { continue }
            var draft = drafts[child] ?? SgroupDraft(sourceIndex: child)
            draft.parentSourceIndices.append(parent)
            draft.parentSourceIndices = uniqueInts(draft.parentSourceIndices)
            drafts[child] = draft
        }
    }

    private static func applySgroupRepeatedTextLine(
        _ fields: [String],
        drafts: inout [Int: SgroupDraft],
        assign: (inout SgroupDraft, String) -> Void
    ) {
        guard let count = parseInt(fields, index: 2), count > 0 else { return }
        for pair in 0..<count {
            let sourceIndex = 3 + pair * 2
            let valueIndex = sourceIndex + 1
            guard valueIndex < fields.count,
                let source = Int(fields[sourceIndex])
            else { continue }
            var draft = drafts[source] ?? SgroupDraft(sourceIndex: source)
            assign(&draft, fields[valueIndex].trimmingCharacters(in: .whitespacesAndNewlines))
            drafts[source] = draft
        }
    }

    private static func applySgroupRepeatedIntLine(
        _ fields: [String],
        drafts: inout [Int: SgroupDraft],
        assign: (inout SgroupDraft, Int) -> Void
    ) {
        guard let count = parseInt(fields, index: 2), count > 0 else { return }
        for pair in 0..<count {
            let sourceIndex = 3 + pair * 2
            let valueIndex = sourceIndex + 1
            guard valueIndex < fields.count,
                let source = Int(fields[sourceIndex]),
                let value = Int(fields[valueIndex])
            else { continue }
            var draft = drafts[source] ?? SgroupDraft(sourceIndex: source)
            assign(&draft, value)
            drafts[source] = draft
        }
    }

    private static func applySDILine(_ fields: [String], drafts: inout [Int: SgroupDraft]) {
        guard let source = parseInt(fields, index: 2),
            let valueCount = parseInt(fields, index: 3),
            valueCount >= 4
        else { return }
        let values = fields.dropFirst(4).compactMap(Double.init)
        guard values.count >= 4 else { return }

        var draft = drafts[source] ?? SgroupDraft(sourceIndex: source)
        draft.brackets.append(
            MoleculeSgroupBracket(
                firstPoint: CGPoint(x: values[0], y: values[1]),
                secondPoint: CGPoint(x: values[2], y: values[3])))
        drafts[source] = draft
    }

    private static func applySingleSgroupTextLine(
        _ line: String,
        drafts: inout [Int: SgroupDraft],
        assign: (inout SgroupDraft, String) -> Void
    ) {
        guard let source = parseFixedInt(line, start: 7, length: 3) else { return }
        let start = line.index(line.startIndex, offsetBy: min(11, line.count))
        let value = String(line[start...]).trimmingCharacters(in: .whitespacesAndNewlines)
        var draft = drafts[source] ?? SgroupDraft(sourceIndex: source)
        assign(&draft, value)
        drafts[source] = draft
    }

    private static func applySDSLine(_ fields: [String], drafts: inout [Int: SgroupDraft]) {
        guard fields.count >= 5, fields[2].uppercased() == "EXP",
            let count = Int(fields[3]), count > 0
        else { return }
        for raw in fields.dropFirst(4).prefix(count) {
            guard let source = Int(raw) else { continue }
            var draft = drafts[source] ?? SgroupDraft(sourceIndex: source)
            draft.expanded = true
            drafts[source] = draft
        }
    }

    private static func applySDTLine(_ line: String, drafts: inout [Int: SgroupDraft]) {
        let fields = line.split(whereSeparator: \.isWhitespace).map(String.init)
        guard let source = parseFixedInt(line, start: 7, length: 3) ?? parseInt(fields, index: 2) else {
            return
        }
        var draft = drafts[source] ?? SgroupDraft(sourceIndex: source)

        if let fixedFieldName = substring(line, start: 11, length: 30)?.trimmingCharacters(
            in: .whitespaces),
            !fixedFieldName.isEmpty
        {
            draft.dataFieldName = fixedFieldName
        } else if let tokenFieldName = parseTokenizedSDTComponent(fields, index: 3) {
            draft.dataFieldName = tokenFieldName
        }

        let unit =
            substring(line, start: 43, length: 20)?.trimmingCharacters(in: .whitespaces)
            ?? parseTokenizedSDTComponent(fields, index: 4)
        let tag =
            substring(line, start: 63, length: 2)?.trimmingCharacters(in: .whitespaces)
            ?? parseTokenizedSDTComponent(fields, index: 5)
        let opStart = min(65, line.count)
        let opIndex = line.index(line.startIndex, offsetBy: opStart)
        let op = String(line[opIndex...]).trimmingCharacters(in: .whitespaces)
        let fallbackOperator = parseTokenizedSDTComponent(fields, index: 6)
        draft.dataUnit = (unit?.isEmpty == false) ? unit : draft.dataUnit
        draft.dataTag = (tag?.isEmpty == false) ? tag : draft.dataTag
        if !op.isEmpty {
            draft.dataOperator = op
        } else if let fallbackOperator, !fallbackOperator.isEmpty {
            draft.dataOperator = fallbackOperator
        }
        drafts[source] = draft
    }

    private static func applySgroupDataLine(_ line: String, drafts: inout [Int: SgroupDraft]) {
        guard let source = parseFixedInt(line, start: 7, length: 3) else { return }
        let start = min(11, line.count)
        let startIndex = line.index(line.startIndex, offsetBy: start)
        let payload = String(line[startIndex...])
        var draft = drafts[source] ?? SgroupDraft(sourceIndex: source)
        draft.dataValueSegments.append(payload)
        drafts[source] = draft
    }

    private static func applyAlias(_ alias: String, atomID: Int, to molecule: inout Molecule) {
        guard let atomIndex = molecule.indexOfAtom(id: atomID) else { return }
        let trimmed = alias.trimmingCharacters(in: .newlines)
        if let rValue = parseRGroupLabel(trimmed) {
            molecule.atoms[atomIndex].rGroupLabel = rValue
            molecule.atoms[atomIndex].aliasLabel = nil
            if molecule.atoms[atomIndex].element.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty {
                molecule.atoms[atomIndex].element = "R#"
            }
            return
        }
        molecule.atoms[atomIndex].aliasLabel = trimmed
    }

    private static func applyAtomValueLine(_ line: String, to molecule: inout Molecule) {
        guard
            let atomID = parseFixedInt(line, start: 3, length: 3)
                ?? parseTrailingPrefixedInt(line, prefix: "V")
        else {
            return
        }
        let commentStart = min(7, line.count)
        let commentIndex = line.index(line.startIndex, offsetBy: commentStart)
        let value = String(line[commentIndex...])
        guard let atomIndex = molecule.indexOfAtom(id: atomID) else { return }
        molecule.atoms[atomIndex].atomValue = value
    }

    private static func applyLegacyAtomListLine(_ line: String, to molecule: inout Molecule) {
        let fields = line.split(whereSeparator: \.isWhitespace).map(String.init)
        guard fields.count >= 4,
            let atomID = Int(fields[0]),
            let count = Int(fields[2]),
            count > 0,
            let atomIndex = molecule.indexOfAtom(id: atomID)
        else { return }

        let negated = fields[1].uppercased() == "T"
        let numbers = fields.dropFirst(3).prefix(count).compactMap(Int.init)
        guard !numbers.isEmpty else { return }
        let entries = numbers.map(symbolForAtomicNumber)

        molecule.atoms[atomIndex].element = "L"
        molecule.atoms[atomIndex].queryType = .anyAtom
        molecule.atoms[atomIndex].atomList = entries
        molecule.atoms[atomIndex].atomListIsNegated = negated
    }

    private static func applyMALSLine(_ fields: [String], to molecule: inout Molecule) {
        guard fields.count >= 6,
            let atomID = Int(fields[2]),
            let count = Int(fields[3]),
            count > 0,
            let atomIndex = molecule.indexOfAtom(id: atomID)
        else {
            return
        }

        let negated = fields[4].uppercased() == "T"
        let entries = Array(fields.dropFirst(5).prefix(count)).map(normalizeElementToken)
        guard entries.count == count else { return }

        molecule.atoms[atomIndex].element = "L"
        molecule.atoms[atomIndex].queryType = .anyAtom
        molecule.atoms[atomIndex].atomList = entries
        molecule.atoms[atomIndex].atomListIsNegated = negated
    }

    private static func applyAtomIntPairs(
        _ fields: [String],
        to molecule: inout Molecule,
        assign: (inout Atom, Int) -> Void
    ) {
        for (atomID, value) in parseAtomValuePairs(fields) {
            guard let atomIndex = molecule.indexOfAtom(id: atomID) else { continue }
            assign(&molecule.atoms[atomIndex], value)
        }
    }

    private static func applySgroupDrafts(_ drafts: [Int: SgroupDraft], to molecule: inout Molecule) {
        guard !drafts.isEmpty else { return }

        let orderedDrafts = drafts.keys.sorted().compactMap { drafts[$0] }
        var sourceToMoleculeIndex: [Int: Int] = [:]

        for draft in orderedDrafts {
            let keyword = draft.keyword?.uppercased()
            var crossingBondIDs = uniqueInts(draft.crossingBondIDs)
            if crossingBondIDs.isEmpty, keyword == "SUP" {
                crossingBondIDs = inferCrossingBondIDs(atomIDs: draft.atomIDs, in: molecule)
            }

            let dataValue = draft.dataValueSegments.isEmpty ? nil : draft.dataValueSegments.joined()

            let sgroup = MoleculeSgroup(
                kind: sgroupKind(for: keyword),
                keyword: keyword,
                atomIDs: uniqueInts(draft.atomIDs),
                crossingBondIDs: crossingBondIDs,
                subscriptText: emptyToNil(draft.subscriptText),
                superscriptText: nil,
                roundBrackets: draft.roundBrackets,
                connectivity: emptyToNil(draft.connectivity),
                dataFieldName: emptyToNil(draft.dataFieldName),
                dataValue: emptyToNil(dataValue),
                dataOperator: emptyToNil(draft.dataOperator),
                dataUnit: emptyToNil(draft.dataUnit),
                dataTag: emptyToNil(draft.dataTag),
                subtype: emptyToNil(draft.subtype),
                parentAtomIDs: uniqueInts(draft.parentAtomIDs),
                componentNumber: draft.componentNumber,
                expanded: draft.expanded,
                brackets: draft.brackets,
                childGroupIndices: [])
            sourceToMoleculeIndex[draft.sourceIndex] = molecule.sgroups.count
            molecule.sgroups.append(sgroup)
        }

        for draft in orderedDrafts {
            guard let childIndex = sourceToMoleculeIndex[draft.sourceIndex] else { continue }
            for parentSource in draft.parentSourceIndices {
                guard let parentIndex = sourceToMoleculeIndex[parentSource],
                    molecule.sgroups.indices.contains(parentIndex)
                else { continue }
                molecule.sgroups[parentIndex].childGroupIndices.append(childIndex)
            }
        }

        for sgroup in molecule.sgroups where sgroup.kind == .data {
            guard let fieldName = sgroup.dataFieldName, !fieldName.isEmpty else { continue }
            molecule.ensureDataField(named: fieldName)
            if let dataValue = sgroup.dataValue, !dataValue.isEmpty {
                molecule.appendDataFieldValue(dataValue, named: fieldName)
            }
        }
    }

    private static func inferCrossingBondIDs(atomIDs: [Int], in molecule: Molecule) -> [Int] {
        let atomIDSet = Set(atomIDs)
        guard !atomIDSet.isEmpty else { return [] }
        return molecule.bonds.compactMap { bond in
            let a1In = atomIDSet.contains(bond.a1)
            let a2In = atomIDSet.contains(bond.a2)
            return a1In == a2In ? nil : bond.id
        }.sorted()
    }

    private static func applyChiralFlag(_ chiralFlag: Int, to molecule: inout Molecule) {
        let stereoIndices = molecule.atoms.indices.filter { molecule.atoms[$0].chirality != .none }
        guard !stereoIndices.isEmpty else { return }

        let kind = chiralFlag == 1 ? "abs" : "or"
        let number = chiralFlag == 1 ? 0 : 1
        for idx in stereoIndices where molecule.atoms[idx].cxStereoGroup == nil {
            molecule.atoms[idx].cxStereoGroup = CDKCxSmilesParser.encodeStereoGroup(
                kind: kind, number: number)
        }
    }

    private static func parseAtomValuePairs(_ fields: [String]) -> [(Int, Int)] {
        guard fields.count >= 4,
            let pairCount = Int(fields[2]),
            pairCount > 0
        else {
            return []
        }

        var pairs: [(Int, Int)] = []
        for pair in 0..<pairCount {
            let atomTokenIndex = 3 + pair * 2
            let valueTokenIndex = atomTokenIndex + 1
            guard valueTokenIndex < fields.count,
                let atomID = Int(fields[atomTokenIndex]),
                let value = Int(fields[valueTokenIndex])
            else {
                continue
            }
            pairs.append((atomID, value))
        }
        return pairs
    }

    private static func parseCountsValue(in countsLine: String, at index: Int) -> Int {
        if let fixed = parseFixedInt(countsLine, start: index * 3, length: 3) {
            return fixed
        }
        let parts = countsLine.split(whereSeparator: \.isWhitespace)
        guard index < parts.count else { return -1 }
        return Int(parts[index]) ?? -1
    }

    private static func parseAtomSymbol(line: String, parts: [String]) -> String {
        if let fixed = substring(line, start: 31, length: 3)?.trimmingCharacters(in: .whitespaces),
            !fixed.isEmpty
        {
            return fixed
        }
        return parts.count > 3 ? parts[3] : "C"
    }

    private static func queryType(for symbol: String) -> AtomQueryType? {
        switch symbol.uppercased() {
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

    private static func decodeChargeAndRadical(_ code: Int) -> (charge: Int, radical: Int?) {
        switch code {
        case 1: return (3, nil)
        case 2: return (2, nil)
        case 3: return (1, nil)
        case 4: return (0, 2)
        case 5: return (-1, nil)
        case 6: return (-2, nil)
        case 7: return (-3, nil)
        default: return (0, nil)
        }
    }

    private static func bondFromTypeCode(_ code: Int, rawLine: String) throws -> (
        BondOrder, BondQueryType?
    ) {
        switch code {
        case 1:
            return (.single, nil)
        case 2:
            return (.double, nil)
        case 3:
            return (.triple, nil)
        case 4:
            return (.aromatic, nil)
        case 5:
            return (.single, .singleOrDouble)
        case 6:
            return (.single, .singleOrAromatic)
        case 7:
            return (.double, .doubleOrAromatic)
        case 8:
            return (.single, .any)
        default:
            throw ChemError.parseFailed("Unsupported V2000 bond type \(code): \(rawLine)")
        }
    }

    private static func stereoFromMolfile(_ code: Int, bondTypeCode: Int) -> BondStereo {
        switch code {
        case 1:
            return .up
        case 3:
            return bondTypeCode == 2 ? .either : .none
        case 4:
            return .either
        case 6:
            return .down
        default:
            return .none
        }
    }

    private static func topologyFromMolfile(_ code: Int) -> BondTopology? {
        switch code {
        case 1:
            return .ring
        case 2:
            return .chain
        default:
            return nil
        }
    }

    private static func normalizeElementToken(_ token: String) -> String {
        let trimmed = token.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { return "C" }
        if trimmed == "*" || trimmed == "R#" { return trimmed }
        if trimmed.uppercased() == "ANY" { return "*" }
        if trimmed.allSatisfy(\.isLetter) {
            guard let first = trimmed.first else { return trimmed }
            let head = String(first).uppercased()
            let tail = String(trimmed.dropFirst()).lowercased()
            return head + tail
        }
        return trimmed
    }

    private static func isAromaticSymbolToken(_ token: String) -> Bool {
        let trimmed = token.trimmingCharacters(in: .whitespacesAndNewlines)
        guard trimmed.allSatisfy(\.isLetter) else { return false }
        let lower = trimmed.lowercased()
        return trimmed == lower && aromaticTokens.contains(lower)
    }

    private static func majorMassNumber(for element: String) -> Int {
        let mass = CDKDescriptorSupport.monoisotopicAtomicMass(forElementSymbol: element)
        guard mass > 0 else { return 0 }
        return Int(mass.rounded())
    }

    private static func parseRGroupLabel(_ label: String) -> Int? {
        let trimmed = label.trimmingCharacters(in: .whitespacesAndNewlines).uppercased()
        guard trimmed.hasPrefix("R") else { return nil }
        if trimmed == "R" { return 0 }
        return Int(trimmed.dropFirst())
    }

    private static func sgroupKind(for keyword: String?) -> MoleculeSgroup.Kind {
        switch keyword?.uppercased() {
        case "SRU":
            return .structureRepeatUnit
        case "DAT":
            return .data
        case "COP", "CRO", "GRA", "MOD", "MON", "MER", "COM", "MIX", "FOR":
            return .polymer
        default:
            return .generic
        }
    }

    private static func isAliasLine(_ line: String) -> Bool {
        guard line.first == "A" else { return false }
        return parseTrailingInt(in: line) != nil
    }

    private static func isLegacyAtomListLine(_ line: String) -> Bool {
        let fields = line.split(whereSeparator: \.isWhitespace)
        guard fields.count >= 4 else { return false }
        return Int(fields[0]) != nil && (fields[1] == "T" || fields[1] == "F")
    }

    private static func parseTrailingInt(in line: String) -> Int? {
        let digits = line.split(whereSeparator: \.isWhitespace).last
        return digits.flatMap { Int(String($0)) }
    }

    private static func parseTrailingPrefixedInt(_ line: String, prefix: String) -> Int? {
        let fields = line.split(whereSeparator: \.isWhitespace).map(String.init)
        guard fields.count >= 2, fields[0] == prefix else { return nil }
        return Int(fields[1])
    }

    private static func symbolForAtomicNumber(_ atomicNumber: Int) -> String {
        guard atomicNumber > 0, atomicNumber < atomicSymbols.count else {
            return "*"
        }
        return atomicSymbols[atomicNumber]
    }

    private static func parseInt(_ parts: [String], index: Int) -> Int? {
        guard index < parts.count else { return nil }
        return Int(parts[index])
    }

    private static func parseTokenizedSDTComponent(_ parts: [String], index: Int) -> String? {
        guard index < parts.count else { return nil }
        let value = parts[index].trimmingCharacters(in: .whitespacesAndNewlines)
        return value.isEmpty ? nil : value
    }

    private static func parseFixedInt(_ s: String, start: Int, length: Int) -> Int? {
        guard let slice = substring(s, start: start, length: length) else { return nil }
        let trimmed = slice.trimmingCharacters(in: .whitespaces)
        guard !trimmed.isEmpty else { return nil }
        return Int(trimmed)
    }

    private static func parseSignedFixedInt(_ s: String, start: Int, length: Int) -> Int? {
        guard let slice = substring(s, start: start, length: length) else { return nil }
        let trimmed = slice.trimmingCharacters(in: .whitespaces)
        guard !trimmed.isEmpty else { return nil }
        return Int(trimmed)
    }

    private static func parseFixedDouble(_ s: String, start: Int, length: Int) -> Double? {
        guard let slice = substring(s, start: start, length: length) else { return nil }
        let trimmed = slice.trimmingCharacters(in: .whitespaces)
        guard !trimmed.isEmpty else { return nil }
        return Double(trimmed)
    }

    private static func substring(_ s: String, start: Int, length: Int) -> String? {
        guard start >= 0, length > 0 else { return nil }
        guard s.count >= start + length else { return nil }

        let lower = s.index(s.startIndex, offsetBy: start)
        let upper = s.index(lower, offsetBy: length)
        return String(s[lower..<upper])
    }

    private static func dropTrailingEmptyLines(_ lines: [String]) -> [String] {
        var out = lines
        while let last = out.last, last.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty {
            out.removeLast()
        }
        return out
    }

    private static func uniqueInts(_ values: [Int]) -> [Int] {
        var seen = Set<Int>()
        return values.filter { seen.insert($0).inserted }
    }

    private static func emptyToNil(_ text: String?) -> String? {
        guard let text else { return nil }
        let trimmed = text.trimmingCharacters(in: .whitespacesAndNewlines)
        return trimmed.isEmpty ? nil : trimmed
    }

    private static let aromaticTokens: Set<String> = ["b", "c", "n", "o", "p", "s", "se", "as"]

    // 0th entry is a placeholder so that array index equals atomic number.
    private static let atomicSymbols: [String] = [
        "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
        "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
        "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
    ]
}
