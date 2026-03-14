import Foundation

public enum CDKMDLV2000Writer {
    public struct Options {
        public var programName: String?

        public init(programName: String? = nil) {
            self.programName = programName
        }
    }

    private struct OutputOrder {
        let atoms: [Atom]
        let bonds: [Bond]
        let atomIndexByID: [Int: Int]
        let bondIndexByID: [Int: Int]
    }

    public static func write(_ molecule: Molecule,
                             options: Options = Options()) throws -> String {
        let order = try outputOrder(for: molecule)
        let atomListCount = order.atoms.reduce(0) { partial, atom in
            partial + ((atom.atomList?.isEmpty == false) ? 1 : 0)
        }
        let chiralFlag = chiralFlag(for: molecule)

        var lines: [String] = []
        let title = molecule.name.trimmingCharacters(in: .whitespacesAndNewlines)
        lines.append(title.isEmpty ? "Molecule" : title)
        lines.append(headerProgramLine(for: molecule, options: options))
        lines.append("")
        lines.append(countsLine(atomCount: order.atoms.count,
                                bondCount: order.bonds.count,
                                atomListCount: atomListCount,
                                chiralFlag: chiralFlag))

        var aliasLines: [String] = []
        var legacyAtomListLines: [String] = []
        var alsLines: [String] = []
        var atomValues: [(Int, String)] = []
        var charges: [(Int, Int)] = []
        var isotopes: [(Int, Int)] = []
        var radicals: [(Int, Int)] = []
        var rgroups: [(Int, Int)] = []

        for (outputIndex, atom) in order.atoms.enumerated() {
            let atomIndex = outputIndex + 1
            let encoding = atomEncoding(for: atom, outputIndex: atomIndex, in: molecule)
            lines.append(encoding.line)

            if let aliasLine = encoding.aliasLine {
                aliasLines.append(contentsOf: aliasLine)
            }
            if let legacyAtomListLine = encoding.legacyAtomListLine {
                legacyAtomListLines.append(legacyAtomListLine)
            }
            if let alsLine = encoding.alsLine {
                alsLines.append(alsLine)
            }
            if let atomValue = atom.atomValue?.trimmingCharacters(in: .newlines), !atomValue.trimmingCharacters(in: .whitespaces).isEmpty {
                atomValues.append((atomIndex, atomValue))
            }
            if atom.charge != 0 {
                charges.append((atomIndex, atom.charge))
            }
            if shouldWriteIsotopeProperty(for: atom) {
                isotopes.append((atomIndex, atom.isotopeMassNumber ?? 0))
            }
            if let radical = serializedRadicalMultiplicity(for: atom), radical > 0 {
                radicals.append((atomIndex, radical))
            }
            if let rGroup = normalizedRGroupLabel(for: atom) {
                rgroups.append((atomIndex, rGroup))
            }
        }

        for bond in order.bonds {
            lines.append(try bondLine(for: bond, order: order))
        }

        for (atomIndex, value) in atomValues {
            lines.append("V  \(formatMDLInt(atomIndex, 3)) \(value)")
        }

        appendPropertyLines(prefix: "M  CHG", pairs: charges, to: &lines)
        appendPropertyLines(prefix: "M  ISO", pairs: isotopes, to: &lines)
        appendPropertyLines(prefix: "M  RAD", pairs: radicals, to: &lines)
        appendPropertyLines(prefix: "M  RGP", pairs: rgroups, to: &lines)

        lines.append(contentsOf: aliasLines)
        lines.append(contentsOf: legacyAtomListLines)
        lines.append(contentsOf: alsLines)
        lines.append(contentsOf: sgroupLines(for: molecule, order: order))
        lines.append("M  END")

        return lines.joined(separator: "\n")
    }

    private struct AtomEncoding {
        let line: String
        let aliasLine: [String]?
        let legacyAtomListLine: String?
        let alsLine: String?
    }

    private static func outputOrder(for molecule: Molecule) throws -> OutputOrder {
        let atoms = molecule.atoms.sorted { $0.id < $1.id }
        let atomIDs = Set(atoms.map(\.id))
        let bonds = molecule.bonds.sorted { $0.id < $1.id }

        for bond in bonds {
            guard atomIDs.contains(bond.a1), atomIDs.contains(bond.a2) else {
                throw ChemError.parseFailed("Bond references unknown atom id while writing Molfile.")
            }
        }

        let atomIndexByID = Dictionary(uniqueKeysWithValues: atoms.enumerated().map { ($1.id, $0 + 1) })
        let bondIndexByID = Dictionary(uniqueKeysWithValues: bonds.enumerated().map { ($1.id, $0 + 1) })
        return OutputOrder(atoms: atoms,
                           bonds: bonds,
                           atomIndexByID: atomIndexByID,
                           bondIndexByID: bondIndexByID)
    }

    private static func atomEncoding(for atom: Atom,
                                     outputIndex: Int,
                                     in molecule: Molecule) -> AtomEncoding {
        let symbolDecision = atomBlockSymbol(for: atom)
        var fields: [Int] = Array(repeating: 0, count: 12)
        fields[0] = inlineMassDifference(for: atom, symbolToken: symbolDecision.symbol)
        fields[1] = atomBlockChargeCode(for: atom)
        fields[2] = atomParity(for: atom)
        fields[5] = atomValenceField(for: atom, in: molecule)
        fields[9] = atom.atomMapNumber ?? atom.atomClass ?? 0

        var line = ""
        line += formatMDLFloat(Double(atom.position.x))
        line += formatMDLFloat(Double(atom.position.y))
        line += formatMDLFloat(atom.zPosition ?? 0)
        line += " "
        line += formatMDLString(symbolDecision.symbol, 3)
        line += formatMDLInt(fields[0], 2)
        for idx in 1..<fields.count {
            line += formatMDLInt(fields[idx], 3)
        }

        let aliasLine: [String]?
        if let alias = symbolDecision.alias {
            aliasLine = ["A\(formatMDLInt(outputIndex, 5))", truncate(alias, to: 70)]
        } else {
            aliasLine = nil
        }

        let legacyAtomListLine = legacyAtomListLine(for: atom, outputIndex: outputIndex)
        let alsLine = malsLine(for: atom, outputIndex: outputIndex)

        return AtomEncoding(line: line,
                            aliasLine: aliasLine,
                            legacyAtomListLine: legacyAtomListLine,
                            alsLine: alsLine)
    }

    private static func atomBlockSymbol(for atom: Atom) -> (symbol: String, alias: String?) {
        if atom.atomList?.isEmpty == false {
            return ("L", nil)
        }

        if let queryType = atom.queryType {
            switch queryType {
            case .anyAtom:
                return ("*", nil)
            case .anyNonHydrogen:
                return ("A", nil)
            case .anyHetero:
                return ("Q", nil)
            }
        }

        if normalizedRGroupLabel(for: atom) != nil {
            return ("R#", nil)
        }

        let canonicalElement = canonicalElementSymbol(atom.element)
        if canonicalElement == "H", atom.isotopeMassNumber == 2 {
            return ("D", nil)
        }
        if canonicalElement == "H", atom.isotopeMassNumber == 3 {
            return ("T", nil)
        }

        if let alias = effectiveAliasLabel(for: atom) {
            if alias.count <= 3 {
                return (alias, nil)
            }
            let fallback = fallbackAtomSymbol(for: atom)
            return (fallback, alias)
        }

        if canonicalElement.isEmpty {
            return ("C", nil)
        }
        if canonicalElement.count <= 3 {
            return (canonicalElement, nil)
        }
        return ("A", canonicalElement)
    }

    private static func effectiveAliasLabel(for atom: Atom) -> String? {
        if let aliasLabel = atom.aliasLabel?.trimmingCharacters(in: .whitespacesAndNewlines),
           !aliasLabel.isEmpty {
            return aliasLabel
        }

        let element = atom.element.trimmingCharacters(in: .whitespacesAndNewlines)
        if element.isEmpty || element == "*" || element == "R#" {
            return nil
        }
        if isPlausibleElementSymbol(element) {
            return nil
        }
        return element
    }

    private static func fallbackAtomSymbol(for atom: Atom) -> String {
        let element = canonicalElementSymbol(atom.element)
        if element.isEmpty {
            return "A"
        }
        if element.count <= 3, isPlausibleElementSymbol(element) {
            return element
        }
        if element == "*" {
            return "*"
        }
        return "A"
    }

    private static func legacyAtomListLine(for atom: Atom, outputIndex: Int) -> String? {
        guard let atomList = atom.atomList, !atomList.isEmpty, atomList.count <= 5 else { return nil }
        let atomicNumbers = atomList.compactMap(atomicNumber(forElementSymbol:))
        guard atomicNumbers.count == atomList.count else { return nil }

        var line = formatMDLInt(outputIndex, 3)
        line += atom.atomListIsNegated ? " T    " : " F    "
        line += formatMDLInt(atomicNumbers.count, 1)
        for atomicNumber in atomicNumbers {
            line += " " + formatMDLInt(atomicNumber, 3)
        }
        return line
    }

    private static func malsLine(for atom: Atom, outputIndex: Int) -> String? {
        guard let atomList = atom.atomList, !atomList.isEmpty else { return nil }
        var line = "M  ALS "
        line += formatMDLInt(outputIndex, 3)
        line += formatMDLInt(atomList.count, 3)
        line += atom.atomListIsNegated ? " T " : " F "
        for symbol in atomList {
            line += formatMDLString(canonicalElementSymbol(symbol), 4)
        }
        return line
    }

    private static func bondLine(for bond: Bond, order: OutputOrder) throws -> String {
        let beginAtomID: Int
        let endAtomID: Int
        switch bond.stereo {
        case .upReversed, .downReversed:
            beginAtomID = bond.a2
            endAtomID = bond.a1
        default:
            beginAtomID = bond.a1
            endAtomID = bond.a2
        }

        guard let beginIndex = order.atomIndexByID[beginAtomID],
              let endIndex = order.atomIndexByID[endAtomID] else {
            throw ChemError.parseFailed("Bond references unknown atom while writing Molfile.")
        }

        var line = ""
        line += formatMDLInt(beginIndex, 3)
        line += formatMDLInt(endIndex, 3)
        line += formatMDLInt(try serializedBondType(for: bond), 3)
        line += formatMDLInt(serializedBondStereo(for: bond), 3)
        line += formatMDLInt(0, 3)
        line += formatMDLInt(serializedBondTopology(for: bond), 3)
        line += formatMDLInt(0, 3)
        return line
    }

    private static func serializedBondType(for bond: Bond) throws -> Int {
        if let queryType = bond.queryType {
            switch queryType {
            case .singleOrDouble:
                return 5
            case .singleOrAromatic:
                return 6
            case .doubleOrAromatic:
                return 7
            case .any:
                return 8
            }
        }

        switch bond.order {
        case .single:
            return 1
        case .double:
            return 2
        case .triple:
            return 3
        case .aromatic:
            return 4
        }
    }

    private static func serializedBondStereo(for bond: Bond) -> Int {
        switch bond.stereo {
        case .up, .upReversed:
            return 1
        case .either:
            return bond.order == .double ? 3 : 4
        case .down, .downReversed:
            return 6
        case .none:
            return 0
        }
    }

    private static func serializedBondTopology(for bond: Bond) -> Int {
        switch bond.topology {
        case .ring:
            return 1
        case .chain:
            return 2
        case .none:
            return 0
        }
    }

    private static func atomBlockChargeCode(for atom: Atom) -> Int {
        switch atom.charge {
        case -3: return 7
        case -2: return 6
        case -1: return 5
        case 0:
            if serializedRadicalMultiplicity(for: atom) == 2 {
                return 4
            }
            return 0
        case 1: return 3
        case 2: return 2
        case 3: return 1
        default:
            return 0
        }
    }

    private static func serializedRadicalMultiplicity(for atom: Atom) -> Int? {
        if let radicalType = atom.radicalType {
            switch radicalType {
            case .monovalent:
                return 2
            case .divalentSinglet:
                return 1
            case .divalentTriplet:
                return 3
            default:
                break
            }
        }
        guard let radical = atom.radical, radical > 0 else { return nil }
        return radical
    }

    private static func inlineMassDifference(for atom: Atom, symbolToken: String) -> Int {
        guard let isotope = atom.isotopeMassNumber else { return 0 }
        if symbolToken.uppercased() == "D" || symbolToken.uppercased() == "T" {
            return 0
        }
        let baseMass = Int(CDKDescriptorSupport.monoisotopicAtomicMass(forElementSymbol: atom.element).rounded())
        guard baseMass > 0 else { return 0 }
        let diff = isotope - baseMass
        return (-3...4).contains(diff) ? diff : 0
    }

    private static func shouldWriteIsotopeProperty(for atom: Atom) -> Bool {
        guard let isotope = atom.isotopeMassNumber, isotope > 0 else { return false }
        let element = canonicalElementSymbol(atom.element)
        if element == "H", isotope == 2 || isotope == 3 {
            return false
        }
        return true
    }

    private static func atomParity(for atom: Atom) -> Int {
        switch atom.chirality {
        case .clockwise:
            return 1
        case .anticlockwise:
            return 2
        case .none:
            return 0
        }
    }

    private static func atomValenceField(for atom: Atom, in molecule: Molecule) -> Int {
        let actualValence: Int?
        if let override = atom.valenceOverride {
            actualValence = override
        } else if let explicitHydrogenCount = atom.explicitHydrogenCount {
            let explicitValence = molecule.bonds(forAtom: atom.id).reduce(0) { partial, bond in
                partial + bondOrderContribution(for: bond.order)
            }
            actualValence = explicitValence + max(0, explicitHydrogenCount)
        } else {
            actualValence = nil
        }

        guard let actualValence else { return 0 }

        if actualValence == 0 {
            return 15
        }
        if !(1...14).contains(actualValence) {
            return 0
        }

        if let implied = defaultValence(for: atom), implied == actualValence {
            return 0
        }
        return actualValence
    }

    private static func appendPropertyLines(prefix: String,
                                            pairs: [(Int, Int)],
                                            to lines: inout [String]) {
        guard !pairs.isEmpty else { return }
        let sorted = pairs.sorted { lhs, rhs in lhs.0 < rhs.0 }

        var start = 0
        while start < sorted.count {
            let end = min(sorted.count, start + 8)
            let chunk = sorted[start..<end]
            var line = prefix + formatMDLInt(chunk.count, 3)
            for (atomIndex, value) in chunk {
                line += formatMDLInt(atomIndex, 4)
                line += formatMDLInt(value, 4)
            }
            lines.append(line)
            start = end
        }
    }

    private static func countsLine(atomCount: Int,
                                   bondCount: Int,
                                   atomListCount: Int,
                                   chiralFlag: Int) -> String {
        formatMDLInt(atomCount, 3)
            + formatMDLInt(bondCount, 3)
            + formatMDLInt(atomListCount, 3)
            + formatMDLInt(0, 3)
            + formatMDLInt(chiralFlag, 3)
            + formatMDLInt(0, 3)
            + "            999 V2000"
    }

    private static func headerProgramLine(for molecule: Molecule,
                                          options: Options) -> String {
        let formatter = DateFormatter()
        formatter.locale = Locale(identifier: "en_US_POSIX")
        formatter.dateFormat = "MMddyyHHmm"
        let programName = normalizedProgramName(options.programName)
        let dimension: String
        switch dimensionality(for: molecule) {
        case 2:
            dimension = "2D"
        case 3:
            dimension = "3D"
        default:
            dimension = ""
        }
        return "  \(programName)\(formatter.string(from: Date()))\(dimension)"
    }

    private static func normalizedProgramName(_ raw: String?) -> String {
        let trimmed = (raw ?? "CDK").trimmingCharacters(in: .whitespacesAndNewlines)
        if trimmed.isEmpty {
            return "        "
        }
        let limited = String(trimmed.prefix(8))
        if limited.count >= 8 {
            return limited
        }
        return limited + String(repeating: " ", count: 8 - limited.count)
    }

    private static func dimensionality(for molecule: Molecule) -> Int {
        for atom in molecule.atoms {
            if let z = atom.zPosition, abs(z) > 0.000001 {
                return 3
            }
            if atom.position != .zero {
                return 2
            }
        }
        return 0
    }

    private static func chiralFlag(for molecule: Molecule) -> Int {
        let stereoAtoms = molecule.atoms.filter { $0.chirality != .none }
        guard !stereoAtoms.isEmpty else { return 0 }

        let groups = stereoAtoms.compactMap { atom -> String? in
            guard let encoded = atom.cxStereoGroup,
                  let decoded = CDKCxSmilesParser.decodeStereoGroup(encoded) else {
                return "abs"
            }
            return decoded.kind
        }

        return groups.allSatisfy { $0 == "abs" } ? 1 : 0
    }

    private static func sgroupLines(for molecule: Molecule, order: OutputOrder) -> [String] {
        let ordered = orderedSgroups(for: molecule)
        guard !ordered.isEmpty else { return [] }

        var lines: [String] = []
        var parentPairs: [(Int, Int)] = []
        var idByOriginalIndex: [Int: Int] = [:]

        for (outputIndex, originalIndex) in ordered.enumerated() {
            idByOriginalIndex[originalIndex] = outputIndex + 1
        }

        for originalIndex in ordered {
            let sgroup = molecule.sgroups[originalIndex]
            let childID = idByOriginalIndex[originalIndex] ?? 0
            for child in sgroup.childGroupIndices {
                if let parentID = idByOriginalIndex[child] {
                    parentPairs.append((parentID, childID))
                }
            }
        }

        var start = 0
        while start < ordered.count {
            let end = min(ordered.count, start + 8)
            var line = "M  STY" + formatMDLInt(end - start, 3)
            for originalIndex in ordered[start..<end] {
                let id = idByOriginalIndex[originalIndex] ?? 0
                line += " " + formatMDLInt(id, 3)
                line += " " + sgroupTypeKey(for: molecule.sgroups[originalIndex])
            }
            lines.append(line)
            start = end
        }

        start = 0
        while start < parentPairs.count {
            let end = min(parentPairs.count, start + 8)
            var line = "M  SPL" + formatMDLInt(end - start, 3)
            for (child, parent) in parentPairs[start..<end] {
                line += " " + formatMDLInt(child, 3)
                line += " " + formatMDLInt(parent, 3)
            }
            lines.append(line)
            start = end
        }

        for originalIndex in ordered {
            let id = idByOriginalIndex[originalIndex] ?? 0
            let sgroup = molecule.sgroups[originalIndex]

            for atomChunk in wrap(sortedOutputIndices(for: sgroup.atomIDs, map: order.atomIndexByID), limit: 15) {
                lines.append(indexListLine(prefix: "M  SAL ", groupID: id, values: atomChunk))
            }

            for bondChunk in wrap(sortedOutputIndices(for: sgroup.crossingBondIDs, map: order.bondIndexByID), limit: 15) {
                lines.append(indexListLine(prefix: "M  SBL ", groupID: id, values: bondChunk))
            }

            if let subscriptText = sgroup.subscriptText, !subscriptText.isEmpty, sgroup.kind != .data {
                lines.append("M  SMT " + formatMDLInt(id, 3) + " " + subscriptText)
            }

            if sgroup.expanded {
                lines.append("M  SDS EXP" + formatMDLInt(1, 3) + " " + formatMDLInt(id, 3))
            }

            for bracket in sgroup.brackets {
                var line = "M  SDI " + formatMDLInt(id, 3) + formatMDLInt(4, 3)
                line += formatMDLFloat(Double(bracket.firstPoint.x))
                line += formatMDLFloat(Double(bracket.firstPoint.y))
                line += formatMDLFloat(Double(bracket.secondPoint.x))
                line += formatMDLFloat(Double(bracket.secondPoint.y))
                lines.append(line)
            }

            if sgroup.roundBrackets {
                lines.append("M  SBT" + formatMDLInt(1, 3) + " " + formatMDLInt(id, 3) + " " + formatMDLInt(1, 3))
            }

            if let connectivity = sgroup.connectivity, !connectivity.isEmpty {
                lines.append("M  SCN" + formatMDLInt(1, 3) + " " + formatMDLInt(id, 3) + " " + connectivity.uppercased())
            }

            if let subtype = sgroup.subtype, !subtype.isEmpty {
                lines.append("M  SST" + formatMDLInt(1, 3) + " " + formatMDLInt(id, 3) + " " + subtype.uppercased())
            }

            for atomChunk in wrap(sortedOutputIndices(for: sgroup.parentAtomIDs, map: order.atomIndexByID), limit: 15) {
                lines.append(indexListLine(prefix: "M  SPA ", groupID: id, values: atomChunk))
            }

            if let componentNumber = sgroup.componentNumber {
                lines.append("M  SNC" + formatMDLInt(1, 3) + " " + formatMDLInt(id, 3) + " " + formatMDLInt(componentNumber, 3))
            }

            if sgroup.kind == .data {
                lines.append(sdtLine(for: sgroup, id: id))
                let payload = (sgroup.dataValue ?? "").replacingOccurrences(of: "\r\n", with: " ")
                    .replacingOccurrences(of: "\r", with: " ")
                    .replacingOccurrences(of: "\n", with: " ")
                if payload.count <= 69 {
                    lines.append("M  SED " + formatMDLInt(id, 3) + " " + payload)
                } else {
                    var remaining = payload[...]
                    while remaining.count > 69 {
                        let end = remaining.index(remaining.startIndex, offsetBy: 69)
                        lines.append("M  SCD " + formatMDLInt(id, 3) + " " + String(remaining[..<end]))
                        remaining = remaining[end...]
                    }
                    lines.append("M  SED " + formatMDLInt(id, 3) + " " + String(remaining))
                }
            }
        }

        return lines
    }

    private static func orderedSgroups(for molecule: Molecule) -> [Int] {
        let included = molecule.sgroups.enumerated().filter { _, sgroup in
            sgroup.kind != .extMulticenter
        }
        guard !included.isEmpty else { return [] }

        var parentByChild: [Int: Int] = [:]
        for (index, sgroup) in included {
            for childIndex in sgroup.childGroupIndices {
                parentByChild[childIndex] = index
            }
        }

        var ordered: [Int] = []
        var visited = Set<Int>()

        func visit(_ index: Int) {
            guard !visited.contains(index), molecule.sgroups.indices.contains(index) else { return }
            visited.insert(index)
            ordered.append(index)
            for child in molecule.sgroups[index].childGroupIndices.sorted() {
                visit(child)
            }
        }

        for root in included.map(\.offset).filter({ parentByChild[$0] == nil }).sorted() {
            visit(root)
        }
        for index in included.map(\.offset).sorted() where !visited.contains(index) {
            visit(index)
        }

        return ordered
    }

    private static func sgroupTypeKey(for sgroup: MoleculeSgroup) -> String {
        if let keyword = sgroup.keyword?.trimmingCharacters(in: .whitespacesAndNewlines),
           !keyword.isEmpty {
            return keyword.uppercased()
        }
        switch sgroup.kind {
        case .structureRepeatUnit:
            return "SRU"
        case .polymer:
            return "MON"
        case .data:
            return "DAT"
        case .generic:
            return "SUP"
        case .extMulticenter:
            return "MUL"
        }
    }

    private static func sdtLine(for sgroup: MoleculeSgroup, id: Int) -> String {
        let fieldName = padRight(truncate(sgroup.dataFieldName ?? "DATA", to: 30), to: 30)
        let formatField = "  "
        let unitField = padRight(truncate(sgroup.dataUnit ?? "", to: 20), to: 20)
        let tagField = padRight(truncate(sgroup.dataTag ?? "", to: 2), to: 2)
        let opField = sgroup.dataOperator ?? ""
        return "M  SDT " + formatMDLInt(id, 3) + " " + fieldName + formatField + unitField + tagField + opField
    }

    private static func indexListLine(prefix: String, groupID: Int, values: [Int]) -> String {
        var line = prefix + formatMDLInt(groupID, 3) + formatMDLInt(values.count, 3)
        for value in values {
            line += " " + formatMDLInt(value, 3)
        }
        return line
    }

    private static func sortedOutputIndices(for ids: [Int], map: [Int: Int]) -> [Int] {
        ids.compactMap { map[$0] }.sorted()
    }

    private static func wrap<T>(_ values: [T], limit: Int) -> [[T]] {
        guard !values.isEmpty else { return [] }
        guard values.count > limit else { return [values] }
        var wrapped: [[T]] = []
        var start = 0
        while start < values.count {
            let end = min(values.count, start + limit)
            wrapped.append(Array(values[start..<end]))
            start = end
        }
        return wrapped
    }

    private static func defaultValence(for atom: Atom) -> Int? {
        switch canonicalElementSymbol(atom.element).uppercased() {
        case "H":
            return 1
        case "B":
            return 3
        case "C":
            return 4
        case "N":
            return atom.charge > 0 ? 4 : 3
        case "O":
            if atom.charge > 0 { return 3 }
            if atom.charge < 0 { return 1 }
            return 2
        case "P":
            return atom.charge > 0 ? 4 : 3
        case "S":
            return atom.charge > 0 ? 3 : 2
        case "F", "CL", "BR", "I":
            return 1
        case "LI", "NA", "K", "RB", "CS":
            return 1
        case "BE", "MG", "CA", "SR", "BA", "ZN":
            return 2
        case "AL":
            return 3
        case "FE":
            if let hydrogenCount = atom.explicitHydrogenCount, hydrogenCount > 0 {
                return hydrogenCount
            }
            return nil
        default:
            return nil
        }
    }

    private static func bondOrderContribution(for order: BondOrder) -> Int {
        switch order {
        case .single:
            return 1
        case .double:
            return 2
        case .triple:
            return 3
        case .aromatic:
            return 1
        }
    }

    private static func normalizedRGroupLabel(for atom: Atom) -> Int? {
        if let rGroup = atom.rGroupLabel {
            return rGroup
        }
        let trimmed = atom.element.trimmingCharacters(in: .whitespacesAndNewlines).uppercased()
        guard trimmed.hasPrefix("R"), trimmed.count > 1 else { return nil }
        return Int(trimmed.dropFirst())
    }

    private static func canonicalElementSymbol(_ raw: String) -> String {
        let trimmed = raw.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { return "" }
        if trimmed == "*" || trimmed == "R#" {
            return trimmed
        }
        guard let first = trimmed.first else { return trimmed }
        if !first.isLetter {
            return trimmed
        }
        return String(first).uppercased() + String(trimmed.dropFirst()).lowercased()
    }

    private static func isPlausibleElementSymbol(_ raw: String) -> Bool {
        let canonical = canonicalElementSymbol(raw)
        guard !canonical.isEmpty else { return false }
        if canonical == "*" || canonical == "R#" {
            return true
        }
        return atomicNumberBySymbol[canonical.uppercased()] != nil
    }

    private static func atomicNumber(forElementSymbol symbol: String) -> Int? {
        atomicNumberBySymbol[canonicalElementSymbol(symbol).uppercased()]
    }

    private static func truncate(_ text: String, to limit: Int) -> String {
        guard text.count > limit else { return text }
        return String(text.prefix(limit))
    }

    private static func padRight(_ text: String, to width: Int) -> String {
        guard text.count < width else { return text }
        return text + String(repeating: " ", count: width - text.count)
    }

    private static func formatMDLInt(_ value: Int, _ width: Int) -> String {
        let rendered = String(value)
        if rendered.count > width {
            return String(repeating: " ", count: max(0, width - 1)) + "0"
        }
        return String(repeating: " ", count: width - rendered.count) + rendered
    }

    private static func formatMDLFloat(_ value: Double) -> String {
        if value.isNaN || value.isInfinite {
            return "    0.0000"
        }
        return String(format: "%10.4f", value == -0 ? 0 : value)
    }

    private static func formatMDLString(_ value: String, _ width: Int) -> String {
        let clipped = truncate(value, to: width)
        if clipped.count >= width {
            return clipped
        }
        return clipped + String(repeating: " ", count: width - clipped.count)
    }

    private static let atomicNumberBySymbol: [String: Int] = {
        var mapping: [String: Int] = [:]
        for (index, symbol) in atomicSymbols.enumerated() where index > 0 {
            mapping[symbol.uppercased()] = index
        }
        return mapping
    }()

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
        "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
    ]
}
