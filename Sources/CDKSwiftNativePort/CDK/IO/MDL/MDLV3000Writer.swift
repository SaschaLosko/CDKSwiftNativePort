import Foundation

/// Swift counterpart for CDK-style MDL V3000 writing used by native import/export flows.
public enum CDKMDLV3000Writer {
    public struct Options {
        public var programName: String?
        public var includeDataFields: Bool
        public var acceptedDataFieldNames: Set<String>?
        public var truncateLongDataFields: Bool

        public init(programName: String? = nil,
                    includeDataFields: Bool = true,
                    acceptedDataFieldNames: Set<String>? = nil,
                    truncateLongDataFields: Bool = false) {
            self.programName = programName
            self.includeDataFields = includeDataFields
            self.acceptedDataFieldNames = acceptedDataFieldNames
            self.truncateLongDataFields = truncateLongDataFields
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
        guard !molecule.atoms.isEmpty else {
            throw ChemError.parseFailed("Cannot write Molfile for empty molecule.")
        }

        let order = try outputOrder(for: molecule)
        let ctabSgroups = orderCtabSgroups(for: molecule)
        let chiralStatus = chiralFlag(for: molecule)

        var lines: [String] = []
        let title = molecule.name.trimmingCharacters(in: .whitespacesAndNewlines)
        lines.append(title.isEmpty ? "Molecule" : title)
        lines.append(headerProgramLine(for: molecule, options: options))
        lines.append("")
        lines.append("  0  0  0     0  0            999 V3000")
        lines.append("M  V30 BEGIN CTAB")
        lines.append("M  V30 COUNTS \(order.atoms.count) \(order.bonds.count) \(ctabSgroups.count) 0 \(chiralStatus == 1 ? 1 : 0)")

        lines.append("M  V30 BEGIN ATOM")
        for (index, atom) in order.atoms.enumerated() {
            lines.append(contentsOf: v30WrappedLines(for: atomLine(atom,
                                                                  outputIndex: index + 1,
                                                                  in: molecule)))
        }
        lines.append("M  V30 END ATOM")

        if !order.bonds.isEmpty {
            lines.append("M  V30 BEGIN BOND")
            for (index, bond) in order.bonds.enumerated() {
                lines.append(contentsOf: v30WrappedLines(for: bondLine(bond,
                                                                      outputIndex: index + 1,
                                                                      in: molecule,
                                                                      order: order)))
            }
            lines.append("M  V30 END BOND")
        }

        let sgroupPayloads = sgroupPayloads(for: ctabSgroups, order: order)
        if !sgroupPayloads.isEmpty {
            lines.append("M  V30 BEGIN SGROUP")
            for payload in sgroupPayloads {
                lines.append(contentsOf: v30WrappedLines(for: payload))
            }
            lines.append("M  V30 END SGROUP")
        }

        let collectionPayloads = collectionPayloads(for: molecule, order: order, chiralStatus: chiralStatus)
        if !collectionPayloads.isEmpty {
            lines.append("M  V30 BEGIN COLLECTION")
            for payload in collectionPayloads {
                lines.append(contentsOf: v30WrappedLines(for: payload))
            }
            lines.append("M  V30 END COLLECTION")
        }

        lines.append("M  V30 END CTAB")
        lines.append("M  END")

        if options.includeDataFields {
            let dataFields = CDKSDFDataFieldSupport.serializeDataFields(for: molecule,
                                                                        acceptedFieldNames: options.acceptedDataFieldNames,
                                                                        truncateLongValues: options.truncateLongDataFields)
            if !dataFields.isEmpty {
                lines.append(contentsOf: dataFields.components(separatedBy: "\n"))
            }
        }

        return lines.joined(separator: "\n")
    }

    private static func outputOrder(for molecule: Molecule) throws -> OutputOrder {
        let sortedAtoms = molecule.atoms.sorted { lhs, rhs in
            let lhsIsHydrogen = lhs.element.trimmingCharacters(in: .whitespacesAndNewlines).uppercased() == "H"
            let rhsIsHydrogen = rhs.element.trimmingCharacters(in: .whitespacesAndNewlines).uppercased() == "H"
            if lhsIsHydrogen != rhsIsHydrogen { return rhsIsHydrogen }
            return lhs.id < rhs.id
        }
        let validAtomIDs = Set(sortedAtoms.map(\.id))
        let sortedBonds = molecule.bonds.sorted { $0.id < $1.id }

        for bond in sortedBonds {
            guard validAtomIDs.contains(bond.a1), validAtomIDs.contains(bond.a2) else {
                throw ChemError.parseFailed("Bond references unknown atom id while writing V3000 Molfile.")
            }
        }

        var atomIndexByID: [Int: Int] = [:]
        for (index, atom) in sortedAtoms.enumerated() {
            atomIndexByID[atom.id] = index + 1
        }
        var bondIndexByID: [Int: Int] = [:]
        for (index, bond) in sortedBonds.enumerated() {
            bondIndexByID[bond.id] = index + 1
        }

        return OutputOrder(atoms: sortedAtoms,
                           bonds: sortedBonds,
                           atomIndexByID: atomIndexByID,
                           bondIndexByID: bondIndexByID)
    }

    private static func headerProgramLine(for molecule: Molecule,
                                          options: Options) -> String {
        let formatter = DateFormatter()
        formatter.locale = Locale(identifier: "en_US_POSIX")
        formatter.dateFormat = "MMddyyHHmm"
        let programName = normalizedProgramName(options.programName)
        let dim = dimensionality(for: molecule)
        let dimSuffix: String
        switch dim {
        case 2:
            dimSuffix = "2D"
        case 3:
            dimSuffix = "3D"
        default:
            dimSuffix = ""
        }
        return "  \(programName)\(formatter.string(from: Date()))\(dimSuffix)"
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

    private static func atomLine(_ atom: Atom, outputIndex: Int, in molecule: Molecule) -> String {
        let symbol = atomSymbol(atom)
        let x = formatCoordinate(Double(atom.position.x))
        let y = formatCoordinate(Double(atom.position.y))
        let z = formatCoordinate(atom.zPosition ?? 0)
        let aamap = atom.atomMapNumber ?? atom.atomClass ?? 0

        var payload = "\(outputIndex) \(symbol) \(x) \(y) \(z) \(aamap)"

        if atom.charge != 0, (-15...15).contains(atom.charge) {
            payload += " CHG=\(atom.charge)"
        }
        if let mass = atom.isotopeMassNumber, mass != 0 {
            payload += " MASS=\(mass)"
        }
        if let radical = atom.radicalType?.rawValue ?? atom.radical, radical > 0, radical < 4 {
            payload += " RAD=\(radical)"
        }
        if let rGroup = normalizedRGroupLabel(atom) {
            payload += " RGROUPS=(1 \(rGroup))"
        }
        if let valence = serializedValence(for: atom, in: molecule) {
            payload += " VAL=\(valence)"
        }
        if let attachmentPoint = atom.attachmentPoint {
            payload += " ATTCHPT=\(attachmentPoint)"
        }
        if atom.chirality != .none {
            payload += atom.chirality == .clockwise ? " CFG=1" : " CFG=2"
        }
        if let atomList = atom.atomList, !atomList.isEmpty {
            payload += " ATOMLIST=(\(atomList.count) \(atomList.joined(separator: " ")))"
            if atom.atomListIsNegated {
                payload += " NOT=TRUE"
            }
            if let hCount = atom.explicitHydrogenCount {
                payload += " HCOUNT=\(hCount)"
            }
        }

        return payload
    }

    private static func serializedValence(for atom: Atom, in molecule: Molecule) -> Int? {
        if let override = atom.valenceOverride {
            return normalizeValence(override)
        }
        guard let hydrogenCount = atom.explicitHydrogenCount else { return nil }
        guard !molecule.bonds(forAtom: atom.id).contains(where: { $0.queryType != nil }) else { return nil }

        let explicitValence = molecule.bonds(forAtom: atom.id).reduce(0) { partial, bond in
            switch bond.order {
            case .single: return partial + 1
            case .double: return partial + 2
            case .triple: return partial + 3
            case .aromatic: return partial + 1
            }
        }

        guard let implicitValence = defaultValence(for: atom) else { return nil }
        if implicitValence - explicitValence == hydrogenCount {
            return nil
        }
        return normalizeValence(explicitValence + hydrogenCount)
    }

    private static func normalizeValence(_ value: Int) -> Int {
        if value <= 0 || value > 14 { return -1 }
        return value
    }

    private static func defaultValence(for atom: Atom) -> Int? {
        switch atom.element.trimmingCharacters(in: .whitespacesAndNewlines).uppercased() {
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
        default:
            return nil
        }
    }

    private static func bondLine(_ bond: Bond,
                                 outputIndex: Int,
                                 in molecule: Molecule,
                                 order: OutputOrder) -> String {
        let typeToken = serializedBondType(for: bond)
        let (beginID, endID) = serializedBondEndpoints(for: bond, order: order)
        var payload = "\(outputIndex) \(typeToken) \(beginID) \(endID)"

        switch bond.stereo {
        case .up, .upReversed:
            payload += " CFG=1"
        case .either:
            payload += " CFG=2"
        case .down, .downReversed:
            payload += " CFG=3"
        case .none:
            break
        }

        if let topology = bond.topology {
            payload += " TOPO=\(topology == .ring ? 1 : 2)"
        }

        if let multicenter = molecule.sgroups.first(where: { $0.kind == .extMulticenter && $0.crossingBondIDs.contains(bond.id) }) {
            let endpoints = multicenter.atomIDs
                .dropFirst()
                .compactMap { order.atomIndexByID[$0] }
                .sorted()
            if !endpoints.isEmpty {
                payload += " ATTACH=ANY ENDPTS=(\(endpoints.count)"
                for endpoint in endpoints {
                    payload += " \(endpoint)"
                }
                payload += ")"
            }
        }

        return payload
    }

    private static func serializedBondEndpoints(for bond: Bond, order: OutputOrder) -> (Int, Int) {
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

        return (order.atomIndexByID[beginAtomID] ?? beginAtomID,
                order.atomIndexByID[endAtomID] ?? endAtomID)
    }

    private static func serializedBondType(for bond: Bond) -> Int {
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

    private static func orderCtabSgroups(for molecule: Molecule) -> [(Int, MoleculeSgroup)] {
        let included = molecule.sgroups.enumerated().filter { $0.element.kind != .extMulticenter }
        guard !included.isEmpty else { return [] }

        var parentByChild: [Int: Int] = [:]
        for (index, sgroup) in included {
            for childIndex in sgroup.childGroupIndices {
                parentByChild[childIndex] = index
            }
        }

        var orderedIndices: [Int] = []
        var visited = Set<Int>()

        func visit(_ index: Int) {
            guard !visited.contains(index),
                  let sgroup = molecule.sgroups.indices.contains(index) ? molecule.sgroups[index] : nil else {
                return
            }
            visited.insert(index)
            orderedIndices.append(index)
            for child in sgroup.childGroupIndices.sorted() {
                visit(child)
            }
        }

        let roots = included.map(\.offset).filter { parentByChild[$0] == nil }.sorted()
        for root in roots {
            visit(root)
        }
        for (index, _) in included where !visited.contains(index) {
            visit(index)
        }

        return orderedIndices.compactMap { index in
            molecule.sgroups.indices.contains(index) ? (index, molecule.sgroups[index]) : nil
        }
    }

    private static func sgroupPayloads(for orderedSgroups: [(Int, MoleculeSgroup)],
                                       order: OutputOrder) -> [String] {
        guard !orderedSgroups.isEmpty else { return [] }

        var payloads: [String] = []
        var writtenIndexByOriginal: [Int: Int] = [:]

        for (writtenIndex, pair) in orderedSgroups.enumerated() {
            let (originalIndex, _) = pair
            writtenIndexByOriginal[originalIndex] = writtenIndex + 1
        }

        var parentByChild: [Int: Int] = [:]
        for (originalIndex, sgroup) in orderedSgroups {
            for child in sgroup.childGroupIndices {
                parentByChild[child] = originalIndex
            }
        }

        for (writtenIndex, pair) in orderedSgroups.enumerated() {
            let (originalIndex, sgroup) = pair
            var payload = "\(writtenIndex + 1) \(sgroupTypeKey(for: sgroup)) 0"

            let atomIDs = sgroup.atomIDs.compactMap { order.atomIndexByID[$0] }.sorted()
            if !atomIDs.isEmpty {
                payload += " ATOMS=(\(atomIDs.count)"
                for atomID in atomIDs {
                    payload += " \(atomID)"
                }
                payload += ")"
            }

            let bondIDs = sgroup.crossingBondIDs.compactMap { order.bondIndexByID[$0] }.sorted()
            if !bondIDs.isEmpty {
                let bondKey = sgroup.kind == .data ? "CBONDS" : "XBONDS"
                payload += " \(bondKey)=(\(bondIDs.count)"
                for bondID in bondIDs {
                    payload += " \(bondID)"
                }
                payload += ")"
            }

            if let parentOriginalIndex = parentByChild[originalIndex],
               let parentWrittenIndex = writtenIndexByOriginal[parentOriginalIndex] {
                payload += " PARENT=\(parentWrittenIndex)"
            }

            if let subtype = sgroup.subtype, !subtype.isEmpty {
                payload += " SUBTYPE=\(quoteIfNeeded(subtype))"
            }
            if let connectivity = sgroup.connectivity, !connectivity.isEmpty {
                payload += " CONNECT=\(quoteIfNeeded(connectivity.uppercased()))"
            }
            if sgroupTypeKey(for: sgroup) == "MUL" {
                if let multiple = sgroup.subscriptText, !multiple.isEmpty {
                    payload += " MULT=\(quoteIfNeeded(multiple))"
                }
            } else if let label = sgroup.subscriptText, !label.isEmpty {
                payload += " LABEL=\(quoteIfNeeded(label))"
            }
            if sgroup.roundBrackets {
                payload += " BRKTYP=PAREN"
            }
            if !sgroup.parentAtomIDs.isEmpty {
                let parentAtomIDs = sgroup.parentAtomIDs.compactMap { order.atomIndexByID[$0] }.sorted()
                if !parentAtomIDs.isEmpty {
                    payload += " PATOMS=(\(parentAtomIDs.count)"
                    for atomID in parentAtomIDs {
                        payload += " \(atomID)"
                    }
                    payload += ")"
                }
            }
            if let componentNumber = sgroup.componentNumber, componentNumber > 0 {
                payload += " COMPNO=\(componentNumber)"
            }
            if sgroup.expanded {
                payload += " ESTATE=E"
            }
            for bracket in sgroup.brackets {
                payload += " BRKXYZ=(9 \(formatCoordinate(Double(bracket.firstPoint.x))) \(formatCoordinate(Double(bracket.firstPoint.y))) 0"
                payload += " \(formatCoordinate(Double(bracket.secondPoint.x))) \(formatCoordinate(Double(bracket.secondPoint.y))) 0 0 0 0)"
            }
            if sgroup.kind == .data {
                if let fieldName = sgroup.dataFieldName, !fieldName.isEmpty {
                    payload += " FIELDNAME=\(quoteIfNeeded(fieldName))"
                }
                if let fieldData = sgroup.dataValue, !fieldData.isEmpty {
                    payload += " FIELDDATA=\(quoteIfNeeded(fieldData))"
                }
            }

            payloads.append(payload)
        }

        return payloads
    }

    private static func sgroupTypeKey(for sgroup: MoleculeSgroup) -> String {
        if let keyword = sgroup.keyword?.trimmingCharacters(in: .whitespacesAndNewlines),
           !keyword.isEmpty {
            return keyword.uppercased()
        }
        switch sgroup.kind {
        case .structureRepeatUnit:
            return "SRU"
        case .extMulticenter:
            return "MUL"
        case .polymer:
            return "MON"
        case .data:
            return "DAT"
        case .generic:
            return "SUP"
        }
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

        if groups.allSatisfy({ $0 == "abs" }) {
            return 1
        }
        return 2
    }

    private static func collectionPayloads(for molecule: Molecule,
                                           order: OutputOrder,
                                           chiralStatus: Int) -> [String] {
        var payloads: [String] = []

        let stereoAtoms = molecule.atoms.filter { $0.chirality != .none }
        if chiralStatus != 1 {
            let grouped = Dictionary(grouping: stereoAtoms) { atom -> Int in
                atom.cxStereoGroup ?? CDKCxSmilesParser.encodeStereoGroup(kind: "abs", number: 0)
            }

            var relativeGroupCount = 0
            var racemicGroupCount = 0
            for key in grouped.keys.sorted() {
                guard let decoded = CDKCxSmilesParser.decodeStereoGroup(key),
                      let atoms = grouped[key], !atoms.isEmpty else {
                    continue
                }
                var payload = "MDLV30/STE"
                switch decoded.kind {
                case "abs":
                    payload += "ABS"
                case "or":
                    racemicGroupCount += 1
                    payload += "RAC\(racemicGroupCount)"
                case "and":
                    relativeGroupCount += 1
                    payload += "REL\(relativeGroupCount)"
                default:
                    continue
                }
                let indices = atoms.compactMap { order.atomIndexByID[$0.id] }.sorted()
                payload += " ATOMS=(\(indices.count)"
                for index in indices {
                    payload += " \(index)"
                }
                payload += ")"
                payloads.append(payload)
            }
        }

        let highlightedAtoms = molecule.highlightedAtomIDs.compactMap { order.atomIndexByID[$0] }
        let highlightedBonds = molecule.highlightedBondIDs.compactMap { order.bondIndexByID[$0] }
        if !highlightedAtoms.isEmpty || !highlightedBonds.isEmpty {
            var payload = "MDLV30/HILITE"
            if !highlightedAtoms.isEmpty {
                payload += " ATOMS=(\(highlightedAtoms.count)"
                for atomID in highlightedAtoms.sorted() {
                    payload += " \(atomID)"
                }
                payload += ")"
            }
            if !highlightedBonds.isEmpty {
                payload += " BONDS=(\(highlightedBonds.count)"
                for bondID in highlightedBonds.sorted() {
                    payload += " \(bondID)"
                }
                payload += ")"
            }
            payloads.append(payload)
        }

        return payloads
    }

    private static func atomSymbol(_ atom: Atom) -> String {
        if atom.atomList != nil {
            return "L"
        }

        switch atom.queryType {
        case .anyAtom:
            return "ANY"
        case .anyNonHydrogen:
            return "A"
        case .anyHetero:
            return "Q"
        case .none:
            break
        }

        if normalizedRGroupLabel(atom) != nil {
            return "R#"
        }

        let trimmed = atom.element.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { return "C" }
        if atom.aromatic, trimmed.allSatisfy(\.isLetter) {
            return trimmed.lowercased()
        }
        return trimmed.count > 3 ? String(trimmed.prefix(3)) : trimmed
    }

    private static func normalizedRGroupLabel(_ atom: Atom) -> Int? {
        if let rGroup = atom.rGroupLabel {
            return rGroup
        }

        let trimmed = atom.element.trimmingCharacters(in: .whitespacesAndNewlines).uppercased()
        guard trimmed.hasPrefix("R"), trimmed.count > 1 else { return nil }
        return Int(trimmed.dropFirst())
    }

    private static func quoteIfNeeded(_ text: String) -> String {
        let trimmed = text.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { return "\"\"" }
        if trimmed.range(of: #"^[A-Za-z0-9_.+\-]+$"#, options: .regularExpression) != nil {
            return trimmed
        }
        let escaped = trimmed.replacingOccurrences(of: "\"", with: "\\\"")
        return "\"\(escaped)\""
    }

    private static func formatCoordinate(_ value: Double) -> String {
        if value == 0 { return "0" }
        var rendered = String(format: "%.4f", value)
        while rendered.contains("."), rendered.last == "0" {
            rendered.removeLast()
        }
        if rendered.last == "." {
            rendered.removeLast()
        }
        if rendered == "-0" {
            return "0"
        }
        return rendered
    }

    private static func v30WrappedLines(for payload: String) -> [String] {
        let prefix = "M  V30 "
        let payloadLimit = 72
        guard payload.count > payloadLimit else {
            return [prefix + payload]
        }

        var remaining = payload[...]
        var lines: [String] = []
        while remaining.count > payloadLimit {
            let end = remaining.index(remaining.startIndex, offsetBy: payloadLimit)
            let searchRange = remaining.startIndex..<end
            let splitIndex = remaining[searchRange].lastIndex(of: " ") ?? end
            let chunk = remaining[..<splitIndex].trimmingCharacters(in: .whitespaces)
            lines.append(prefix + chunk + " -")
            remaining = remaining[splitIndex...].trimmingCharacters(in: .whitespaces)[...]
        }
        if !remaining.isEmpty {
            lines.append(prefix + String(remaining))
        }
        return lines
    }
}
