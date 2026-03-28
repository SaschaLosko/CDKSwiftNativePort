import Foundation

/// Writer for MDL RGfiles using the advanced R-group logic model from upstream CDK.
public enum CDKRGFileWriter {
    public static func canWrite(_ molecule: Molecule) -> Bool {
        guard let query = try? CDKRGroupQueryManipulator.toRGroupQuery(molecule) else { return false }
        return canWrite(query)
    }

    public static func canWrite(_ query: CDKRGroupQuery) -> Bool {
        !query.rootStructure.atoms.isEmpty
            && !query.rGroupDefinitions.isEmpty
            && query.areRootAtomsDefined()
            && query.areSubstituentsDefined()
    }

    public static func write(_ molecule: Molecule) throws -> String {
        try write(CDKRGroupQueryManipulator.toRGroupQuery(molecule))
    }

    public static func write(_ query: CDKRGroupQuery) throws -> String {
        guard canWrite(query) else {
            throw ChemError.unsupported("RGfile export requires a root structure with R-group query atoms and matching substituent definitions.")
        }

        let timestamp = currentTimestamp()
        var lines: [String] = [
            "$MDL  REV  1   \(timestamp)",
            "$MOL",
            "$HDR",
            "  Rgroup query file (RGFile)",
            "  CDK    \(timestamp)2D",
            "",
            "$END HDR",
            "$CTAB"
        ]

        var rootLines = try ctabBody(for: query.rootStructure)
        guard let mendIndex = rootLines.lastIndex(where: { $0.trimmingCharacters(in: .whitespaces) == "M  END" }) else {
            throw ChemError.parseFailed("Unable to locate M  END while building RGfile root CTAB.")
        }

        let insertionIndex = mendIndex
        rootLines.insert(contentsOf: logLines(for: query), at: insertionIndex)
        rootLines.insert(contentsOf: aalLines(for: query.rootStructure), at: insertionIndex + logLines(for: query).count)
        lines.append(contentsOf: rootLines)
        lines.append("$END CTAB")

        for rGroupNumber in query.rGroupDefinitions.keys.sorted() {
            guard let definition = query.rGroupDefinitions[rGroupNumber],
                  !definition.rGroups.isEmpty else { continue }

            lines.append("$RGP")
            lines.append(formatMDLInt(rGroupNumber, width: 4))

            for rGroup in definition.rGroups {
                lines.append("$CTAB")
                lines.append(contentsOf: try ctabBody(for: rGroup.group))
                lines.append("$END CTAB")
            }

            lines.append("$END RGP")
        }

        lines.append("$END MOL")
        return lines.joined(separator: "\n") + "\n"
    }

    private static func currentTimestamp() -> String {
        let formatter = DateFormatter()
        formatter.locale = Locale(identifier: "en_US_POSIX")
        formatter.dateFormat = "MMddyyHHmm"
        return formatter.string(from: Date())
    }

    private static func ctabBody(for molecule: Molecule) throws -> [String] {
        let molfile = try CDKMDLV2000Writer.write(molecule)
        let lines = molfile
            .replacingOccurrences(of: "\r\n", with: "\n")
            .replacingOccurrences(of: "\r", with: "\n")
            .components(separatedBy: "\n")
        guard lines.count >= 4 else {
            throw ChemError.parseFailed("Unexpected V2000 molfile while creating RGfile.")
        }
        return Array(lines.dropFirst(3)).filter { !($0.isEmpty && $0 == lines.last) }
    }

    private static func logLines(for query: CDKRGroupQuery) -> [String] {
        query.rGroupDefinitions.keys.sorted().compactMap { rGroupNumber in
            guard let definition = query.rGroupDefinitions[rGroupNumber] else { return nil }
            let restH = definition.restH ? 1 : 0
            return "M  LOG"
                + formatMDLInt(1, width: 3)
                + formatMDLInt(rGroupNumber, width: 4)
                + formatMDLInt(definition.requiredRGroupNumber, width: 4)
                + formatMDLInt(restH, width: 4)
                + "   "
                + definition.occurrence
        }
    }

    private static func aalLines(for molecule: Molecule) -> [String] {
        let atomOrder = Dictionary(uniqueKeysWithValues: molecule.atoms.enumerated().map { ($0.element.id, $0.offset + 1) })
        var lines: [String] = []

        for atom in molecule.atoms where atom.rGroupLabel != nil {
            let neighbors = molecule.neighbors(of: atom.id)
            guard neighbors.count > 1,
                  let rootPosition = atomOrder[atom.id] else { continue }

            let naturalOrder = neighbors.sorted { lhs, rhs in
                (atomOrder[lhs] ?? .max) < (atomOrder[rhs] ?? .max)
            }
            let orderedPartners = sanitizeLigandOrder(atom.ligandOrderingAtomIDs,
                                                      neighbors: naturalOrder)
            guard orderedPartners != naturalOrder else { continue }

            var payload = "M  AAL"
                + formatMDLInt(rootPosition, width: 4)
                + formatMDLInt(orderedPartners.count, width: 3)
            for (index, partnerID) in orderedPartners.enumerated() {
                guard let partnerPosition = atomOrder[partnerID] else { continue }
                payload += formatMDLInt(partnerPosition, width: 4)
                payload += formatMDLInt(index + 1, width: 4)
            }
            lines.append(payload)
        }

        return lines
    }

    private static func sanitizeLigandOrder(_ requested: [Int]?, neighbors: [Int]) -> [Int] {
        guard let requested, !requested.isEmpty else { return neighbors }
        var ordered: [Int] = []
        var seen = Set<Int>()

        for atomID in requested where neighbors.contains(atomID) && seen.insert(atomID).inserted {
            ordered.append(atomID)
        }
        for atomID in neighbors where seen.insert(atomID).inserted {
            ordered.append(atomID)
        }
        return ordered
    }

    private static func formatMDLInt(_ value: Int, width: Int) -> String {
        let text = String(value)
        if text.count >= width { return text }
        return String(repeating: " ", count: width - text.count) + text
    }
}
