import Foundation

/// CDK-style MDL RXN V3000 reader.
public enum CDKRXNV3000Reader {
    public static func canParseReactionBlock(_ lines: [String]) -> Bool {
        lines.contains { line in
            let trimmed = line.trimmingCharacters(in: .whitespacesAndNewlines).uppercased()
            return trimmed == "$RXN V3000" || trimmed.hasPrefix("M  V30 COUNTS")
        }
    }

    public static func readReaction(text: String) throws -> CDKReaction {
        let reactions = try readReactions(text: text)
        guard let first = reactions.first else {
            throw ChemError.parseFailed("RXN V3000 file did not contain any reaction blocks.")
        }
        return first
    }

    public static func readReactions(text: String) throws -> [CDKReaction] {
        let normalized = text
            .replacingOccurrences(of: "\r\n", with: "\n")
            .replacingOccurrences(of: "\r", with: "\n")
        let lines = normalized.components(separatedBy: "\n")
        let rxnStarts = lines.enumerated().compactMap { index, line in
            line.trimmingCharacters(in: .whitespacesAndNewlines).uppercased().hasPrefix("$RXN") ? index : nil
        }
        guard !rxnStarts.isEmpty else {
            throw ChemError.parseFailed("RXN V3000 file is missing a $RXN header.")
        }

        var reactions: [CDKReaction] = []
        for (index, start) in rxnStarts.enumerated() {
            let end = index + 1 < rxnStarts.count ? rxnStarts[index + 1] : lines.count
            let block = Array(lines[start..<end])
            guard canParseReactionBlock(block) else {
                continue
            }
            reactions.append(try parseReactionBlock(block, reactionIndex: index + 1))
        }
        return reactions
    }

    static func parseReactionBlock(_ lines: [String], reactionIndex: Int) throws -> CDKReaction {
        guard lines.count >= 5 else {
            throw ChemError.parseFailed("RXN V3000 block is truncated.")
        }

        let reactionName = lines[1].trimmingCharacters(in: .whitespacesAndNewlines)
        let countsIndex = try countsLineIndex(in: lines)
        let counts = parseCounts(line: lines[countsIndex])

        var lineIndex = countsIndex + 1
        let reactants = try parseSection(named: "REACTANT",
                                         lines: lines,
                                         lineIndex: &lineIndex,
                                         expectedCount: counts.reactants,
                                         reactionName: reactionName,
                                         reactionIndex: reactionIndex)
        let products = try parseSection(named: "PRODUCT",
                                        lines: lines,
                                        lineIndex: &lineIndex,
                                        expectedCount: counts.products,
                                        reactionName: reactionName,
                                        reactionIndex: reactionIndex)
        let agents = try parseSection(named: "AGENT",
                                      lines: lines,
                                      lineIndex: &lineIndex,
                                      expectedCount: counts.agents,
                                      reactionName: reactionName,
                                      reactionIndex: reactionIndex)

        return CDKReaction(reactantParticipants: reactants.map { CDKReactionParticipant(molecule: $0, role: .reactant) },
                           agentParticipants: agents.map { CDKReactionParticipant(molecule: $0, role: .agent) },
                           productParticipants: products.map { CDKReactionParticipant(molecule: $0, role: .product) },
                           name: reactionName.isEmpty ? nil : reactionName)
    }

    private static func countsLineIndex(in lines: [String]) throws -> Int {
        if let index = lines.firstIndex(where: { $0.trimmingCharacters(in: .whitespacesAndNewlines).uppercased().hasPrefix("M  V30 COUNTS") }) {
            return index
        }
        throw ChemError.parseFailed("RXN V3000 block is missing an M  V30 COUNTS line.")
    }

    private static func parseCounts(line: String) -> (reactants: Int, products: Int, agents: Int) {
        let parts = line.split(whereSeparator: \.isWhitespace).map(String.init)
        guard let countsIndex = parts.firstIndex(where: { $0.uppercased() == "COUNTS" }) else {
            return (0, 0, 0)
        }
        let numeric = parts.suffix(from: parts.index(after: countsIndex)).compactMap(Int.init)
        guard numeric.count >= 2 else {
            return (0, 0, 0)
        }
        return (max(0, numeric[0]), max(0, numeric[1]), max(0, numeric.count > 2 ? numeric[2] : 0))
    }

    private static func parseSection(named sectionName: String,
                                     lines: [String],
                                     lineIndex: inout Int,
                                     expectedCount: Int,
                                     reactionName: String,
                                     reactionIndex: Int) throws -> [Molecule] {
        let beginToken = "M  V30 BEGIN \(sectionName)"
        let endToken = "M  V30 END \(sectionName)"

        guard let beginIndex = lines[lineIndex...].firstIndex(where: {
            $0.trimmingCharacters(in: .whitespacesAndNewlines).uppercased() == beginToken
        }) else {
            return []
        }

        lineIndex = beginIndex + 1
        var molecules: [Molecule] = []

        while lineIndex < lines.count {
            let trimmed = lines[lineIndex].trimmingCharacters(in: .whitespacesAndNewlines).uppercased()
            if trimmed == endToken {
                lineIndex += 1
                break
            }
            if trimmed == "M  V30 BEGIN CTAB" {
                let ctabLines = try captureCTAB(lines: lines, lineIndex: &lineIndex)
                var molecule = try parseMolecule(fromCTABLines: ctabLines,
                                                 name: roleName(sectionName,
                                                                ordinal: molecules.count + 1,
                                                                reactionName: reactionName,
                                                                reactionIndex: reactionIndex))
                if molecule.name.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty || molecule.name == "Molecule" {
                    molecule.name = roleName(sectionName,
                                             ordinal: molecules.count + 1,
                                             reactionName: reactionName,
                                             reactionIndex: reactionIndex)
                }
                molecules.append(molecule)
                continue
            }
            lineIndex += 1
        }

        if expectedCount > 0, molecules.count > expectedCount {
            return Array(molecules.prefix(expectedCount))
        }
        return molecules
    }

    private static func captureCTAB(lines: [String], lineIndex: inout Int) throws -> [String] {
        guard lineIndex < lines.count else {
            throw ChemError.parseFailed("RXN V3000 CTAB block is truncated.")
        }

        var ctabLines: [String] = []
        var depth = 0
        while lineIndex < lines.count {
            let line = lines[lineIndex]
            let trimmed = line.trimmingCharacters(in: .whitespacesAndNewlines).uppercased()
            if trimmed == "M  V30 BEGIN CTAB" {
                depth += 1
            } else if trimmed == "M  V30 END CTAB" {
                depth -= 1
            }
            ctabLines.append(line)
            lineIndex += 1
            if depth == 0 {
                return ctabLines
            }
        }
        throw ChemError.parseFailed("RXN V3000 CTAB block is unterminated.")
    }

    private static func parseMolecule(fromCTABLines ctabLines: [String], name: String) throws -> Molecule {
        let molfileLines = [
            name,
            "  CDKSwiftNativePort",
            "",
            "  0  0  0     0  0            999 V3000"
        ] + ctabLines + ["M  END"]
        return try CDKMDLReader.read(lines: molfileLines)
    }

    private static func roleName(_ sectionName: String,
                                 ordinal: Int,
                                 reactionName: String,
                                 reactionIndex: Int) -> String {
        let normalizedSection = sectionName.capitalized
        let base = reactionName.isEmpty ? "RXN \(reactionIndex)" : reactionName
        return "\(base) \(normalizedSection) \(ordinal)"
    }
}
