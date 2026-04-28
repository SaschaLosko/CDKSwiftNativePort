import Foundation

/// CDK-style MDL RXN reader (V2000 reaction blocks).
public enum CDKRXNReader {
    public enum Mode: String, Sendable {
        case relaxed
        case strict
    }

    public struct Options: Sendable {
        public var mode: Mode

        public init(mode: Mode = .relaxed) {
            self.mode = mode
        }
    }

    public static func readReaction(text: String,
                                    options: Options = Options()) throws -> CDKReaction {
        let reactions = try readReactions(text: text, options: options)
        guard let first = reactions.first else {
            throw ChemError.parseFailed("RXN file did not contain any reaction blocks.")
        }
        return first
    }

    public static func readReactions(text: String,
                                     options: Options = Options()) throws -> [CDKReaction] {
        let normalized = text
            .replacingOccurrences(of: "\r\n", with: "\n")
            .replacingOccurrences(of: "\r", with: "\n")
        let lines = normalized.components(separatedBy: "\n")

        let rxnStarts = lines.enumerated().compactMap { idx, line in
            line.trimmingCharacters(in: .whitespacesAndNewlines).uppercased().hasPrefix("$RXN") ? idx : nil
        }

        guard !rxnStarts.isEmpty else {
            throw ChemError.parseFailed("RXN file is missing $RXN header.")
        }

        var reactions: [CDKReaction] = []
        for (idx, start) in rxnStarts.enumerated() {
            let end = idx + 1 < rxnStarts.count ? rxnStarts[idx + 1] : lines.count
            let block = Array(lines[start..<end])
            let blockText = trimTrailingEmptyLines(block).joined(separator: "\n")
            let parsed: CDKReaction
            if CDKRXNV3000Reader.canParseReactionBlock(block) {
                parsed = try CDKRXNV3000Reader.parseReactionBlock(block, reactionIndex: idx + 1)
            } else {
                parsed = try parseReactionBlock(block, reactionIndex: idx + 1, options: options)
            }
            reactions.append(CDKRInChIReferenceSourceCache.annotating(parsed, sourceText: blockText))
        }
        return reactions
    }

    public static func read(text: String,
                            options: Options = Options()) throws -> [Molecule] {
        let reactions = try readReactions(text: text, options: options)
        let molecules = reactions.flatMap { $0.reactants + $0.products + $0.agents }

        guard !molecules.isEmpty else {
            throw ChemError.parseFailed("RXN file did not contain any molecule blocks.")
        }
        return molecules
    }

    private static func parseReactionBlock(_ lines: [String],
                                           reactionIndex: Int,
                                           options: Options) throws -> CDKReaction {
        guard lines.count >= 5 else {
            throw ChemError.parseFailed("RXN block is truncated.")
        }

        let reactionName = lines[1].trimmingCharacters(in: .whitespacesAndNewlines)
        let counts = parseCounts(line: lines[4])
        let reactantCount = counts.reactants
        let productCount = counts.products
        let declaredAgentCount = counts.agents
        let declaredCoreTotal = reactantCount + productCount

        if options.mode == .strict, counts.hasAgentExtension, declaredAgentCount > 0 {
            throw ChemError.parseFailed("RXN files uses agent count extension. This is not supported in mode STRICT")
        }

        var parsedParticipants: [(molecule: Molecule, stoichiometry: Double?)] = []
        var lineIndex = 5
        while lineIndex < lines.count {
            let trimmed = lines[lineIndex].trimmingCharacters(in: .whitespacesAndNewlines)
            guard trimmed.hasPrefix("$MOL") else {
                lineIndex += 1
                continue
            }

            let stoichiometry = parseStoichiometry(fromMolHeader: trimmed)
            lineIndex += 1
            var molBlock: [String] = []
            while lineIndex < lines.count {
                let candidate = lines[lineIndex]
                let candidateTrimmed = candidate.trimmingCharacters(in: .whitespacesAndNewlines)
                if candidateTrimmed.hasPrefix("$MOL")
                    || candidateTrimmed == "$RXN"
                    || candidateTrimmed == "$RFMT" {
                    break
                }
                molBlock.append(candidate)
                lineIndex += 1
            }

            let trimmedMol = trimTrailingEmptyLines(molBlock)
            if trimmedMol.isEmpty { continue }

            var molecule = try CDKMDLReader.read(lines: trimmedMol)
            let ordinal = parsedParticipants.count + 1
            if molecule.name.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty || molecule.name == "Molecule" {
                molecule.name = roleName(for: ordinal,
                                         reactantCount: reactantCount,
                                         productCount: productCount,
                                         reactionName: reactionName,
                                         reactionIndex: reactionIndex)
            }
            parsedParticipants.append((molecule: molecule, stoichiometry: stoichiometry))
        }

        if options.mode == .strict,
           !counts.hasAgentExtension,
           parsedParticipants.count > declaredCoreTotal {
            throw ChemError.parseFailed(
                "Agents are not supported in mode STRICT. Found \(parsedParticipants.count) molecular entities, but there are only \(declaredCoreTotal) molecular entities declared on the counts line."
            )
        }

        let effectiveAgentCount: Int
        if counts.hasAgentExtension {
            if options.mode == .relaxed,
               declaredAgentCount == 0,
               parsedParticipants.count > declaredCoreTotal {
                effectiveAgentCount = parsedParticipants.count - declaredCoreTotal
            } else {
                effectiveAgentCount = declaredAgentCount
            }
        } else if options.mode == .relaxed, parsedParticipants.count > declaredCoreTotal {
            effectiveAgentCount = parsedParticipants.count - declaredCoreTotal
        } else {
            effectiveAgentCount = 0
        }

        let expectedTotal = reactantCount + productCount + effectiveAgentCount
        if expectedTotal > 0, parsedParticipants.count > expectedTotal {
            parsedParticipants = Array(parsedParticipants.prefix(expectedTotal))
        }

        let reactantUpperBound = min(reactantCount, parsedParticipants.count)
        let productUpperBound = min(reactantUpperBound + productCount, parsedParticipants.count)
        let agentUpperBound = min(productUpperBound + effectiveAgentCount, parsedParticipants.count)

        let reactants = Array(parsedParticipants[0..<reactantUpperBound]).map { participant in
            CDKReactionParticipant(molecule: participant.molecule, role: .reactant, stoichiometry: participant.stoichiometry)
        }
        let products = Array(parsedParticipants[reactantUpperBound..<productUpperBound]).map { participant in
            CDKReactionParticipant(molecule: participant.molecule, role: .product, stoichiometry: participant.stoichiometry)
        }
        let agents = Array(parsedParticipants[productUpperBound..<agentUpperBound]).map { participant in
            CDKReactionParticipant(molecule: participant.molecule, role: .agent, stoichiometry: participant.stoichiometry)
        }

        if reactants.isEmpty, products.isEmpty, agents.isEmpty, !parsedParticipants.isEmpty {
            // Fallback for malformed counts lines: retain all parsed molecules as reactants.
            return CDKReaction(reactantParticipants: parsedParticipants.map { participant in
                CDKReactionParticipant(molecule: participant.molecule,
                                       role: .reactant,
                                       stoichiometry: participant.stoichiometry)
            }, agentParticipants: [], productParticipants: [], name: reactionName.isEmpty ? nil : reactionName)
        }
        return CDKReaction(reactantParticipants: reactants,
                           agentParticipants: agents,
                           productParticipants: products,
                           name: reactionName.isEmpty ? nil : reactionName)
    }

    private static func parseCounts(line: String) -> (reactants: Int, products: Int, agents: Int, hasAgentExtension: Bool) {
        let fixedReactants = fixedInt(line, start: 0, length: 3)
        let fixedProducts = fixedInt(line, start: 3, length: 3)
        let fixedAgents = fixedInt(line, start: 6, length: 3)

        if let fixedReactants, let fixedProducts {
            return (max(0, fixedReactants),
                    max(0, fixedProducts),
                    max(0, fixedAgents ?? 0),
                    fixedAgents != nil)
        }

        let parts = line.split(whereSeparator: \.isWhitespace).compactMap { Int($0) }
        if parts.count >= 2 {
            return (max(0, parts[0]),
                    max(0, parts[1]),
                    max(0, parts.count > 2 ? parts[2] : 0),
                    parts.count > 2)
        }

        return (0, 0, 0, false)
    }

    private static func fixedInt(_ line: String, start: Int, length: Int) -> Int? {
        guard start < line.count else { return nil }
        let lower = line.index(line.startIndex, offsetBy: start)
        let upper = line.index(lower,
                               offsetBy: min(length, line.distance(from: lower, to: line.endIndex)),
                               limitedBy: line.endIndex) ?? line.endIndex
        return Int(line[lower..<upper].trimmingCharacters(in: .whitespacesAndNewlines))
    }

    private static func roleName(for ordinal: Int,
                                 reactantCount: Int,
                                 productCount: Int,
                                 reactionName: String,
                                 reactionIndex: Int) -> String {
        let prefix: String
        let localIndex: Int
        if ordinal <= reactantCount {
            prefix = "Reactant"
            localIndex = ordinal
        } else if ordinal <= reactantCount + productCount {
            prefix = "Product"
            localIndex = ordinal - reactantCount
        } else {
            prefix = "Agent"
            localIndex = ordinal - reactantCount - productCount
        }

        let base = reactionName.isEmpty ? "RXN \(reactionIndex)" : reactionName
        return "\(base) \(prefix) \(localIndex)"
    }

    private static func trimTrailingEmptyLines(_ lines: [String]) -> [String] {
        var result = lines
        while let last = result.last,
              last.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty {
            result.removeLast()
        }
        return result
    }

    private static func parseStoichiometry(fromMolHeader line: String) -> Double? {
        let parts = line.split(whereSeparator: \.isWhitespace).map(String.init)
        guard parts.count >= 2, parts[0].uppercased() == "$MOL" else { return nil }

        if let direct = parseStoichiometryValue(parts[1]) {
            return direct
        }

        if parts.count >= 3,
           parts[1].uppercased() == "STOICH" {
            return parseStoichiometryValue(parts[2])
        }

        return nil
    }

    private static func parseStoichiometryValue(_ raw: String) -> Double? {
        let normalized = raw.replacingOccurrences(of: ",", with: ".")
        guard let value = Double(normalized), value.isFinite else { return nil }
        return value
    }
}
