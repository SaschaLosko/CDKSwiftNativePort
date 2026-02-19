import Foundation

/// CDK-style MDL RXN reader (V2000 reaction blocks).
public enum CDKRXNReader {
    public static func read(text: String) throws -> [Molecule] {
        let normalized = text
            .replacingOccurrences(of: "\r\n", with: "\n")
            .replacingOccurrences(of: "\r", with: "\n")
        let lines = normalized.components(separatedBy: "\n")

        let rxnStarts = lines.enumerated().compactMap { idx, line in
            line.trimmingCharacters(in: .whitespacesAndNewlines) == "$RXN" ? idx : nil
        }

        guard !rxnStarts.isEmpty else {
            throw ChemError.parseFailed("RXN file is missing $RXN header.")
        }

        var molecules: [Molecule] = []
        for (idx, start) in rxnStarts.enumerated() {
            let end = idx + 1 < rxnStarts.count ? rxnStarts[idx + 1] : lines.count
            let block = Array(lines[start..<end])
            molecules.append(contentsOf: try parseReactionBlock(block, reactionIndex: idx + 1))
        }

        guard !molecules.isEmpty else {
            throw ChemError.parseFailed("RXN file did not contain any molecule blocks.")
        }
        return molecules
    }

    private static func parseReactionBlock(_ lines: [String], reactionIndex: Int) throws -> [Molecule] {
        guard lines.count >= 5 else {
            throw ChemError.parseFailed("RXN block is truncated.")
        }

        let reactionName = lines[1].trimmingCharacters(in: .whitespacesAndNewlines)
        let counts = parseCounts(line: lines[4])
        let reactantCount = counts.reactants
        let productCount = counts.products
        let agentCount = counts.agents
        let expectedTotal = reactantCount + productCount + agentCount

        var molecules: [Molecule] = []
        var lineIndex = 5
        while lineIndex < lines.count {
            let trimmed = lines[lineIndex].trimmingCharacters(in: .whitespacesAndNewlines)
            guard trimmed == "$MOL" else {
                lineIndex += 1
                continue
            }

            lineIndex += 1
            var molBlock: [String] = []
            while lineIndex < lines.count {
                let candidate = lines[lineIndex]
                let candidateTrimmed = candidate.trimmingCharacters(in: .whitespacesAndNewlines)
                if candidateTrimmed == "$MOL" || candidateTrimmed == "$RXN" || candidateTrimmed == "$RFMT" {
                    break
                }
                molBlock.append(candidate)
                lineIndex += 1
            }

            let trimmedMol = trimTrailingEmptyLines(molBlock)
            if trimmedMol.isEmpty { continue }

            var molecule = try CDKMDLReader.read(lines: trimmedMol)
            let ordinal = molecules.count + 1
            if molecule.name.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty || molecule.name == "Molecule" {
                molecule.name = roleName(for: ordinal,
                                         reactantCount: reactantCount,
                                         productCount: productCount,
                                         reactionName: reactionName,
                                         reactionIndex: reactionIndex)
            }
            molecules.append(molecule)
        }

        if expectedTotal > 0, molecules.count > expectedTotal {
            molecules = Array(molecules.prefix(expectedTotal))
        }

        return molecules
    }

    private static func parseCounts(line: String) -> (reactants: Int, products: Int, agents: Int) {
        let fixedReactants = fixedInt(line, start: 0, length: 3)
        let fixedProducts = fixedInt(line, start: 3, length: 3)
        let fixedAgents = fixedInt(line, start: 6, length: 3)

        if let fixedReactants, let fixedProducts {
            return (max(0, fixedReactants),
                    max(0, fixedProducts),
                    max(0, fixedAgents ?? 0))
        }

        let parts = line.split(whereSeparator: \.isWhitespace).compactMap { Int($0) }
        if parts.count >= 2 {
            return (max(0, parts[0]),
                    max(0, parts[1]),
                    max(0, parts.count > 2 ? parts[2] : 0))
        }

        return (0, 0, 0)
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
}
