import Foundation

/// CDK-style MDL RDF reader (extracts embedded $RXN records).
public enum CDKRDFReader {
    public static func readReaction(text: String) throws -> CDKReaction {
        let reactions = try readReactions(text: text)
        guard let first = reactions.first else {
            throw ChemError.parseFailed("RDF file did not contain embedded reaction blocks.")
        }
        return first
    }

    public static func readReactions(text: String) throws -> [CDKReaction] {
        let normalized = text
            .replacingOccurrences(of: "\r\n", with: "\n")
            .replacingOccurrences(of: "\r", with: "\n")
        let lines = normalized.components(separatedBy: "\n")

        let rxnStarts = lines.enumerated().compactMap { idx, line in
            line.trimmingCharacters(in: .whitespacesAndNewlines).uppercased().hasPrefix("$RXN") ? idx : nil
        }

        guard !rxnStarts.isEmpty else {
            throw ChemError.parseFailed("RDF file does not contain embedded $RXN blocks.")
        }

        var reactions: [CDKReaction] = []
        for (index, start) in rxnStarts.enumerated() {
            let end = (index + 1 < rxnStarts.count) ? rxnStarts[index + 1] : lines.count
            let blockLines = Array(lines[start..<end])
            let blockText = blockLines.joined(separator: "\n")
            let parsed = try CDKRXNReader.readReactions(text: blockText)
            reactions.append(contentsOf: parsed)
        }
        return reactions
    }

    public static func read(text: String) throws -> [Molecule] {
        let reactions = try readReactions(text: text)
        let molecules = reactions.flatMap { $0.reactants + $0.products + $0.agents }

        guard !molecules.isEmpty else {
            throw ChemError.parseFailed("No molecules could be extracted from RDF reaction blocks.")
        }
        return molecules
    }
}
