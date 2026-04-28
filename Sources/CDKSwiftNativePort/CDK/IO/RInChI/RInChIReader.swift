import Foundation

/// Line-oriented RInChI reader.
/// Each non-empty line may contain:
/// - `RInChI=...`
/// - `RInChI=... <whitespace>name`
public enum CDKRInChIReader {
    public static func read(text: String) throws -> [CDKReaction] {
        let records = try parseRecords(text: text)
        guard !records.isEmpty else {
            throw ChemError.emptyInput
        }
        return records.map(\.reaction)
    }

    static func parseRecords(text: String) throws -> [(reaction: CDKReaction, rinchi: String)] {
        var records: [(reaction: CDKReaction, rinchi: String)] = []

        for rawLine in text.components(separatedBy: .newlines) {
            let line = rawLine.trimmingCharacters(in: .whitespacesAndNewlines)
            guard !line.isEmpty, !line.hasPrefix("#"), !line.hasPrefix("//") else {
                continue
            }
            guard line.hasPrefix("RInChI=") else {
                continue
            }

            let (rinchiToken, optionalName) = split(line: line)
            let parser = CDKRInChIToReaction(rinchi: rinchiToken)
            guard parser.getStatus() != .error,
                  var reaction = parser.getReaction() else {
                let message = parser.getMessages().joined(separator: " ")
                throw ChemError.parseFailed(message.isEmpty ? "Unable to parse RInChI reaction." : message)
            }

            if let name = optionalName, !name.isEmpty {
                reaction.name = name
            }
            reaction = CDKRInChIReferenceSourceCache.annotating(reaction, sourceText: rinchiToken)
            records.append((reaction: reaction, rinchi: rinchiToken))
        }

        return records
    }

    private static func split(line: String) -> (rinchi: String, name: String?) {
        let separators = CharacterSet.whitespacesAndNewlines
        guard let firstSeparatorRange = line.rangeOfCharacter(from: separators) else {
            return (line, nil)
        }

        let token = String(line[..<firstSeparatorRange.lowerBound])
        let trailing = String(line[firstSeparatorRange.upperBound...])
            .trimmingCharacters(in: separators)
        return (token, trailing.isEmpty ? nil : trailing)
    }
}
