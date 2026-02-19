import Foundation

public enum CDKMDLReader {

    public static func read(text: String) throws -> Molecule {
        let normalized = text
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
        return try CDKMDLV2000Reader.read(lines: trimmed)
    }

    private static func dropTrailingEmptyLines(_ lines: [String]) -> [String] {
        var out = lines
        while let last = out.last, last.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty {
            out.removeLast()
        }
        return out
    }
}
