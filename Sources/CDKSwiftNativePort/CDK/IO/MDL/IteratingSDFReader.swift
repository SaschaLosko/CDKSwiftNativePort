import Foundation

public enum CDKIteratingSDFReader {

    public static func readFile(url: URL) throws -> [Molecule] {
        let data = try Data(contentsOf: url)
        guard let text = String(data: data, encoding: .utf8)
                ?? String(data: data, encoding: .isoLatin1) else {
            throw ChemError.parseFailed("File is not valid text (UTF-8/Latin-1).")
        }
        return try read(text: text)
    }

    public static func read(text: String) throws -> [Molecule] {
        let normalized = text
            .replacingOccurrences(of: "\r\n", with: "\n")
            .replacingOccurrences(of: "\r", with: "\n")

        if normalized.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty {
            throw ChemError.emptyInput
        }

        var records: [[String]] = []
        var current: [String] = []

        for line in normalized.components(separatedBy: "\n") {
            if line.trimmingCharacters(in: .whitespacesAndNewlines) == "$$$$" {
                appendRecordIfMeaningful(current, to: &records)
                current.removeAll(keepingCapacity: true)
            } else {
                current.append(line)
            }
        }
        appendRecordIfMeaningful(current, to: &records)

        if records.isEmpty {
            throw ChemError.parseFailed("No SD records found.")
        }

        var molecules: [Molecule] = []
        molecules.reserveCapacity(records.count)

        for record in records {
            if let parsed = try? CDKMDLReader.read(lines: record) {
                molecules.append(parsed)
            }
        }

        if molecules.isEmpty {
            throw ChemError.parseFailed("No valid MDL records found.")
        }

        return molecules
    }

    private static func appendRecordIfMeaningful(_ record: [String], to records: inout [[String]]) {
        let joined = record.joined(separator: "\n")
        guard !joined.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else { return }
        records.append(record)
    }
}
