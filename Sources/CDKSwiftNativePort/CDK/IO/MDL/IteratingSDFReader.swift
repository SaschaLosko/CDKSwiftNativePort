import Foundation

public struct CDKIteratingSDFReaderOptions {
    public var skipInvalidRecords: Bool
    public var forceReadAs3DCoordinates: Bool

    public init(skipInvalidRecords: Bool = true,
                forceReadAs3DCoordinates: Bool = false) {
        self.skipInvalidRecords = skipInvalidRecords
        self.forceReadAs3DCoordinates = forceReadAs3DCoordinates
    }
}

public enum CDKIteratingSDFReader {

    public static func readFile(url: URL,
                                options: CDKIteratingSDFReaderOptions = CDKIteratingSDFReaderOptions()) throws -> [Molecule] {
        let text = try CDKFileAccess.decodeText(from: url, preferredEncodings: [.utf8, .isoLatin1])
        return try read(text: text, options: options)
    }

    public static func read(text: String,
                            options: CDKIteratingSDFReaderOptions = CDKIteratingSDFReaderOptions()) throws -> [Molecule] {
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
        var firstError: Error?

        for record in records {
            do {
                var parsed = try CDKMDLReader.read(lines: record)
                if options.forceReadAs3DCoordinates {
                    parsed = forceReadAs3DCoordinates(in: parsed)
                }
                molecules.append(parsed)
            } catch {
                firstError = firstError ?? error
                if !options.skipInvalidRecords {
                    throw error
                }
            }
        }

        if molecules.isEmpty {
            if let firstError {
                throw firstError
            }
            throw ChemError.parseFailed("No valid MDL records found.")
        }

        return molecules
    }

    private static func forceReadAs3DCoordinates(in molecule: Molecule) -> Molecule {
        var updated = molecule
        for index in updated.atoms.indices where updated.atoms[index].zPosition == nil {
            updated.atoms[index].zPosition = 0
        }
        return updated
    }

    private static func appendRecordIfMeaningful(_ record: [String], to records: inout [[String]]) {
        let joined = record.joined(separator: "\n")
        guard !joined.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else { return }
        records.append(record)
    }
}
