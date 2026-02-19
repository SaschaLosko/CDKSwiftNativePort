import Foundation

/// CDK-style MDL RDF reader (extracts embedded $RXN records).
public enum CDKRDFReader {
    public static func read(text: String) throws -> [Molecule] {
        let normalized = text
            .replacingOccurrences(of: "\r\n", with: "\n")
            .replacingOccurrences(of: "\r", with: "\n")
        let lines = normalized.components(separatedBy: "\n")

        let rxnStarts = lines.enumerated().compactMap { idx, line in
            line.trimmingCharacters(in: .whitespacesAndNewlines) == "$RXN" ? idx : nil
        }

        guard !rxnStarts.isEmpty else {
            throw ChemError.parseFailed("RDF file does not contain embedded $RXN blocks.")
        }

        var molecules: [Molecule] = []
        for (index, start) in rxnStarts.enumerated() {
            let end = (index + 1 < rxnStarts.count) ? rxnStarts[index + 1] : lines.count
            let blockLines = Array(lines[start..<end])
            let blockText = blockLines.joined(separator: "\n")
            let parsed = try CDKRXNReader.read(text: blockText)
            molecules.append(contentsOf: parsed)
        }

        guard !molecules.isEmpty else {
            throw ChemError.parseFailed("No molecules could be extracted from RDF reaction blocks.")
        }
        return molecules
    }
}
