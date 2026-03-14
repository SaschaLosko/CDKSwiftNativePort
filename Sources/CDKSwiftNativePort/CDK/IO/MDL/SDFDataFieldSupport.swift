import Foundation

enum CDKSDFDataFieldSupport {
    static let maxDataLineLength = 200

    static func serializeDataFields(for molecule: Molecule,
                                    acceptedFieldNames: Set<String>? = nil,
                                    truncateLongValues: Bool = false) -> String {
        var lines: [String] = []

        for fieldName in molecule.orderedDataFieldNames {
            if let acceptedFieldNames, !acceptedFieldNames.contains(fieldName) {
                continue
            }

            let serializedFieldName = sanitizeHeaderFieldName(fieldName)
            guard !serializedFieldName.isEmpty else { continue }

            lines.append("> <\(serializedFieldName)>")

            let values = molecule.dataFieldValues(named: fieldName)
            if values.isEmpty {
                lines.append("")
                continue
            }

            for value in values {
                lines.append(contentsOf: normalizedValueLines(value, truncateLongValues: truncateLongValues))
            }
            lines.append("")
        }

        return lines.joined(separator: "\n")
    }

    private static func sanitizeHeaderFieldName(_ fieldName: String) -> String {
        fieldName.replacingOccurrences(of: "[-<>.=% ]",
                                       with: "_",
                                       options: .regularExpression)
    }

    private static func normalizedValueLines(_ value: String,
                                             truncateLongValues: Bool) -> [String] {
        let normalized = value
            .replacingOccurrences(of: "\r\n", with: "\n")
            .replacingOccurrences(of: "\r", with: "\n")

        return normalized.components(separatedBy: "\n").map { line in
            guard truncateLongValues, line.count > maxDataLineLength else { return line }
            return String(line.prefix(maxDataLineLength))
        }
    }
}
