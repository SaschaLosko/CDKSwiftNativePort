import Foundation

/// CDK-style SDFile writer built on the V2000 mol writer.
public enum CDKSDFWriter {
    public static func write(_ molecules: [Molecule]) throws -> String {
        guard !molecules.isEmpty else { throw ChemError.emptyInput }

        var chunks: [String] = []
        chunks.reserveCapacity(molecules.count)
        for molecule in molecules {
            let mol = try CDKMDLV2000Writer.write(molecule)
                .trimmingCharacters(in: .whitespacesAndNewlines)
            let dataFields = serializedDataFields(for: molecule)
            let chunk = dataFields.isEmpty
                ? "\(mol)\n$$$$"
                : "\(mol)\n\(dataFields)\n$$$$"
            chunks.append(chunk)
        }

        return chunks.joined(separator: "\n") + "\n"
    }

    private static func serializedDataFields(for molecule: Molecule) -> String {
        var lines: [String] = []

        for fieldName in molecule.orderedDataFieldNames {
            lines.append(">  <\(fieldName)>")
            let values = molecule.dataFieldValues(named: fieldName)
            if values.isEmpty {
                lines.append("")
                continue
            }

            for value in values {
                let normalized = value
                    .replacingOccurrences(of: "\r\n", with: "\n")
                    .replacingOccurrences(of: "\r", with: "\n")
                lines.append(contentsOf: normalized.components(separatedBy: "\n"))
            }
            lines.append("")
        }

        return lines.joined(separator: "\n")
    }
}
