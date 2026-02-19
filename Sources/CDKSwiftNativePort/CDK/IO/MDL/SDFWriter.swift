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
            chunks.append("\(mol)\n$$$$")
        }

        return chunks.joined(separator: "\n") + "\n"
    }
}
