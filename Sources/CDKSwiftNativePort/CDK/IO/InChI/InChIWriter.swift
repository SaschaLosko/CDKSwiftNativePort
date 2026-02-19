import Foundation

/// Line-oriented InChI writer mirroring CDK's text-based writer behavior.
public enum CDKInChIWriter {
    public static func write(_ molecules: [Molecule]) throws -> String {
        guard !molecules.isEmpty else { throw ChemError.emptyInput }

        var lines: [String] = []
        lines.reserveCapacity(molecules.count)

        for molecule in molecules {
            let generator = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)
            let inchi: String
            do {
                inchi = try generator.getInchi().trimmingCharacters(in: .whitespacesAndNewlines)
            } catch {
                let message = generator.getMessage().trimmingCharacters(in: .whitespacesAndNewlines)
                if message.isEmpty {
                    throw error
                }
                throw ChemError.unsupported("InChI export failed: \(message)")
            }

            guard !inchi.isEmpty else {
                throw ChemError.parseFailed("InChI generator produced an empty line.")
            }
            let name = normalizedName(molecule.name)
            if name.isEmpty {
                lines.append(inchi)
            } else {
                lines.append("\(inchi) \(name)")
            }
        }

        return lines.joined(separator: "\n") + "\n"
    }

    private static func normalizedName(_ raw: String) -> String {
        raw.replacingOccurrences(of: "\r", with: " ")
            .replacingOccurrences(of: "\n", with: " ")
            .trimmingCharacters(in: .whitespacesAndNewlines)
    }
}
