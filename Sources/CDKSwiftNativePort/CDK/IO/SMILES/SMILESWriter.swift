import Foundation

/// CDK-style line-oriented SMILES writer.
public enum CDKSMILESWriter {
    public static func write(_ molecules: [Molecule],
                             flavor: CDKSmiFlavor = [.useAromaticSymbols, .strict]) throws -> String {
        guard !molecules.isEmpty else { throw ChemError.emptyInput }

        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor: flavor)
        var lines: [String] = []
        lines.reserveCapacity(molecules.count)

        for molecule in molecules {
            let smiles = generator.create(molecule).trimmingCharacters(in: .whitespacesAndNewlines)
            guard !smiles.isEmpty else {
                throw ChemError.parseFailed("SMILES generator produced an empty line.")
            }
            let name = normalizedName(molecule.name)
            if name.isEmpty {
                lines.append(smiles)
            } else {
                lines.append("\(smiles) \(name)")
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
