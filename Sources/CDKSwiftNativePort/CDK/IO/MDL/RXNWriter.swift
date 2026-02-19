import Foundation

/// CDK-style MDL RXN writer (V2000 mol blocks).
public enum CDKRXNWriter {
    public static func write(reactants: [Molecule],
                             products: [Molecule] = [],
                             agents: [Molecule] = [],
                             reactionName: String = "CDKSwiftNativePort Reaction") throws -> String {
        let total = reactants.count + products.count + agents.count
        guard total > 0 else { throw ChemError.emptyInput }

        for molecule in reactants + products + agents where molecule.atoms.isEmpty {
            throw ChemError.parseFailed("Cannot write RXN with an empty molecule block.")
        }

        var lines: [String] = []
        lines.append("$RXN")
        lines.append(normalizedName(reactionName, fallback: "CDKSwiftNativePort Reaction"))
        lines.append("  CDKSwiftNativePort")
        lines.append("")
        lines.append(String(format: "%3d%3d%3d",
                            min(999, reactants.count),
                            min(999, products.count),
                            min(999, agents.count)))

        for molecule in reactants + products + agents {
            lines.append("$MOL")
            let mol = try CDKMDLV2000Writer.write(molecule)
                .trimmingCharacters(in: .whitespacesAndNewlines)
            lines.append(mol)
        }

        return lines.joined(separator: "\n") + "\n"
    }

    private static func normalizedName(_ raw: String, fallback: String) -> String {
        let cleaned = raw.replacingOccurrences(of: "\r", with: " ")
            .replacingOccurrences(of: "\n", with: " ")
            .trimmingCharacters(in: .whitespacesAndNewlines)
        return cleaned.isEmpty ? fallback : cleaned
    }
}
