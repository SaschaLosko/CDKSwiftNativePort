import Foundation

/// CDK-style XYZ writer for one or more molecule blocks.
public enum CDKXYZWriter {
    public static func write(_ molecules: [Molecule]) throws -> String {
        guard !molecules.isEmpty else { throw ChemError.emptyInput }

        let blocks = try molecules.map(writeBlock(for:))
        return blocks.joined(separator: "\n")
    }

    private static func writeBlock(for molecule: Molecule) throws -> String {
        guard !molecule.atoms.isEmpty else {
            throw ChemError.parseFailed("Cannot write XYZ for an empty molecule.")
        }

        let atoms = molecule.atoms.sorted { $0.id < $1.id }
        var lines: [String] = []
        lines.reserveCapacity(atoms.count + 2)
        lines.append(String(atoms.count))
        lines.append(normalizedName(molecule.name))

        for atom in atoms {
            let symbol = CDKDescriptorSupport.canonicalElementSymbol(atom.element)
            let x = fmt(atom.position.x)
            let y = fmt(atom.position.y)
            lines.append("\(symbol) \(x) \(y) 0.00000")
        }

        return lines.joined(separator: "\n") + "\n"
    }

    private static func normalizedName(_ raw: String) -> String {
        let cleaned = raw.replacingOccurrences(of: "\r", with: " ")
            .replacingOccurrences(of: "\n", with: " ")
            .trimmingCharacters(in: .whitespacesAndNewlines)
        return cleaned.isEmpty ? "XYZ Molecule" : cleaned
    }

    private static func fmt(_ value: CGFloat) -> String {
        String(format: "%.5f", Double(value))
    }
}
