import Foundation

/// Swift counterpart of CDK's line-oriented SMILES reader behavior.
/// Accepts one SMILES entry per non-empty line and lets `SmilesParser`
/// consume any trailing CXSMILES/title suffix on the full line.
public enum CDKSMILESReader {
    public static func read(text: String,
                            flavor: CDKSmiFlavor = .cdkDefault) throws -> [Molecule] {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: flavor)
        var molecules: [Molecule] = []

        for rawLine in text.components(separatedBy: .newlines) {
            let line = rawLine.trimmingCharacters(in: .whitespacesAndNewlines)
            guard !line.isEmpty, !line.hasPrefix("#"), !line.hasPrefix("//") else {
                continue
            }

            let molecule = try parser.parseSmiles(line)
            molecules.append(molecule)
        }

        guard !molecules.isEmpty else {
            throw ChemError.emptyInput
        }
        return molecules
    }
}
