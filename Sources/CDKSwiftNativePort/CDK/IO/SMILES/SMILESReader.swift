import Foundation

/// Swift counterpart of CDK's line-oriented SMILES reader behavior.
/// Accepts one SMILES entry per non-empty line, with an optional trailing name.
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

            let (smilesToken, optionalName) = split(line: line)
            var molecule = try parser.parseSmiles(smilesToken)
            if let name = optionalName, !name.isEmpty {
                molecule.name = name
            }
            molecules.append(molecule)
        }

        guard !molecules.isEmpty else {
            throw ChemError.emptyInput
        }
        return molecules
    }

    private static func split(line: String) -> (smiles: String, name: String?) {
        let separators = CharacterSet.whitespacesAndNewlines
        guard let firstSeparatorRange = line.rangeOfCharacter(from: separators) else {
            return (line, nil)
        }

        let smiles = String(line[..<firstSeparatorRange.lowerBound])
        let trailing = String(line[firstSeparatorRange.upperBound...])
            .trimmingCharacters(in: separators)
        return (smiles, trailing.isEmpty ? nil : trailing)
    }
}
