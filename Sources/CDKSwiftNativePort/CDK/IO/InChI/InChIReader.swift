import Foundation

/// Line-oriented InChI reader using the CDK-style InChI parser bridge.
/// Each non-empty line may contain:
/// - `InChI=...`
/// - `InChI=... <whitespace>name`
public enum CDKInChIReader {
    public static func read(text: String) throws -> [Molecule] {
        var molecules: [Molecule] = []

        for rawLine in text.components(separatedBy: .newlines) {
            let line = rawLine.trimmingCharacters(in: .whitespacesAndNewlines)
            guard !line.isEmpty, !line.hasPrefix("#"), !line.hasPrefix("//") else {
                continue
            }

            let (inchiToken, optionalName) = split(line: line)
            let parser = CDKInChIGeneratorFactory.shared.getInChIToStructure(inchiToken)
            var molecule = try parser.getAtomContainer()
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

    private static func split(line: String) -> (inchi: String, name: String?) {
        let separators = CharacterSet.whitespacesAndNewlines
        guard let firstSeparatorRange = line.rangeOfCharacter(from: separators) else {
            return (line, nil)
        }

        let token = String(line[..<firstSeparatorRange.lowerBound])
        let trailing = String(line[firstSeparatorRange.upperBound...])
            .trimmingCharacters(in: separators)
        return (token, trailing.isEmpty ? nil : trailing)
    }
}
