import Foundation

/// Convenience facade that bundles common molecule identifiers produced by CDK-style generators.
public struct CDKMoleculeIdentifiers: Equatable {
    public let smiles: String
    public let isoSmiles: String
    public let inchi: String
    public let inchiKey: String

    public init(smiles: String, isoSmiles: String, inchi: String, inchiKey: String) {
        self.smiles = smiles
        self.isoSmiles = isoSmiles
        self.inchi = inchi
        self.inchiKey = inchiKey
    }
}

public enum CDKMoleculeIdentifierService {
    public static func compute(for molecule: Molecule,
                               smilesFlavor: CDKSmiFlavor = [.useAromaticSymbols, .strict],
                               isoSmilesFlavor: CDKSmiFlavor = [.useAromaticSymbols, .isomeric, .strict]) -> CDKMoleculeIdentifiers {
        let smilesGenerator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor: smilesFlavor)
        let isoSmilesGenerator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor: isoSmilesFlavor)

        let smiles = smilesGenerator.create(molecule)
        let isoSmiles = isoSmilesGenerator.create(molecule)

        let inchiGenerator = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)
        let inchi = (try? inchiGenerator.getInchi()) ?? unavailableText(from: inchiGenerator.getMessage())
        let inchiKey = (try? inchiGenerator.getInchiKey()) ?? unavailableText(from: inchiGenerator.getMessage())

        return CDKMoleculeIdentifiers(smiles: smiles,
                                      isoSmiles: isoSmiles,
                                      inchi: inchi,
                                      inchiKey: inchiKey)
    }

    public static func unavailableText(from message: String) -> String {
        let trimmed = message.trimmingCharacters(in: .whitespacesAndNewlines)
        if trimmed.isEmpty {
            return "Unavailable"
        }
        return "Unavailable (\(trimmed))"
    }
}
