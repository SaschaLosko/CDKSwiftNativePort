import Foundation

/// Convenience facade that bundles common molecule identifiers produced by CDK-style generators.
public struct CDKMoleculeIdentifiers: Equatable {
    public let smiles: String
    public let isoSmiles: String
    public let inchi: String
    public let inchiKey: String
    public let fixedHInchi: String
    public let fixedHInchiKey: String

    public init(
        smiles: String,
        isoSmiles: String,
        inchi: String,
        inchiKey: String,
        fixedHInchi: String = "Unavailable",
        fixedHInchiKey: String = "Unavailable"
    ) {
        self.smiles = smiles
        self.isoSmiles = isoSmiles
        self.inchi = inchi
        self.inchiKey = inchiKey
        self.fixedHInchi = fixedHInchi
        self.fixedHInchiKey = fixedHInchiKey
    }
}

public enum CDKMoleculeIdentifierService {
    public static func compute(
        for molecule: Molecule,
        smilesFlavor: CDKSmiFlavor = [.useAromaticSymbols, .strict],
        isoSmilesFlavor: CDKSmiFlavor = [.useAromaticSymbols, .isomeric, .strict],
        recalculateInChI: Bool = false
    ) -> CDKMoleculeIdentifiers {
        let smilesGenerator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor: smilesFlavor)
        let isoSmilesGenerator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor: isoSmilesFlavor)

        let smiles = smilesGenerator.create(molecule)
        let isoSmiles = isoSmilesGenerator.create(molecule)

        let standardIdentifiers: (inchi: String, inchiKey: String)
        if recalculateInChI {
            standardIdentifiers = generatedInChI(for: molecule, mode: .standard)
        } else {
            let inchiGenerator = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)
            standardIdentifiers = (
                (try? inchiGenerator.getInchi()) ?? unavailableText(from: inchiGenerator.getMessage()),
                (try? inchiGenerator.getInchiKey()) ?? unavailableText(from: inchiGenerator.getMessage())
            )
        }
        let fixedHIdentifiers = generatedInChI(for: molecule, mode: .fixedH)

        return CDKMoleculeIdentifiers(
            smiles: smiles,
            isoSmiles: isoSmiles,
            inchi: standardIdentifiers.inchi,
            inchiKey: standardIdentifiers.inchiKey,
            fixedHInchi: fixedHIdentifiers.inchi,
            fixedHInchiKey: fixedHIdentifiers.inchiKey)
    }

    public static func unavailableText(from message: String) -> String {
        let trimmed = message.trimmingCharacters(in: .whitespacesAndNewlines)
        if trimmed.isEmpty {
            return "Unavailable"
        }
        return "Unavailable (\(trimmed))"
    }

    private static func generatedInChI(
        for molecule: Molecule,
        mode: CDKInChIOfficialLibraryGenerator.Mode
    ) -> (inchi: String, inchiKey: String) {
        do {
            let result = try CDKInChIOfficialLibraryGenerator.generate(for: molecule, mode: mode)
            return (result.inchi, result.inchiKey)
        } catch {
            let message = (error as? LocalizedError)?.errorDescription ?? error.localizedDescription
            let unavailable = unavailableText(from: message)
            return (unavailable, unavailable)
        }
    }
}
