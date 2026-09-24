import Foundation

/// Compatibility entry points only. No ChemDraw encoder is included in this package.
/// Applications must supply their own format extension.
public enum CDKChemDrawWriter {
    public static func cdxml(molecules: [Molecule]) throws -> String { throw unavailable() }
    public static func cdx(molecules: [Molecule]) throws -> Data { throw unavailable() }
    public static func cdxml(reactions: [CDKReaction]) throws -> String { throw unavailable() }
    public static func cdx(reactions: [CDKReaction]) throws -> Data { throw unavailable() }
    private static func unavailable() -> ChemError {
        .unsupported("ChemDraw export requires an external format extension.")
    }
}
