import Foundation

/// Lightweight Swift counterpart of CDK's `SmiFlavor` options.
public struct CDKSmiFlavor: OptionSet {
    public let rawValue: Int

    public init(rawValue: Int) {
        self.rawValue = rawValue
    }

    /// Parse and preserve aromatic atom/bond notation.
    public static let useAromaticSymbols = CDKSmiFlavor(rawValue: 1 << 0)
    /// Parse isomeric information (chiral tags, directional bond hints).
    public static let isomeric = CDKSmiFlavor(rawValue: 1 << 1)
    /// Fail on malformed branch/ring syntax instead of recovering silently.
    public static let strict = CDKSmiFlavor(rawValue: 1 << 2)
    /// Parse and apply CXSMILES trailing layers when present.
    public static let cxsmiles = CDKSmiFlavor(rawValue: 1 << 3)

    public static let cdkDefault: CDKSmiFlavor = [.useAromaticSymbols, .isomeric, .strict, .cxsmiles]
}
