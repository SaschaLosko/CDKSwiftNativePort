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
    public static let cxAtomLabel = CDKSmiFlavor(rawValue: 1 << 4)
    public static let cxAtomValue = CDKSmiFlavor(rawValue: 1 << 5)
    public static let cxCoordinates = CDKSmiFlavor(rawValue: 1 << 6)
    public static let cxMulticenter = CDKSmiFlavor(rawValue: 1 << 7)
    public static let cxPolymer = CDKSmiFlavor(rawValue: 1 << 8)
    public static let cxRadical = CDKSmiFlavor(rawValue: 1 << 9)
    public static let cxFragmentGroup = CDKSmiFlavor(rawValue: 1 << 10)
    public static let cxLigandOrder = CDKSmiFlavor(rawValue: 1 << 11)
    public static let cxEnhancedStereo = CDKSmiFlavor(rawValue: 1 << 12)
    public static let cxDataSgroups = CDKSmiFlavor(rawValue: 1 << 13)
    public static let cxRGroups = CDKSmiFlavor(rawValue: 1 << 14)

    public static let cxAll: CDKSmiFlavor = [
        .cxAtomLabel,
        .cxAtomValue,
        .cxCoordinates,
        .cxMulticenter,
        .cxPolymer,
        .cxRadical,
        .cxFragmentGroup,
        .cxLigandOrder,
        .cxEnhancedStereo,
        .cxDataSgroups,
        .cxRGroups
    ]

    public static let cdkDefault: CDKSmiFlavor = [.useAromaticSymbols, .isomeric, .strict, .cxsmiles]
}
