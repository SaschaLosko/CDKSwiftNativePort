import Foundation
#if canImport(CoreGraphics)
import CoreGraphics
#endif

public enum CDKHighlightStyle: String, CaseIterable, Hashable {
    case none
    case colored
    case outerGlow
    case outerGlowWhiteEdge
}

public struct RenderStyle: Hashable {
    public var showCarbons: Bool = false
    public var showImplicitHydrogens: Bool = true
    public var showExplicitHydrogens: Bool = false
    public var showAtomIDs: Bool = false
    public var showAtomMapNumbers: Bool = false
    public var atomColoringMode: CDKAtomColoringMode = .monochrome
    public var highlightStyle: CDKHighlightStyle = .colored
    public var colorBondsByAtom: Bool = false
    public var aromaticDisplayMode: CDKAromaticDisplayMode = .innerLine
    public var bondWidth: CGFloat = 1.95
    public var fontSize: CGFloat = 21.0
    public var padding: CGFloat = 24.0
    public var stereoAttenuation: CGFloat = 1.0
    public var hashedStereoAttenuation: CGFloat = 1.0

    public init(showCarbons: Bool = false,
                showImplicitHydrogens: Bool = true,
                showExplicitHydrogens: Bool = false,
                showAtomIDs: Bool = false,
                showAtomMapNumbers: Bool = false,
                atomColoringMode: CDKAtomColoringMode = .monochrome,
                highlightStyle: CDKHighlightStyle = .colored,
                colorBondsByAtom: Bool = false,
                aromaticDisplayMode: CDKAromaticDisplayMode = .innerLine,
                bondWidth: CGFloat = 1.95,
                fontSize: CGFloat = 21.0,
                padding: CGFloat = 24.0,
                stereoAttenuation: CGFloat = 1.0,
                hashedStereoAttenuation: CGFloat = 1.0) {
        self.showCarbons = showCarbons
        self.showImplicitHydrogens = showImplicitHydrogens
        self.showExplicitHydrogens = showExplicitHydrogens
        self.showAtomIDs = showAtomIDs
        self.showAtomMapNumbers = showAtomMapNumbers
        self.atomColoringMode = atomColoringMode
        self.highlightStyle = highlightStyle
        self.colorBondsByAtom = colorBondsByAtom
        self.aromaticDisplayMode = aromaticDisplayMode
        self.bondWidth = bondWidth
        self.fontSize = fontSize
        self.padding = padding
        self.stereoAttenuation = stereoAttenuation
        self.hashedStereoAttenuation = hashedStereoAttenuation
    }
}

// Keep the historical symbol name available without linking SwiftUI into
// headless consumers such as Spotlight and Quick Look extensions.
public enum CDKStandardGenerator {}
