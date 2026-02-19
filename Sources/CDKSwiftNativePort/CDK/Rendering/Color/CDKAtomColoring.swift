import CoreGraphics
import Foundation

public enum CDKAtomColoringMode: String, CaseIterable, Hashable {
    case monochrome
    case cdk2D
}

public enum CDKAromaticDisplayMode: String, CaseIterable, Hashable {
    case innerLine
    case circle
}

public struct CDKRenderColor: Hashable {
    public var red: CGFloat
    public var green: CGFloat
    public var blue: CGFloat
    public var alpha: CGFloat

    public init(red: CGFloat, green: CGFloat, blue: CGFloat, alpha: CGFloat = 1.0) {
        self.red = red.clamped(to: 0...1)
        self.green = green.clamped(to: 0...1)
        self.blue = blue.clamped(to: 0...1)
        self.alpha = alpha.clamped(to: 0...1)
    }

    public func withAlpha(_ alpha: CGFloat) -> CDKRenderColor {
        CDKRenderColor(red: red, green: green, blue: blue, alpha: alpha)
    }

    public func mixed(with other: CDKRenderColor, ratio: CGFloat = 0.5) -> CDKRenderColor {
        let t = ratio.clamped(to: 0...1)
        return CDKRenderColor(red: red * (1 - t) + other.red * t,
                              green: green * (1 - t) + other.green * t,
                              blue: blue * (1 - t) + other.blue * t,
                              alpha: alpha * (1 - t) + other.alpha * t)
    }

    public func svgHexRGB() -> String {
        let r = Int((red * 255).rounded())
        let g = Int((green * 255).rounded())
        let b = Int((blue * 255).rounded())
        return String(format: "#%02X%02X%02X", r, g, b)
    }

    public static let ink = CDKRenderColor(red: 0.07, green: 0.07, blue: 0.09)
    public static let aromaticInk = CDKRenderColor(red: 0.14, green: 0.14, blue: 0.16)
    public static let grid = CDKRenderColor(red: 0.50, green: 0.50, blue: 0.54)
}

public enum CDKRenderingStyleResolver {
    public static func atomColor(for atom: Atom, style: RenderStyle) -> CDKRenderColor {
        switch style.atomColoringMode {
        case .monochrome:
            return atom.aromatic ? .aromaticInk : .ink
        case .cdk2D:
            return cdk2DPaletteColor(symbol: atom.element)
        }
    }

    public static func bondColor(for bond: Bond, molecule: Molecule, style: RenderStyle) -> CDKRenderColor {
        guard style.atomColoringMode == .cdk2D, style.colorBondsByAtom else {
            return bond.order == .aromatic ? .aromaticInk : .ink
        }
        guard let a1 = molecule.atoms.first(where: { $0.id == bond.a1 }),
              let a2 = molecule.atoms.first(where: { $0.id == bond.a2 }) else {
            return bond.order == .aromatic ? .aromaticInk : .ink
        }
        let c1 = atomColor(for: a1, style: style)
        let c2 = atomColor(for: a2, style: style)
        return c1.mixed(with: c2, ratio: 0.5)
    }

    public static func aromaticRingColor(atomIDs: [Int], molecule: Molecule, style: RenderStyle) -> CDKRenderColor {
        guard style.atomColoringMode == .cdk2D, style.colorBondsByAtom else {
            return .aromaticInk
        }
        let colors = atomIDs.compactMap { atomID in
            molecule.atoms.first(where: { $0.id == atomID }).map { atomColor(for: $0, style: style) }
        }
        guard !colors.isEmpty else { return .aromaticInk }
        let sum = colors.reduce((r: CGFloat(0), g: CGFloat(0), b: CGFloat(0), a: CGFloat(0))) { acc, color in
            (r: acc.r + color.red, g: acc.g + color.green, b: acc.b + color.blue, a: acc.a + color.alpha)
        }
        let count = CGFloat(colors.count)
        return CDKRenderColor(red: sum.r / count,
                              green: sum.g / count,
                              blue: sum.b / count,
                              alpha: sum.a / count)
    }

    private static func cdk2DPaletteColor(symbol: String) -> CDKRenderColor {
        switch symbol.uppercased() {
        case "H":
            return CDKRenderColor(red: 0.98, green: 0.98, blue: 0.98)
        case "C":
            return CDKRenderColor(red: 0.20, green: 0.20, blue: 0.22)
        case "N":
            return CDKRenderColor(red: 0.20, green: 0.31, blue: 0.94)
        case "O":
            return CDKRenderColor(red: 0.90, green: 0.14, blue: 0.14)
        case "F", "CL":
            return CDKRenderColor(red: 0.16, green: 0.64, blue: 0.22)
        case "BR":
            return CDKRenderColor(red: 0.66, green: 0.16, blue: 0.12)
        case "I":
            return CDKRenderColor(red: 0.50, green: 0.15, blue: 0.64)
        case "S":
            return CDKRenderColor(red: 0.92, green: 0.74, blue: 0.14)
        case "P":
            return CDKRenderColor(red: 0.95, green: 0.50, blue: 0.12)
        case "B":
            return CDKRenderColor(red: 0.90, green: 0.52, blue: 0.42)
        case "SI":
            return CDKRenderColor(red: 0.68, green: 0.61, blue: 0.49)
        case "NA", "K", "LI", "RB", "CS":
            return CDKRenderColor(red: 0.45, green: 0.27, blue: 0.82)
        case "MG", "CA", "SR", "BA":
            return CDKRenderColor(red: 0.18, green: 0.64, blue: 0.18)
        case "FE":
            return CDKRenderColor(red: 0.88, green: 0.40, blue: 0.20)
        case "CU":
            return CDKRenderColor(red: 0.80, green: 0.50, blue: 0.20)
        case "ZN":
            return CDKRenderColor(red: 0.50, green: 0.50, blue: 0.70)
        default:
            return .ink
        }
    }
}

private extension Comparable {
    func clamped(to range: ClosedRange<Self>) -> Self {
        min(max(self, range.lowerBound), range.upperBound)
    }
}
