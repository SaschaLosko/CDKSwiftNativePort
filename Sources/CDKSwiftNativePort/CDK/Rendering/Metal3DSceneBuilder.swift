#if canImport(CoreGraphics)
import CoreGraphics
#endif
import Foundation

public struct CDKPoint3D: Hashable, Codable, Sendable {
    public var x: Double
    public var y: Double
    public var z: Double

    public init(x: Double, y: Double, z: Double) {
        self.x = x
        self.y = y
        self.z = z
    }
}

public struct CDK3DBoundingBox: Hashable, Sendable {
    public let min: CDKPoint3D
    public let max: CDKPoint3D

    public var center: CDKPoint3D {
        CDKPoint3D(
            x: (min.x + max.x) * 0.5,
            y: (min.y + max.y) * 0.5,
            z: (min.z + max.z) * 0.5)
    }

    public var size: CDKPoint3D {
        CDKPoint3D(
            x: max.x - min.x,
            y: max.y - min.y,
            z: max.z - min.z)
    }

    public var maximumExtent: Double {
        Swift.max(size.x, Swift.max(size.y, size.z))
    }

    public init(min: CDKPoint3D, max: CDKPoint3D) {
        self.min = min
        self.max = max
    }
}

public enum CDK3DAtomColorPalette: String, CaseIterable, Hashable, Sendable {
    case jmol
    case cpk
    case okabeIto
    case viridis
    case cividis
    case magma
    case inferno
    case plasma
    case colorBrewerSet2
    case colorBrewerDark2
    case cartoSafe
    case cartoVivid
    case matplotlibTab10
    case matplotlibTab20

    public var displayName: String {
        switch self {
        case .jmol:
            "Jmol"
        case .cpk:
            "CPK"
        case .okabeIto:
            "Okabe-Ito"
        case .viridis:
            "Viridis"
        case .cividis:
            "Cividis"
        case .magma:
            "Magma"
        case .inferno:
            "Inferno"
        case .plasma:
            "Plasma"
        case .colorBrewerSet2:
            "ColorBrewer Set2"
        case .colorBrewerDark2:
            "ColorBrewer Dark2"
        case .cartoSafe:
            "CARTO Safe"
        case .cartoVivid:
            "CARTO Vivid"
        case .matplotlibTab10:
            "Matplotlib Tab10"
        case .matplotlibTab20:
            "Matplotlib Tab20"
        }
    }
}

public enum CDK3DRepresentationMode: String, CaseIterable, Hashable, Sendable {
    case ballAndStick
    case spaceFilling

    public var displayName: String {
        switch self {
        case .ballAndStick:
            "Ball and Stick"
        case .spaceFilling:
            "Space Filling"
        }
    }
}

public struct CDKRenderer3DModel: Hashable, Sendable {
    public var atomRadiusScale: Double
    public var bondRadius: Double
    public var minimumAtomRadius: Double
    public var maximumAtomRadius: Double
    public var atomColoringMode: CDKAtomColoringMode
    public var atomColorPalette: CDK3DAtomColorPalette
    public var colorBondsByAtom: Bool
    public var representationMode: CDK3DRepresentationMode

    public init(
        atomRadiusScale: Double = 0.34,
        bondRadius: Double = 0.10,
        minimumAtomRadius: Double = 0.22,
        maximumAtomRadius: Double = 0.58,
        atomColoringMode: CDKAtomColoringMode = .cdk2D,
        atomColorPalette: CDK3DAtomColorPalette = .jmol,
        colorBondsByAtom: Bool = true,
        representationMode: CDK3DRepresentationMode = .ballAndStick
    ) {
        self.atomRadiusScale = atomRadiusScale
        self.bondRadius = bondRadius
        self.minimumAtomRadius = minimumAtomRadius
        self.maximumAtomRadius = maximumAtomRadius
        self.atomColoringMode = atomColoringMode
        self.atomColorPalette = atomColorPalette
        self.colorBondsByAtom = colorBondsByAtom
        self.representationMode = representationMode
    }
}

public struct CDKMetal3DScene: Hashable, Sendable {
    public struct AtomSphere: Hashable, Sendable {
        public let id: Int
        public let element: String
        public let center: CDKPoint3D
        public let radius: Double
        public let color: CDKRenderColor

        public init(id: Int, element: String, center: CDKPoint3D, radius: Double, color: CDKRenderColor) {
            self.id = id
            self.element = element
            self.center = center
            self.radius = radius
            self.color = color
        }
    }

    public struct BondCylinder: Hashable, Sendable {
        public let id: Int
        public let fromAtomID: Int
        public let toAtomID: Int
        public let from: CDKPoint3D
        public let to: CDKPoint3D
        public let radius: Double
        public let order: BondOrder
        public let color: CDKRenderColor

        public init(
            id: Int,
            fromAtomID: Int,
            toAtomID: Int,
            from: CDKPoint3D,
            to: CDKPoint3D,
            radius: Double,
            order: BondOrder,
            color: CDKRenderColor
        ) {
            self.id = id
            self.fromAtomID = fromAtomID
            self.toAtomID = toAtomID
            self.from = from
            self.to = to
            self.radius = radius
            self.order = order
            self.color = color
        }
    }

    public let atoms: [AtomSphere]
    public let bonds: [BondCylinder]
    public let boundingBox: CDK3DBoundingBox?
    public let hasExplicit3DCoordinates: Bool

    public init(
        atoms: [AtomSphere],
        bonds: [BondCylinder],
        boundingBox: CDK3DBoundingBox?,
        hasExplicit3DCoordinates: Bool
    ) {
        self.atoms = atoms
        self.bonds = bonds
        self.boundingBox = boundingBox
        self.hasExplicit3DCoordinates = hasExplicit3DCoordinates
    }
}

public enum CDKMetal3DSceneBuilder {
    public static func build(
        molecule: Molecule,
        rendererModel: CDKRenderer3DModel = CDKRenderer3DModel()
    ) -> CDKMetal3DScene {
        let hasExplicit3DCoordinates = CDKModelBuilder3D.hasMeaningful3DCoordinates(molecule)
        let sceneMolecule = hasExplicit3DCoordinates
            ? molecule
            : CDKModelBuilder3D.generate3DCoordinates(molecule: molecule)
        let style = RenderStyle(
            atomColoringMode: rendererModel.atomColoringMode,
            colorBondsByAtom: rendererModel.colorBondsByAtom)
        let atomPoints = Dictionary(uniqueKeysWithValues: sceneMolecule.atoms.map { atom in
            (atom.id, point3D(for: atom))
        })

        let atoms = sceneMolecule.atoms.map { atom in
            CDKMetal3DScene.AtomSphere(
                id: atom.id,
                element: atom.element,
                center: atomPoints[atom.id] ?? point3D(for: atom),
                radius: radius(for: atom, rendererModel: rendererModel),
                color: atomColor(for: atom, rendererModel: rendererModel, style: style))
        }

        let bonds: [CDKMetal3DScene.BondCylinder]
        switch rendererModel.representationMode {
        case .ballAndStick:
            bonds = sceneMolecule.bonds.compactMap { bond -> CDKMetal3DScene.BondCylinder? in
                guard let from = atomPoints[bond.a1],
                      let to = atomPoints[bond.a2] else {
                    return nil
                }
                return CDKMetal3DScene.BondCylinder(
                    id: bond.id,
                    fromAtomID: bond.a1,
                    toAtomID: bond.a2,
                    from: from,
                    to: to,
                    radius: max(0.025, rendererModel.bondRadius),
                    order: bond.order,
                    color: bondColor(for: bond, molecule: sceneMolecule, rendererModel: rendererModel, style: style))
            }
        case .spaceFilling:
            bonds = []
        }

        return CDKMetal3DScene(
            atoms: atoms,
            bonds: bonds,
            boundingBox: boundingBox(points: Array(atomPoints.values)),
            hasExplicit3DCoordinates: hasExplicit3DCoordinates)
    }

    private static func point3D(for atom: Atom) -> CDKPoint3D {
        CDKPoint3D(
            x: Double(atom.position.x),
            y: Double(atom.position.y),
            z: atom.zPosition ?? 0.0)
    }

    private static func radius(for atom: Atom, rendererModel: CDKRenderer3DModel) -> Double {
        let normalized = atom.element.trimmingCharacters(in: .whitespacesAndNewlines).uppercased()
        if rendererModel.representationMode == .spaceFilling {
            return vanDerWaalsRadii[normalized] ?? 1.70
        }
        let covalentRadius: Double =
            switch normalized {
            case "H": 0.31
            case "B": 0.84
            case "C": 0.76
            case "N": 0.71
            case "O": 0.66
            case "F": 0.57
            case "P": 1.07
            case "S": 1.05
            case "CL": 1.02
            case "BR": 1.20
            case "I": 1.39
            default: 0.80
            }
        let scaled = covalentRadius * rendererModel.atomRadiusScale
        return min(rendererModel.maximumAtomRadius, max(rendererModel.minimumAtomRadius, scaled))
    }

    private static func boundingBox(points: [CDKPoint3D]) -> CDK3DBoundingBox? {
        guard var minX = points.first?.x,
              var minY = points.first?.y,
              var minZ = points.first?.z,
              var maxX = points.first?.x,
              var maxY = points.first?.y,
              var maxZ = points.first?.z else {
            return nil
        }
        for point in points.dropFirst() {
            minX = min(minX, point.x)
            minY = min(minY, point.y)
            minZ = min(minZ, point.z)
            maxX = max(maxX, point.x)
            maxY = max(maxY, point.y)
            maxZ = max(maxZ, point.z)
        }
        return CDK3DBoundingBox(
            min: CDKPoint3D(x: minX, y: minY, z: minZ),
            max: CDKPoint3D(x: maxX, y: maxY, z: maxZ))
    }

    private static func atomColor(
        for atom: Atom,
        rendererModel: CDKRenderer3DModel,
        style: RenderStyle
    ) -> CDKRenderColor {
        switch rendererModel.atomColoringMode {
        case .monochrome, .atomMapHighlight:
            return CDKRenderingStyleResolver.atomColor(for: atom, style: style)
        case .cdk2D:
            return paletteColor(for: atom.element, palette: rendererModel.atomColorPalette)
        }
    }

    private static func bondColor(
        for bond: Bond,
        molecule: Molecule,
        rendererModel: CDKRenderer3DModel,
        style: RenderStyle
    ) -> CDKRenderColor {
        guard rendererModel.colorBondsByAtom,
              rendererModel.atomColoringMode != .monochrome,
              let a1 = molecule.atoms.first(where: { $0.id == bond.a1 }),
              let a2 = molecule.atoms.first(where: { $0.id == bond.a2 }) else {
            return bond.order == .aromatic ? .aromaticInk : .ink
        }
        let c1 = atomColor(for: a1, rendererModel: rendererModel, style: style)
        let c2 = atomColor(for: a2, rendererModel: rendererModel, style: style)
        return c1.mixed(with: c2, ratio: 0.5)
    }

    private static func paletteColor(for element: String, palette: CDK3DAtomColorPalette) -> CDKRenderColor {
        let symbol = element.trimmingCharacters(in: .whitespacesAndNewlines).uppercased()
        switch palette {
        case .jmol:
            return color(hex: jmolPalette[symbol] ?? 0xB31FBA)
        case .cpk:
            return color(hex: cpkPalette[symbol] ?? 0xFF1493)
        case .okabeIto:
            return categoricalPaletteColor(for: symbol, palette: okabeItoPalette, fallback: 0x000000)
        case .viridis:
            return categoricalPaletteColor(for: symbol, palette: viridisPalette, fallback: 0x440154)
        case .cividis:
            return categoricalPaletteColor(for: symbol, palette: cividisPalette, fallback: 0x00224E)
        case .magma:
            return categoricalPaletteColor(for: symbol, palette: magmaPalette, fallback: 0x000004)
        case .inferno:
            return categoricalPaletteColor(for: symbol, palette: infernoPalette, fallback: 0x000004)
        case .plasma:
            return categoricalPaletteColor(for: symbol, palette: plasmaPalette, fallback: 0x0D0887)
        case .colorBrewerSet2:
            return categoricalPaletteColor(for: symbol, palette: colorBrewerSet2Palette, fallback: 0xB3B3B3)
        case .colorBrewerDark2:
            return categoricalPaletteColor(for: symbol, palette: colorBrewerDark2Palette, fallback: 0x666666)
        case .cartoSafe:
            return categoricalPaletteColor(for: symbol, palette: cartoSafePalette, fallback: 0x888888)
        case .cartoVivid:
            return categoricalPaletteColor(for: symbol, palette: cartoVividPalette, fallback: 0xA5AA99)
        case .matplotlibTab10:
            return categoricalPaletteColor(for: symbol, palette: matplotlibTab10Palette, fallback: 0x7F7F7F)
        case .matplotlibTab20:
            return categoricalPaletteColor(for: symbol, palette: matplotlibTab20Palette, fallback: 0x7F7F7F)
        }
    }

    private static func categoricalPaletteColor(for symbol: String, palette: [UInt32], fallback: UInt32) -> CDKRenderColor {
        guard !palette.isEmpty,
              let elementIndex = elementOrder.firstIndex(of: symbol) else {
            return color(hex: fallback)
        }
        return color(hex: palette[elementIndex % palette.count])
    }

    private static func color(hex: UInt32) -> CDKRenderColor {
        CDKRenderColor(
            red: CGFloat((hex >> 16) & 0xFF) / 255.0,
            green: CGFloat((hex >> 8) & 0xFF) / 255.0,
            blue: CGFloat(hex & 0xFF) / 255.0)
    }

    private static let cpkPalette: [String: UInt32] = [
        "H": 0xFFFFFF,
        "HE": 0xFFC0CB,
        "LI": 0xB22222,
        "B": 0x00FF00,
        "C": 0xC8C8C8,
        "N": 0x8F8FFF,
        "O": 0xF00000,
        "F": 0xDAA520,
        "NA": 0x0000FF,
        "MG": 0x228B22,
        "AL": 0x808090,
        "SI": 0xDAA520,
        "P": 0xFFA500,
        "S": 0xFFC832,
        "CL": 0x00FF00,
        "CA": 0x808090,
        "TI": 0x808090,
        "CR": 0x808090,
        "MN": 0x808090,
        "FE": 0xFFA500,
        "NI": 0xA52A2A,
        "CU": 0xA52A2A,
        "ZN": 0xA52A2A,
        "BR": 0xA52A2A,
        "AG": 0x808090,
        "I": 0xA020F0,
        "BA": 0xFFA500,
        "AU": 0xDAA520,
    ]

    private static let jmolPalette: [String: UInt32] = [
        "H": 0xFFFFFF,
        "HE": 0xD9FFFF,
        "LI": 0xCC80FF,
        "BE": 0xC2FF00,
        "B": 0xFFB5B5,
        "C": 0x909090,
        "N": 0x3050F8,
        "O": 0xFF0D0D,
        "F": 0x90E050,
        "NE": 0xB3E3F5,
        "NA": 0xAB5CF2,
        "MG": 0x8AFF00,
        "AL": 0xBFA6A6,
        "SI": 0xF0C8A0,
        "P": 0xFF8000,
        "S": 0xFFFF30,
        "CL": 0x1FF01F,
        "AR": 0x80D1E3,
        "K": 0x8F40D4,
        "CA": 0x3DFF00,
        "SC": 0xE6E6E6,
        "TI": 0xBFC2C7,
        "V": 0xA6A6AB,
        "CR": 0x8A99C7,
        "MN": 0x9C7AC7,
        "FE": 0xE06633,
        "CO": 0xF090A0,
        "NI": 0x50D050,
        "CU": 0xC88033,
        "ZN": 0x7D80B0,
        "GA": 0xC28F8F,
        "GE": 0x668F8F,
        "AS": 0xBD80E3,
        "SE": 0xFFA100,
        "BR": 0xA62929,
        "KR": 0x5CB8D1,
        "RB": 0x702EB0,
        "SR": 0x00FF00,
        "Y": 0x94FFFF,
        "ZR": 0x94E0E0,
        "NB": 0x73C2C9,
        "MO": 0x54B5B5,
        "TC": 0x3B9E9E,
        "RU": 0x248F8F,
        "RH": 0x0A7D8C,
        "PD": 0x006985,
        "AG": 0xC0C0C0,
        "CD": 0xFFD98F,
        "IN": 0xA67573,
        "SN": 0x668080,
        "SB": 0x9E63B5,
        "TE": 0xD47A00,
        "I": 0x940094,
        "XE": 0x429EB0,
        "CS": 0x57178F,
        "BA": 0x00C900,
        "LA": 0x70D4FF,
        "CE": 0xFFFFC7,
        "PR": 0xD9FFC7,
        "ND": 0xC7FFC7,
        "PM": 0xA3FFC7,
        "SM": 0x8FFFC7,
        "EU": 0x61FFC7,
        "GD": 0x45FFC7,
        "TB": 0x30FFC7,
        "DY": 0x1FFFC7,
        "HO": 0x00FF9C,
        "ER": 0x00E675,
        "TM": 0x00D452,
        "YB": 0x00BF38,
        "LU": 0x00AB24,
        "HF": 0x4DC2FF,
        "TA": 0x4DA6FF,
        "W": 0x2194D6,
        "RE": 0x267DAB,
        "OS": 0x266696,
        "IR": 0x175487,
        "PT": 0xD0D0E0,
        "AU": 0xFFD123,
        "HG": 0xB8B8D0,
        "TL": 0xA6544D,
        "PB": 0x575961,
        "BI": 0x9E4FB5,
        "PO": 0xAB5C00,
        "AT": 0x754F45,
        "RN": 0x428296,
        "FR": 0x420066,
        "RA": 0x007D00,
        "AC": 0x70ABFA,
        "TH": 0x00BAFF,
        "PA": 0x00A1FF,
        "U": 0x008FFF,
        "NP": 0x0080FF,
        "PU": 0x006BFF,
        "AM": 0x545CF2,
        "CM": 0x785CE3,
        "BK": 0x8A4FE3,
        "CF": 0xA136D4,
        "ES": 0xB31FD4,
        "FM": 0xB31FBA,
        "MD": 0xB30DA6,
        "NO": 0xBD0D87,
        "LR": 0xC70066,
        "RF": 0xCC0059,
        "DB": 0xD1004F,
        "SG": 0xD90045,
        "BH": 0xE00038,
        "HS": 0xE6002E,
        "MT": 0xEB0026,
    ]

    private static let elementOrder = [
        "H", "HE", "LI", "BE", "B", "C", "N", "O", "F", "NE",
        "NA", "MG", "AL", "SI", "P", "S", "CL", "AR", "K", "CA",
        "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN",
        "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y", "ZR",
        "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN",
        "SB", "TE", "I", "XE", "CS", "BA", "LA", "CE", "PR", "ND",
        "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB",
        "LU", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG",
        "TL", "PB", "BI", "PO", "AT", "RN", "FR", "RA", "AC", "TH",
        "PA", "U", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM",
        "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS",
        "RG", "CN", "NH", "FL", "MC", "LV", "TS", "OG",
    ]

    private static let okabeItoPalette: [UInt32] = [
        0xE69F00, 0x56B4E9, 0x009E73, 0xF0E442,
        0x0072B2, 0xD55E00, 0xCC79A7, 0x000000,
    ]

    private static let viridisPalette: [UInt32] = [
        0x440154, 0x482878, 0x3E4A89, 0x31688E, 0x26828E,
        0x1F9E89, 0x35B779, 0x6CCD5A, 0xB5DE2B, 0xFDE725,
    ]

    private static let cividisPalette: [UInt32] = [
        0x00224E, 0x123570, 0x3C4A6C, 0x575D6D, 0x707173,
        0x8A8678, 0xA59C74, 0xC2B369, 0xE1CC55, 0xFEE838,
    ]

    private static let magmaPalette: [UInt32] = [
        0x000004, 0x180F3D, 0x451077, 0x721F81, 0x9E2F7F,
        0xCD4071, 0xF1605D, 0xFD9467, 0xFECA8D, 0xFCFDBF,
    ]

    private static let infernoPalette: [UInt32] = [
        0x000004, 0x1B0C41, 0x4C0C6B, 0x781C6D, 0xA52C60,
        0xCF4446, 0xED6925, 0xFB9906, 0xF7D13D, 0xFCFFA4,
    ]

    private static let plasmaPalette: [UInt32] = [
        0x0D0887, 0x46039F, 0x7401A8, 0x9C179E, 0xBD3786,
        0xD8576B, 0xED7953, 0xFA9E3B, 0xFDCA26, 0xF0F921,
    ]

    private static let colorBrewerSet2Palette: [UInt32] = [
        0x66C2A5, 0xFC8D62, 0x8DA0CB, 0xE78AC3,
        0xA6D854, 0xFFD92F, 0xE5C494, 0xB3B3B3,
    ]

    private static let colorBrewerDark2Palette: [UInt32] = [
        0x1B9E77, 0xD95F02, 0x7570B3, 0xE7298A,
        0x66A61E, 0xE6AB02, 0xA6761D, 0x666666,
    ]

    private static let cartoSafePalette: [UInt32] = [
        0x88CCEE, 0xCC6677, 0xDDCC77, 0x117733, 0x332288,
        0xAA4499, 0x44AA99, 0x999933, 0x882255, 0x661100,
    ]

    private static let cartoVividPalette: [UInt32] = [
        0xE58606, 0x5D69B1, 0x52BCA3, 0x99C945, 0xCC61B0,
        0x24796C, 0xDAA51B, 0x2F8AC4, 0x764E9F, 0xED645A,
    ]

    private static let matplotlibTab10Palette: [UInt32] = [
        0x1F77B4, 0xFF7F0E, 0x2CA02C, 0xD62728, 0x9467BD,
        0x8C564B, 0xE377C2, 0x7F7F7F, 0xBCBD22, 0x17BECF,
    ]

    private static let matplotlibTab20Palette: [UInt32] = [
        0x1F77B4, 0xAEC7E8, 0xFF7F0E, 0xFFBB78, 0x2CA02C,
        0x98DF8A, 0xD62728, 0xFF9896, 0x9467BD, 0xC5B0D5,
        0x8C564B, 0xC49C94, 0xE377C2, 0xF7B6D2, 0x7F7F7F,
        0xC7C7C7, 0xBCBD22, 0xDBDB8D, 0x17BECF, 0x9EDAE5,
    ]

    private static let vanDerWaalsRadii: [String: Double] = [
        "H": 1.20,
        "HE": 1.40,
        "LI": 2.20,
        "BE": 1.90,
        "B": 1.80,
        "C": 1.70,
        "N": 1.60,
        "O": 1.55,
        "F": 1.50,
        "NE": 1.54,
        "NA": 2.40,
        "MG": 2.20,
        "AL": 2.10,
        "SI": 2.10,
        "P": 1.95,
        "S": 1.80,
        "CL": 1.80,
        "AR": 1.88,
        "K": 2.80,
        "CA": 2.40,
        "SC": 2.30,
        "TI": 2.15,
        "V": 2.05,
        "CR": 2.05,
        "MN": 2.05,
        "FE": 2.05,
        "ZN": 2.10,
        "GA": 2.10,
        "GE": 2.10,
        "AS": 2.05,
        "SE": 1.90,
        "BR": 1.90,
        "KR": 2.02,
        "RB": 2.90,
        "SR": 2.55,
        "Y": 2.40,
        "ZR": 2.30,
        "NB": 2.15,
        "MO": 2.10,
        "TC": 2.05,
        "RU": 2.05,
        "PD": 2.05,
        "AG": 2.10,
        "CD": 2.20,
        "IN": 2.20,
        "SN": 2.25,
        "SB": 2.20,
        "TE": 2.10,
        "I": 2.10,
        "XE": 2.16,
        "CS": 3.00,
        "BA": 2.70,
        "LA": 2.50,
        "CE": 2.48,
        "PR": 2.47,
        "ND": 2.45,
        "PM": 2.43,
        "SM": 2.42,
        "EU": 2.40,
        "GD": 2.38,
        "TB": 2.37,
        "DY": 2.35,
        "HO": 2.33,
        "ER": 2.32,
        "TM": 2.30,
        "YB": 2.28,
        "LU": 2.27,
        "HF": 2.25,
        "TA": 2.20,
        "W": 2.10,
        "RE": 2.05,
        "PT": 2.05,
        "AU": 2.10,
        "HG": 2.05,
        "TL": 2.20,
        "PB": 2.30,
        "BI": 2.30,
        "TH": 2.40,
        "U": 2.30,
    ]
}
