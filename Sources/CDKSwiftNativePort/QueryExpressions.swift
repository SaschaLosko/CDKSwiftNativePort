import Foundation

/// Atom predicates use graph properties, not depiction labels. Boolean expressions preserve
/// unknown stereo results under negation; recursive queries are anchored at their first atom.
public indirect enum CDKAtomExpression: Hashable, Codable, Sendable {
    case any, none
    case element(Int), aromatic(Bool), isotope(Int), charge(Int)
    case degree(Int), heavyDegree(Int), totalDegree(Int), valence(Int), hydrogens(Int), implicitHydrogens(Int)
    case hasImplicitHydrogen, inRing, ringCount(Int), smallestRing(Int), ringBonds(Int)
    case and(CDKAtomExpression, CDKAtomExpression)
    case or(CDKAtomExpression, CDKAtomExpression)
    case not(CDKAtomExpression)
    case recursive(CDKQueryPattern)
}

public indirect enum CDKBondExpression: Hashable, Codable, Sendable {
    case any, none, order(BondOrder), inRing
    case and(CDKBondExpression, CDKBondExpression)
    case or(CDKBondExpression, CDKBondExpression)
    case not(CDKBondExpression)
}

/// Explicit, ordered coordination stereo. Configuration numbers follow OpenSMILES/CDK.
/// A focus ID in the carrier list denotes an implicit ligand, as in CDK stereo elements.
public struct CDKQueryAtomStereo: Hashable, Codable, Sendable {
    public enum Geometry: String, Codable, Sendable {
        case tetrahedral, squarePlanar, trigonalBipyramidal, octahedral
        var carrierCount: Int {
            switch self {
            case .tetrahedral, .squarePlanar: 4;
            case .trigonalBipyramidal: 5;
            case .octahedral: 6
            }
        }
        var configurationCount: Int {
            switch self {
            case .tetrahedral: 2;
            case .squarePlanar: 3;
            case .trigonalBipyramidal: 20;
            case .octahedral: 30
            }
        }
    }
    public var focus: Int
    public var carriers: [Int]
    public var geometry: Geometry
    public var configuration: Int
    public var allowsUnspecified: Bool
    /// Nonzero groups permit joint inversion of their tetrahedral centres, not independent inversion.
    public var relativeGroup: Int?

    public init(
        focus: Int, carriers: [Int], geometry: Geometry, configuration: Int,
        allowsUnspecified: Bool = false, relativeGroup: Int? = nil
    ) {
        self.focus = focus; self.carriers = carriers; self.geometry = geometry
        self.configuration = configuration; self.allowsUnspecified = allowsUnspecified
        self.relativeGroup = relativeGroup
    }
}

/// Query-only metadata is kept separate from authoritative molecule serialization.
public struct CDKQueryPattern: Hashable, Codable, Sendable {
    public var molecule: Molecule
    public var atoms: [Int: CDKAtomExpression]
    public var bonds: [Int: CDKBondExpression]
    public var stereo: [CDKQueryAtomStereo]
    public init(
        molecule: Molecule, atoms: [Int: CDKAtomExpression], bonds: [Int: CDKBondExpression],
        stereo: [CDKQueryAtomStereo] = []
    ) {
        self.molecule = molecule; self.atoms = atoms; self.bonds = bonds; self.stereo = stereo
    }
}

public enum CDKQueryError: Error, Equatable, Sendable {
    case invalid(String), unsupported(String), limit(String)
}

enum QueryTruth {
    case yes, no, unknown
    var negated: Self {
        switch self {
        case .yes: .no;
        case .no: .yes;
        case .unknown: .unknown
        }
    }
    func and(_ other: Self) -> Self {
        if self == .no || other == .no { return .no }
        return self == .yes && other == .yes ? .yes : .unknown
    }
    func or(_ other: Self) -> Self {
        if self == .yes || other == .yes { return .yes }
        return self == .no && other == .no ? .no : .unknown
    }
    init(_ value: Bool) { self = value ? .yes : .no }
}

enum QueryElements {
    static let symbols =
        ["*"]
        + "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
        .split(separator: " ").map(String.init)
    static func number(_ symbol: String) -> Int? { symbols.firstIndex(of: symbol) }
}
