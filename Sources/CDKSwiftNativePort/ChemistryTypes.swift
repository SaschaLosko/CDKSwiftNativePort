import Foundation

#if canImport(CoreGraphics)
    import CoreGraphics
#endif

public enum ChemFormat: String, CaseIterable, Identifiable, Sendable {
    case sdf = "SDF / MOL"
    case smiles = "SMILES"
    case inchi = "InChI"

    public var id: String { rawValue }
}

public enum BondOrder: Int, Codable, Hashable, Sendable {
    case single = 1
    case double = 2
    case triple = 3
    case aromatic = 4

    public var displayName: String {
        switch self {
        case .single: "Single"
        case .double: "Double"
        case .triple: "Triple"
        case .aromatic: "Aromatic"
        }
    }
}

extension BondOrder {
    /// Contribution used for simple valence / implicit-hydrogen estimation.
    var valenceContribution: Double {
        switch self {
        case .single: 1.0
        case .double: 2.0
        case .triple: 3.0
        case .aromatic: 1.5
        }
    }
}

public enum BondStereo: Codable, Hashable, Sendable {
    case none
    case up
    case down
    case either
    // Used for inferred wedges when the stereocenter is at atom a2.
    case upReversed
    case downReversed

    public static func fromMolfile(_ code: Int) -> BondStereo {
        switch code {
        case 1: .up
        case 4: .either
        case 6: .down
        default: .none
        }
    }
}

public enum DoubleBondStereo: String, Codable, Hashable, Sendable {
    case cis
    case trans
}

/// Reacting-center status encoded by MDL V2000/V3000 bond records.
public enum CDKMDLReactingCenterStatus: Int, Codable, CaseIterable, Hashable, Sendable {
    case notReactingCenter = -1
    case unmarked = 0
    case genericReactingCenter = 1
    case noChange = 2
    case bondMadeOrBroken = 4
    case genericReactingCenterAndBondMadeOrBroken = 5
    case bondOrderChanges = 8
    case genericReactingCenterAndBondOrderChanges = 9
    case bondMadeOrBrokenAndBondOrderChanges = 12
    case genericReactingCenterAndBondMadeOrBrokenAndBondOrderChanges = 13
}

public enum AtomChirality: Codable, Hashable, Sendable {
    case none
    case clockwise
    case anticlockwise
}

public enum AtomHybridization: String, Codable, Hashable, Sendable {
    case sp1
    case sp2
    case sp3
    case planar3
    case s
    case other
}

public enum AtomQueryType: String, Codable, Hashable, Sendable {
    case anyAtom
    case anyNonHydrogen
    case anyHetero
}

public enum BondQueryType: String, Codable, Hashable, Sendable {
    case singleOrDouble
    case singleOrAromatic
    case doubleOrAromatic
    case any
}

public enum BondTopology: String, Codable, Hashable, Sendable {
    case ring
    case chain
}

public struct CxCoordinate: Hashable, Codable, Sendable {
    public var x: Double
    public var y: Double
    public var z: Double

    public init(x: Double, y: Double, z: Double = 0) {
        self.x = x
        self.y = y
        self.z = z
    }
}

public enum CxRadicalType: Int, CaseIterable, Codable, Hashable, Sendable {
    case monovalent = 1
    case divalent = 2
    case divalentSinglet = 3
    case divalentTriplet = 4
    case trivalent = 5
    case trivalentDoublet = 6
    case trivalentQuartet = 7

    public var electronCount: Int {
        switch self {
        case .monovalent:
            return 1
        case .divalent, .divalentSinglet, .divalentTriplet:
            return 2
        case .trivalent, .trivalentDoublet, .trivalentQuartet:
            return 3
        }
    }
}

public struct MoleculeSgroupBracket: Hashable, Codable, Sendable {
    public var firstPoint: CGPoint
    public var secondPoint: CGPoint

    public init(firstPoint: CGPoint, secondPoint: CGPoint) {
        self.firstPoint = firstPoint
        self.secondPoint = secondPoint
    }
}

public struct MoleculeSgroup: Hashable, Codable, Sendable {
    public enum Kind: String, Codable, Hashable, Sendable {
        case structureRepeatUnit
        case extMulticenter
        case polymer
        case data
        case generic
    }

    public let kind: Kind
    public var keyword: String?
    public var atomIDs: [Int]
    public var crossingBondIDs: [Int]
    public var subscriptText: String?
    public var superscriptText: String?
    public var roundBrackets: Bool
    public var connectivity: String?
    public var dataFieldName: String?
    public var dataValue: String?
    public var dataOperator: String?
    public var dataUnit: String?
    public var dataTag: String?
    public var subtype: String?
    public var parentAtomIDs: [Int]
    public var componentNumber: Int?
    public var expanded: Bool
    public var brackets: [MoleculeSgroupBracket]
    public var childGroupIndices: [Int]

    public init(
        kind: Kind,
        keyword: String? = nil,
        atomIDs: [Int],
        crossingBondIDs: [Int] = [],
        subscriptText: String? = nil,
        superscriptText: String? = nil,
        roundBrackets: Bool = false,
        connectivity: String? = nil,
        dataFieldName: String? = nil,
        dataValue: String? = nil,
        dataOperator: String? = nil,
        dataUnit: String? = nil,
        dataTag: String? = nil,
        subtype: String? = nil,
        parentAtomIDs: [Int] = [],
        componentNumber: Int? = nil,
        expanded: Bool = false,
        brackets: [MoleculeSgroupBracket] = [],
        childGroupIndices: [Int] = []
    ) {
        self.kind = kind
        self.keyword = keyword
        self.atomIDs = atomIDs
        self.crossingBondIDs = crossingBondIDs
        self.subscriptText = subscriptText
        self.superscriptText = superscriptText
        self.roundBrackets = roundBrackets
        self.connectivity = connectivity
        self.dataFieldName = dataFieldName
        self.dataValue = dataValue
        self.dataOperator = dataOperator
        self.dataUnit = dataUnit
        self.dataTag = dataTag
        self.subtype = subtype
        self.parentAtomIDs = parentAtomIDs
        self.componentNumber = componentNumber
        self.expanded = expanded
        self.brackets = brackets
        self.childGroupIndices = childGroupIndices
    }
}
