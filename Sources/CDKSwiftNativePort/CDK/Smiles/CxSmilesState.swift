import Foundation

public struct CDKCxSmilesState: Equatable, Hashable, Codable, Sendable {
    public enum WedgeDisplay: String, Equatable, Hashable, Codable, Sendable {
        case wedgeBegin
        case wedgedHashBegin
    }

    public struct BondDisplayEntry: Equatable, Hashable, Codable, Sendable {
        public var atomIndex: Int
        public var bondIndex: Int
        public var display: WedgeDisplay

        public init(atomIndex: Int, bondIndex: Int, display: WedgeDisplay) {
            self.atomIndex = atomIndex
            self.bondIndex = bondIndex
            self.display = display
        }
    }

    public struct LinkNode: Equatable, Hashable, Codable, Sendable {
        public var atomIndex: Int
        public var lowerBound: Int
        public var upperBound: Int
        public var bondIndices: [Int]
        /// CX LN trailing indices refer to neighboring atoms, not bonds.
        public var outerAtomIndices: [Int]? = nil

        public init(
            atomIndex: Int,
            lowerBound: Int,
            upperBound: Int,
            bondIndices: [Int] = []
        ) {
            self.atomIndex = atomIndex
            self.lowerBound = lowerBound
            self.upperBound = upperBound
            self.bondIndices = bondIndices
        }
    }

    public struct SgroupDefinition: Equatable, Hashable, Codable, Sendable {
        public enum Kind: String, Codable, Hashable, Sendable {
            case polymer
            case data
        }

        public var kind: Kind
        public var keyword: String?
        public var atomIndices: [Int]
        public var bondIndices: [Int]
        public var subscriptText: String?
        public var superscriptText: String?
        public var dataFieldName: String?
        public var dataValue: String?
        public var dataOperator: String?
        public var dataUnit: String?
        public var dataTag: String?
        public var childIndices: [Int]

        public init(
            kind: Kind,
            keyword: String? = nil,
            atomIndices: [Int] = [],
            bondIndices: [Int] = [],
            subscriptText: String? = nil,
            superscriptText: String? = nil,
            dataFieldName: String? = nil,
            dataValue: String? = nil,
            dataOperator: String? = nil,
            dataUnit: String? = nil,
            dataTag: String? = nil,
            childIndices: [Int] = []
        ) {
            self.kind = kind
            self.keyword = keyword
            self.atomIndices = atomIndices
            self.bondIndices = bondIndices
            self.subscriptText = subscriptText
            self.superscriptText = superscriptText
            self.dataFieldName = dataFieldName
            self.dataValue = dataValue
            self.dataOperator = dataOperator
            self.dataUnit = dataUnit
            self.dataTag = dataTag
            self.childIndices = childIndices
        }
    }

    public var atomLabels: [Int: String] = [:]
    public var atomValues: [Int: String] = [:]
    public var atomCoordinates: [CxCoordinate] = []
    public var has3DCoordinates: Bool = false
    public var fragmentGroups: [[Int]] = []
    public var racemic: Bool = false
    public var racemicFragments: [Int] = []
    public var positionalVariations: [Int: [Int]] = [:]
    public var ligandOrdering: [Int: [Int]] = [:]
    public var stereoGroups: [Int: Int] = [:]
    public var atomRadicals: [Int: CxRadicalType] = [:]
    public var sgroups: [SgroupDefinition] = []
    public var atomHighlights: [Int] = []
    public var bondHighlights: [Int] = []
    public var bondDisplays: [BondDisplayEntry] = []
    public var coordinateBonds: [Int: [Int]]? = nil
    public var linkNodes: [LinkNode] = []
    public var rGroupDefinitions: [String: [String]] = [:]
    public var rGroupOrder: [String] = []
    /// Optional for compatibility with older serialized CX state.
    public var rGroupLogic: [Int: MoleculeRGroupLogic]? = nil

    public init(
        atomLabels: [Int: String] = [:],
        atomValues: [Int: String] = [:],
        atomCoordinates: [CxCoordinate] = [],
        has3DCoordinates: Bool = false,
        fragmentGroups: [[Int]] = [],
        racemic: Bool = false,
        racemicFragments: [Int] = [],
        positionalVariations: [Int: [Int]] = [:],
        ligandOrdering: [Int: [Int]] = [:],
        stereoGroups: [Int: Int] = [:],
        atomRadicals: [Int: CxRadicalType] = [:],
        sgroups: [SgroupDefinition] = [],
        atomHighlights: [Int] = [],
        bondHighlights: [Int] = [],
        bondDisplays: [BondDisplayEntry] = [],
        coordinateBonds: [Int: [Int]]? = nil,
        linkNodes: [LinkNode] = [],
        rGroupDefinitions: [String: [String]] = [:],
        rGroupOrder: [String] = []
    ) {
        self.atomLabels = atomLabels
        self.atomValues = atomValues
        self.atomCoordinates = atomCoordinates
        self.has3DCoordinates = has3DCoordinates
        self.fragmentGroups = fragmentGroups
        self.racemic = racemic
        self.racemicFragments = racemicFragments
        self.positionalVariations = positionalVariations
        self.ligandOrdering = ligandOrdering
        self.stereoGroups = stereoGroups
        self.atomRadicals = atomRadicals
        self.sgroups = sgroups
        self.atomHighlights = atomHighlights
        self.bondHighlights = bondHighlights
        self.bondDisplays = bondDisplays
        self.coordinateBonds = coordinateBonds
        self.linkNodes = linkNodes
        self.rGroupDefinitions = rGroupDefinitions
        self.rGroupOrder = rGroupOrder
    }
}
