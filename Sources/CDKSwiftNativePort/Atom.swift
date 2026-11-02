import Foundation

#if canImport(CoreGraphics)
    import CoreGraphics
#endif

public struct Atom: Identifiable, Hashable, Codable, Sendable {
    public let id: Int
    public var externalID: String? = nil
    public var element: String
    public var position: CGPoint
    public var zPosition: Double? = nil
    public var charge: Int = 0
    public var isotopeMassNumber: Int? = nil
    public var aromatic: Bool = false
    public var chirality: AtomChirality = .none
    // Bracket-specified hydrogen count (e.g. [nH], [CH2], [nH0+]).
    // `nil` means unspecified and should be inferred heuristically.
    public var explicitHydrogenCount: Int? = nil
    public var queryType: AtomQueryType? = nil
    public var atomList: [String]? = nil
    public var atomListIsNegated: Bool = false
    public var radical: Int? = nil
    public var radicalType: CxRadicalType? = nil
    public var atomValue: String? = nil
    public var rGroupLabel: Int? = nil
    public var rGroupMembership: String? = nil
    public var componentGroupID: Int? = nil
    public var reactionRole: CDKReactionRole? = nil
    public var substitutionCount: Int? = nil
    public var unsaturated: Int? = nil
    public var ringBondCount: Int? = nil
    public var attachmentPoint: Int? = nil
    public var valenceOverride: Int? = nil
    public var cxStereoGroup: Int? = nil
    public var ligandOrderingAtomIDs: [Int]? = nil
    public var atomClass: Int? = nil
    public var atomMapNumber: Int? = nil
    public var aliasLabel: String? = nil
    public var properties: [String: String] = [:]
    public var atomTypeName: String? = nil
    public var maximumBondOrder: BondOrder? = nil
    public var bondOrderSum: Double? = nil
    public var valency: Int? = nil
    public var formalNeighbourCount: Int? = nil
    public var hybridization: AtomHybridization? = nil

    public var symbolToDraw: String {
        if let aliasLabel {
            let trimmed = aliasLabel.trimmingCharacters(in: .whitespacesAndNewlines)
            if !trimmed.isEmpty {
                return trimmed
            }
        }
        if let rGroupLabel {
            let normalized = element.trimmingCharacters(in: .whitespacesAndNewlines).uppercased()
            if rGroupLabel == 0, normalized == "R" {
                return "R"
            }
            if normalized.isEmpty || normalized == "*" || normalized == "R" || normalized == "R#" {
                return "R\(rGroupLabel)"
            }
        }
        if let isotopeMassNumber {
            return "\(isotopeMassNumber)\(element)"
        }
        return element
    }

    public init(
        id: Int,
        externalID: String? = nil,
        element: String,
        position: CGPoint,
        zPosition: Double? = nil,
        charge: Int = 0,
        isotopeMassNumber: Int? = nil,
        aromatic: Bool = false,
        chirality: AtomChirality = .none,
        explicitHydrogenCount: Int? = nil,
        queryType: AtomQueryType? = nil,
        atomList: [String]? = nil,
        atomListIsNegated: Bool = false,
        radical: Int? = nil,
        radicalType: CxRadicalType? = nil,
        atomValue: String? = nil,
        rGroupLabel: Int? = nil,
        rGroupMembership: String? = nil,
        componentGroupID: Int? = nil,
        reactionRole: CDKReactionRole? = nil,
        substitutionCount: Int? = nil,
        unsaturated: Int? = nil,
        ringBondCount: Int? = nil,
        attachmentPoint: Int? = nil,
        valenceOverride: Int? = nil,
        cxStereoGroup: Int? = nil,
        ligandOrderingAtomIDs: [Int]? = nil,
        atomClass: Int? = nil,
        atomMapNumber: Int? = nil,
        aliasLabel: String? = nil,
        properties: [String: String] = [:],
        atomTypeName: String? = nil,
        maximumBondOrder: BondOrder? = nil,
        bondOrderSum: Double? = nil,
        valency: Int? = nil,
        formalNeighbourCount: Int? = nil,
        hybridization: AtomHybridization? = nil
    ) {
        self.id = id
        self.externalID = externalID
        self.element = element
        self.position = position
        self.zPosition = zPosition
        self.charge = charge
        self.isotopeMassNumber = isotopeMassNumber
        self.aromatic = aromatic
        self.chirality = chirality
        self.explicitHydrogenCount = explicitHydrogenCount
        self.queryType = queryType
        self.atomList = atomList
        self.atomListIsNegated = atomListIsNegated
        self.radical = radical
        self.radicalType = radicalType
        self.atomValue = atomValue
        self.rGroupLabel = rGroupLabel
        self.rGroupMembership = rGroupMembership
        self.componentGroupID = componentGroupID
        self.reactionRole = reactionRole
        self.substitutionCount = substitutionCount
        self.unsaturated = unsaturated
        self.ringBondCount = ringBondCount
        self.attachmentPoint = attachmentPoint
        self.valenceOverride = valenceOverride
        self.cxStereoGroup = cxStereoGroup
        self.ligandOrderingAtomIDs = ligandOrderingAtomIDs
        self.atomClass = atomClass
        self.atomMapNumber = atomMapNumber
        self.aliasLabel = aliasLabel
        self.properties = properties
        self.atomTypeName = atomTypeName
        self.maximumBondOrder = maximumBondOrder
        self.bondOrderSum = bondOrderSum
        self.valency = valency
        self.formalNeighbourCount = formalNeighbourCount
        self.hybridization = hybridization
    }

    public func hash(into hasher: inout Hasher) {
        hasher.combine(id)
        hasher.combine(externalID)
        hasher.combine(element)
        hasher.combine(position.x)
        hasher.combine(position.y)
        hasher.combine(zPosition)
        hasher.combine(charge)
        hasher.combine(isotopeMassNumber)
        hasher.combine(aromatic)
        hasher.combine(chirality)
        hasher.combine(explicitHydrogenCount)
        hasher.combine(queryType)
        hasher.combine(atomList)
        hasher.combine(atomListIsNegated)
        hasher.combine(radical)
        hasher.combine(radicalType)
        hasher.combine(atomValue)
        hasher.combine(rGroupLabel)
        hasher.combine(rGroupMembership)
        hasher.combine(componentGroupID)
        hasher.combine(reactionRole)
        hasher.combine(substitutionCount)
        hasher.combine(unsaturated)
        hasher.combine(ringBondCount)
        hasher.combine(attachmentPoint)
        hasher.combine(valenceOverride)
        hasher.combine(cxStereoGroup)
        hasher.combine(ligandOrderingAtomIDs)
        hasher.combine(atomClass)
        hasher.combine(atomMapNumber)
        hasher.combine(aliasLabel)
        for key in properties.keys.sorted() {
            hasher.combine(key)
            hasher.combine(properties[key])
        }
        hasher.combine(atomTypeName)
        hasher.combine(maximumBondOrder)
        hasher.combine(bondOrderSum)
        hasher.combine(valency)
        hasher.combine(formalNeighbourCount)
        hasher.combine(hybridization)
    }
}
