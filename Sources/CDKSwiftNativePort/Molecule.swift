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

    public init(kind: Kind,
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
                childGroupIndices: [Int] = []) {
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
            if normalized.isEmpty || normalized == "*" || normalized == "R" || normalized == "R#" {
                return "R\(rGroupLabel)"
            }
        }
        if let isotopeMassNumber {
            return "\(isotopeMassNumber)\(element)"
        }
        return element
    }

    public init(id: Int,
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
                hybridization: AtomHybridization? = nil) {
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

public struct Bond: Identifiable, Hashable, Codable, Sendable {
    public let id: Int
    public var externalID: String? = nil
    public let a1: Int
    public let a2: Int
    public var order: BondOrder
    public var stereo: BondStereo = .none
    public var doubleBondStereo: DoubleBondStereo? = nil
    public var stereoReferenceAtomIDs: [Int]? = nil
    public var queryType: BondQueryType? = nil
    public var topology: BondTopology? = nil

    public init(id: Int,
                externalID: String? = nil,
                a1: Int,
                a2: Int,
                order: BondOrder,
                stereo: BondStereo = .none,
                doubleBondStereo: DoubleBondStereo? = nil,
                stereoReferenceAtomIDs: [Int]? = nil,
                queryType: BondQueryType? = nil,
                topology: BondTopology? = nil) {
        self.id = id
        self.externalID = externalID
        self.a1 = a1
        self.a2 = a2
        self.order = order
        self.stereo = stereo
        self.doubleBondStereo = doubleBondStereo
        self.stereoReferenceAtomIDs = stereoReferenceAtomIDs
        self.queryType = queryType
        self.topology = topology
    }
}

public struct Molecule: Hashable, Codable, Sendable {
    public var name: String = "Untitled"
    public var externalID: String? = nil
    public var atoms: [Atom] = []
    public var bonds: [Bond] = []
    public var sgroups: [MoleculeSgroup] = []
    public var highlightedAtomIDs: [Int] = []
    public var highlightedBondIDs: [Int] = []
    public var cxState: CDKCxSmilesState? = nil
    public var dataFields: [String: [String]] = [:]
    public var dataFieldOrder: [String] = []
    public var rGroupLogicDefinitions: [Int: MoleculeRGroupLogic] = [:]

    public var atomCount: Int { atoms.count }
    public var bondCount: Int { bonds.count }
    public var orderedDataFieldNames: [String] {
        Molecule.normalizedDataFieldOrder(preferredOrder: dataFieldOrder, availableFields: dataFields)
    }

    private enum CodingKeys: String, CodingKey {
        case name
        case externalID
        case atoms
        case bonds
        case sgroups
        case highlightedAtomIDs
        case highlightedBondIDs
        case cxState
        case dataFields
        case dataFieldOrder
        case rGroupLogicDefinitions
    }

    public func indexOfAtom(id: Int) -> Int? {
        atoms.firstIndex(where: { $0.id == id })
    }

    public func bonds(forAtom atomID: Int) -> [Bond] {
        bonds.filter { $0.a1 == atomID || $0.a2 == atomID }
    }

    public func atom(id atomID: Int) -> Atom? {
        atoms.first(where: { $0.id == atomID })
    }

    public func bond(between a: Int, and b: Int) -> Bond? {
        bonds.first(where: { ($0.a1 == a && $0.a2 == b) || ($0.a1 == b && $0.a2 == a) })
    }

    public func dataFieldValues(named fieldName: String) -> [String] {
        let normalized = Molecule.normalizedDataFieldName(fieldName)
        guard !normalized.isEmpty else { return [] }
        return dataFields[normalized] ?? []
    }

    public func neighbors(of atomID: Int) -> [Int] {
        bonds.compactMap { b in
            if b.a1 == atomID { return b.a2 }
            if b.a2 == atomID { return b.a1 }
            return nil
        }
    }

    /// Returns simple cycles as ordered atom id lists (without repeating the first atom at the end).
    /// This is intentionally bounded and lightweight for small molecules used in depiction.
    public func simpleCycles(maxSize: Int = 8) -> [[Int]] {
        guard atoms.count >= 3 else { return [] }

        let adjacency = adjacencyMap()
        let atomIDs = atoms.map(\.id).sorted()
        var uniqueCycles = Set<[Int]>()

        func canonical(_ cycle: [Int]) -> [Int] {
            guard !cycle.isEmpty else { return cycle }
            let n = cycle.count

            func rotations(of seq: [Int]) -> [[Int]] {
                (0..<n).map { i in
                    Array(seq[i..<n]) + Array(seq[0..<i])
                }
            }

            let forward = rotations(of: cycle)
            let reversed = rotations(of: Array(cycle.reversed()))
            return (forward + reversed).min(by: isLexicographicallySmaller) ?? cycle
        }

        func dfs(start: Int, current: Int, path: [Int], visited: Set<Int>) {
            if path.count > maxSize { return }
            let nextNodes = adjacency[current] ?? []

            for next in nextNodes {
                if next == start {
                    if path.count >= 3 {
                        uniqueCycles.insert(canonical(path))
                    }
                    continue
                }
                if visited.contains(next) { continue }
                if path.count >= maxSize { continue }
                // Prune search: only keep cycles where start is the minimum atom id.
                if next < start { continue }

                var nextVisited = visited
                nextVisited.insert(next)
                dfs(start: start, current: next, path: path + [next], visited: nextVisited)
            }
        }

        for start in atomIDs {
            dfs(start: start, current: start, path: [start], visited: [start])
        }

        return Array(uniqueCycles)
    }

    /// Rings that should be rendered with aromatic styling.
    /// Includes explicit aromatic rings and common alternating single/double rings.
    public func aromaticDisplayRings() -> [[Int]] {
        let atomByID = Dictionary(uniqueKeysWithValues: atoms.map { ($0.id, $0) })
        let cycles = simpleCycles(maxSize: 8)

        return cycles.filter { ring in
            guard ring.count >= 5 && ring.count <= 7 else { return false }
            let ringBonds = bonds(inCycle: ring)
            guard ringBonds.count == ring.count else { return false }

            let allAromaticAtoms = ring.allSatisfy { atomByID[$0]?.aromatic ?? false }
            let allAromaticBonds = ringBonds.allSatisfy { $0.order == .aromatic }
            if allAromaticAtoms || allAromaticBonds { return true }

            return isAlternatingSingleDouble(ringBonds)
        }
    }

    public func aromaticDisplayBondIDs() -> Set<Int> {
        var ids = Set(bonds.filter { $0.order == .aromatic }.map(\.id))
        for ring in aromaticDisplayRings() {
            for (a, b) in cycleEdges(ring) {
                if let bond = bond(between: a, and: b) {
                    ids.insert(bond.id)
                }
            }
        }
        return ids
    }

    /// Heuristic implicit hydrogen count for depiction labels.
    public func implicitHydrogenCount(for atomID: Int) -> Int {
        guard let atom = atom(id: atomID) else { return 0 }
        guard atom.element.uppercased() != "H" else { return 0 }
        if let explicit = atom.explicitHydrogenCount {
            return max(0, explicit)
        }

        let targetValence = preferredValence(for: atom)
        guard targetValence > 0 else { return 0 }

        let bondOrderSum = bonds(forAtom: atomID)
            .reduce(0.0) { $0 + $1.order.valenceContribution }

        return max(0, Int(round(targetValence - bondOrderSum)))
    }

    public mutating func assignWedgeHashFromChiralCenters() {
        let degreeByAtom = Dictionary(uniqueKeysWithValues: atoms.map { ($0.id, neighbors(of: $0.id).count) })
        let positionsByAtom = Dictionary(uniqueKeysWithValues: atoms.map { ($0.id, $0.position) })

        for atom in atoms where atom.chirality != .none {
            let candidateIndices = bonds.indices.filter { idx in
                let b = bonds[idx]
                guard b.order == .single, b.stereo == .none else { return false }
                return b.a1 == atom.id || b.a2 == atom.id
            }
            guard !candidateIndices.isEmpty else { continue }

            let picked = candidateIndices.min { lhs, rhs in
                let l = bondStereoPriority(for: bonds[lhs], around: atom.id, degreeByAtom: degreeByAtom)
                let r = bondStereoPriority(for: bonds[rhs], around: atom.id, degreeByAtom: degreeByAtom)
                if l != r { return l < r }
                let lc = bondStereoClearance(for: bonds[lhs], around: atom.id, positionsByAtom: positionsByAtom)
                let rc = bondStereoClearance(for: bonds[rhs], around: atom.id, positionsByAtom: positionsByAtom)
                if abs(lc - rc) > 0.0001 { return lc > rc }
                return bonds[lhs].id < bonds[rhs].id
            }
            guard let picked else { continue }

            let fromA1 = bonds[picked].a1 == atom.id
            switch (atom.chirality, fromA1) {
            case (.clockwise, true):
                bonds[picked].stereo = .up
            case (.clockwise, false):
                bonds[picked].stereo = .upReversed
            case (.anticlockwise, true):
                bonds[picked].stereo = .down
            case (.anticlockwise, false):
                bonds[picked].stereo = .downReversed
            case (.none, _):
                break
            }
        }
    }

    public func boundingBox() -> CGRect? {
        guard !atoms.isEmpty else { return nil }
        var minX = atoms[0].position.x
        var minY = atoms[0].position.y
        var maxX = atoms[0].position.x
        var maxY = atoms[0].position.y
        for a in atoms {
            minX = min(minX, a.position.x)
            minY = min(minY, a.position.y)
            maxX = max(maxX, a.position.x)
            maxY = max(maxY, a.position.y)
        }
        return CGRect(x: minX, y: minY, width: max(0.0001, maxX - minX), height: max(0.0001, maxY - minY))
    }

    private func adjacencyMap() -> [Int: [Int]] {
        var map: [Int: Set<Int>] = [:]
        for atom in atoms {
            map[atom.id] = []
        }
        for bond in bonds {
            map[bond.a1, default: []].insert(bond.a2)
            map[bond.a2, default: []].insert(bond.a1)
        }
        return map.mapValues { Array($0).sorted() }
    }

    private func cycleEdges(_ ring: [Int]) -> [(Int, Int)] {
        guard ring.count >= 2 else { return [] }
        return (0..<ring.count).map { i in
            let a = ring[i]
            let b = ring[(i + 1) % ring.count]
            return (a, b)
        }
    }

    private func bonds(inCycle ring: [Int]) -> [Bond] {
        cycleEdges(ring).compactMap { pair in
            bond(between: pair.0, and: pair.1)
        }
    }

    private func isAlternatingSingleDouble(_ cycleBonds: [Bond]) -> Bool {
        let orders = cycleBonds.map(\.order)
        guard !orders.isEmpty, orders.count % 2 == 0 else { return false }
        guard orders.allSatisfy({ $0 == .single || $0 == .double }) else { return false }

        func matchesPattern(start: BondOrder) -> Bool {
            for (idx, order) in orders.enumerated() {
                let expected: BondOrder = (idx % 2 == 0) ? start : (start == .single ? .double : .single)
                if order != expected { return false }
            }
            return true
        }

        return matchesPattern(start: .single) || matchesPattern(start: .double)
    }

    private func preferredValence(for atom: Atom) -> Double {
        switch atom.element.uppercased() {
        case "C":
            return atom.aromatic ? 3.0 : 4.0
        case "N":
            if atom.aromatic { return 3.0 }
            return atom.charge > 0 ? 4.0 : 3.0
        case "O":
            if atom.charge > 0 { return 3.0 }
            if atom.charge < 0 { return 1.0 }
            return 2.0
        case "S":
            return atom.charge > 0 ? 3.0 : 2.0
        case "P":
            return atom.charge > 0 ? 4.0 : 3.0
        case "B":
            return atom.charge < 0 ? 4.0 : 3.0
        case "F", "CL", "BR", "I":
            return 1.0
        default:
            return 0.0
        }
    }

    private func bondStereoPriority(for bond: Bond, around atomID: Int, degreeByAtom: [Int: Int]) -> Int {
        let neighborID = (bond.a1 == atomID) ? bond.a2 : bond.a1
        let neighborDegree = degreeByAtom[neighborID] ?? 0
        // Prefer terminal substituents; then stable ordering by atom id.
        return (neighborDegree == 1 ? 0 : 10) + neighborID
    }

    private func bondStereoClearance(for bond: Bond,
                                     around atomID: Int,
                                     positionsByAtom: [Int: CGPoint]) -> CGFloat {
        let neighborID = (bond.a1 == atomID) ? bond.a2 : bond.a1
        guard let center = positionsByAtom[atomID], let neighbor = positionsByAtom[neighborID] else { return 0 }
        let dx = neighbor.x - center.x
        let dy = neighbor.y - center.y
        let len = hypot(dx, dy)
        guard len > 0.0001 else { return 0 }

        let ux = dx / len
        let uy = dy / len
        let tip = CGPoint(x: neighbor.x + ux * 0.35, y: neighbor.y + uy * 0.35)

        var minDistance = CGFloat.greatestFiniteMagnitude
        var crowdPenalty: CGFloat = 0
        for atom in atoms where atom.id != atomID && atom.id != neighborID {
            guard let p = positionsByAtom[atom.id] else { continue }
            let d = tip.distance(to: p)
            minDistance = min(minDistance, d)
            crowdPenalty += 1 / max(0.15, d)
        }
        if minDistance == .greatestFiniteMagnitude { return 0 }
        return minDistance - crowdPenalty * 0.14
    }

    private func isLexicographicallySmaller(_ lhs: [Int], _ rhs: [Int]) -> Bool {
        for (a, b) in zip(lhs, rhs) {
            if a != b { return a < b }
        }
        return lhs.count < rhs.count
    }

    public init(name: String = "Untitled",
                externalID: String? = nil,
                atoms: [Atom] = [],
                bonds: [Bond] = []) {
        self.name = name
        self.externalID = externalID
        self.atoms = atoms
        self.bonds = bonds
        self.dataFields = [:]
        self.dataFieldOrder = []
        self.rGroupLogicDefinitions = [:]
    }

    public init(from decoder: Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        name = try container.decodeIfPresent(String.self, forKey: .name) ?? "Untitled"
        externalID = try container.decodeIfPresent(String.self, forKey: .externalID)
        atoms = try container.decodeIfPresent([Atom].self, forKey: .atoms) ?? []
        bonds = try container.decodeIfPresent([Bond].self, forKey: .bonds) ?? []
        sgroups = try container.decodeIfPresent([MoleculeSgroup].self, forKey: .sgroups) ?? []
        highlightedAtomIDs = try container.decodeIfPresent([Int].self, forKey: .highlightedAtomIDs) ?? []
        highlightedBondIDs = try container.decodeIfPresent([Int].self, forKey: .highlightedBondIDs) ?? []
        cxState = try container.decodeIfPresent(CDKCxSmilesState.self, forKey: .cxState)
        dataFields = try container.decodeIfPresent([String: [String]].self, forKey: .dataFields) ?? [:]
        dataFieldOrder = try container.decodeIfPresent([String].self, forKey: .dataFieldOrder) ?? []
        dataFieldOrder = Molecule.normalizedDataFieldOrder(preferredOrder: dataFieldOrder, availableFields: dataFields)
        rGroupLogicDefinitions = try container.decodeIfPresent([Int: MoleculeRGroupLogic].self,
                                                               forKey: .rGroupLogicDefinitions) ?? [:]
    }

    public func encode(to encoder: Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)
        try container.encode(name, forKey: .name)
        try container.encodeIfPresent(externalID, forKey: .externalID)
        try container.encode(atoms, forKey: .atoms)
        try container.encode(bonds, forKey: .bonds)
        try container.encode(sgroups, forKey: .sgroups)
        try container.encode(highlightedAtomIDs, forKey: .highlightedAtomIDs)
        try container.encode(highlightedBondIDs, forKey: .highlightedBondIDs)
        try container.encodeIfPresent(cxState, forKey: .cxState)
        try container.encode(dataFields, forKey: .dataFields)
        try container.encode(orderedDataFieldNames, forKey: .dataFieldOrder)
        try container.encode(rGroupLogicDefinitions, forKey: .rGroupLogicDefinitions)
    }

    public mutating func ensureDataField(named fieldName: String) {
        let normalized = Molecule.normalizedDataFieldName(fieldName)
        guard !normalized.isEmpty else { return }

        if dataFields[normalized] == nil {
            dataFields[normalized] = []
        }
        if !dataFieldOrder.contains(normalized) {
            dataFieldOrder.append(normalized)
        }
    }

    public mutating func appendDataFieldValue(_ value: String, named fieldName: String) {
        let normalized = Molecule.normalizedDataFieldName(fieldName)
        guard !normalized.isEmpty else { return }

        ensureDataField(named: normalized)
        dataFields[normalized, default: []].append(value)
    }

    private static func normalizedDataFieldName(_ fieldName: String) -> String {
        fieldName.trimmingCharacters(in: .whitespacesAndNewlines)
    }

    private static func normalizedDataFieldOrder(preferredOrder: [String],
                                                 availableFields: [String: [String]]) -> [String] {
        var orderedNames: [String] = []
        orderedNames.reserveCapacity(availableFields.count)
        var seen = Set<String>()

        for rawName in preferredOrder {
            let normalized = normalizedDataFieldName(rawName)
            guard !normalized.isEmpty,
                  availableFields[normalized] != nil,
                  seen.insert(normalized).inserted else {
                continue
            }
            orderedNames.append(normalized)
        }

        for fieldName in availableFields.keys.sorted() where seen.insert(fieldName).inserted {
            orderedNames.append(fieldName)
        }

        return orderedNames
    }
}

public enum Depiction2DGenerator {
    public static func generate(for molecule: Molecule) -> Molecule {
        guard let collapsed = CDKHydrogenLayoutCollapse.collapse(in: molecule) else {
            return CDKStructureDiagramGenerator.apply(to: molecule)
        }

        let laidOutScaffold = CDKStructureDiagramGenerator.apply(to: collapsed.scaffold)
        return CDKHydrogenLayoutCollapse.restore(original: molecule,
                                                 from: laidOutScaffold,
                                                 hydrogenAnchorByID: collapsed.hydrogenAnchorByID)
    }
}

private struct CDKHydrogenLayoutCollapseResult {
    let scaffold: Molecule
    let hydrogenAnchorByID: [Int: Int]
}

private enum CDKHydrogenLayoutCollapse {
    static func collapse(in molecule: Molecule) -> CDKHydrogenLayoutCollapseResult? {
        let hydrogenAnchorByID = Dictionary(uniqueKeysWithValues: molecule.atoms.compactMap { atom -> (Int, Int)? in
            guard isSuppressibleHydrogen(atomID: atom.id, in: molecule),
                  let bond = molecule.bonds(forAtom: atom.id).first else {
                return nil
            }
            let anchorID = (bond.a1 == atom.id) ? bond.a2 : bond.a1
            return (atom.id, anchorID)
        })

        guard !hydrogenAnchorByID.isEmpty else { return nil }

        let suppressedHydrogenIDs = Set(hydrogenAnchorByID.keys)
        var removedHydrogenCountByNeighbor: [Int: Int] = [:]
        for anchorID in hydrogenAnchorByID.values {
            removedHydrogenCountByNeighbor[anchorID, default: 0] += 1
        }

        var scaffold = molecule
        scaffold.atoms = molecule.atoms.compactMap { atom in
            guard !suppressedHydrogenIDs.contains(atom.id) else { return nil }
            var updated = atom
            if let removed = removedHydrogenCountByNeighbor[atom.id], removed > 0 {
                updated.explicitHydrogenCount = max(0, updated.explicitHydrogenCount ?? 0) + removed
            }
            return updated
        }
        scaffold.bonds = molecule.bonds.filter { bond in
            !suppressedHydrogenIDs.contains(bond.a1) && !suppressedHydrogenIDs.contains(bond.a2)
        }

        return CDKHydrogenLayoutCollapseResult(scaffold: scaffold,
                                               hydrogenAnchorByID: hydrogenAnchorByID)
    }

    static func restore(original: Molecule,
                        from scaffold: Molecule,
                        hydrogenAnchorByID: [Int: Int]) -> Molecule {
        var restored = original
        let scaffoldPositionByAtomID = Dictionary(uniqueKeysWithValues: scaffold.atoms.map { ($0.id, $0.position) })

        for index in restored.atoms.indices {
            if let position = scaffoldPositionByAtomID[restored.atoms[index].id] {
                restored.atoms[index].position = position
            }
        }

        placeSuppressedHydrogens(in: &restored, hydrogenAnchorByID: hydrogenAnchorByID)
        return restored
    }

    private static func placeSuppressedHydrogens(in molecule: inout Molecule,
                                                 hydrogenAnchorByID: [Int: Int]) {
        let groupedByAnchor = Dictionary(grouping: hydrogenAnchorByID.keys, by: { hydrogenAnchorByID[$0] ?? -1 })

        for (anchorID, hydrogenIDs) in groupedByAnchor {
            guard let anchorIndex = molecule.atoms.firstIndex(where: { $0.id == anchorID }) else { continue }
            let anchor = molecule.atoms[anchorIndex]
            let anchorPoint = anchor.position

            let placedNeighborIDs = molecule.neighbors(of: anchorID).filter { neighborID in
                !hydrogenIDs.contains(neighborID)
            }
            let placedNeighborPoints = placedNeighborIDs.compactMap { molecule.atom(id: $0)?.position }
            let heavyNeighborPoints = placedNeighborIDs.compactMap { neighborID -> CGPoint? in
                guard let atom = molecule.atom(id: neighborID),
                      !isHydrogen(atom) else {
                    return nil
                }
                return atom.position
            }

            let baseAngle = preferredHydrogenBaseAngle(anchor: anchorPoint,
                                                       neighborPoints: heavyNeighborPoints.isEmpty ? placedNeighborPoints : heavyNeighborPoints)
            let step = CGFloat.pi / 3.0
            let centeredOffsets = hydrogenOffsets(count: hydrogenIDs.count, step: step)
            let referenceNeighbors = heavyNeighborPoints.isEmpty ? placedNeighborPoints : heavyNeighborPoints
            let meanNeighborDistance = referenceNeighbors.isEmpty
                ? CGFloat(1.4)
                : referenceNeighbors.reduce(CGFloat.zero) { partial, point in
                    partial + anchorPoint.distance(to: point)
                } / CGFloat(referenceNeighbors.count)
            let hydrogenBondLength = max(0.72, min(1.15, meanNeighborDistance * 0.78))

            for (offset, hydrogenID) in zip(centeredOffsets, hydrogenIDs.sorted()) {
                guard let hydrogenIndex = molecule.atoms.firstIndex(where: { $0.id == hydrogenID }) else { continue }
                let angle = baseAngle + offset
                molecule.atoms[hydrogenIndex].position = CGPoint(x: anchorPoint.x + cos(angle) * hydrogenBondLength,
                                                                 y: anchorPoint.y + sin(angle) * hydrogenBondLength)
            }
        }
    }

    private static func preferredHydrogenBaseAngle(anchor: CGPoint,
                                                   neighborPoints: [CGPoint]) -> CGFloat {
        guard !neighborPoints.isEmpty else { return 0 }

        let unitVectors = neighborPoints.compactMap { point -> CGVector? in
            let dx = point.x - anchor.x
            let dy = point.y - anchor.y
            let length = hypot(dx, dy)
            guard length > 0.0001 else { return nil }
            return CGVector(dx: dx / length, dy: dy / length)
        }

        guard !unitVectors.isEmpty else { return 0 }

        let summed = unitVectors.reduce(CGVector.zero) { partial, vector in
            CGVector(dx: partial.dx + vector.dx, dy: partial.dy + vector.dy)
        }
        let magnitude = hypot(summed.dx, summed.dy)
        if magnitude > 0.2 {
            return atan2(-summed.dy, -summed.dx)
        }

        let sortedAngles = unitVectors
            .map { atan2($0.dy, $0.dx) }
            .sorted()
        guard sortedAngles.count >= 2 else {
            return sortedAngles[0] + .pi
        }

        var widestGap: CGFloat = -.greatestFiniteMagnitude
        var bestMidpoint = sortedAngles[0] + .pi
        for index in sortedAngles.indices {
            let start = sortedAngles[index]
            let end = sortedAngles[(index + 1) % sortedAngles.count] + (index + 1 == sortedAngles.count ? (.pi * 2.0) : 0)
            let gap = end - start
            if gap > widestGap {
                widestGap = gap
                bestMidpoint = start + (gap * 0.5)
            }
        }

        return normalizedAngle(bestMidpoint)
    }

    private static func hydrogenOffsets(count: Int, step: CGFloat) -> [CGFloat] {
        guard count > 1 else { return [0] }
        let center = CGFloat(count - 1) * 0.5
        return (0..<count).map { index in
            (CGFloat(index) - center) * step
        }
    }

    private static func normalizedAngle(_ angle: CGFloat) -> CGFloat {
        var normalized = angle.truncatingRemainder(dividingBy: .pi * 2.0)
        if normalized <= -.pi {
            normalized += .pi * 2.0
        } else if normalized > .pi {
            normalized -= .pi * 2.0
        }
        return normalized
    }

    private static func isSuppressibleHydrogen(atomID: Int, in molecule: Molecule) -> Bool {
        guard let atom = molecule.atom(id: atomID), isHydrogen(atom) else { return false }
        guard atom.charge == 0,
              atom.isotopeMassNumber == nil,
              atom.radical == nil,
              atom.queryType == nil,
              atom.atomList == nil else {
            return false
        }

        let attachedBonds = molecule.bonds(forAtom: atom.id)
        guard attachedBonds.count == 1, let bond = attachedBonds.first else { return false }
        guard bond.order == .single, bond.stereo == .none, bond.queryType == nil else { return false }

        let neighborID = (bond.a1 == atom.id) ? bond.a2 : bond.a1
        guard let neighbor = molecule.atom(id: neighborID), !isHydrogen(neighbor) else { return false }
        return true
    }

    private static func isHydrogen(_ atom: Atom) -> Bool {
        atom.element.trimmingCharacters(in: .whitespacesAndNewlines).uppercased() == "H"
    }
}

public enum ChemError: LocalizedError {
    case emptyInput
    case unsupported(String)
    case parseFailed(String)

    public var errorDescription: String? {
        switch self {
        case .emptyInput:
            return "No input provided."
        case .unsupported(let s):
            return "Unsupported: \(s)"
        case .parseFailed(let s):
            return "Could not parse: \(s)"
        }
    }
}
