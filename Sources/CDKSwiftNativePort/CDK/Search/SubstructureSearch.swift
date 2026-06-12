import Foundation

public struct CDKSubstructureSearchOptions: Hashable, Sendable {
    public var matchFormalCharge: Bool
    public var matchIsotopeMass: Bool
    public var matchAromaticity: Bool
    public var strictBondOrder: Bool
    public var matchBondTopology: Bool
    public var maximumMatches: Int?

    public init(
        matchFormalCharge: Bool = true,
        matchIsotopeMass: Bool = true,
        matchAromaticity: Bool = true,
        strictBondOrder: Bool = false,
        matchBondTopology: Bool = true,
        maximumMatches: Int? = nil
    ) {
        self.matchFormalCharge = matchFormalCharge
        self.matchIsotopeMass = matchIsotopeMass
        self.matchAromaticity = matchAromaticity
        self.strictBondOrder = strictBondOrder
        self.matchBondTopology = matchBondTopology
        self.maximumMatches = maximumMatches
    }

    public static let `default` = CDKSubstructureSearchOptions()
}

public struct CDKSubstructureMatch: Hashable, Codable, Sendable {
    public let atomIDByQueryAtomID: [Int: Int]
    public let bondIDByQueryBondID: [Int: Int]

    public var queryAtomIDs: [Int] {
        atomIDByQueryAtomID.keys.sorted()
    }

    public var targetAtomIDs: [Int] {
        queryAtomIDs.compactMap { atomIDByQueryAtomID[$0] }
    }

    public init(
        atomIDByQueryAtomID: [Int: Int],
        bondIDByQueryBondID: [Int: Int]
    ) {
        self.atomIDByQueryAtomID = atomIDByQueryAtomID
        self.bondIDByQueryBondID = bondIDByQueryBondID
    }
}

public struct CDKSubstructurePattern: Sendable {
    public let query: Molecule
    public let options: CDKSubstructureSearchOptions

    public init(
        query: Molecule,
        options: CDKSubstructureSearchOptions = .default
    ) {
        self.query = query
        self.options = options
    }

    public func matches(_ target: Molecule) -> Bool {
        firstMatch(in: target) != nil
    }

    public func firstMatch(in target: Molecule) -> CDKSubstructureMatch? {
        var limitedOptions = options
        limitedOptions.maximumMatches = 1
        return CDKSubstructureSearch.matches(query: query, target: target, options: limitedOptions).first
    }

    public func matchAll(in target: Molecule) -> [CDKSubstructureMatch] {
        CDKSubstructureSearch.matches(query: query, target: target, options: options)
    }
}

public enum CDKSubstructureSearch {
    public static func contains(
        query: Molecule,
        in target: Molecule,
        options: CDKSubstructureSearchOptions = .default
    ) -> Bool {
        CDKSubstructurePattern(query: query, options: options).matches(target)
    }

    public static func firstMatch(
        query: Molecule,
        in target: Molecule,
        options: CDKSubstructureSearchOptions = .default
    ) -> CDKSubstructureMatch? {
        CDKSubstructurePattern(query: query, options: options).firstMatch(in: target)
    }

    public static func matches(
        query: Molecule,
        target: Molecule,
        options: CDKSubstructureSearchOptions = .default
    ) -> [CDKSubstructureMatch] {
        guard !query.atoms.isEmpty else {
            return [CDKSubstructureMatch(atomIDByQueryAtomID: [:], bondIDByQueryBondID: [:])]
        }
        guard query.atoms.count <= target.atoms.count,
            query.bonds.count <= target.bonds.count
        else {
            return []
        }

        var matcher = SubstructureMatcher(query: query, target: target, options: options)
        return matcher.matchAll()
    }
}

public extension Molecule {
    func containsSubstructure(
        _ query: Molecule,
        options: CDKSubstructureSearchOptions = .default
    ) -> Bool {
        CDKSubstructureSearch.contains(query: query, in: self, options: options)
    }

    func substructureMatches(
        of query: Molecule,
        options: CDKSubstructureSearchOptions = .default
    ) -> [CDKSubstructureMatch] {
        CDKSubstructureSearch.matches(query: query, target: self, options: options)
    }
}

private struct SubstructureMatcher {
    private let query: Molecule
    private let target: Molecule
    private let options: CDKSubstructureSearchOptions
    private let queryAtoms: [Atom]
    private let targetAtoms: [Atom]
    private let queryBonds: [Bond]
    private let targetBonds: [Bond]
    private let queryAtomIndexByID: [Int: Int]
    private let targetAtomIndexByID: [Int: Int]
    private let queryAdjacency: [[QueryEdge]]
    private let targetAdjacency: [[TargetEdge]]
    private let queryAromaticBondIDs: Set<Int>
    private let targetAromaticBondIDs: Set<Int>
    private let targetRingBondIDs: Set<Int>

    private var mapping: [Int?]
    private var reverseMapping: [Int?]
    private var matches: [CDKSubstructureMatch] = []

    init(query: Molecule, target: Molecule, options: CDKSubstructureSearchOptions) {
        let queryAtomIndexByID = Dictionary(
            uniqueKeysWithValues: query.atoms.enumerated().map { ($0.element.id, $0.offset) })
        let targetAtomIndexByID = Dictionary(
            uniqueKeysWithValues: target.atoms.enumerated().map { ($0.element.id, $0.offset) })
        let queryAdjacency = Self.makeQueryAdjacency(
            atoms: query.atoms,
            bonds: query.bonds,
            atomIndexByID: queryAtomIndexByID)
        let targetAdjacency = Self.makeTargetAdjacency(
            atoms: target.atoms,
            bonds: target.bonds,
            atomIndexByID: targetAtomIndexByID)

        self.query = query
        self.target = target
        self.options = options
        self.queryAtoms = query.atoms
        self.targetAtoms = target.atoms
        self.queryBonds = query.bonds
        self.targetBonds = target.bonds
        self.queryAtomIndexByID = queryAtomIndexByID
        self.targetAtomIndexByID = targetAtomIndexByID
        self.queryAdjacency = queryAdjacency
        self.targetAdjacency = targetAdjacency
        self.queryAromaticBondIDs = Set(query.bonds.filter { $0.order == .aromatic }.map(\.id))
            .union(query.aromaticDisplayBondIDs())
        self.targetAromaticBondIDs = Set(target.bonds.filter { $0.order == .aromatic }.map(\.id))
            .union(target.aromaticDisplayBondIDs())
        self.targetRingBondIDs = Set(
            target.bonds.filter { CDKDescriptorSupport.isBondInRing($0, in: target) }.map(\.id))
        self.mapping = Array(repeating: nil, count: query.atoms.count)
        self.reverseMapping = Array(repeating: nil, count: target.atoms.count)
    }

    mutating func matchAll() -> [CDKSubstructureMatch] {
        guard candidateCountsAreFeasible() else { return [] }
        search(depth: 0)
        return matches
    }

    private func candidateCountsAreFeasible() -> Bool {
        for queryIndex in queryAtoms.indices {
            if targetAtoms.indices.allSatisfy({ !atomMatches(queryAtoms[queryIndex], targetAtoms[$0]) }) {
                return false
            }
        }
        return true
    }

    private mutating func search(depth: Int) {
        if let maximumMatches = options.maximumMatches, matches.count >= maximumMatches {
            return
        }

        guard depth < queryAtoms.count else {
            appendCurrentMatch()
            return
        }

        guard let queryIndex = selectNextQueryAtom() else {
            appendCurrentMatch()
            return
        }

        for targetIndex in targetAtoms.indices where reverseMapping[targetIndex] == nil {
            guard feasible(queryIndex: queryIndex, targetIndex: targetIndex) else { continue }
            mapping[queryIndex] = targetIndex
            reverseMapping[targetIndex] = queryIndex
            search(depth: depth + 1)
            reverseMapping[targetIndex] = nil
            mapping[queryIndex] = nil
        }
    }

    private func selectNextQueryAtom() -> Int? {
        var best: (index: Int, mappedNeighbors: Int, candidateCount: Int, degree: Int)?

        for queryIndex in queryAtoms.indices where mapping[queryIndex] == nil {
            let mappedNeighbors = queryAdjacency[queryIndex].filter { mapping[$0.neighborIndex] != nil }.count
            let candidateCount = targetAtoms.indices.reduce(0) { partial, targetIndex in
                guard reverseMapping[targetIndex] == nil else { return partial }
                return partial + (atomMatches(queryAtoms[queryIndex], targetAtoms[targetIndex]) ? 1 : 0)
            }
            let degree = queryAdjacency[queryIndex].count
            let candidate = (
                index: queryIndex,
                mappedNeighbors: mappedNeighbors,
                candidateCount: candidateCount,
                degree: degree
            )
            if let current = best {
                if candidate.mappedNeighbors > current.mappedNeighbors
                    || candidate.mappedNeighbors == current.mappedNeighbors
                        && candidate.candidateCount < current.candidateCount
                    || candidate.mappedNeighbors == current.mappedNeighbors
                        && candidate.candidateCount == current.candidateCount && candidate.degree > current.degree
                {
                    best = candidate
                }
            } else {
                best = candidate
            }
        }

        return best?.index
    }

    private func feasible(queryIndex: Int, targetIndex: Int) -> Bool {
        let queryAtom = queryAtoms[queryIndex]
        let targetAtom = targetAtoms[targetIndex]
        guard atomMatches(queryAtom, targetAtom) else { return false }

        let queryDegree = queryAdjacency[queryIndex].count
        let targetDegree = targetAdjacency[targetIndex].count
        guard targetDegree >= queryDegree else { return false }

        var unmappedQueryNeighborCount = 0
        var availableTargetNeighborCount = 0

        for edge in queryAdjacency[queryIndex] {
            if let mappedTarget = mapping[edge.neighborIndex] {
                guard let targetEdge = targetEdgeBetween(targetIndex, mappedTarget),
                    bondMatches(edge.bond, targetEdge.bond)
                else {
                    return false
                }
            } else {
                unmappedQueryNeighborCount += 1
            }
        }

        for edge in targetAdjacency[targetIndex] where reverseMapping[edge.neighborIndex] == nil {
            availableTargetNeighborCount += 1
        }
        guard availableTargetNeighborCount >= unmappedQueryNeighborCount else { return false }

        return true
    }

    private func atomMatches(_ queryAtom: Atom, _ targetAtom: Atom) -> Bool {
        let targetElement = canonicalElement(targetAtom)
        switch queryAtom.queryType {
        case .anyAtom:
            break
        case .anyNonHydrogen:
            if targetElement == "H" { return false }
        case .anyHetero:
            if targetElement == "H" || targetElement == "C" { return false }
        case nil:
            if let atomList = queryAtom.atomList, !atomList.isEmpty {
                let normalizedList = Set(atomList.map { CDKDescriptorSupport.canonicalElementSymbol($0).uppercased() })
                let contains = normalizedList.contains(targetElement)
                if queryAtom.atomListIsNegated ? contains : !contains {
                    return false
                }
            } else {
                let queryElement = canonicalElement(queryAtom)
                if queryElement != "*" && !queryElement.isEmpty && queryElement != targetElement {
                    return false
                }
            }
        }

        if options.matchFormalCharge, queryAtom.charge != targetAtom.charge {
            return false
        }
        if options.matchIsotopeMass, queryAtom.isotopeMassNumber != targetAtom.isotopeMassNumber {
            return false
        }
        if let explicitHydrogens = queryAtom.explicitHydrogenCount {
            let targetHydrogens = CDKDescriptorSupport.totalHydrogenCount(on: targetAtom.id, in: target)
            if explicitHydrogens != targetHydrogens {
                return false
            }
        }
        if queryAtom.aromatic, options.matchAromaticity, !targetAtom.aromatic {
            return false
        }

        return true
    }

    private func bondMatches(_ queryBond: Bond, _ targetBond: Bond) -> Bool {
        if options.matchBondTopology, let topology = queryBond.topology {
            let targetIsRing = targetRingBondIDs.contains(targetBond.id)
            switch topology {
            case .ring where !targetIsRing:
                return false
            case .chain where targetIsRing:
                return false
            default:
                break
            }
        }

        if let queryType = queryBond.queryType {
            switch queryType {
            case .any:
                return true
            case .singleOrDouble:
                return targetBond.order == .single || targetBond.order == .double
            case .singleOrAromatic:
                return targetBond.order == .single || isAromatic(targetBond, inQuery: false)
            case .doubleOrAromatic:
                return targetBond.order == .double || isAromatic(targetBond, inQuery: false)
            }
        }

        let queryIsAromatic = isAromatic(queryBond, inQuery: true)
        let targetIsAromatic = isAromatic(targetBond, inQuery: false)
        if options.matchAromaticity, queryIsAromatic || targetIsAromatic {
            if options.strictBondOrder {
                return queryIsAromatic == targetIsAromatic
                    && (queryBond.order == targetBond.order || queryIsAromatic && targetIsAromatic)
            }
            return queryIsAromatic && targetIsAromatic || queryBond.order == targetBond.order
        }

        return queryBond.order == targetBond.order
    }

    private func isAromatic(_ bond: Bond, inQuery: Bool) -> Bool {
        let aromaticIDs = inQuery ? queryAromaticBondIDs : targetAromaticBondIDs
        return aromaticIDs.contains(bond.id)
    }

    private func targetEdgeBetween(_ first: Int, _ second: Int) -> TargetEdge? {
        targetAdjacency[first].first { $0.neighborIndex == second }
    }

    private mutating func appendCurrentMatch() {
        var atomMap: [Int: Int] = [:]
        atomMap.reserveCapacity(queryAtoms.count)
        for queryIndex in queryAtoms.indices {
            guard let targetIndex = mapping[queryIndex] else { return }
            atomMap[queryAtoms[queryIndex].id] = targetAtoms[targetIndex].id
        }

        var bondMap: [Int: Int] = [:]
        bondMap.reserveCapacity(queryBonds.count)
        for queryBond in queryBonds {
            guard let queryFirst = queryAtomIndexByID[queryBond.a1],
                let querySecond = queryAtomIndexByID[queryBond.a2],
                let targetFirst = mapping[queryFirst],
                let targetSecond = mapping[querySecond],
                let targetEdge = targetEdgeBetween(targetFirst, targetSecond)
            else {
                return
            }
            bondMap[queryBond.id] = targetEdge.bond.id
        }

        let match = CDKSubstructureMatch(
            atomIDByQueryAtomID: atomMap,
            bondIDByQueryBondID: bondMap)
        matches.append(match)
    }

    private func canonicalElement(_ atom: Atom) -> String {
        CDKDescriptorSupport.canonicalElementSymbol(atom.element).uppercased()
    }

    private static func makeQueryAdjacency(
        atoms: [Atom],
        bonds: [Bond],
        atomIndexByID: [Int: Int]
    ) -> [[QueryEdge]] {
        var adjacency = Array(repeating: [QueryEdge](), count: atoms.count)
        for bond in bonds {
            guard let first = atomIndexByID[bond.a1], let second = atomIndexByID[bond.a2] else { continue }
            adjacency[first].append(QueryEdge(neighborIndex: second, bond: bond))
            adjacency[second].append(QueryEdge(neighborIndex: first, bond: bond))
        }
        return adjacency
    }

    private static func makeTargetAdjacency(
        atoms: [Atom],
        bonds: [Bond],
        atomIndexByID: [Int: Int]
    ) -> [[TargetEdge]] {
        var adjacency = Array(repeating: [TargetEdge](), count: atoms.count)
        for bond in bonds {
            guard let first = atomIndexByID[bond.a1], let second = atomIndexByID[bond.a2] else { continue }
            adjacency[first].append(TargetEdge(neighborIndex: second, bond: bond))
            adjacency[second].append(TargetEdge(neighborIndex: first, bond: bond))
        }
        return adjacency
    }

    private struct QueryEdge {
        let neighborIndex: Int
        let bond: Bond
    }

    private struct TargetEdge {
        let neighborIndex: Int
        let bond: Bond
    }
}
