import Foundation

public struct CDKMaximumCommonSubstructureMatch: Hashable, Codable, Sendable {
    public let atomIDByFirstAtomID: [Int: Int]
    public let bondIDByFirstBondID: [Int: Int]

    public var atomCount: Int {
        atomIDByFirstAtomID.count
    }

    public var bondCount: Int {
        bondIDByFirstBondID.count
    }

    public init(
        atomIDByFirstAtomID: [Int: Int],
        bondIDByFirstBondID: [Int: Int]
    ) {
        self.atomIDByFirstAtomID = atomIDByFirstAtomID
        self.bondIDByFirstBondID = bondIDByFirstBondID
    }
}

public enum UniversalIsomorphismTester {
    public static func maximumCommonInducedSubgraph(
        query: Molecule,
        target: Molecule,
        options: CDKSubstructureSearchOptions = .default
    ) -> CDKMaximumCommonSubstructureMatch {
        var matcher = MCISMatcher(first: query, second: target, options: options)
        return matcher.maximumCommonInducedSubgraph()
    }

    public static func maximumCommonInducedSubgraphSize(
        query: Molecule,
        target: Molecule,
        options: CDKSubstructureSearchOptions = .default
    ) -> Int {
        maximumCommonInducedSubgraph(query: query, target: target, options: options).atomCount
    }

    public static func getOverlaps(
        _ first: Molecule,
        _ second: Molecule,
        options: CDKSubstructureSearchOptions = .default
    ) -> [CDKMaximumCommonSubstructureMatch] {
        let overlap = maximumCommonInducedSubgraph(query: first, target: second, options: options)
        return overlap.atomCount == 0 ? [] : [overlap]
    }
}

private struct MCISMatcher {
    private let first: Molecule
    private let second: Molecule
    private let options: CDKSubstructureSearchOptions
    private let firstAtoms: [Atom]
    private let secondAtoms: [Atom]
    private let firstBonds: [Bond]
    private let firstAtomIndexByID: [Int: Int]
    private let secondAtomIndexByID: [Int: Int]
    private let firstBondByAtomPair: [AtomPair: Bond]
    private let secondBondByAtomPair: [AtomPair: Bond]
    private let firstAromaticBondIDs: Set<Int>
    private let secondAromaticBondIDs: Set<Int>
    private let secondRingBondIDs: Set<Int>
    private let searchOrder: [Int]
    private let maximumPossibleAtomCount: Int

    private var secondIndexByFirstIndex: [Int?]
    private var usedSecondAtomIndices: Set<Int> = []
    private var bestSecondIndexByFirstIndex: [Int?]
    private var bestAtomCount = 0

    init(first: Molecule, second: Molecule, options: CDKSubstructureSearchOptions) {
        let firstAtomIndexByID = Dictionary(
            uniqueKeysWithValues: first.atoms.enumerated().map { ($0.element.id, $0.offset) })
        let secondAtomIndexByID = Dictionary(
            uniqueKeysWithValues: second.atoms.enumerated().map { ($0.element.id, $0.offset) })
        let firstBondByAtomPair = Self.makeBondLookup(first.bonds, atomIndexByID: firstAtomIndexByID)
        let secondBondByAtomPair = Self.makeBondLookup(second.bonds, atomIndexByID: secondAtomIndexByID)
        let firstDegreeByIndex = Self.makeDegreeLookup(first.atoms, bonds: first.bonds, atomIndexByID: firstAtomIndexByID)

        self.first = first
        self.second = second
        self.options = options
        self.firstAtoms = first.atoms
        self.secondAtoms = second.atoms
        self.firstBonds = first.bonds
        self.firstAtomIndexByID = firstAtomIndexByID
        self.secondAtomIndexByID = secondAtomIndexByID
        self.firstBondByAtomPair = firstBondByAtomPair
        self.secondBondByAtomPair = secondBondByAtomPair
        self.firstAromaticBondIDs = Set(first.bonds.filter { $0.order == .aromatic }.map(\.id))
            .union(first.aromaticDisplayBondIDs())
        self.secondAromaticBondIDs = Set(second.bonds.filter { $0.order == .aromatic }.map(\.id))
            .union(second.aromaticDisplayBondIDs())
        self.secondRingBondIDs = Set(
            second.bonds.filter { CDKDescriptorSupport.isBondInRing($0, in: second) }.map(\.id))
        self.searchOrder = first.atoms.indices.sorted {
            let leftDegree = firstDegreeByIndex[$0, default: 0]
            let rightDegree = firstDegreeByIndex[$1, default: 0]
            if leftDegree != rightDegree {
                return leftDegree > rightDegree
            }
            return first.atoms[$0].id < first.atoms[$1].id
        }
        self.maximumPossibleAtomCount = min(first.atoms.count, second.atoms.count)
        self.secondIndexByFirstIndex = Array(repeating: nil, count: first.atoms.count)
        self.bestSecondIndexByFirstIndex = Array(repeating: nil, count: first.atoms.count)
    }

    mutating func maximumCommonInducedSubgraph() -> CDKMaximumCommonSubstructureMatch {
        guard !firstAtoms.isEmpty, !secondAtoms.isEmpty else {
            return CDKMaximumCommonSubstructureMatch(atomIDByFirstAtomID: [:], bondIDByFirstBondID: [:])
        }

        if firstAtoms.count <= secondAtoms.count,
            findCompleteInducedMatch(orderPosition: 0)
        {
            bestSecondIndexByFirstIndex = secondIndexByFirstIndex
            bestAtomCount = firstAtoms.count
            return makeBestMatch()
        }

        search(orderPosition: 0, mappedAtomCount: 0)
        return makeBestMatch()
    }

    private mutating func findCompleteInducedMatch(orderPosition: Int) -> Bool {
        guard orderPosition < searchOrder.count else {
            return true
        }

        let firstIndex = searchOrder[orderPosition]
        for secondIndex in secondAtoms.indices where !usedSecondAtomIndices.contains(secondIndex) {
            guard feasible(firstIndex: firstIndex, secondIndex: secondIndex) else { continue }
            secondIndexByFirstIndex[firstIndex] = secondIndex
            usedSecondAtomIndices.insert(secondIndex)
            if findCompleteInducedMatch(orderPosition: orderPosition + 1) {
                return true
            }
            usedSecondAtomIndices.remove(secondIndex)
            secondIndexByFirstIndex[firstIndex] = nil
        }
        return false
    }

    private mutating func search(orderPosition: Int, mappedAtomCount: Int) {
        if mappedAtomCount > bestAtomCount {
            bestAtomCount = mappedAtomCount
            bestSecondIndexByFirstIndex = secondIndexByFirstIndex
            if bestAtomCount == maximumPossibleAtomCount {
                return
            }
        }

        guard orderPosition < searchOrder.count else { return }
        let remainingAtoms = searchOrder.count - orderPosition
        guard mappedAtomCount + remainingAtoms > bestAtomCount else { return }
        guard mappedAtomCount + (secondAtoms.count - usedSecondAtomIndices.count) > bestAtomCount else { return }

        let firstIndex = searchOrder[orderPosition]
        for secondIndex in secondAtoms.indices where !usedSecondAtomIndices.contains(secondIndex) {
            guard feasible(firstIndex: firstIndex, secondIndex: secondIndex) else { continue }
            secondIndexByFirstIndex[firstIndex] = secondIndex
            usedSecondAtomIndices.insert(secondIndex)
            search(orderPosition: orderPosition + 1, mappedAtomCount: mappedAtomCount + 1)
            usedSecondAtomIndices.remove(secondIndex)
            secondIndexByFirstIndex[firstIndex] = nil
            if bestAtomCount == maximumPossibleAtomCount {
                return
            }
        }

        search(orderPosition: orderPosition + 1, mappedAtomCount: mappedAtomCount)
    }

    private func feasible(firstIndex: Int, secondIndex: Int) -> Bool {
        guard atomMatches(firstAtoms[firstIndex], secondAtoms[secondIndex]) else { return false }

        for mappedFirstIndex in firstAtoms.indices {
            guard let mappedSecondIndex = secondIndexByFirstIndex[mappedFirstIndex] else { continue }
            let firstBond = firstBondByAtomPair[AtomPair(firstIndex, mappedFirstIndex)]
            let secondBond = secondBondByAtomPair[AtomPair(secondIndex, mappedSecondIndex)]

            switch (firstBond, secondBond) {
            case (.none, .none):
                continue
            case let (.some(queryBond), .some(targetBond)):
                guard bondMatches(queryBond, targetBond) else { return false }
            default:
                return false
            }
        }

        return true
    }

    private func makeBestMatch() -> CDKMaximumCommonSubstructureMatch {
        var atomMap: [Int: Int] = [:]
        for firstIndex in firstAtoms.indices {
            guard let secondIndex = bestSecondIndexByFirstIndex[firstIndex] else { continue }
            atomMap[firstAtoms[firstIndex].id] = secondAtoms[secondIndex].id
        }

        var bondMap: [Int: Int] = [:]
        for firstBond in firstBonds {
            guard let firstIndex = firstAtomIndexByID[firstBond.a1],
                let secondIndex = firstAtomIndexByID[firstBond.a2],
                let mappedFirst = bestSecondIndexByFirstIndex[firstIndex],
                let mappedSecond = bestSecondIndexByFirstIndex[secondIndex],
                let secondBond = secondBondByAtomPair[AtomPair(mappedFirst, mappedSecond)],
                bondMatches(firstBond, secondBond)
            else {
                continue
            }
            bondMap[firstBond.id] = secondBond.id
        }

        return CDKMaximumCommonSubstructureMatch(
            atomIDByFirstAtomID: atomMap,
            bondIDByFirstBondID: bondMap)
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
            let targetHydrogens = CDKDescriptorSupport.totalHydrogenCount(on: targetAtom.id, in: second)
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
            let targetIsRing = secondRingBondIDs.contains(targetBond.id)
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
                return targetBond.order == .single || isAromatic(targetBond, inFirst: false)
            case .doubleOrAromatic:
                return targetBond.order == .double || isAromatic(targetBond, inFirst: false)
            }
        }

        let queryIsAromatic = isAromatic(queryBond, inFirst: true)
        let targetIsAromatic = isAromatic(targetBond, inFirst: false)
        if options.matchAromaticity, queryIsAromatic || targetIsAromatic {
            if options.strictBondOrder {
                return queryIsAromatic == targetIsAromatic
                    && (queryBond.order == targetBond.order || queryIsAromatic && targetIsAromatic)
            }
            return queryIsAromatic && targetIsAromatic || queryBond.order == targetBond.order
        }

        return queryBond.order == targetBond.order
    }

    private func isAromatic(_ bond: Bond, inFirst: Bool) -> Bool {
        let aromaticIDs = inFirst ? firstAromaticBondIDs : secondAromaticBondIDs
        return aromaticIDs.contains(bond.id)
    }

    private func canonicalElement(_ atom: Atom) -> String {
        CDKDescriptorSupport.canonicalElementSymbol(atom.element).uppercased()
    }

    private static func makeBondLookup(_ bonds: [Bond], atomIndexByID: [Int: Int]) -> [AtomPair: Bond] {
        var bondByAtomPair: [AtomPair: Bond] = [:]
        for bond in bonds {
            guard let first = atomIndexByID[bond.a1], let second = atomIndexByID[bond.a2] else {
                continue
            }
            bondByAtomPair[AtomPair(first, second)] = bond
        }
        return bondByAtomPair
    }

    private static func makeDegreeLookup(
        _ atoms: [Atom],
        bonds: [Bond],
        atomIndexByID: [Int: Int]
    ) -> [Int: Int] {
        var degreeByIndex = Dictionary(uniqueKeysWithValues: atoms.indices.map { ($0, 0) })
        for bond in bonds {
            guard let first = atomIndexByID[bond.a1], let second = atomIndexByID[bond.a2] else { continue }
            degreeByIndex[first, default: 0] += 1
            degreeByIndex[second, default: 0] += 1
        }
        return degreeByIndex
    }
}

private struct AtomPair: Hashable {
    let first: Int
    let second: Int

    init(_ first: Int, _ second: Int) {
        self.first = min(first, second)
        self.second = max(first, second)
    }
}
