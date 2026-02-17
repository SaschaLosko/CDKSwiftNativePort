import Foundation
import CoreGraphics

public enum ChemFormat: String, CaseIterable, Identifiable {
    case sdf = "SDF / MOL"
    case smiles = "SMILES"
    case inchi = "InChI"

    public var id: String { rawValue }
}

public enum BondOrder: Int, Codable, Hashable {
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

public enum BondStereo: Codable, Hashable {
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

public enum AtomChirality: Codable, Hashable {
    case none
    case clockwise
    case anticlockwise
}

public struct Atom: Identifiable, Hashable, Codable {
    public let id: Int
    public var element: String
    public var position: CGPoint
    public var charge: Int = 0
    public var isotopeMassNumber: Int? = nil
    public var aromatic: Bool = false
    public var chirality: AtomChirality = .none
    // Bracket-specified hydrogen count (e.g. [nH], [CH2], [nH0+]).
    // `nil` means unspecified and should be inferred heuristically.
    public var explicitHydrogenCount: Int? = nil

    public var symbolToDraw: String {
        if let isotopeMassNumber {
            return "\(isotopeMassNumber)\(element)"
        }
        return element
    }

    public init(id: Int,
                element: String,
                position: CGPoint,
                charge: Int = 0,
                isotopeMassNumber: Int? = nil,
                aromatic: Bool = false,
                chirality: AtomChirality = .none,
                explicitHydrogenCount: Int? = nil) {
        self.id = id
        self.element = element
        self.position = position
        self.charge = charge
        self.isotopeMassNumber = isotopeMassNumber
        self.aromatic = aromatic
        self.chirality = chirality
        self.explicitHydrogenCount = explicitHydrogenCount
    }
}

public struct Bond: Identifiable, Hashable, Codable {
    public let id: Int
    public let a1: Int
    public let a2: Int
    public var order: BondOrder
    public var stereo: BondStereo = .none

    public init(id: Int, a1: Int, a2: Int, order: BondOrder, stereo: BondStereo = .none) {
        self.id = id
        self.a1 = a1
        self.a2 = a2
        self.order = order
        self.stereo = stereo
    }
}

public struct Molecule: Hashable, Codable {
    public var name: String = "Untitled"
    public var atoms: [Atom] = []
    public var bonds: [Bond] = []

    public var atomCount: Int { atoms.count }
    public var bondCount: Int { bonds.count }

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
            return 3.0
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

    public init(name: String = "Untitled", atoms: [Atom] = [], bonds: [Bond] = []) {
        self.name = name
        self.atoms = atoms
        self.bonds = bonds
    }
}

public enum Depiction2DGenerator {
    public static func generate(for molecule: Molecule) -> Molecule {
        CDKStructureDiagramGenerator.apply(to: molecule)
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
