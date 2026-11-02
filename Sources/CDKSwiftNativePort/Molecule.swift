import Foundation

#if canImport(CoreGraphics)
    import CoreGraphics
#endif

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
    /// `true` when coordinates were synthesized for depiction rather than supplied by the source.
    /// Identity generators must not infer otherwise unspecified stereochemistry from these coordinates.
    public var coordinatesAreGenerated: Bool? = nil

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
        case coordinatesAreGenerated
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

    /// Returns the size of the smallest ring containing an atom, or `nil` when
    /// no ring up to `maxSize` contains it.
    public func smallestRingSize(
        containingAtomID atomID: Int,
        maxSize: Int? = nil
    ) -> Int? {
        guard atoms.contains(where: { $0.id == atomID }) else { return nil }
        let limit = max(3, min(maxSize ?? atoms.count, atoms.count))
        return simpleCycles(maxSize: limit)
            .lazy
            .filter { $0.contains(atomID) }
            .map(\.count)
            .min()
    }

    /// Returns the size of the smallest ring containing a bond, or `nil` when
    /// no ring up to `maxSize` contains it.
    public func smallestRingSize(
        containingBondID bondID: Int,
        maxSize: Int? = nil
    ) -> Int? {
        guard let bond = bonds.first(where: { $0.id == bondID }) else { return nil }
        let limit = max(3, min(maxSize ?? atoms.count, atoms.count))
        return simpleCycles(maxSize: limit)
            .lazy
            .filter { cycle in
                cycle.indices.contains { index in
                    let first = cycle[index]
                    let second = cycle[(index + 1) % cycle.count]
                    return (first == bond.a1 && second == bond.a2)
                        || (first == bond.a2 && second == bond.a1)
                }
            }
            .map(\.count)
            .min()
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
            .reduce(0.0) { partial, bond in
                partial + (bond.coordinateBondReferenceAtomID == nil ? bond.order.valenceContribution : 0)
            }

        return max(0, Int(round(targetValence - bondOrderSum)))
    }

    public mutating func assignWedgeHashFromChiralCenters() {
        let degreeByAtom = Dictionary(
            uniqueKeysWithValues: atoms.map { ($0.id, neighbors(of: $0.id).count) })
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
                let lc = bondStereoClearance(
                    for: bonds[lhs], around: atom.id, positionsByAtom: positionsByAtom)
                let rc = bondStereoClearance(
                    for: bonds[rhs], around: atom.id, positionsByAtom: positionsByAtom)
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
        return CGRect(
            x: minX, y: minY, width: max(0.0001, maxX - minX), height: max(0.0001, maxY - minY))
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

    private func bondStereoPriority(for bond: Bond, around atomID: Int, degreeByAtom: [Int: Int])
        -> Int
    {
        let neighborID = (bond.a1 == atomID) ? bond.a2 : bond.a1
        let neighborDegree = degreeByAtom[neighborID] ?? 0
        // Prefer terminal substituents; then stable ordering by atom id.
        return (neighborDegree == 1 ? 0 : 10) + neighborID
    }

    private func bondStereoClearance(
        for bond: Bond,
        around atomID: Int,
        positionsByAtom: [Int: CGPoint]
    ) -> CGFloat {
        let neighborID = (bond.a1 == atomID) ? bond.a2 : bond.a1
        guard let center = positionsByAtom[atomID], let neighbor = positionsByAtom[neighborID] else {
            return 0
        }
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

    public init(
        name: String = "Untitled",
        externalID: String? = nil,
        atoms: [Atom] = [],
        bonds: [Bond] = []
    ) {
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
        highlightedAtomIDs =
            try container.decodeIfPresent([Int].self, forKey: .highlightedAtomIDs) ?? []
        highlightedBondIDs =
            try container.decodeIfPresent([Int].self, forKey: .highlightedBondIDs) ?? []
        cxState = try container.decodeIfPresent(CDKCxSmilesState.self, forKey: .cxState)
        coordinatesAreGenerated = try container.decodeIfPresent(Bool.self, forKey: .coordinatesAreGenerated)
        dataFields = try container.decodeIfPresent([String: [String]].self, forKey: .dataFields) ?? [:]
        dataFieldOrder = try container.decodeIfPresent([String].self, forKey: .dataFieldOrder) ?? []
        dataFieldOrder = Molecule.normalizedDataFieldOrder(
            preferredOrder: dataFieldOrder, availableFields: dataFields)
        rGroupLogicDefinitions =
            try container.decodeIfPresent(
                [Int: MoleculeRGroupLogic].self,
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
        try container.encodeIfPresent(coordinatesAreGenerated, forKey: .coordinatesAreGenerated)
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

    private static func normalizedDataFieldOrder(
        preferredOrder: [String],
        availableFields: [String: [String]]
    ) -> [String] {
        var orderedNames: [String] = []
        orderedNames.reserveCapacity(availableFields.count)
        var seen = Set<String>()

        for rawName in preferredOrder {
            let normalized = normalizedDataFieldName(rawName)
            guard !normalized.isEmpty,
                availableFields[normalized] != nil,
                seen.insert(normalized).inserted
            else {
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
