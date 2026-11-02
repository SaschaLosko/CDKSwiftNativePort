import Foundation

final class QueryGraphFacts {
    let molecule: Molecule
    let atoms: [Int: Atom]
    let edges: [Int: [Int: Bond]]
    let aromaticBonds: Set<Int>
    let aromaticAtoms: Set<Int>
    let implicitHydrogens: [Int: Int]
    private var ringBasis: [CDKDescriptorRingSupport.Cycle]?
    private var localizedValences: [Int: Int]?

    init(_ molecule: Molecule, budget: MarkushMatchBudget) throws {
        try budget.check()
        try Self.validate(molecule)
        self.molecule = molecule
        atoms = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0) })
        var edges: [Int: [Int: Bond]] = [:]
        for bond in molecule.bonds {
            edges[bond.a1, default: [:]][bond.a2] = bond; edges[bond.a2, default: [:]][bond.a1] = bond
        }
        self.edges = edges
        let aromaticBonds = CDKDescriptorRingSupport.aromaticBondIDs(in: molecule, perceiveAromaticity: true)
        self.aromaticBonds = aromaticBonds
        aromaticAtoms = Set(molecule.atoms.filter(\.aromatic).map(\.id))
            .union(molecule.bonds.filter { aromaticBonds.contains($0.id) }.flatMap { [$0.a1, $0.a2] })
        implicitHydrogens = Dictionary(
            uniqueKeysWithValues: molecule.atoms.map {
                ($0.id, CDKDescriptorSupport.implicitHydrogenCount(on: $0.id, in: molecule))
            })
        try budget.check()
    }
    static func validate(_ molecule: Molecule) throws {
        let ids = Set(molecule.atoms.map(\.id))
        guard !ids.isEmpty, ids.count == molecule.atoms.count,
            Set(molecule.bonds.map(\.id)).count == molecule.bonds.count
        else {
            throw CDKQueryError.invalid("Empty graph or duplicate IDs.")
        }
        var edges: Set<[Int]> = []
        for bond in molecule.bonds {
            guard ids.contains(bond.a1), ids.contains(bond.a2), bond.a1 != bond.a2,
                edges.insert([min(bond.a1, bond.a2), max(bond.a1, bond.a2)]).inserted
            else {
                throw CDKQueryError.invalid("Invalid endpoints, self bond or duplicate edge.")
            }
        }
    }
    func rings(_ budget: MarkushMatchBudget) throws -> [CDKDescriptorRingSupport.Cycle] {
        if let ringBasis { return ringBasis }
        try budget.check()
        let rings = CDKDescriptorRingSupport.smallestRingBasis(in: molecule)
        try budget.check(); ringBasis = rings
        return rings
    }
    func bond(_ expression: CDKBondExpression, _ bond: Bond, budget: MarkushMatchBudget, depth: Int = 0) throws -> Bool
    {
        try budget.check()
        guard depth < 128 else { throw CDKQueryError.limit("Bond expression depth limit.") }
        switch expression {
        case .any: return true
        case .none: return false
        case .order(let order): return (aromaticBonds.contains(bond.id) ? .aromatic : bond.order) == order
        case .inRing: return try rings(budget).contains { $0.bondIDs.contains(bond.id) }
        case .and(let a, let b):
            return try self.bond(a, bond, budget: budget, depth: depth + 1)
                && self.bond(b, bond, budget: budget, depth: depth + 1)
        case .or(let a, let b):
            return try self.bond(a, bond, budget: budget, depth: depth + 1)
                || self.bond(b, bond, budget: budget, depth: depth + 1)
        case .not(let value): return try !self.bond(value, bond, budget: budget, depth: depth + 1)
        }
    }
    func primitive(_ expression: CDKAtomExpression, at id: Int, budget: MarkushMatchBudget) throws -> QueryTruth {
        let atom = atoms[id]!
        let adjacent = edges[id, default: [:]]
        let implicitH = implicitHydrogens[id, default: 0]
        let h = implicitH + adjacent.keys.filter { ["H", "D", "T"].contains(atoms[$0]!.element) }.count
        switch expression {
        case .any: return .yes
        case .none: return .no
        case .element(let value):
            return .init(QueryElements.number(["D", "T"].contains(atom.element) ? "H" : atom.element) == value)
        case .aromatic(let value): return .init(aromaticAtoms.contains(id) == value)
        case .isotope(let value):
            return .init((atom.isotopeMassNumber ?? (atom.element == "D" ? 2 : atom.element == "T" ? 3 : 0)) == value)
        case .charge(let value): return .init(atom.charge == value)
        case .degree(let value): return .init(adjacent.count == value)
        case .heavyDegree(let value):
            return .init(adjacent.keys.filter { !["H", "D", "T"].contains(atoms[$0]!.element) }.count == value)
        case .totalDegree(let value): return .init(adjacent.count + implicitH == value)
        case .hydrogens(let value): return .init(h == value)
        case .implicitHydrogens(let value): return .init(implicitH == value)
        case .hasImplicitHydrogen: return .init(implicitH > 0)
        case .valence(let value):
            if localizedValences == nil {
                try budget.check()
                let localized: Molecule
                do { localized = try CDKKekulization.kekulized(molecule) } catch {
                    throw CDKQueryError.unsupported("Valence query requires a valid localized bond assignment.")
                }
                localizedValences = Dictionary(
                    uniqueKeysWithValues: localized.atoms.map { atom in
                        (
                            atom.id,
                            localized.bonds(forAtom: atom.id).reduce(implicitHydrogens[atom.id, default: 0]) {
                                $0 + $1.order.rawValue
                            }
                        )
                    })
                try budget.check()
            }
            return .init(localizedValences?[id] == value)
        case .inRing: return .init(try rings(budget).contains { $0.atoms.contains(id) })
        case .ringCount(let value): return .init(try rings(budget).filter { $0.atoms.contains(id) }.count == value)
        case .smallestRing(let value):
            return .init(try rings(budget).filter { $0.atoms.contains(id) }.map(\.size).min() == value)
        case .ringBonds(let value):
            let ringBonds = Set(try rings(budget).flatMap(\.bondIDs))
            return .init(adjacent.values.filter { ringBonds.contains($0.id) }.count == value)
        case .and, .or, .not, .recursive: preconditionFailure("Logical predicates are evaluated by the query engine.")
        }
    }
}
