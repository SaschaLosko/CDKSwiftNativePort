import Foundation

public struct CDKMarkushMatchOptions: Sendable {
    public var maximumConfigurations: Int
    /// Bounds both the input graphs and each expanded member, before hydrogen suppression.
    public var maximumAtoms: Int
    /// Shared bound on enumeration checkpoints and graph-search steps.
    public var maximumSearchSteps: Int
    /// Cooperative monotonic deadline; this is not a hard real-time guarantee.
    public var timeoutSeconds: TimeInterval

    public init(
        maximumConfigurations: Int = 1_000, maximumAtoms: Int = 2_000,
        maximumSearchSteps: Int = 1_000_000, timeoutSeconds: TimeInterval = 30
    ) {
        self.maximumConfigurations = maximumConfigurations
        self.maximumAtoms = maximumAtoms
        self.maximumSearchSteps = maximumSearchSteps
        self.timeoutSeconds = timeoutSeconds
    }
}

public struct CDKMarkushMatchResult: Sendable {
    public enum Verdict: String, Codable, Sendable { case match, noMatch, unresolved }
    public enum Reason: String, Codable, Sendable {
        case invalidInput, unsupported, configurationLimit, atomLimit, searchLimit, timedOut, cancelled
        case unspecifiedStereochemistry
    }
    public struct Witness: Sendable {
        public let configuration: CDKMarkushConfiguration
        /// Member atom ID -> compound atom ID. Suppressed ordinary hydrogen vertices are omitted.
        public let compoundAtomIDByMemberAtomID: [Int: Int]
    }
    public let verdict: Verdict
    public let reason: Reason?
    public let detail: String?
    public let witness: Witness?
    public let configurationsChecked: Int
    /// True only after complete enumeration and comparison of every permitted assignment.
    public let isExhaustive: Bool
    public let searchSteps: Int
    public let variableNumberBySourceLabel: [String: Int]
}

/// Closed, finite compound-to-family membership, not substructure search.
/// Charges, isotopes, hydrogen counts and the entire graph must agree.
/// Declared tetrahedral and alkene stereo constrain membership; unspecified family stereo is free.
/// Missing compound stereo required by a member is unresolved. Unsupported chemistry is never a negative.
public enum CDKMarkushMatcher {
    public static func match(
        _ compound: Molecule, in family: Molecule, options: CDKMarkushMatchOptions = .init()
    ) -> CDKMarkushMatchResult {
        match(compound, in: family, budget: MarkushMatchBudget(options))
    }

    static func match(_ compound: Molecule, in family: Molecule, budget: MarkushMatchBudget) -> CDKMarkushMatchResult {
        let options = budget.options
        var checked = 0
        var labels: [String: Int] = [:]
        var witness: CDKMarkushMatchResult.Witness?
        var uncertainStereo = false
        func result(
            _ verdict: CDKMarkushMatchResult.Verdict, _ reason: CDKMarkushMatchResult.Reason? = nil,
            _ detail: String? = nil, exhaustive: Bool = false
        ) -> CDKMarkushMatchResult {
            CDKMarkushMatchResult(
                verdict: verdict, reason: reason, detail: detail, witness: witness,
                configurationsChecked: checked, isExhaustive: exhaustive, searchSteps: budget.steps,
                variableNumberBySourceLabel: labels)
        }
        do {
            guard options.maximumConfigurations > 0, options.maximumAtoms > 0,
                options.maximumSearchSteps > 0, options.timeoutSeconds.isFinite, options.timeoutSeconds > 0
            else { throw MatchFailure(.invalidInput, "Limits must be positive and finite.") }
            try budget.check()
            guard compound.atoms.count <= options.maximumAtoms, family.atoms.count <= options.maximumAtoms else {
                throw MatchFailure(.atomLimit, "An input exceeds the atom limit.")
            }
            _ = try MarkushStereoPlan(family)
            let resolution = try CDKMarkushVariableLabels.resolve(family)
            labels = resolution.numberBySourceLabel
            let target = try MatchingGraph(compound)
            try CDKMarkushEnumerator.enumerate(
                resolution.molecule, maximumConfigurations: options.maximumConfigurations,
                maximumAtoms: options.maximumAtoms, checkpoint: budget.check
            ) { configuration in
                let member = try MatchingGraph(configuration.molecule)
                checked += 1
                var comparison = FiniteGraphComparison(member: member, target: target, budget: budget)
                if let mapping = try comparison.findMatch() {
                    witness = .init(configuration: configuration, compoundAtomIDByMemberAtomID: mapping)
                    return false
                }
                uncertainStereo = uncertainStereo || comparison.uncertainStereo
                return true
            }
            if witness != nil { return result(.match) }
            if uncertainStereo {
                return result(
                    .unresolved, .unspecifiedStereochemistry,
                    "The compound omits stereochemistry required by a structurally compatible member.", exhaustive: true
                )
            }
            return result(.noMatch, exhaustive: true)
        } catch is CancellationError {
            return result(.unresolved, .cancelled, "Matching was cancelled.")
        } catch let failure as MatchFailure {
            return result(.unresolved, failure.reason, failure.detail)
        } catch let error as CDKRGroupEnumerationError {
            switch error {
            case .atomLimitExceeded:
                return result(.unresolved, .atomLimit, "A generated member exceeds the atom limit.")
            case .configurationLimitExceeded:
                return result(.unresolved, .configurationLimit, "The finite family was not exhausted.")
            case .invalidQuery(let detail): return result(.unresolved, .invalidInput, detail)
            case .unsupportedFeature(let detail): return result(.unresolved, .unsupported, detail)
            }
        } catch { return result(.unresolved, .invalidInput, String(describing: error)) }
    }
}

struct MatchFailure: Error {
    let reason: CDKMarkushMatchResult.Reason
    let detail: String
    init(_ reason: CDKMarkushMatchResult.Reason, _ detail: String) {
        self.reason = reason
        self.detail = detail
    }
}

final class MarkushMatchBudget {
    let options: CDKMarkushMatchOptions
    let start = ProcessInfo.processInfo.systemUptime
    var steps = 0
    init(_ options: CDKMarkushMatchOptions) { self.options = options }
    func check() throws {
        try Task.checkCancellation()
        guard ProcessInfo.processInfo.systemUptime - start < options.timeoutSeconds else {
            throw MatchFailure(.timedOut, "The cooperative matching deadline was reached.")
        }
        guard steps < options.maximumSearchSteps else {
            throw MatchFailure(.searchLimit, "The search-step budget was exhausted.")
        }
        steps += 1
    }
}

private struct MatchingGraph {
    struct AtomKey: Hashable {
        let element: String
        let isotope: Int?
        let charge: Int
        let hydrogens: Int
        let aromatic: Bool
    }
    struct Alkene {
        let left: Int
        let right: Int
        let leftLigand: Int
        let rightLigand: Int
        let cis: Bool
    }
    let molecule: Molecule
    let keys: [Int: AtomKey]
    let edges: [Int: [Int: BondOrder]]
    let tetrahedra: [Int: [Int]]
    let alkenes: [Alkene]

    init(_ input: Molecule) throws {
        let source = try MarkushStereoPlan(input, concrete: true).source
        try RGroupEnumerationPlan.validateGraph(source)
        guard !source.atoms.isEmpty, !CDKMarkushStructure.requiresSpecifiedStructure(source) else {
            throw MatchFailure(.unsupported, "A concrete compound or member is required.")
        }
        guard source.atoms.allSatisfy({ ($0.radical ?? 0) == 0 && $0.radicalType == nil }),
            source.bonds.allSatisfy({ $0.coordinateBondReferenceAtomID == nil })
        else { throw MatchFailure(.unsupported, "Radicals and coordinate bonds are not yet supported.") }
        // Suppress only neutral, unlabelled, terminal protium vertices. Isotopic H stays in the graph.
        var removed: Set<Int> = []
        for atom in source.atoms where atom.element == "H" && atom.isotopeMassNumber == nil && atom.charge == 0 {
            let bonds = source.bonds(forAtom: atom.id)
            if bonds.count == 1, bonds[0].order == .single,
                let neighbor = source.atom(id: bonds[0].a1 == atom.id ? bonds[0].a2 : bonds[0].a1),
                neighbor.element != "H", neighbor.element != "D", neighbor.element != "T"
            {
                removed.insert(atom.id)
            }
        }
        var graph = source
        graph.atoms.removeAll { removed.contains($0.id) }
        graph.bonds.removeAll { removed.contains($0.a1) || removed.contains($0.a2) }
        for i in graph.atoms.indices {
            let atom = graph.atoms[i]
            if atom.element == "D" || atom.element == "T" {
                graph.atoms[i].element = "H"
                graph.atoms[i].isotopeMassNumber = atom.element == "D" ? 2 : 3
            }
            graph.atoms[i].ligandOrderingAtomIDs = atom.ligandOrderingAtomIDs?.map {
                removed.contains($0) ? atom.id : $0
            }
        }
        let aromaticBonds = Set(graph.bonds.filter { $0.order == .aromatic }.map(\.id))
            .union(graph.aromaticDisplayBondIDs())
        let aromaticAtoms = Set(graph.bonds.filter { aromaticBonds.contains($0.id) }.flatMap { [$0.a1, $0.a2] })
        var keys: [Int: AtomKey] = [:]
        var edges: [Int: [Int: BondOrder]] = Dictionary(uniqueKeysWithValues: graph.atoms.map { ($0.id, [:]) })
        for atom in graph.atoms {
            let suppressed = source.neighbors(of: atom.id).filter { removed.contains($0) }.count
            keys[atom.id] = AtomKey(
                element: atom.element, isotope: atom.isotopeMassNumber, charge: atom.charge,
                hydrogens: CDKDescriptorSupport.implicitHydrogenCount(on: atom.id, in: source) + suppressed,
                aromatic: aromaticAtoms.contains(atom.id))
        }
        for bond in graph.bonds {
            guard edges[bond.a1]?[bond.a2] == nil else {
                throw MatchFailure(.invalidInput, "Parallel edges are not supported.")
            }
            let order: BondOrder = aromaticBonds.contains(bond.id) ? .aromatic : bond.order
            edges[bond.a1]?[bond.a2] = order
            edges[bond.a2]?[bond.a1] = order
        }
        var tetrahedra: [Int: [Int]] = [:]
        for atom in graph.atoms where atom.chirality != .none {
            guard let order = atom.ligandOrderingAtomIDs, order.count == 4, Set(order).count == 4,
                Set(order.filter { $0 != atom.id }) == Set(edges[atom.id, default: [:]].keys),
                edges[atom.id, default: [:]].values.allSatisfy({ $0 == .single }),
                !order.contains(atom.id) || keys[atom.id]?.hydrogens == 1
            else { throw MatchFailure(.unsupported, "Tetrahedral stereo requires four valid ordered ligands.") }
            tetrahedra[atom.id] = order
        }
        var alkenes: [Alkene] = []
        var interpretedDirectionalBonds: Set<Int> = []
        func direction(_ bond: Bond, awayFrom id: Int) -> Int? {
            let sign: Int
            switch bond.stereo {
            case .up, .downReversed: sign = 1
            case .down, .upReversed: sign = -1
            default: return nil
            }
            return bond.a1 == id ? sign : -sign
        }
        guard source.bonds.allSatisfy({ $0.doubleBondStereo == nil || $0.order == .double }) else {
            throw MatchFailure(.invalidInput, "Alkene stereo is attached to a non-double bond.")
        }
        for bond in source.bonds where bond.order == .double {
            var value: Alkene?
            if let stereo = bond.doubleBondStereo {
                guard let refs = bond.stereoReferenceAtomIDs, refs.count == 4,
                    Set(refs[1...2]) == Set([bond.a1, bond.a2]),
                    source.bond(between: refs[0], and: refs[1]) != nil,
                    source.bond(between: refs[2], and: refs[3]) != nil,
                    refs[0] != refs[2], refs[3] != refs[1]
                else { throw MatchFailure(.invalidInput, "Invalid double-bond stereo references.") }
                value = Alkene(
                    left: refs[1], right: refs[2],
                    leftLigand: removed.contains(refs[0]) ? refs[1] : refs[0],
                    rightLigand: removed.contains(refs[3]) ? refs[2] : refs[3], cis: stereo == .cis)
            } else {
                func isAlkeneDirection(_ side: Bond, center: Int) -> Bool {
                    guard side.order == .single, direction(side, awayFrom: center) != nil else { return false }
                    // The SMILES parser marks a complete directional pair on the double bond.
                    // Generated tetrahedral wedges alone must not invent alkene stereochemistry.
                    return bond.stereo == .either || (tetrahedra[side.a1] == nil && tetrahedra[side.a2] == nil)
                }
                let left = source.bonds(forAtom: bond.a1).filter { isAlkeneDirection($0, center: bond.a1) }
                let right = source.bonds(forAtom: bond.a2).filter { isAlkeneDirection($0, center: bond.a2) }
                if !left.isEmpty && !right.isEmpty {
                    for (side, center) in [(left, bond.a1), (right, bond.a2)] {
                        if side.count > 2
                            || (side.count == 2
                                && direction(side[0], awayFrom: center) == direction(side[1], awayFrom: center))
                        {
                            throw MatchFailure(.invalidInput, "Conflicting alkene directions.")
                        }
                        interpretedDirectionalBonds.formUnion(side.map(\.id))
                    }
                    let a = left[0].a1 == bond.a1 ? left[0].a2 : left[0].a1
                    let b = right[0].a1 == bond.a2 ? right[0].a2 : right[0].a1
                    value = Alkene(
                        left: bond.a1, right: bond.a2,
                        leftLigand: removed.contains(a) ? bond.a1 : a,
                        rightLigand: removed.contains(b) ? bond.a2 : b,
                        cis: direction(left[0], awayFrom: bond.a1) == direction(right[0], awayFrom: bond.a2))
                } else if bond.stereo != .none || !left.isEmpty || !right.isEmpty {
                    throw MatchFailure(.unsupported, "Incomplete or unknown alkene stereochemistry.")
                }
            }
            if let value {
                for center in [value.left, value.right] {
                    let neighbors = edges[center, default: [:]]
                    let coordination = neighbors.count + (keys[center]?.hydrogens ?? 0)
                    let ordinaryCarbon = keys[center]?.element == "C" && coordination == 3
                    let neutralImine = keys[center]?.element == "N" && keys[center]?.charge == 0 && coordination == 2
                    guard ordinaryCarbon || neutralImine,
                        keys[center]?.aromatic == false,
                        neighbors.values.filter({ $0 == .double }).count == 1
                    else {
                        throw MatchFailure(.unsupported, "Only ordinary alkene and neutral imine stereo is supported.")
                    }
                }
                alkenes.append(value)
            }
        }
        for bond in source.bonds where bond.order != .double && bond.stereo != .none {
            guard
                interpretedDirectionalBonds.contains(bond.id)
                    || (bond.stereo != .either && (tetrahedra[bond.a1] != nil || tetrahedra[bond.a2] != nil))
            else { throw MatchFailure(.unsupported, "Wedge-only or unknown stereo requires explicit ordered stereo.") }
        }
        molecule = graph
        self.keys = keys
        self.edges = edges
        self.tetrahedra = tetrahedra
        self.alkenes = alkenes
    }
}

private struct FiniteGraphComparison {
    let member: MatchingGraph
    let target: MatchingGraph
    let budget: MarkushMatchBudget
    var uncertainStereo = false
    var mapping: [Int: Int] = [:]
    var used: Set<Int> = []

    mutating func findMatch() throws -> [Int: Int]? {
        try budget.check()
        guard member.keys.count == target.keys.count,
            member.molecule.bonds.count == target.molecule.bonds.count,
            Dictionary(grouping: member.keys.values, by: { $0 }).mapValues(\.count)
                == Dictionary(grouping: target.keys.values, by: { $0 }).mapValues(\.count)
        else { return nil }
        return try search()
    }

    mutating func search() throws -> [Int: Int]? {
        try budget.check()
        if mapping.count == member.keys.count {
            switch stereoCompatibility() {
            case .match: return mapping
            case .unresolved: uncertainStereo = true
            case .noMatch: break
            }
            return nil
        }
        let id = member.keys.keys.filter { mapping[$0] == nil }.sorted { a, b in
            let connectedA = member.edges[a, default: [:]].keys.filter { mapping[$0] != nil }.count
            let connectedB = member.edges[b, default: [:]].keys.filter { mapping[$0] != nil }.count
            if connectedA != connectedB { return connectedA > connectedB }
            let degreeA = member.edges[a, default: [:]].count
            let degreeB = member.edges[b, default: [:]].count
            return degreeA == degreeB ? a < b : degreeA > degreeB
        }[0]
        for candidate in target.keys.keys.sorted() where !used.contains(candidate) {
            try budget.check()
            guard member.keys[id] == target.keys[candidate],
                member.edges[id]?.count == target.edges[candidate]?.count,
                mapping.allSatisfy({ member.edges[id]?[$0.key] == target.edges[candidate]?[$0.value] })
            else { continue }
            mapping[id] = candidate
            used.insert(candidate)
            if let found = try search() { return found }
            mapping[id] = nil
            used.remove(candidate)
        }
        return nil
    }

    func stereoCompatibility() -> CDKMarkushMatchResult.Verdict {
        var missing = false
        for (center, order) in member.tetrahedra {
            let mappedCenter = mapping[center]!
            guard let targetOrder = target.tetrahedra[mappedCenter] else { missing = true; continue }
            let mapped = order.map { mapping[$0]! }
            let permutation = mapped.compactMap { targetOrder.firstIndex(of: $0) }
            guard permutation.count == 4 else { return .noMatch }
            var inversions = 0
            for i in 0..<4 { for j in (i + 1)..<4 where permutation[i] > permutation[j] { inversions += 1 } }
            let sameWinding =
                member.molecule.atom(id: center)!.chirality == target.molecule.atom(id: mappedCenter)!.chirality
            if sameWinding != inversions.isMultiple(of: 2) { return .noMatch }
        }
        for alkene in member.alkenes {
            let left = mapping[alkene.left]!
            let right = mapping[alkene.right]!
            guard let other = target.alkenes.first(where: { Set([$0.left, $0.right]) == Set([left, right]) }) else {
                missing = true
                continue
            }
            let targetLeft = other.left == left ? other.leftLigand : other.rightLigand
            let targetRight = other.left == left ? other.rightLigand : other.leftLigand
            let flipLeft = mapping[alkene.leftLigand]! != targetLeft
            let flipRight = mapping[alkene.rightLigand]! != targetRight
            if (alkene.cis == other.cis) != (flipLeft == flipRight) { return .noMatch }
        }
        return missing ? .unresolved : .match
    }
}
