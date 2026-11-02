import Foundation

public struct CDKQueryMatchOptions: Sendable {
    public var maximumAtoms: Int
    public var maximumSearchSteps: Int
    public var timeoutSeconds: TimeInterval
    public init(maximumAtoms: Int = 2_000, maximumSearchSteps: Int = 1_000_000, timeoutSeconds: TimeInterval = 30) {
        self.maximumAtoms = maximumAtoms; self.maximumSearchSteps = maximumSearchSteps
        self.timeoutSeconds = timeoutSeconds
    }
    var budgetOptions: CDKMarkushMatchOptions {
        .init(maximumAtoms: maximumAtoms, maximumSearchSteps: maximumSearchSteps, timeoutSeconds: timeoutSeconds)
    }
}

public struct CDKQueryMatchResult: Sendable {
    public let verdict: CDKMarkushMatchResult.Verdict
    public let reason: CDKMarkushMatchResult.Reason?
    public let detail: String?
    /// Query atom IDs to target atom IDs, including resolved variable-attachment placeholders.
    public let atomMapping: [Int: Int]?
    public let isExhaustive: Bool
    public let searchSteps: Int
}

public enum CDKQueryMatcher {
    public enum Mode: Sendable { case substructure, wholeGraph }

    /// Matches a query graph. This is not closed finite R-group family membership.
    public static func match(
        _ query: CDKQueryPattern, in molecule: Molecule,
        mode: Mode = .substructure, targetStereo: [CDKQueryAtomStereo] = [],
        options: CDKQueryMatchOptions = .init()
    ) -> CDKQueryMatchResult {
        let budget = MarkushMatchBudget(options.budgetOptions)
        func result(
            _ verdict: CDKMarkushMatchResult.Verdict, _ reason: CDKMarkushMatchResult.Reason? = nil,
            _ detail: String? = nil, mapping: [Int: Int]? = nil, exhaustive: Bool = false
        ) -> CDKQueryMatchResult {
            .init(
                verdict: verdict, reason: reason, detail: detail, atomMapping: mapping,
                isExhaustive: exhaustive, searchSteps: budget.steps)
        }
        do {
            guard options.maximumAtoms > 0, options.maximumSearchSteps > 0,
                options.timeoutSeconds.isFinite, options.timeoutSeconds > 0
            else {
                throw CDKQueryError.invalid("Invalid matching limits.")
            }
            guard molecule.atoms.count <= options.maximumAtoms, query.molecule.atoms.count <= options.maximumAtoms
            else {
                return result(.unresolved, .atomLimit, "Input exceeds the atom limit.")
            }
            guard
                molecule.atoms.allSatisfy({
                    (QueryElements.number($0.element) ?? 0) > 0 || ["D", "T"].contains($0.element)
                }),
                molecule.atoms.allSatisfy({
                    $0.queryType == nil && $0.atomList == nil && $0.rGroupLabel == nil && !$0.isUnresolvedPseudoAtom
                        && $0.rGroupMembership == nil && $0.attachmentPoint == nil
                        && $0.substitutionCount == nil && $0.ringBondCount == nil && $0.unsaturated == nil
                }),
                molecule.bonds.allSatisfy({
                    $0.queryType == nil && $0.topology == nil && $0.coordinateBondReferenceAtomID == nil
                }),
                molecule.cxState?.linkNodes.isEmpty != false,
                molecule.rGroupLogicDefinitions.isEmpty,
                !molecule.sgroups.contains(where: {
                    [.structureRepeatUnit, .polymer, .extMulticenter].contains($0.kind)
                })
            else {
                throw CDKQueryError.unsupported("The target must have a specified molecular graph.")
            }
            let facts = try QueryGraphFacts(molecule, budget: budget)
            let stereo = try QueryStereoSupport.elements(molecule, explicit: targetStereo)
            var engine = try QuerySearch(query: query, target: facts, stereo: stereo, mode: mode, budget: budget)
            if let found = try engine.search() { return result(.match, mapping: found) }
            return engine.uncertain
                ? result(
                    .unresolved, .unspecifiedStereochemistry,
                    "At least one compatible mapping has unresolved stereo.", exhaustive: true)
                : result(.noMatch, exhaustive: true)
        } catch is CancellationError { return result(.unresolved, .cancelled) } catch let failure as MatchFailure {
            return result(.unresolved, failure.reason, failure.detail)
        } catch let failure as CDKQueryError {
            switch failure {
            case .invalid(let detail): return result(.unresolved, .invalidInput, detail)
            case .unsupported(let detail): return result(.unresolved, .unsupported, detail)
            case .limit(let detail): return result(.unresolved, .searchLimit, detail)
            }
        } catch { return result(.unresolved, .invalidInput, String(describing: error)) }
    }
}

struct QuerySearch {
    let query: CDKQueryPattern
    let target: QueryGraphFacts
    let targetStereo: [CDKQueryAtomStereo]
    let queryStereo: [CDKQueryAtomStereo]
    let mode: CDKQueryMatcher.Mode
    let budget: MarkushMatchBudget
    let depth: Int
    let positions: QueryVariableAttachments
    let ids: [Int]
    var mapping: [Int: Int] = [:]
    var used: Set<Int> = []
    var uncertain = false

    init(
        query: CDKQueryPattern, target: QueryGraphFacts, stereo: [CDKQueryAtomStereo],
        mode: CDKQueryMatcher.Mode, budget: MarkushMatchBudget, depth: Int = 0
    ) throws {
        try budget.check()
        guard depth <= 16, query.molecule.atoms.count <= budget.options.maximumAtoms else {
            throw CDKQueryError.limit("Recursive query depth or atom limit.")
        }
        try QueryGraphFacts.validate(query.molecule)
        guard Set(query.atoms.keys) == Set(query.molecule.atoms.map(\.id)),
            Set(query.bonds.keys) == Set(query.molecule.bonds.map(\.id)),
            query.molecule.rGroupLogicDefinitions.isEmpty,
            !query.molecule.sgroups.contains(where: { [.structureRepeatUnit, .polymer].contains($0.kind) })
        else {
            throw CDKQueryError.invalid("Incomplete expressions or unexpanded family semantics.")
        }
        self.query = query; self.target = target; targetStereo = stereo; self.mode = mode
        self.budget = budget; self.depth = depth
        queryStereo = try QueryStereoSupport.elements(query.molecule, explicit: query.stereo)
        let positions = try QueryVariableAttachments(query)
        self.positions = positions
        ids = query.molecule.atoms.map(\.id).filter { !positions.placeholders.contains($0) }.sorted()
        guard !ids.isEmpty else { throw CDKQueryError.invalid("A query requires non-placeholder atoms.") }
    }

    mutating func search() throws -> [Int: Int]? {
        try budget.check()
        if ids.count > target.atoms.count || (mode == .wholeGraph && ids.count != target.atoms.count) { return nil }
        if mapping.count == ids.count {
            return try positions.resolve(mapping, target: target, budget: budget) { resolved in
                var status = QueryTruth.yes
                for id in ids { status = status.and(try expression(query.atoms[id]!, at: resolved[id]!, level: 0)) }
                status = status.and(
                    try QueryStereoSupport.compare(
                        query: queryStereo, target: targetStereo, mapping: resolved,
                        targetMolecule: target.molecule))
                status = status.and(
                    try QueryStereoSupport.doubleBonds(query.molecule, target: target.molecule, mapping: resolved))
                if mode == .wholeGraph {
                    let mappedBonds = query.molecule.bonds.compactMap { bond -> Int? in
                        guard let a = resolved[bond.a1], let b = resolved[bond.a2] else { return nil }
                        return target.edges[a]?[b]?.id
                    }
                    if Set(mappedBonds).count != target.molecule.bonds.count { return false }
                }
                if status == .unknown { uncertain = true }
                return status == .yes
            }
        }
        let id = ids.filter { mapping[$0] == nil }.max { a, b in
            query.molecule.neighbors(of: a).filter { mapping[$0] != nil }.count
                < query.molecule.neighbors(of: b).filter { mapping[$0] != nil }.count
        }!
        for candidate in target.atoms.keys.sorted() where !used.contains(candidate) {
            try budget.check()
            if try expression(query.atoms[id]!, at: candidate, level: 0) == .no { continue }
            var compatible = true
            for edge in query.molecule.bonds(forAtom: id) where !positions.variableBondIDs.contains(edge.id) {
                let other = edge.a1 == id ? edge.a2 : edge.a1
                if let assigned = mapping[other] {
                    guard let bond = target.edges[candidate]?[assigned],
                        try target.bond(query.bonds[edge.id]!, bond, budget: budget)
                    else { compatible = false; break }
                }
            }
            if !compatible { continue }
            mapping[id] = candidate; used.insert(candidate)
            if let found = try search() { return found }
            mapping[id] = nil; used.remove(candidate)
        }
        return nil
    }

    func expression(_ value: CDKAtomExpression, at id: Int, level: Int) throws -> QueryTruth {
        try budget.check()
        guard level < 128 else { throw CDKQueryError.limit("Atom expression depth limit.") }
        switch value {
        case .and(let a, let b):
            let lhs = try expression(a, at: id, level: level + 1)
            return lhs == .no ? .no : lhs.and(try expression(b, at: id, level: level + 1))
        case .or(let a, let b):
            let lhs = try expression(a, at: id, level: level + 1)
            return lhs == .yes ? .yes : lhs.or(try expression(b, at: id, level: level + 1))
        case .not(let inner): return try expression(inner, at: id, level: level + 1).negated
        case .recursive(let pattern):
            var nested = try QuerySearch(
                query: pattern, target: target, stereo: targetStereo,
                mode: .substructure, budget: budget, depth: depth + 1)
            guard let root = pattern.molecule.atoms.first?.id, nested.ids.contains(root) else {
                throw CDKQueryError.invalid("Recursive root cannot be an attachment placeholder.")
            }
            nested.mapping[root] = id; nested.used.insert(id)
            if try nested.search() != nil { return .yes }
            return nested.uncertain ? .unknown : .no
        default: return try target.primitive(value, at: id, budget: budget)
        }
    }
}
