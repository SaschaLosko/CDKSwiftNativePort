import Foundation

public enum CDKMarkushRelation: String, Codable, Sendable {
    case overlaps, containedIn, equivalent, disjoint
}

public struct CDKMarkushRelationResult: Sendable {
    public enum Verdict: String, Codable, Sendable { case holds, doesNotHold, unresolved }
    public enum Direction: String, Codable, Sendable { case leftToRight, rightToLeft }
    public let verdict: Verdict
    public let reason: CDKMarkushMatchResult.Reason?
    public let detail: String?
    /// Shared member for overlap; an uncovered source member for failed containment.
    public let witness: CDKMarkushConfiguration?
    public let membershipWitness: CDKMarkushMatchResult.Witness?
    public let direction: Direction
    public let sourceConfigurationsChecked: Int
    public let searchSteps: Int
    public let isExhaustive: Bool
}

/// Conservative relations between closed finite families. Every traversal shares one deadline
/// and search-step budget. A limit or unsupported feature cannot establish a negative result.
public enum CDKMarkushRelations {
    public static func compare(
        _ left: Molecule, _ right: Molecule, relation: CDKMarkushRelation,
        options: CDKMarkushMatchOptions = .init()
    ) -> CDKMarkushRelationResult {
        let budget = MarkushMatchBudget(options)
        var checked = 0
        var direction = CDKMarkushRelationResult.Direction.leftToRight
        var witness: CDKMarkushConfiguration?
        var membership: CDKMarkushMatchResult.Witness?
        func result(
            _ verdict: CDKMarkushRelationResult.Verdict,
            _ reason: CDKMarkushMatchResult.Reason? = nil, _ detail: String? = nil,
            exhaustive: Bool = false
        ) -> CDKMarkushRelationResult {
            .init(
                verdict: verdict, reason: reason, detail: detail, witness: witness,
                membershipWitness: membership, direction: direction,
                sourceConfigurationsChecked: checked, searchSteps: budget.steps, isExhaustive: exhaustive)
        }
        // Returns a decisive overlap/counterexample, or nil after exhaustion.
        func traverse(_ source: Molecule, _ target: Molecule, overlap: Bool) throws -> Bool? {
            var uncertainty: CDKMarkushMatchResult?
            var decisive = false
            try CDKMarkushEnumerator.enumerate(
                source, maximumConfigurations: options.maximumConfigurations,
                maximumAtoms: options.maximumAtoms, checkpoint: budget.check
            ) { member in
                checked += 1
                let match = CDKMarkushMatcher.match(member.molecule, in: target, budget: budget)
                if (overlap && match.verdict == .match) || (!overlap && match.verdict == .noMatch) {
                    witness = member
                    membership = match.witness
                    decisive = true
                    return false
                }
                if match.verdict == .unresolved {
                    if [.searchLimit, .timedOut, .cancelled, .invalidInput].contains(match.reason) {
                        throw MatchFailure(match.reason!, match.detail ?? "Comparison was interrupted.")
                    }
                    uncertainty = uncertainty ?? match
                }
                return true
            }
            if decisive { return true }
            if let uncertainty {
                throw MatchFailure(
                    uncertainty.reason ?? .unsupported, uncertainty.detail ?? "A member could not be compared.")
            }
            return nil
        }
        do {
            guard options.maximumConfigurations > 0, options.maximumAtoms > 0,
                options.maximumSearchSteps > 0, options.timeoutSeconds.isFinite, options.timeoutSeconds > 0
            else { throw MatchFailure(.invalidInput, "Limits must be positive and finite.") }
            guard left.atoms.count <= options.maximumAtoms, right.atoms.count <= options.maximumAtoms else {
                throw MatchFailure(.atomLimit, "An input exceeds the atom limit.")
            }
            try budget.check()
            let overlap = relation == .overlaps || relation == .disjoint
            let found: Bool?
            do { found = try traverse(left, right, overlap: overlap) } catch let failure as MatchFailure
                where overlap && failure.reason == .unspecifiedStereochemistry
            {
                // A specified member in the other direction may resolve an existential query.
                direction = .rightToLeft
                found = try traverse(right, left, overlap: true)
            }
            if found == true {
                return result(relation == .overlaps ? .holds : .doesNotHold)
            }
            if relation == .equivalent {
                direction = .rightToLeft
                if try traverse(right, left, overlap: false) == true { return result(.doesNotHold) }
            }
            return result(relation == .overlaps ? .doesNotHold : .holds, exhaustive: true)
        } catch is CancellationError {
            return result(.unresolved, .cancelled, "Comparison was cancelled.")
        } catch let failure as MatchFailure {
            return result(.unresolved, failure.reason, failure.detail)
        } catch let error as CDKRGroupEnumerationError {
            switch error {
            case .atomLimitExceeded:
                return result(.unresolved, .atomLimit, "A generated member exceeds the atom limit.")
            case .configurationLimitExceeded:
                return result(.unresolved, .configurationLimit, "A finite family was not exhausted.")
            case .invalidQuery(let detail): return result(.unresolved, .invalidInput, detail)
            case .unsupportedFeature(let detail): return result(.unresolved, .unsupported, detail)
            }
        } catch { return result(.unresolved, .invalidInput, String(describing: error)) }
    }
}
