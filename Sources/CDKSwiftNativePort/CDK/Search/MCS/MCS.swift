import Foundation

/// Portable exact common-subgraph search, with bounded best-so-far results.
/// Default scoring is atoms then common bonds; disconnected scoring is bonds then atoms.
/// This matches specified graphs or explicit query predicates, not Markush families.
public struct CDKMCS: Sendable {
    public let query: Molecule
    public let options: CDKMCSOptions
    let pattern: CDKQueryPattern?

    public init(query: Molecule, options: CDKMCSOptions = .init()) {
        self.query = query
        self.options = options
        pattern = nil
    }

    /// Top-level SMARTS stereo is not checked; recursive predicates retain their own stereo rules.
    public init(query: CDKQueryPattern, options: CDKMCSOptions = .init()) {
        self.query = query.molecule
        self.options = options
        pattern = query
    }

    public func search(in target: Molecule) throws -> CDKMCSResult {
        try options.validate()
        guard query.atoms.count <= options.maximumAtoms, target.atoms.count <= options.maximumAtoms else {
            throw CDKMCSError.invalidInput("MCS input exceeds the configured atom limit.")
        }
        if pattern != nil
            && (options.matchStereochemistry || options.atomMatching != .symbol
                || options.ignoreBondOrder || options.matchFormalCharge || options.matchIsotopes
                || options.matchRingTopology)
        {
            throw CDKMCSError.invalidInput(
                "Query patterns own their matching predicates; molecule matching overrides are not applicable.")
        }
        let budget = MarkushMatchBudget(
            .init(
                maximumAtoms: options.maximumAtoms, maximumSearchSteps: options.maximumSearchSteps,
                timeoutSeconds: options.timeoutSeconds))
        var engine = try MCSEngine(query: query, target: target, pattern: pattern, options: options, budget: budget)
        var termination = CDKMCSTermination.completed
        do {
            try engine.prepare()
            try engine.visit(0)
        } catch is CancellationError {
            termination = .cancelled
        } catch let failure as MatchFailure {
            switch failure.reason {
            case .timedOut: termination = .timedOut
            case .searchLimit: termination = .searchLimit
            default: throw CDKMCSError.unsupported(failure.detail)
            }
        } catch let failure as CDKQueryError {
            switch failure {
            case .invalid(let detail): throw CDKMCSError.invalidInput(detail)
            case .unsupported(let detail): throw CDKMCSError.unsupported(detail)
            case .limit: termination = .searchLimit
            }
        }
        return CDKMCSResult(
            mappings: engine.best, mappedAtomCount: engine.bestAtoms, mappedBondCount: engine.bestBonds,
            termination: termination, mappingsTruncated: engine.truncated, searchSteps: budget.steps)
    }

    /// Require proof of maximum score; a mapping limit does not prevent that proof.
    public func match(in target: Molecule) throws -> CDKSubstructureMatch? {
        let result = try search(in: target)
        guard result.isMaximumProven else { throw CDKMCSError.incomplete(result.termination) }
        return result.mappings.first
    }
}

public enum CDKMCSTermination: String, Codable, Sendable {
    case completed, timedOut, searchLimit, cancelled
}

public enum CDKMCSError: Error, Equatable, Sendable {
    case invalidInput(String), unsupported(String), incomplete(CDKMCSTermination)
}

public struct CDKMCSResult: Sendable {
    public let mappings: [CDKSubstructureMatch]
    public let mappedAtomCount: Int
    public let mappedBondCount: Int
    public let termination: CDKMCSTermination
    /// More equal-best mappings exist than can be retained. Search still proves the best score if completed.
    public let mappingsTruncated: Bool
    public let searchSteps: Int
    public var isMaximumProven: Bool { termination == .completed }
    public var allMaximumMappingsReturned: Bool { isMaximumProven && !mappingsTruncated }
}

public struct CDKMCSOptions: Sendable {
    public enum AtomMatching: Sendable { case symbol, element, any }
    public let atomMatching: AtomMatching
    public let ignoreBondOrder: Bool
    public let matchFormalCharge: Bool
    public let matchIsotopes: Bool
    public let matchRingTopology: Bool
    public let completeRings: Bool
    public let matchStereochemistry: Bool
    /// Nil selects connected MCS. A positive value is the minimum common bonds per disconnected fragment.
    public let minimumBondsPerFragment: Int?
    public let maximumMappings: Int
    public let maximumAtoms: Int
    public let maximumSearchSteps: Int
    public let timeoutSeconds: TimeInterval

    public init(
        atomMatching: AtomMatching = .symbol, ignoreBondOrder: Bool = false,
        matchFormalCharge: Bool = false, matchIsotopes: Bool = false, matchRingTopology: Bool = false,
        completeRings: Bool = false, matchStereochemistry: Bool = false, minimumBondsPerFragment: Int? = nil,
        maximumMappings: Int = 1_000, maximumAtoms: Int = 256,
        maximumSearchSteps: Int = 1_000_000, timeoutSeconds: TimeInterval = 30
    ) {
        self.atomMatching = atomMatching; self.ignoreBondOrder = ignoreBondOrder
        self.matchFormalCharge = matchFormalCharge; self.matchIsotopes = matchIsotopes
        self.matchRingTopology = matchRingTopology; self.completeRings = completeRings
        self.matchStereochemistry = matchStereochemistry; self.minimumBondsPerFragment = minimumBondsPerFragment
        self.maximumMappings = maximumMappings; self.maximumAtoms = maximumAtoms
        self.maximumSearchSteps = maximumSearchSteps; self.timeoutSeconds = timeoutSeconds
    }

    func validate() throws {
        guard (1...256).contains(maximumAtoms), maximumSearchSteps > 0,
            (1...100_000).contains(maximumMappings), timeoutSeconds.isFinite, timeoutSeconds > 0,
            minimumBondsPerFragment.map({ $0 > 0 }) ?? true
        else { throw CDKMCSError.invalidInput("Invalid MCS resource or fragment limits.") }
    }
}
