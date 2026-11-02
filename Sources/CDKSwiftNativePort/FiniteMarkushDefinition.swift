import Foundation

/// A reviewed, explicit interpretation of a source. Text is retained as evidence,
/// not automatically interpreted as chemistry. Unknown conditions block compilation.
public struct CDKFiniteMarkushDefinition: Codable, Equatable, Sendable {
    public struct Variable: Codable, Equatable, Sendable {
        public var label: String
        public var sourceText: String
        /// Each SMILES starts at its attachment atom. Nested labelled variables are allowed.
        public var alternatives: [String]
        public var excludedAlternativeIndices: [Int]
        public var occurrence: String
        public var requiresLabel: String?
        public var restH: Bool
        /// Finite classes supplement the explicit alternatives; exclusions index only explicit alternatives.
        public var chemicalClasses: [CDKFiniteChemicalClass]?
        /// Definitions local to one explicit alternative. They shadow enclosing definitions.
        public var localDefinitions: [AlternativeScope]?

        public struct AlternativeScope: Codable, Equatable, Sendable {
            public var alternativeIndex: Int
            public var variables: [Variable]

            public init(alternativeIndex: Int, variables: [Variable]) {
                self.alternativeIndex = alternativeIndex
                self.variables = variables
            }
        }

        public init(
            label: String, sourceText: String, alternatives: [String],
            excludedAlternativeIndices: [Int] = [], occurrence: String = ">0",
            requiresLabel: String? = nil, restH: Bool = false,
            chemicalClasses: [CDKFiniteChemicalClass]? = nil, localDefinitions: [AlternativeScope]? = nil
        ) {
            self.label = label; self.sourceText = sourceText; self.alternatives = alternatives
            self.excludedAlternativeIndices = excludedAlternativeIndices; self.occurrence = occurrence
            self.requiresLabel = requiresLabel; self.restH = restH
            self.chemicalClasses = chemicalClasses; self.localDefinitions = localDefinitions
        }
    }
    public var sourceIdentifier: String
    public var backbone: String
    public var reviewed: Bool
    public var unresolvedConditions: [String]
    public var variables: [Variable]

    public init(
        sourceIdentifier: String, backbone: String, reviewed: Bool = false,
        unresolvedConditions: [String] = [], variables: [Variable]
    ) {
        self.sourceIdentifier = sourceIdentifier; self.backbone = backbone
        self.reviewed = reviewed; self.unresolvedConditions = unresolvedConditions; self.variables = variables
    }

    /// Compile finite definitions with nested variables into the existing portable graph model.
    /// Source text and excluded alternatives remain in this sidecar, suitable for independent review.
    public func compile(maximumAlternatives: Int = 1_000, maximumAtoms: Int = 2_000) throws -> Molecule {
        try compileWithEvidence(maximumAlternatives: maximumAlternatives, maximumAtoms: maximumAtoms).molecule
    }

    public func compileWithEvidence(
        maximumAlternatives: Int = 1_000, maximumAtoms: Int = 2_000
    ) throws -> CDKFiniteDefinitionCompilation {
        var compiler = FiniteDefinitionCompiler(maximumAlternatives: maximumAlternatives, maximumAtoms: maximumAtoms)
        return try compiler.compile(self)
    }
}
