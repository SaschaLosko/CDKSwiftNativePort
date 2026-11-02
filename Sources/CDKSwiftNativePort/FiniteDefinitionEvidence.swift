import Foundation

public struct CDKFiniteDefinitionCompilation: Sendable {
    public let molecule: Molecule
    public let evidence: CDKFiniteDefinitionEvidence
}

/// Records how a compiled alternative was derived. These are selections from a reviewed
/// interpretation, not assertions that a complete compound was explicitly disclosed.
public struct CDKFiniteDefinitionEvidence: Codable, Equatable, Sendable {
    public struct Choice: Codable, Equatable, Sendable {
        public let definitionPath: String
        public let variableLabel: String
        public let sourceText: String
        public let sourceAlternativeIndex: Int?
        public let chemicalClass: CDKFiniteChemicalClass?
        public let fragmentSMILES: String?
        /// Atom ID in this choice's local scaffold, not the outer compound.
        public var siteAtomID: Int?
        public let omitted: Bool
    }
    public struct Variable: Codable, Equatable, Sendable {
        public let rGroupNumber: Int
        public let definitionPath: String
        public let label: String
        public let sourceText: String
        /// Indexed by the compiled alternative index, after exclusions and nested expansion.
        public let alternatives: [[Choice]]
    }
    public let sourceIdentifier: String
    public let variables: [Variable]

    /// Use only with assignments produced from this compilation's molecule.
    /// Missing or inconsistent evidence is an error, never an invented explanation.
    public func choices(for assignments: [CDKRGroupAssignment]) throws -> [Choice] {
        try assignments.flatMap { assignment -> [Choice] in
            guard let variable = variables.first(where: { $0.rGroupNumber == assignment.rGroupNumber }) else {
                throw CDKRGroupEnumerationError.invalidQuery("Assignment has no definition evidence.")
            }
            guard let index = assignment.alternativeIndex else {
                return [
                    .init(
                        definitionPath: variable.definitionPath, variableLabel: variable.label,
                        sourceText: variable.sourceText, sourceAlternativeIndex: nil, chemicalClass: nil,
                        fragmentSMILES: nil, siteAtomID: assignment.siteAtomID, omitted: true)
                ]
            }
            guard variable.alternatives.indices.contains(index), !variable.alternatives[index].isEmpty else {
                throw CDKRGroupEnumerationError.invalidQuery("Assignment has an unknown compiled alternative.")
            }
            var choices = variable.alternatives[index]
            choices[0].siteAtomID = assignment.siteAtomID
            return choices
        }
    }
}
