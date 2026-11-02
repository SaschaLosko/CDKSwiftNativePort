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

        public init(
            label: String, sourceText: String, alternatives: [String],
            excludedAlternativeIndices: [Int] = [], occurrence: String = ">0",
            requiresLabel: String? = nil, restH: Bool = false
        ) {
            self.label = label; self.sourceText = sourceText; self.alternatives = alternatives
            self.excludedAlternativeIndices = excludedAlternativeIndices; self.occurrence = occurrence
            self.requiresLabel = requiresLabel; self.restH = restH
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
        guard reviewed, !sourceIdentifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
            unresolvedConditions.isEmpty
        else {
            throw CDKRGroupEnumerationError.unsupportedFeature("Source definitions or conditions still require review.")
        }
        guard maximumAlternatives > 0, maximumAtoms > 0,
            Set(variables.map(\.label)).count == variables.count,
            variables.allSatisfy({
                !$0.sourceText.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty
                    && !$0.label.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty
            })
        else {
            throw CDKRGroupEnumerationError.invalidQuery("Definition labels, provenance and limits must be valid.")
        }
        let definitions = Dictionary(uniqueKeysWithValues: variables.map { ($0.label, $0) })
        let parser = CDKSmilesParser(flavor: .cdkDefault, generateCoordinates: false)
        var produced = 0
        var retainedAtoms = 0
        func query(for input: Molecule, path: [String]) throws -> CDKRGroupQuery {
            try Task.checkCancellation()
            guard input.atoms.count <= maximumAtoms, path.count <= 16,
                input.atoms.allSatisfy({ $0.rGroupMembership == nil })
            else {
                throw CDKRGroupEnumerationError.unsupportedFeature(
                    "Nested definition depth, graph size or embedded definitions are unsupported.")
            }
            let resolution = try CDKMarkushVariableLabels.resolve(input)
            var lists: [Int: CDKRGroupList] = [:]
            for (label, number) in resolution.numberBySourceLabel.sorted(by: { $0.key < $1.key }) {
                guard !path.contains(label), let definition = definitions[label] else {
                    throw CDKRGroupEnumerationError.invalidQuery("Missing or cyclic definition for \(label).")
                }
                guard definition.excludedAlternativeIndices.allSatisfy(definition.alternatives.indices.contains) else {
                    throw CDKRGroupEnumerationError.invalidQuery("An exclusion refers to an unknown alternative.")
                }
                var groups: [CDKRGroup] = []
                for (index, text) in definition.alternatives.enumerated()
                where !definition.excludedAlternativeIndices.contains(index) {
                    let fragment = try parser.parseSmiles(text)
                    guard let attachmentID = fragment.atoms.first?.id else {
                        throw CDKRGroupEnumerationError.invalidQuery("An alternative is empty.")
                    }
                    let nested = try query(for: fragment, path: path + [label])
                    func retain(_ member: Molecule) throws {
                        guard member.atom(id: attachmentID) != nil, produced < maximumAlternatives else {
                            throw CDKRGroupEnumerationError.unsupportedFeature(
                                "Alternative expansion exceeds its limit or replaces the attachment atom.")
                        }
                        guard member.atoms.count <= maximumAtoms - retainedAtoms else {
                            throw CDKRGroupEnumerationError.atomLimitExceeded(maximumAtoms)
                        }
                        produced += 1
                        retainedAtoms += member.atoms.count
                        groups.append(.init(group: member, firstAttachmentPointAtomID: attachmentID))
                    }
                    if nested.rGroupDefinitions.isEmpty {
                        try retain(fragment)
                    } else {
                        try nested.enumerateConfigurations(
                            maximumConfigurations: maximumAlternatives,
                            maximumAtoms: maximumAtoms
                        ) { member in
                            try retain(member)
                            return true
                        }
                    }
                }
                guard !groups.isEmpty else {
                    throw CDKRGroupEnumerationError.invalidQuery("A definition has no permitted alternatives.")
                }
                let required: Int
                if let label = definition.requiresLabel {
                    guard let n = resolution.numberBySourceLabel[label], n != 0 else {
                        throw CDKRGroupEnumerationError.unsupportedFeature(
                            "A dependency must refer to a nonzero variable on the same scaffold.")
                    }
                    required = n
                } else {
                    required = 0
                }
                lists[number] = try .init(
                    rGroupNumber: number, restH: definition.restH,
                    occurrence: definition.occurrence, requiredRGroupNumber: required, rGroups: groups)
            }
            return .init(rootStructure: resolution.molecule, rGroupDefinitions: lists)
        }
        let input = try parser.parseSmiles(backbone)
        let compiledQuery = try query(for: input, path: [])
        guard input.atoms.count <= maximumAtoms - retainedAtoms else {
            throw CDKRGroupEnumerationError.atomLimitExceeded(maximumAtoms)
        }
        let compiled = CDKRGroupQueryManipulator.toFlatMolecule(compiledQuery)
        guard compiled.atoms.count <= maximumAtoms else {
            throw CDKRGroupEnumerationError.unsupportedFeature("Compiled definitions exceed the atom limit.")
        }
        return compiled
    }
}
