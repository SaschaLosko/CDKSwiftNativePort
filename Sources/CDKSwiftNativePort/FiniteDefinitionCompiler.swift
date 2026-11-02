import Foundation

/// Lexical scopes are attached to explicit source alternatives. A referenced enclosing
/// definition keeps its declaration scope; a caller cannot accidentally redefine its variables.
struct FiniteDefinitionCompiler {
    typealias Definition = CDKFiniteMarkushDefinition
    typealias Evidence = CDKFiniteDefinitionEvidence
    struct Binding {
        let variable: Definition.Variable
        let path: String
        let scope: String
    }
    struct Scope {
        let parent: String?
        let bindings: [String: Binding]
    }
    struct Query {
        let query: CDKRGroupQuery
        let evidence: Evidence
    }
    let maximumAlternatives: Int
    let maximumAtoms: Int
    private var scopes: [String: Scope] = [:]
    private var produced = 0
    private var retainedAtoms = 0
    private var declaredVariables = 0
    private var declaredAlternatives = 0
    private let parser = CDKSmilesParser(flavor: .cdkDefault, generateCoordinates: false)

    init(maximumAlternatives: Int, maximumAtoms: Int) {
        self.maximumAlternatives = maximumAlternatives; self.maximumAtoms = maximumAtoms
    }

    mutating func compile(_ definition: Definition) throws -> CDKFiniteDefinitionCompilation {
        guard definition.reviewed, !definition.sourceIdentifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
            definition.unresolvedConditions.isEmpty
        else {
            throw CDKRGroupEnumerationError.unsupportedFeature("Source definitions or conditions still require review.")
        }
        guard maximumAlternatives > 0, maximumAtoms > 0 else {
            throw CDKRGroupEnumerationError.invalidQuery("Definition limits must be positive.")
        }
        try register(definition.variables, scope: "/variables", parent: nil, depth: 0)
        let input = try parse(definition.backbone)
        let result = try query(input, scope: "/variables", path: [])
        guard input.atoms.count <= maximumAtoms - retainedAtoms else {
            throw CDKRGroupEnumerationError.atomLimitExceeded(maximumAtoms)
        }
        let molecule = CDKRGroupQueryManipulator.toFlatMolecule(result.query)
        guard molecule.atoms.count <= maximumAtoms else {
            throw CDKRGroupEnumerationError.atomLimitExceeded(maximumAtoms)
        }
        return .init(
            molecule: molecule,
            evidence: .init(sourceIdentifier: definition.sourceIdentifier, variables: result.evidence.variables))
    }

    private mutating func register(_ variables: [Definition.Variable], scope: String, parent: String?, depth: Int)
        throws
    {
        try Task.checkCancellation()
        guard depth <= 16, variables.count <= 128 - declaredVariables,
            Set(variables.map(\.label)).count == variables.count
        else {
            throw CDKRGroupEnumerationError.invalidQuery("Duplicate labels or excessive definition scope.")
        }
        declaredVariables += variables.count
        var bindings: [String: Binding] = [:]
        for (index, variable) in variables.enumerated() {
            guard !variable.label.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
                !variable.sourceText.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
                variable.alternatives.count <= 1_000 - declaredAlternatives
            else {
                throw CDKRGroupEnumerationError.invalidQuery("Definition labels, provenance and limits must be valid.")
            }
            declaredAlternatives += variable.alternatives.count
            let classes = variable.chemicalClasses ?? []
            guard classes.count <= 1_000 - declaredAlternatives else {
                throw CDKRGroupEnumerationError.invalidQuery("Too many chemical classes.")
            }
            declaredAlternatives += classes.count
            for chemicalClass in classes { _ = try chemicalClass.fragments() }
            let locals = variable.localDefinitions ?? []
            guard locals.count <= variable.alternatives.count,
                Set(locals.map(\.alternativeIndex)).count == locals.count,
                locals.allSatisfy({ variable.alternatives.indices.contains($0.alternativeIndex) }),
                variable.excludedAlternativeIndices.allSatisfy(variable.alternatives.indices.contains)
            else {
                throw CDKRGroupEnumerationError.invalidQuery(
                    "An exclusion or local scope refers to an unknown alternative.")
            }
            let path = "\(scope)/\(index)"
            bindings[variable.label] = Binding(variable: variable, path: path, scope: scope)
            for local in locals {
                try register(
                    local.variables, scope: "\(path)/alternatives/\(local.alternativeIndex)",
                    parent: scope, depth: depth + 1)
            }
        }
        scopes[scope] = Scope(parent: parent, bindings: bindings)
    }

    private func binding(_ label: String, in scope: String) -> Binding? {
        guard let current = scopes[scope] else { return nil }
        return current.bindings[label] ?? current.parent.flatMap { binding(label, in: $0) }
    }

    private func parse(_ text: String) throws -> Molecule {
        try Task.checkCancellation()
        guard text.utf8.count <= 32_768 else {
            throw CDKRGroupEnumerationError.unsupportedFeature("A definition fragment exceeds the text limit.")
        }
        return try parser.parseSmiles(text)
    }

    private mutating func query(_ input: Molecule, scope: String, path: [String]) throws -> Query {
        try Task.checkCancellation()
        guard input.atoms.count <= maximumAtoms, path.count <= 16,
            input.atoms.allSatisfy({ $0.rGroupMembership == nil })
        else {
            throw CDKRGroupEnumerationError.unsupportedFeature(
                "Nested definition depth, graph size or embedded definitions are unsupported.")
        }
        let resolution = try CDKMarkushVariableLabels.resolve(input)
        var lists: [Int: CDKRGroupList] = [:]
        var evidence: [Evidence.Variable] = []
        for (label, number) in resolution.numberBySourceLabel.sorted(by: { $0.key < $1.key }) {
            guard let binding = binding(label, in: scope), !path.contains(binding.path) else {
                throw CDKRGroupEnumerationError.invalidQuery("Missing or cyclic definition for \(label).")
            }
            let definition = binding.variable
            var groups: [CDKRGroup] = []
            var choices: [[Evidence.Choice]] = []
            for (index, text) in definition.alternatives.enumerated()
            where !definition.excludedAlternativeIndices.contains(index) {
                let childScope = "\(binding.path)/alternatives/\(index)"
                let fragment = try parse(text)
                let nested = try query(
                    fragment, scope: scopes[childScope] == nil ? binding.scope : childScope,
                    path: path + [binding.path])
                let choice = Evidence.Choice(
                    definitionPath: binding.path, variableLabel: label,
                    sourceText: definition.sourceText, sourceAlternativeIndex: index, chemicalClass: nil,
                    fragmentSMILES: text, siteAtomID: nil, omitted: false)
                if nested.query.rGroupDefinitions.isEmpty {
                    try retain(fragment, attachment: fragment.atoms.first?.id, groups: &groups)
                    choices.append([choice])
                } else {
                    try nested.query.enumerateConfigurationsWithAssignments(
                        maximumConfigurations: maximumAlternatives, maximumAtoms: maximumAtoms
                    ) { member in
                        try retain(member.molecule, attachment: fragment.atoms.first?.id, groups: &groups)
                        choices.append([choice] + (try nested.evidence.choices(for: member.assignments)))
                        return true
                    }
                }
            }
            for chemicalClass in definition.chemicalClasses ?? [] {
                try validateClassAttachments(resolution.molecule, number: number)
                for text in try chemicalClass.fragments() {
                    let fragment = try parse(text)
                    try retain(fragment, attachment: fragment.atoms.first?.id, groups: &groups)
                    choices.append([
                        .init(
                            definitionPath: binding.path, variableLabel: label,
                            sourceText: definition.sourceText, sourceAlternativeIndex: nil,
                            chemicalClass: chemicalClass,
                            fragmentSMILES: text, siteAtomID: nil, omitted: false)
                    ])
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
            evidence.append(
                .init(
                    rGroupNumber: number, definitionPath: binding.path, label: label,
                    sourceText: definition.sourceText, alternatives: choices))
        }
        return Query(
            query: .init(rootStructure: resolution.molecule, rGroupDefinitions: lists),
            evidence: .init(sourceIdentifier: "", variables: evidence))
    }

    private mutating func retain(_ molecule: Molecule, attachment: Int?, groups: inout [CDKRGroup]) throws {
        guard let attachment, molecule.atom(id: attachment) != nil, produced < maximumAlternatives else {
            throw CDKRGroupEnumerationError.unsupportedFeature(
                "Alternative expansion exceeds its limit or replaces the attachment atom.")
        }
        guard molecule.atoms.count <= maximumAtoms - retainedAtoms else {
            throw CDKRGroupEnumerationError.atomLimitExceeded(maximumAtoms)
        }
        produced += 1; retainedAtoms += molecule.atoms.count
        groups.append(.init(group: molecule, firstAttachmentPointAtomID: attachment))
    }

    private func validateClassAttachments(_ molecule: Molecule, number: Int) throws {
        for atom in molecule.atoms where atom.rGroupLabel == number {
            let bonds = molecule.bonds.filter { $0.a1 == atom.id || $0.a2 == atom.id }
            guard bonds.count == 1, bonds[0].order == .single else {
                throw CDKRGroupEnumerationError.unsupportedFeature(
                    "Alkyl and alkoxy classes require one single-bond attachment at each site.")
            }
        }
    }
}
