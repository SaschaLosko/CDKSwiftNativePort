import Foundation

// Inspired by CDK 2.13 VarAttachFilter, Copyright 2025 NextMove Software Ltd.
// SPDX-License-Identifier: LGPL-2.1-or-later
// All compatible choices are searched, so an early overlapping bond choice cannot hide a match.
struct QueryVariableAttachments {
    struct VariableBond {
        let bond: Bond
        let beginnings: [Int]
        let endings: [Int]
        let expression: CDKBondExpression
    }
    let placeholders: Set<Int>
    let variableBondIDs: Set<Int>
    let variables: [VariableBond]
    let fixedBonds: [Bond]

    init(_ pattern: CDKQueryPattern) throws {
        let molecule = pattern.molecule
        var endpoints: [Int: [Int]] = [:]
        for group in molecule.sgroups where group.kind == .extMulticenter {
            guard group.keyword == "m", let placeholder = group.atomIDs.first,
                let atom = molecule.atom(id: placeholder), atom.element == "*", atom.aliasLabel == nil,
                molecule.bonds(forAtom: placeholder).count == 1,
                group.atomIDs.count > 1, Set(group.atomIDs).count == group.atomIDs.count,
                group.atomIDs.allSatisfy({ molecule.atom(id: $0) != nil }), endpoints[placeholder] == nil,
                pattern.atoms[placeholder] == .any
            else {
                throw CDKQueryError.unsupported(
                    "Variable attachments require an unconstrained explicit dummy endpoint.")
            }
            endpoints[placeholder] = Array(group.atomIDs.dropFirst())
        }
        let placeholders = Set(endpoints.keys)
        self.placeholders = placeholders
        guard endpoints.values.allSatisfy({ Set($0).isDisjoint(with: placeholders) }) else {
            throw CDKQueryError.invalid("A variable endpoint refers to another placeholder.")
        }
        for (index, targets) in molecule.cxState?.positionalVariations ?? [:] {
            guard molecule.atoms.indices.contains(index), targets.allSatisfy(molecule.atoms.indices.contains),
                endpoints[molecule.atoms[index].id] == targets.map({ molecule.atoms[$0].id })
            else {
                throw CDKQueryError.invalid("Inconsistent variable attachment metadata.")
            }
        }
        variables = molecule.bonds.filter { placeholders.contains($0.a1) || placeholders.contains($0.a2) }.map {
            .init(
                bond: $0, beginnings: endpoints[$0.a1] ?? [$0.a1], endings: endpoints[$0.a2] ?? [$0.a2],
                expression: pattern.bonds[$0.id]!)
        }
        let variableBondIDs = Set(variables.map { $0.bond.id })
        self.variableBondIDs = variableBondIDs
        fixedBonds = molecule.bonds.filter { !variableBondIDs.contains($0.id) }
    }

    func resolve(
        _ mapping: [Int: Int], target: QueryGraphFacts, budget: MarkushMatchBudget,
        accept: ([Int: Int]) throws -> Bool
    ) throws -> [Int: Int]? {
        // Reserve fixed edges as well as chosen variable edges. Two query bonds cannot
        // consume the same target edge, even when their source endpoint sets overlap.
        var occupied: Set<Int> = []
        for bond in fixedBonds {
            if let a = mapping[bond.a1], let b = mapping[bond.a2], let edge = target.edges[a]?[b] {
                occupied.insert(edge.id)
            }
        }
        func visit(_ index: Int, _ current: [Int: Int]) throws -> [Int: Int]? {
            try budget.check()
            if index == variables.count { return try accept(current) ? current : nil }
            let variable = variables[index]
            for begin in variable.beginnings {
                for end in variable.endings {
                    try budget.check()
                    guard let a = current[begin], let b = current[end], a != b,
                        let bond = target.edges[a]?[b], !occupied.contains(bond.id),
                        try target.bond(variable.expression, bond, budget: budget)
                    else { continue }
                    var resolved = current
                    if placeholders.contains(variable.bond.a1) { resolved[variable.bond.a1] = a }
                    if placeholders.contains(variable.bond.a2) { resolved[variable.bond.a2] = b }
                    occupied.insert(bond.id)
                    if let found = try visit(index + 1, resolved) { return found }
                    occupied.remove(bond.id)
                }
            }
            return nil
        }
        return try visit(0, mapping)
    }
}
