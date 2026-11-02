import Foundation

/// Explicit, reversible label resolution for executable finite definitions.
/// Numeric R labels retain their number. Other labels receive unused numbers in sorted order.
/// Only declared pseudo atoms are variables; ordinary elements such as B and Y remain elements.
public enum CDKMarkushVariableLabels {
    public struct Resolution: Sendable {
        public let molecule: Molecule
        public let numberBySourceLabel: [String: Int]
    }

    public static func resolve(_ source: Molecule) throws -> Resolution {
        func number(_ label: String) -> Int? {
            guard label == "R" || label.range(of: #"^R[0-9]+$"#, options: .regularExpression) != nil else {
                return nil
            }
            return label == "R" ? 0 : Int(label.dropFirst())
        }
        func label(_ atom: Atom) -> String? {
            if let n = atom.rGroupLabel { return number(atom.element) == n ? atom.element : "R\(n)" }
            guard atom.isUnresolvedPseudoAtom else { return nil }
            return atom.aliasLabel ?? atom.element
        }
        guard source.atoms.allSatisfy({ ($0.rGroupLabel ?? 0) >= 0 }),
            source.rGroupLogicDefinitions.keys.allSatisfy({ $0 >= 0 })
        else { throw CDKRGroupEnumerationError.invalidQuery("Variable numbers must be nonnegative.") }
        let positionDummies = Set(source.sgroups.filter { $0.kind == .extMulticenter }.compactMap { $0.atomIDs.first })
        var labels = Set(source.atoms.compactMap(\.rGroupMembership))
        labels.formUnion(
            source.atoms.filter { $0.rGroupMembership == nil && !positionDummies.contains($0.id) }.compactMap(label))
        var used = Set(source.atoms.compactMap(\.rGroupLabel)).union(source.rGroupLogicDefinitions.keys)
        var mapping: [String: Int] = [:]
        for value in labels.sorted() {
            guard !value.isEmpty, value != "*", value == value.trimmingCharacters(in: .whitespacesAndNewlines) else {
                throw CDKRGroupEnumerationError.invalidQuery("An unnamed or malformed variable cannot be bound.")
            }
            if let n = number(value) {
                guard !mapping.values.contains(n) else {
                    throw CDKRGroupEnumerationError.invalidQuery("Ambiguous numeric aliases for R\(n).")
                }
                mapping[value] = n
                used.insert(n)
            } else if value.range(of: #"^R[0-9]+$"#, options: .regularExpression) != nil {
                throw CDKRGroupEnumerationError.invalidQuery("Variable number is out of range.")
            }
        }
        var next = 1
        for value in labels.sorted() where mapping[value] == nil {
            while used.contains(next) { next += 1 }
            mapping[value] = next
            used.insert(next)
        }
        var result = source
        for index in result.atoms.indices {
            let atom = source.atoms[index]
            if let membership = atom.rGroupMembership {
                result.atoms[index].rGroupMembership = mapping[membership].map { "R\($0)" }
            } else if let original = label(atom), let n = mapping[original] {
                result.atoms[index].rGroupLabel = n
                result.atoms[index].element = "R\(n)"
            }
        }
        return Resolution(molecule: result, numberBySourceLabel: mapping)
    }
}
