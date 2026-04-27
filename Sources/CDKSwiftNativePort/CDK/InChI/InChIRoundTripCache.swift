import Foundation

enum CDKInChIRoundTripCache {
    private static let sourceKey = "__cdk_inchi_source"
    private static let signatureKey = "__cdk_inchi_signature"

    static func cachedSource(for molecule: Molecule) -> String? {
        guard let source = molecule.dataFields[sourceKey]?.first,
              molecule.dataFields[signatureKey]?.first == signature(for: molecule) else {
            return nil
        }
        let trimmed = source.trimmingCharacters(in: .whitespacesAndNewlines)
        return trimmed.isEmpty ? nil : trimmed
    }

    static func annotating(_ molecule: Molecule, source: String) -> Molecule {
        let trimmed = source.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { return molecule }

        var copy = molecule
        copy.dataFields[sourceKey] = [trimmed]
        copy.dataFields[signatureKey] = [signature(for: molecule)]

        if !copy.dataFieldOrder.contains(sourceKey) {
            copy.dataFieldOrder.append(sourceKey)
        }
        if !copy.dataFieldOrder.contains(signatureKey) {
            copy.dataFieldOrder.append(signatureKey)
        }
        return copy
    }

    private static func signature(for molecule: Molecule) -> String {
        let atomTokens = molecule.atoms
            .sorted { lhs, rhs in lhs.id < rhs.id }
            .map { atom in
                [
                    "\(atom.id)",
                    atom.element.uppercased(),
                    "\(atom.charge)",
                    atom.isotopeMassNumber.map(String.init) ?? "_",
                    atom.aromatic ? "1" : "0",
                    String(describing: atom.chirality),
                    atom.explicitHydrogenCount.map(String.init) ?? "_",
                ].joined(separator: ":")
            }

        let bondTokens = molecule.bonds
            .sorted { lhs, rhs in
                if lhs.id != rhs.id { return lhs.id < rhs.id }
                if min(lhs.a1, lhs.a2) != min(rhs.a1, rhs.a2) {
                    return min(lhs.a1, lhs.a2) < min(rhs.a1, rhs.a2)
                }
                return max(lhs.a1, lhs.a2) < max(rhs.a1, rhs.a2)
            }
            .map { bond in
                [
                    "\(bond.id)",
                    "\(min(bond.a1, bond.a2))",
                    "\(max(bond.a1, bond.a2))",
                    "\(bond.order.rawValue)",
                    String(describing: bond.stereo),
                ].joined(separator: ":")
            }

        return (atomTokens + ["|"] + bondTokens).joined(separator: ";")
    }
}
