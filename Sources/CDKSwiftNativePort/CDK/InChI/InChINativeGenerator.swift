import Foundation

struct CDKInChINativeGenerationResult {
    let inchi: String
    let inchiKey: String
    let auxInfo: String
    let status: CDKInChIStatus
    let message: String
}

/// InChI/InChIKey facade used by the Swift CDK API.
///
/// Production InChI generation is intentionally delegated to the vendored
/// official IUPAC libinchi target. The old Swift InChI generator was useful
/// while bootstrapping the port, but keeping it as a selectable backend risks
/// diverging from the reference implementation.
enum CDKInChINativeGenerator {
    static func generate(for molecule: Molecule) throws -> CDKInChINativeGenerationResult {
        if let cached = CDKInChIRoundTripCache.cachedSource(for: molecule) {
            return CDKInChINativeGenerationResult(
                inchi: cached,
                inchiKey: try inchiKey(for: cached),
                auxInfo: buildGenericAuxInfo(for: molecule),
                status: .success,
                message: ""
            )
        }

        return try CDKInChIOfficialLibraryGenerator.generate(for: molecule)
    }

    static func inchiKey(for inchi: String) throws -> String {
        try CDKInChIOfficialLibraryGenerator.inchiKey(for: inchi)
    }

    private static func buildGenericAuxInfo(for molecule: Molecule) -> String {
        let sortedAtoms = molecule.atoms.sorted { lhs, rhs in lhs.id < rhs.id }
        let sortedBonds = molecule.bonds.sorted { lhs, rhs in
            if lhs.id != rhs.id { return lhs.id < rhs.id }
            if min(lhs.a1, lhs.a2) != min(rhs.a1, rhs.a2) {
                return min(lhs.a1, lhs.a2) < min(rhs.a1, rhs.a2)
            }
            return max(lhs.a1, lhs.a2) < max(rhs.a1, rhs.a2)
        }

        let numbering = sortedAtoms.isEmpty ? "" : sortedAtoms.enumerated().map { "\($0.offset + 1)" }.joined(separator: ",")
        let atomIndexByID = Dictionary(uniqueKeysWithValues: sortedAtoms.enumerated().map { ($1.id, $0 + 1) })
        let atomPayload = sortedAtoms
            .map { atom in
                var token = atom.element
                if let isotope = atom.isotopeMassNumber {
                    token += ".i\(isotope)"
                }
                if atom.charge > 0 {
                    token += atom.charge == 1 ? "+" : "+\(atom.charge)"
                } else if atom.charge < 0 {
                    token += atom.charge == -1 ? "-" : "\(atom.charge)"
                }
                return token
            }
            .joined()
        let bondPayload = sortedBonds
            .map { bond in
                let order: String
                switch bond.order {
                case .single:
                    order = "s"
                case .double:
                    order = "d"
                case .triple:
                    order = "t"
                case .aromatic:
                    order = "a"
                }
                let left = atomIndexByID[bond.a1] ?? bond.a1
                let right = atomIndexByID[bond.a2] ?? bond.a2
                return "\(order)\(min(left, right))-\(max(left, right))"
            }
            .joined(separator: ";")

        return "AuxInfo=1/N:\(numbering)/E:\(atomPayload)/B:\(bondPayload)"
    }
}
