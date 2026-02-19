import Foundation

/// CDK-style CML writer for atomArray/bondArray output.
public enum CDKCMLWriter {
    public static func write(_ molecules: [Molecule]) throws -> String {
        guard !molecules.isEmpty else { throw ChemError.emptyInput }
        for molecule in molecules where molecule.atoms.isEmpty {
            throw ChemError.parseFailed("Cannot write CML for an empty molecule.")
        }

        var lines: [String] = []
        lines.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
        lines.append("<cml xmlns=\"http://www.xml-cml.org/schema\">")

        for (idx, molecule) in molecules.enumerated() {
            let name = normalizedName(molecule.name, fallback: "Molecule \(idx + 1)")
            lines.append("  <molecule id=\"m\(idx + 1)\" title=\"\(xmlEsc(name))\">")
            lines.append("    <atomArray>")

            let atoms = molecule.atoms.sorted { $0.id < $1.id }
            for atom in atoms {
                var attrs = [
                    "id=\"a\(atom.id)\"",
                    "elementType=\"\(xmlEsc(CDKDescriptorSupport.canonicalElementSymbol(atom.element)))\"",
                    "x2=\"\(fmt(atom.position.x))\"",
                    "y2=\"\(fmt(atom.position.y))\""
                ]
                if atom.charge != 0 {
                    attrs.append("formalCharge=\"\(atom.charge)\"")
                }
                if let isotope = atom.isotopeMassNumber {
                    attrs.append("isotopeNumber=\"\(isotope)\"")
                }
                if let h = atom.explicitHydrogenCount {
                    attrs.append("hydrogenCount=\"\(h)\"")
                }
                if atom.aromatic {
                    attrs.append("aromatic=\"true\"")
                }
                lines.append("      <atom \(attrs.joined(separator: " ")) />")
            }
            lines.append("    </atomArray>")

            let bonds = molecule.bonds.sorted { lhs, rhs in
                if lhs.a1 != rhs.a1 { return lhs.a1 < rhs.a1 }
                if lhs.a2 != rhs.a2 { return lhs.a2 < rhs.a2 }
                return lhs.id < rhs.id
            }
            lines.append("    <bondArray>")
            for (idxBond, bond) in bonds.enumerated() {
                let order = cmlOrder(for: bond.order)
                lines.append("      <bond id=\"b\(idxBond + 1)\" atomRefs2=\"a\(bond.a1) a\(bond.a2)\" order=\"\(order)\" />")
            }
            lines.append("    </bondArray>")
            lines.append("  </molecule>")
        }

        lines.append("</cml>")
        return lines.joined(separator: "\n") + "\n"
    }

    private static func cmlOrder(for order: BondOrder) -> String {
        switch order {
        case .single:
            return "1"
        case .double:
            return "2"
        case .triple:
            return "3"
        case .aromatic:
            return "A"
        }
    }

    private static func normalizedName(_ raw: String, fallback: String) -> String {
        let cleaned = raw.replacingOccurrences(of: "\r", with: " ")
            .replacingOccurrences(of: "\n", with: " ")
            .trimmingCharacters(in: .whitespacesAndNewlines)
        return cleaned.isEmpty ? fallback : cleaned
    }

    private static func fmt(_ value: CGFloat) -> String {
        String(format: "%.5f", Double(value))
    }

    private static func xmlEsc(_ s: String) -> String {
        var result = s
        result = result.replacingOccurrences(of: "&", with: "&amp;")
        result = result.replacingOccurrences(of: "\"", with: "&quot;")
        result = result.replacingOccurrences(of: "<", with: "&lt;")
        result = result.replacingOccurrences(of: ">", with: "&gt;")
        return result
    }
}
