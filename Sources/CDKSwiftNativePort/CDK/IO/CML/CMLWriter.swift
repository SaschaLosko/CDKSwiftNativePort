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
            let has3D = molecule.atoms.contains { $0.zPosition != nil }
            let aromaticBondIDs = molecule.aromaticDisplayBondIDs()

            lines.append("  <molecule id=\"m\(idx + 1)\" title=\"\(xmlEsc(name))\">")
            lines.append("    <atomArray>")

            let atoms = molecule.atoms.sorted { $0.id < $1.id }
            for atom in atoms {
                let pseudoLabel = pseudoLabel(for: atom)
                let elementType = pseudoLabel == nil
                    ? normalizedElementType(for: atom)
                    : "Du"

                var attrs = [
                    "id=\"a\(atom.id)\"",
                    "elementType=\"\(xmlEsc(elementType))\"",
                    "x2=\"\(fmt(Double(atom.position.x)))\"",
                    "y2=\"\(fmt(Double(atom.position.y)))\""
                ]
                if has3D {
                    attrs.append("x3=\"\(fmt(Double(atom.position.x)))\"")
                    attrs.append("y3=\"\(fmt(Double(atom.position.y)))\"")
                    attrs.append("z3=\"\(fmt(atom.zPosition ?? 0))\"")
                }
                if let pseudoLabel {
                    attrs.append("title=\"\(xmlEsc(pseudoLabel))\"")
                }
                if atom.charge != 0 {
                    attrs.append("formalCharge=\"\(atom.charge)\"")
                }
                if pseudoLabel == nil, let isotope = atom.isotopeMassNumber {
                    attrs.append("isotopeNumber=\"\(isotope)\"")
                }

                let totalHydrogenCount = hydrogenCountForCML(atomID: atom.id, in: molecule)
                if totalHydrogenCount > 0 || atom.explicitHydrogenCount != nil {
                    attrs.append("hydrogenCount=\"\(totalHydrogenCount)\"")
                }
                if atom.aromatic {
                    attrs.append("aromatic=\"true\"")
                }

                var childLines: [String] = []
                if atom.aromatic {
                    childLines.append("        <scalar dictRef=\"cdk:aromaticAtom\" />")
                }

                if childLines.isEmpty {
                    lines.append("      <atom \(attrs.joined(separator: " ")) />")
                } else {
                    lines.append("      <atom \(attrs.joined(separator: " "))>")
                    lines.append(contentsOf: childLines)
                    lines.append("      </atom>")
                }
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
                var childLines: [String] = []
                if aromaticBondIDs.contains(bond.id) {
                    childLines.append("        <bondType dictRef=\"cdk:aromaticBond\" />")
                }
                if let stereo = cmlBondStereo(for: bond.stereo) {
                    childLines.append("        <bondStereo dictRef=\"\(stereo.dictRef)\">\(stereo.content)</bondStereo>")
                }

                let attrs = [
                    "id=\"b\(idxBond + 1)\"",
                    "atomRefs2=\"a\(bond.a1) a\(bond.a2)\"",
                    "order=\"\(order)\""
                ]

                if childLines.isEmpty {
                    lines.append("      <bond \(attrs.joined(separator: " ")) />")
                } else {
                    lines.append("      <bond \(attrs.joined(separator: " "))>")
                    lines.append(contentsOf: childLines)
                    lines.append("      </bond>")
                }
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
            return "S"
        case .double:
            return "D"
        case .triple:
            return "T"
        case .aromatic:
            return "A"
        }
    }

    private static func cmlBondStereo(for stereo: BondStereo) -> (dictRef: String, content: String)? {
        switch stereo {
        case .up, .upReversed:
            return ("cml:W", "W")
        case .down, .downReversed:
            return ("cml:H", "H")
        default:
            return nil
        }
    }

    private static func hydrogenCountForCML(atomID: Int, in molecule: Molecule) -> Int {
        CDKDescriptorSupport.totalHydrogenCount(on: atomID, in: molecule)
    }

    private static func pseudoLabel(for atom: Atom) -> String? {
        if let alias = atom.aliasLabel?.trimmingCharacters(in: .whitespacesAndNewlines),
           !alias.isEmpty {
            return alias
        }

        let trimmed = atom.element.trimmingCharacters(in: .whitespacesAndNewlines)
        if trimmed.isEmpty || trimmed == "*" {
            return nil
        }
        let upper = trimmed.uppercased()
        if upper == "R" || upper == "R#" || upper == "DU" || upper == "DUMMY" {
            return "R"
        }
        return isPlausibleElementSymbol(trimmed) ? nil : trimmed
    }

    private static func normalizedElementType(for atom: Atom) -> String {
        let canonical = CDKDescriptorSupport.canonicalElementSymbol(atom.element)
        if canonical.isEmpty {
            return "C"
        }
        return canonical
    }

    private static func isPlausibleElementSymbol(_ raw: String) -> Bool {
        let canonical = CDKDescriptorSupport.canonicalElementSymbol(raw)
        if canonical.isEmpty {
            return false
        }
        let upper = canonical.uppercased()
        if upper == "D" || upper == "T" {
            return true
        }
        return CDKDescriptorSupport.averageAtomicMass(forElementSymbol: canonical) > 0
    }

    private static func normalizedName(_ raw: String, fallback: String) -> String {
        let cleaned = raw.replacingOccurrences(of: "\r", with: " ")
            .replacingOccurrences(of: "\n", with: " ")
            .trimmingCharacters(in: .whitespacesAndNewlines)
        return cleaned.isEmpty ? fallback : cleaned
    }

    private static func fmt(_ value: Double) -> String {
        String(format: "%.5f", value)
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
