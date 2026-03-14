import Foundation

public enum CDKCMLReactionWriter {
    public static func write(_ reaction: CDKReaction) throws -> String {
        try write([reaction])
    }

    public static func write(_ reactions: [CDKReaction]) throws -> String {
        guard !reactions.isEmpty else { throw ChemError.emptyInput }

        var lines: [String] = []
        lines.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
        lines.append("<cml xmlns=\"http://www.xml-cml.org/schema\">")

        if reactions.count == 1 {
            lines.append(contentsOf: reactionLines(reactions[0], indent: "  ", defaultIndex: 1))
        } else {
            lines.append("  <reactionList>")
            for (index, reaction) in reactions.enumerated() {
                lines.append(contentsOf: reactionLines(reaction, indent: "    ", defaultIndex: index + 1))
            }
            lines.append("  </reactionList>")
        }

        lines.append("</cml>")
        return lines.joined(separator: "\n") + "\n"
    }

    private static func reactionLines(_ reaction: CDKReaction,
                                      indent: String,
                                      defaultIndex: Int) -> [String] {
        let reactionID = sanitizedIdentifier(reaction.id)
            ?? sanitizedIdentifier(reaction.name)
            ?? "reaction\(defaultIndex)"
        let reactionName = normalizedName(reaction.name)

        var attrs = ["id=\"\(xmlEsc(reactionID))\""]
        if let reactionName, reactionName != reactionID {
            attrs.append("title=\"\(xmlEsc(reactionName))\"")
        }

        var lines: [String] = []
        lines.append("\(indent)<reaction \(attrs.joined(separator: " "))>")

        for key in reaction.properties.keys.sorted() {
            guard let value = reaction.properties[key] else { continue }
            lines.append("\(indent)  <scalar dictRef=\"cdk:reactionProperty\" title=\"\(xmlEsc(key))\" dataType=\"xsd:string\">\(xmlEsc(value))</scalar>")
        }

        lines.append(contentsOf: participantListLines(listTag: "reactantList",
                                                      participantTag: "reactant",
                                                      participants: reaction.reactants,
                                                      indent: indent + "  ",
                                                      rolePrefix: "reactant"))
        lines.append(contentsOf: participantListLines(listTag: "productList",
                                                      participantTag: "product",
                                                      participants: reaction.products,
                                                      indent: indent + "  ",
                                                      rolePrefix: "product"))

        if !reaction.agents.isEmpty {
            lines.append(contentsOf: participantListLines(listTag: "substanceList",
                                                          participantTag: "substance",
                                                          participants: reaction.agents,
                                                          indent: indent + "  ",
                                                          rolePrefix: "agent"))
        }

        lines.append("\(indent)</reaction>")
        return lines
    }

    private static func participantListLines(listTag: String,
                                             participantTag: String,
                                             participants: [Molecule],
                                             indent: String,
                                             rolePrefix: String) -> [String] {
        guard !participants.isEmpty else { return [] }

        var lines: [String] = []
        lines.append("\(indent)<\(listTag)>")
        for (index, molecule) in participants.enumerated() {
            lines.append("\(indent)  <\(participantTag)>")
            lines.append(contentsOf: moleculeLines(molecule,
                                                   indent: indent + "    ",
                                                   defaultIdentifier: "\(rolePrefix)\(index + 1)"))
            lines.append("\(indent)  </\(participantTag)>")
        }
        lines.append("\(indent)</\(listTag)>")
        return lines
    }

    private static func moleculeLines(_ molecule: Molecule,
                                      indent: String,
                                      defaultIdentifier: String) -> [String] {
        let has3D = molecule.atoms.contains { $0.zPosition != nil }
        let aromaticBondIDs = molecule.aromaticDisplayBondIDs()
        let moleculeID = sanitizedIdentifier(molecule.externalID)
            ?? sanitizedIdentifier(molecule.name)
            ?? defaultIdentifier
        let moleculeName = normalizedName(molecule.name)
        let formulaValues = molecule.dataFieldValues(named: "Formula")

        var attrs = ["id=\"\(xmlEsc(moleculeID))\""]
        if let moleculeName, moleculeName != moleculeID {
            attrs.append("title=\"\(xmlEsc(moleculeName))\"")
        }

        if molecule.atoms.isEmpty && formulaValues.isEmpty {
            return ["\(indent)<molecule \(attrs.joined(separator: " ")) />"]
        }

        var lines: [String] = []
        lines.append("\(indent)<molecule \(attrs.joined(separator: " "))>")

        if molecule.atoms.isEmpty {
            for formula in formulaValues {
                lines.append("\(indent)  <formula concise=\"\(xmlEsc(formula))\" />")
            }
            lines.append("\(indent)</molecule>")
            return lines
        }

        lines.append("\(indent)  <atomArray>")
        let atoms = molecule.atoms.sorted { $0.id < $1.id }
        for atom in atoms {
            let pseudoLabel = pseudoLabel(for: atom)
            let elementType = pseudoLabel == nil ? normalizedElementType(for: atom) : "Du"

            var atomAttrs = [
                "id=\"\(xmlEsc(atom.externalID ?? "a\(atom.id)"))\"",
                "elementType=\"\(xmlEsc(elementType))\"",
                "x2=\"\(fmt(Double(atom.position.x)))\"",
                "y2=\"\(fmt(Double(atom.position.y)))\""
            ]
            if has3D {
                atomAttrs.append("x3=\"\(fmt(Double(atom.position.x)))\"")
                atomAttrs.append("y3=\"\(fmt(Double(atom.position.y)))\"")
                atomAttrs.append("z3=\"\(fmt(atom.zPosition ?? 0))\"")
            }
            if let pseudoLabel {
                atomAttrs.append("title=\"\(xmlEsc(pseudoLabel))\"")
            }
            if atom.charge != 0 {
                atomAttrs.append("formalCharge=\"\(atom.charge)\"")
            }
            if pseudoLabel == nil, let isotope = atom.isotopeMassNumber {
                atomAttrs.append("isotopeNumber=\"\(isotope)\"")
            }

            let totalHydrogenCount = CDKDescriptorSupport.totalHydrogenCount(on: atom.id, in: molecule)
            if totalHydrogenCount > 0 || atom.explicitHydrogenCount != nil {
                atomAttrs.append("hydrogenCount=\"\(totalHydrogenCount)\"")
            }
            if atom.aromatic {
                atomAttrs.append("aromatic=\"true\"")
            }

            var childLines: [String] = []
            if atom.aromatic {
                childLines.append("\(indent)      <scalar dictRef=\"cdk:aromaticAtom\" />")
            }

            if childLines.isEmpty {
                lines.append("\(indent)    <atom \(atomAttrs.joined(separator: " ")) />")
            } else {
                lines.append("\(indent)    <atom \(atomAttrs.joined(separator: " "))>")
                lines.append(contentsOf: childLines)
                lines.append("\(indent)    </atom>")
            }
        }
        lines.append("\(indent)  </atomArray>")

        let bonds = molecule.bonds.sorted { lhs, rhs in
            if lhs.a1 != rhs.a1 { return lhs.a1 < rhs.a1 }
            if lhs.a2 != rhs.a2 { return lhs.a2 < rhs.a2 }
            return lhs.id < rhs.id
        }
        lines.append("\(indent)  <bondArray>")
        for (index, bond) in bonds.enumerated() {
            let order = cmlOrder(for: bond.order)
            var childLines: [String] = []
            if aromaticBondIDs.contains(bond.id) {
                childLines.append("\(indent)      <bondType dictRef=\"cdk:aromaticBond\" />")
            }
            if let stereo = cmlBondStereo(for: bond.stereo) {
                childLines.append("\(indent)      <bondStereo dictRef=\"\(stereo.dictRef)\">\(stereo.content)</bondStereo>")
            }

            let attrs = [
                "id=\"\(xmlEsc(bond.externalID ?? "b\(index + 1)"))\"",
                "atomRefs2=\"\(xmlEsc(atomReference(for: bond.a1, in: molecule))) \(xmlEsc(atomReference(for: bond.a2, in: molecule)))\"",
                "order=\"\(order)\""
            ]

            if childLines.isEmpty {
                lines.append("\(indent)    <bond \(attrs.joined(separator: " ")) />")
            } else {
                lines.append("\(indent)    <bond \(attrs.joined(separator: " "))>")
                lines.append(contentsOf: childLines)
                lines.append("\(indent)    </bond>")
            }
        }
        lines.append("\(indent)  </bondArray>")
        lines.append("\(indent)</molecule>")
        return lines
    }

    private static func atomReference(for atomID: Int, in molecule: Molecule) -> String {
        molecule.atom(id: atomID)?.externalID ?? "a\(atomID)"
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
        let canonical = CDKDescriptorSupport.canonicalElementSymbol(trimmed)
        if canonical.isEmpty {
            return trimmed
        }
        let canonicalUpper = canonical.uppercased()
        if canonicalUpper == "D" || canonicalUpper == "T" {
            return nil
        }
        return CDKDescriptorSupport.averageAtomicMass(forElementSymbol: canonical) > 0 ? nil : trimmed
    }

    private static func normalizedElementType(for atom: Atom) -> String {
        let canonical = CDKDescriptorSupport.canonicalElementSymbol(atom.element)
        return canonical.isEmpty ? "C" : canonical
    }

    private static func sanitizedIdentifier(_ raw: String?) -> String? {
        guard let raw else { return nil }
        let trimmed = raw.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { return nil }
        let filtered = trimmed.unicodeScalars.map { scalar -> Character in
            if CharacterSet.alphanumerics.contains(scalar) || scalar == "_" || scalar == "-" || scalar == "." || scalar == "+" {
                return Character(scalar)
            }
            return "_"
        }
        let result = String(filtered)
        return result.isEmpty ? nil : result
    }

    private static func normalizedName(_ raw: String?) -> String? {
        guard let raw else { return nil }
        let cleaned = raw.replacingOccurrences(of: "\r", with: " ")
            .replacingOccurrences(of: "\n", with: " ")
            .trimmingCharacters(in: .whitespacesAndNewlines)
        return cleaned.isEmpty ? nil : cleaned
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
