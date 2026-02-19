import Foundation

/// CDK-style Tripos MOL2 writer.
public enum CDKMol2Writer {
    public static func write(_ molecules: [Molecule]) throws -> String {
        guard !molecules.isEmpty else { throw ChemError.emptyInput }
        let blocks = try molecules.enumerated().map { (idx, molecule) in
            try writeBlock(molecule, index: idx + 1)
        }
        return blocks.joined(separator: "\n")
    }

    private static func writeBlock(_ molecule: Molecule, index: Int) throws -> String {
        guard !molecule.atoms.isEmpty else {
            throw ChemError.parseFailed("Cannot write MOL2 for an empty molecule.")
        }

        let atoms = molecule.atoms.sorted { $0.id < $1.id }
        let atomSerialByID = Dictionary(uniqueKeysWithValues: atoms.enumerated().map { (idx, atom) in
            (atom.id, idx + 1)
        })

        let bonds = molecule.bonds.sorted { lhs, rhs in
            if lhs.a1 != rhs.a1 { return lhs.a1 < rhs.a1 }
            if lhs.a2 != rhs.a2 { return lhs.a2 < rhs.a2 }
            return lhs.id < rhs.id
        }

        let name = normalizedName(molecule.name, fallback: "MOL2 Molecule \(index)")

        var lines: [String] = []
        lines.append("@<TRIPOS>MOLECULE")
        lines.append(name)
        lines.append(String(format: "%5d %5d 0 0 0", atoms.count, bonds.count))
        lines.append("SMALL")
        lines.append("USER_CHARGES")
        lines.append("")
        lines.append("@<TRIPOS>ATOM")

        for atom in atoms {
            let serial = atomSerialByID[atom.id] ?? atom.id
            let atomName = "\(CDKDescriptorSupport.canonicalElementSymbol(atom.element))\(serial)"
            let type = atomType(for: atom, in: molecule)
            lines.append(String(format: "%7d %-8@ %10.4f %10.4f %10.4f %-8@ %3d MOL %8.4f",
                                serial,
                                atomName as NSString,
                                Double(atom.position.x),
                                Double(atom.position.y),
                                0.0,
                                type as NSString,
                                1,
                                Double(atom.charge)))
        }

        lines.append("@<TRIPOS>BOND")
        for (idx, bond) in bonds.enumerated() {
            guard let a1 = atomSerialByID[bond.a1], let a2 = atomSerialByID[bond.a2] else {
                throw ChemError.parseFailed("Cannot write MOL2: bond references an unknown atom ID.")
            }
            lines.append(String(format: "%6d %4d %4d %@",
                                idx + 1,
                                a1,
                                a2,
                                bondType(for: bond) as NSString))
        }

        return lines.joined(separator: "\n") + "\n"
    }

    private static func atomType(for atom: Atom, in molecule: Molecule) -> String {
        let symbol = CDKDescriptorSupport.canonicalElementSymbol(atom.element)
        let lower = symbol.lowercased()
        if atom.aromatic || molecule.bonds(forAtom: atom.id).contains(where: { $0.order == .aromatic }) {
            return "\(lower).ar"
        }
        if symbol.uppercased() == "N" {
            let hasDouble = molecule.bonds(forAtom: atom.id).contains(where: { $0.order == .double })
            return hasDouble ? "N.2" : "N.3"
        }
        if symbol.uppercased() == "O" {
            let hasDouble = molecule.bonds(forAtom: atom.id).contains(where: { $0.order == .double })
            return hasDouble ? "O.2" : "O.3"
        }
        if symbol.uppercased() == "C" {
            let hasDouble = molecule.bonds(forAtom: atom.id).contains(where: { $0.order == .double })
            return hasDouble ? "C.2" : "C.3"
        }
        return "\(symbol).3"
    }

    private static func bondType(for bond: Bond) -> String {
        switch bond.order {
        case .single:
            return "1"
        case .double:
            return "2"
        case .triple:
            return "3"
        case .aromatic:
            return "ar"
        }
    }

    private static func normalizedName(_ raw: String, fallback: String) -> String {
        let cleaned = raw.replacingOccurrences(of: "\r", with: " ")
            .replacingOccurrences(of: "\n", with: " ")
            .trimmingCharacters(in: .whitespacesAndNewlines)
        return cleaned.isEmpty ? fallback : cleaned
    }
}
