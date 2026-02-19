import Foundation

/// Swift counterpart for CDK-style MDL V2000 writing used by native InChI generation.
public enum CDKMDLV2000Writer {

    public static func write(_ molecule: Molecule) throws -> String {
        guard !molecule.atoms.isEmpty else {
            throw ChemError.parseFailed("Cannot write Molfile for empty molecule.")
        }

        if molecule.atoms.contains(where: { $0.queryType != nil || ($0.atomList != nil && !($0.atomList?.isEmpty ?? true)) }) {
            throw ChemError.unsupported("MDL writer does not support query atoms for InChI generation.")
        }
        if molecule.bonds.contains(where: { $0.queryType != nil }) {
            throw ChemError.unsupported("MDL writer does not support query bonds for InChI generation.")
        }

        let atoms = molecule.atoms.sorted { $0.id < $1.id }
        let atomIndexByID = Dictionary(uniqueKeysWithValues: atoms.enumerated().map { (idx, atom) in
            (atom.id, idx + 1)
        })

        let bonds = molecule.bonds.sorted { lhs, rhs in
            if lhs.a1 != rhs.a1 { return lhs.a1 < rhs.a1 }
            if lhs.a2 != rhs.a2 { return lhs.a2 < rhs.a2 }
            return lhs.id < rhs.id
        }

        for bond in bonds {
            guard atomIndexByID[bond.a1] != nil, atomIndexByID[bond.a2] != nil else {
                throw ChemError.parseFailed("Bond references unknown atom id while writing Molfile.")
            }
        }

        var lines: [String] = []
        let title = molecule.name.trimmingCharacters(in: .whitespacesAndNewlines)
        lines.append(title.isEmpty ? "Molecule" : title)
        lines.append("CDKSwiftNativePort")
        lines.append("")

        lines.append(String(format: "%3d%3d  0  0  0  0            999 V2000", atoms.count, bonds.count))

        for atom in atoms {
            let symbol = normalizedElementSymbol(atom.element)
            let x = atom.position.x
            let y = atom.position.y
            let z: Double = 0
            // Write charge via M  CHG for broader compatibility.
            lines.append(String(format: "%10.4f%10.4f%10.4f %-3@  0  0  0  0  0  0  0  0  0  0  0  0", x, y, z, symbol))
        }

        for bond in bonds {
            let a1 = atomIndexByID[bond.a1] ?? 0
            let a2 = atomIndexByID[bond.a2] ?? 0
            let order = molfileBondOrderCode(bond.order)
            let stereo = molfileStereoCode(bond.stereo)
            lines.append(String(format: "%3d%3d%3d%3d  0  0  0", a1, a2, order, stereo))
        }

        appendAtomPropertyLines(prefix: "M  CHG", valuesByAtomID: Dictionary(uniqueKeysWithValues: atoms.compactMap { atom in
            atom.charge == 0 ? nil : (atom.id, atom.charge)
        }), atomIndexByID: atomIndexByID, to: &lines)

        appendAtomPropertyLines(prefix: "M  ISO", valuesByAtomID: Dictionary(uniqueKeysWithValues: atoms.compactMap { atom in
            guard let mass = atom.isotopeMassNumber else { return nil }
            return (atom.id, mass)
        }), atomIndexByID: atomIndexByID, to: &lines)

        appendAtomPropertyLines(prefix: "M  RAD", valuesByAtomID: Dictionary(uniqueKeysWithValues: atoms.compactMap { atom in
            guard let radical = atom.radical, radical > 0 else { return nil }
            return (atom.id, radical)
        }), atomIndexByID: atomIndexByID, to: &lines)

        lines.append("M  END")
        return lines.joined(separator: "\n")
    }

    private static func appendAtomPropertyLines(prefix: String,
                                                valuesByAtomID: [Int: Int],
                                                atomIndexByID: [Int: Int],
                                                to lines: inout [String]) {
        let sorted = valuesByAtomID.keys.sorted().compactMap { atomID -> (Int, Int)? in
            guard let index = atomIndexByID[atomID], let value = valuesByAtomID[atomID] else { return nil }
            return (index, value)
        }

        guard !sorted.isEmpty else { return }

        let chunkSize = 8
        var start = 0
        while start < sorted.count {
            let end = min(sorted.count, start + chunkSize)
            let chunk = sorted[start..<end]
            var line = "\(prefix)\(String(format: "%3d", chunk.count))"
            for (atomIndex, value) in chunk {
                line += String(format: "%4d%4d", atomIndex, value)
            }
            lines.append(line)
            start = end
        }
    }

    private static func normalizedElementSymbol(_ raw: String) -> String {
        let trimmed = raw.trimmingCharacters(in: .whitespacesAndNewlines)
        if trimmed == "*" { return "C" }
        guard let first = trimmed.first else { return "C" }
        let head = String(first).uppercased()
        let tail = String(trimmed.dropFirst()).lowercased()
        return (head + tail).prefix(3).padding(toLength: 3, withPad: " ", startingAt: 0)
    }

    private static func molfileBondOrderCode(_ order: BondOrder) -> Int {
        switch order {
        case .single:
            return 1
        case .double:
            return 2
        case .triple:
            return 3
        case .aromatic:
            return 4
        }
    }

    private static func molfileStereoCode(_ stereo: BondStereo) -> Int {
        switch stereo {
        case .up, .upReversed:
            return 1
        case .either:
            return 4
        case .down, .downReversed:
            return 6
        case .none:
            return 0
        }
    }
}
