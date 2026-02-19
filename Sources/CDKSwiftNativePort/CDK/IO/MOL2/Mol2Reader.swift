import Foundation
import CoreGraphics

/// CDK-style TRIPOS MOL2 reader.
public enum CDKMol2Reader {
    public static func read(text: String) throws -> [Molecule] {
        let normalized = text
            .replacingOccurrences(of: "\r\n", with: "\n")
            .replacingOccurrences(of: "\r", with: "\n")
        let lines = normalized.components(separatedBy: "\n")

        let marker = "@<TRIPOS>MOLECULE"
        var starts: [Int] = []
        for (idx, line) in lines.enumerated() where line.trimmingCharacters(in: .whitespacesAndNewlines).uppercased() == marker {
            starts.append(idx)
        }

        guard !starts.isEmpty else {
            throw ChemError.parseFailed("MOL2 file is missing @<TRIPOS>MOLECULE section.")
        }

        var molecules: [Molecule] = []
        for (blockIndex, start) in starts.enumerated() {
            let end = (blockIndex + 1 < starts.count) ? starts[blockIndex + 1] : lines.count
            let blockLines = Array(lines[start..<end])
            if let molecule = try parseMoleculeBlock(blockLines, index: blockIndex + 1) {
                molecules.append(molecule)
            }
        }

        guard !molecules.isEmpty else {
            throw ChemError.parseFailed("No MOL2 molecule blocks could be parsed.")
        }
        return molecules
    }

    private static func parseMoleculeBlock(_ lines: [String], index: Int) throws -> Molecule? {
        guard lines.count >= 2 else { return nil }

        var section = ""
        var moleculeName = lines[1].trimmingCharacters(in: .whitespacesAndNewlines)
        if moleculeName.isEmpty {
            moleculeName = "MOL2 Molecule \(index)"
        }

        struct PendingAtom {
            let sourceID: Int
            let atom: Atom
        }

        var pendingAtoms: [PendingAtom] = []
        var sourceToAtomID: [Int: Int] = [:]
        var bonds: [Bond] = []

        var nextAtomID = 1
        var nextBondID = 1

        for line in lines.dropFirst(2) {
            let trimmed = line.trimmingCharacters(in: .whitespacesAndNewlines)
            if trimmed.isEmpty { continue }

            if trimmed.uppercased().hasPrefix("@<TRIPOS>") {
                section = trimmed.uppercased()
                continue
            }

            if section == "@<TRIPOS>ATOM" {
                let parts = trimmed.split(whereSeparator: \.isWhitespace).map(String.init)
                guard parts.count >= 6,
                      let sourceID = Int(parts[0]),
                      let x = Double(parts[2]),
                      let y = Double(parts[3]) else {
                    continue
                }

                let atomName = parts[1]
                let atomType = parts[5]
                let symbol = elementSymbol(fromType: atomType, fallbackName: atomName)
                let aromatic = atomType.lowercased().contains(".ar")

                let charge: Int
                if parts.count >= 9, let parsedCharge = Double(parts[8]) {
                    charge = Int(parsedCharge.rounded())
                } else {
                    charge = 0
                }

                let atom = Atom(id: nextAtomID,
                                element: symbol,
                                position: CGPoint(x: x, y: y),
                                charge: charge,
                                aromatic: aromatic)
                pendingAtoms.append(PendingAtom(sourceID: sourceID, atom: atom))
                sourceToAtomID[sourceID] = nextAtomID
                nextAtomID += 1
                continue
            }

            if section == "@<TRIPOS>BOND" {
                let parts = trimmed.split(whereSeparator: \.isWhitespace).map(String.init)
                guard parts.count >= 4,
                      let sourceA = Int(parts[1]),
                      let sourceB = Int(parts[2]),
                      let a1 = sourceToAtomID[sourceA],
                      let a2 = sourceToAtomID[sourceB] else {
                    continue
                }

                let typeToken = parts[3].lowercased()
                let order: BondOrder
                switch typeToken {
                case "2", "d", "dbl", "double":
                    order = .double
                case "3", "t", "tpl", "triple":
                    order = .triple
                case "ar", "aromatic":
                    order = .aromatic
                default:
                    order = .single
                }

                bonds.append(Bond(id: nextBondID,
                                  a1: a1,
                                  a2: a2,
                                  order: order,
                                  stereo: .none))
                nextBondID += 1
            }
        }

        guard !pendingAtoms.isEmpty else { return nil }
        var atoms = pendingAtoms.map(\.atom)

        if bonds.contains(where: { $0.order == .aromatic }) {
            var aromaticAtomIDs = Set<Int>()
            for bond in bonds where bond.order == .aromatic {
                aromaticAtomIDs.insert(bond.a1)
                aromaticAtomIDs.insert(bond.a2)
            }
            for idx in atoms.indices where aromaticAtomIDs.contains(atoms[idx].id) {
                atoms[idx].aromatic = true
            }
        }

        if bonds.isEmpty {
            bonds = CDKBondPerception.inferSingleBonds(for: atoms)
        }

        var molecule = Molecule(name: moleculeName, atoms: atoms, bonds: bonds)
        if let box = molecule.boundingBox(), box.width <= 0.0001 && box.height <= 0.0001 {
            molecule = Depiction2DGenerator.generate(for: molecule)
        }

        return molecule
    }

    private static func elementSymbol(fromType atomType: String, fallbackName: String) -> String {
        let raw = atomType.split(separator: ".", maxSplits: 1).first.map(String.init) ?? atomType
        let letters = raw.filter { $0.isLetter }
        if !letters.isEmpty {
            return normalizeCandidateSymbol(letters)
        }
        return normalizeCandidateSymbol(fallbackName.filter { $0.isLetter })
    }

    private static func normalizeCandidateSymbol(_ raw: String) -> String {
        let uppercase = raw.uppercased()
        if uppercase.isEmpty { return "C" }

        let twoLetterElements: Set<String> = [
            "CL", "BR", "SI", "NA", "MG", "AL", "CA", "FE", "ZN", "CU", "SE", "LI", "BE", "NE", "AR"
        ]
        if uppercase.count >= 2 {
            let two = String(uppercase.prefix(2))
            if twoLetterElements.contains(two) {
                return CDKDescriptorSupport.canonicalElementSymbol(two)
            }
        }
        return CDKDescriptorSupport.canonicalElementSymbol(String(uppercase.prefix(1)))
    }
}
