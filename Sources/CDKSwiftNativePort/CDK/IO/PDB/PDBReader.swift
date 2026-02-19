import Foundation
import CoreGraphics

/// CDK-style PDB reader (ATOM/HETATM + CONECT).
public enum CDKPDBReader {
    public static func read(text: String) throws -> [Molecule] {
        let normalized = text
            .replacingOccurrences(of: "\r\n", with: "\n")
            .replacingOccurrences(of: "\r", with: "\n")
        let lines = normalized.components(separatedBy: "\n")

        var atomRecords: [(serial: Int, atom: Atom)] = []
        var serialToAtomID: [Int: Int] = [:]
        var bondPairs = Set<String>()
        var bonds: [Bond] = []

        var nextAtomID = 1
        var titleParts: [String] = []
        var name = ""

        var modelSeen = false
        var inFirstModel = true

        for line in lines {
            let recordName = field(in: line, start: 0, length: 6).trimmingCharacters(in: .whitespacesAndNewlines)

            if recordName == "MODEL" {
                if !modelSeen {
                    modelSeen = true
                    inFirstModel = true
                } else {
                    inFirstModel = false
                }
                continue
            }
            if recordName == "ENDMDL" {
                if modelSeen { break }
                continue
            }
            if modelSeen && !inFirstModel { continue }

            if recordName == "HEADER", name.isEmpty {
                let headerName = field(in: line, start: 10, length: 40).trimmingCharacters(in: .whitespacesAndNewlines)
                if !headerName.isEmpty {
                    name = headerName
                }
                continue
            }

            if recordName == "TITLE" {
                let part = field(in: line, start: 10, length: 70).trimmingCharacters(in: .whitespacesAndNewlines)
                if !part.isEmpty { titleParts.append(part) }
                continue
            }

            if recordName == "ATOM" || recordName == "HETATM" {
                guard let serial = Int(field(in: line, start: 6, length: 5).trimmingCharacters(in: .whitespacesAndNewlines)),
                      let x = Double(field(in: line, start: 30, length: 8).trimmingCharacters(in: .whitespacesAndNewlines)),
                      let y = Double(field(in: line, start: 38, length: 8).trimmingCharacters(in: .whitespacesAndNewlines)) else {
                    continue
                }

                let explicitElement = field(in: line, start: 76, length: 2).trimmingCharacters(in: .whitespacesAndNewlines)
                let atomName = field(in: line, start: 12, length: 4)
                let symbol = normalizeElementSymbol(explicitElement, atomName: atomName)

                let atom = Atom(id: nextAtomID,
                                element: symbol,
                                position: CGPoint(x: x, y: y))
                atomRecords.append((serial: serial, atom: atom))
                serialToAtomID[serial] = nextAtomID
                nextAtomID += 1
                continue
            }

            if recordName == "CONECT" {
                let parts = line.split(whereSeparator: \.isWhitespace).map(String.init)
                guard parts.count >= 3, let sourceSerial = Int(parts[1]), let sourceID = serialToAtomID[sourceSerial] else {
                    continue
                }

                for targetToken in parts.dropFirst(2) {
                    guard let targetSerial = Int(targetToken),
                          let targetID = serialToAtomID[targetSerial],
                          targetID != sourceID else {
                        continue
                    }
                    let key = canonicalPair(sourceID, targetID)
                    if bondPairs.contains(key) { continue }
                    bondPairs.insert(key)
                    bonds.append(Bond(id: bonds.count + 1,
                                      a1: sourceID,
                                      a2: targetID,
                                      order: .single,
                                      stereo: .none))
                }
            }
        }

        var atoms = atomRecords.map(\.atom)
        guard !atoms.isEmpty else {
            throw ChemError.parseFailed("PDB file does not contain ATOM/HETATM coordinates.")
        }

        if bonds.isEmpty {
            bonds = CDKBondPerception.inferSingleBonds(for: atoms)
        }

        if !titleParts.isEmpty {
            name = titleParts.joined(separator: " ")
        }
        if name.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty {
            name = "PDB Molecule"
        }

        var molecule = Molecule(name: name, atoms: atoms, bonds: bonds)
        if let box = molecule.boundingBox(), box.width <= 0.0001 && box.height <= 0.0001 {
            molecule = Depiction2DGenerator.generate(for: molecule)
            atoms = molecule.atoms
        }

        return [Molecule(name: molecule.name, atoms: atoms, bonds: molecule.bonds)]
    }

    private static func normalizeElementSymbol(_ explicit: String, atomName: String) -> String {
        let explicitTrimmed = explicit.trimmingCharacters(in: .whitespacesAndNewlines)
        if !explicitTrimmed.isEmpty {
            return CDKDescriptorSupport.canonicalElementSymbol(explicitTrimmed)
        }

        let lettersOnly = atomName
            .filter { $0.isLetter }
            .uppercased()

        guard !lettersOnly.isEmpty else { return "C" }
        let twoLetterElements: Set<String> = [
            "HE", "LI", "BE", "NE", "NA", "MG", "AL", "SI", "CL", "AR",
            "CA", "SC", "TI", "CR", "MN", "FE", "CO", "NI", "CU", "ZN",
            "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "ZR", "NB",
            "MO", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE",
            "XE", "CS", "BA", "LA", "CE", "PR", "ND", "SM", "EU", "GD",
            "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "RE",
            "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI"
        ]

        if lettersOnly.count >= 2 {
            let pair = String(lettersOnly.prefix(2))
            if twoLetterElements.contains(pair) {
                return CDKDescriptorSupport.canonicalElementSymbol(pair)
            }
        }

        return CDKDescriptorSupport.canonicalElementSymbol(String(lettersOnly.prefix(1)))
    }

    private static func field(in line: String, start: Int, length: Int) -> String {
        guard start < line.count else { return "" }
        let lower = line.index(line.startIndex, offsetBy: max(0, start))
        let upper = line.index(lower, offsetBy: min(length, line.distance(from: lower, to: line.endIndex)), limitedBy: line.endIndex) ?? line.endIndex
        return String(line[lower..<upper])
    }

    private static func canonicalPair(_ a: Int, _ b: Int) -> String {
        if a < b { return "\(a)-\(b)" }
        return "\(b)-\(a)"
    }
}
