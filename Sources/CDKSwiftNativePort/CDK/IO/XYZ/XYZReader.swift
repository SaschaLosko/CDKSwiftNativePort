import Foundation
import CoreGraphics

/// CDK-style XYZ reader for one or more XYZ blocks in a single text stream.
public enum CDKXYZReader {
    public static func read(text: String) throws -> [Molecule] {
        let normalized = text
            .replacingOccurrences(of: "\r\n", with: "\n")
            .replacingOccurrences(of: "\r", with: "\n")
        let lines = normalized.components(separatedBy: "\n")

        var molecules: [Molecule] = []
        var lineIndex = 0

        while lineIndex < lines.count {
            while lineIndex < lines.count,
                  lines[lineIndex].trimmingCharacters(in: .whitespacesAndNewlines).isEmpty {
                lineIndex += 1
            }
            guard lineIndex < lines.count else { break }

            let header = lines[lineIndex].trimmingCharacters(in: .whitespacesAndNewlines)
            if let atomCount = Int(header), atomCount > 0 {
                let comment = lineIndex + 1 < lines.count ? lines[lineIndex + 1] : ""
                let startAtoms = lineIndex + 2
                var atomLines: [String] = []
                atomLines.reserveCapacity(atomCount)

                var scan = startAtoms
                while scan < lines.count && atomLines.count < atomCount {
                    let candidate = lines[scan].trimmingCharacters(in: .whitespacesAndNewlines)
                    if !candidate.isEmpty {
                        atomLines.append(lines[scan])
                    }
                    scan += 1
                }

                guard atomLines.count == atomCount else {
                    throw ChemError.parseFailed("XYZ atom count does not match coordinate rows.")
                }

                let moleculeName = comment.trimmingCharacters(in: .whitespacesAndNewlines)
                molecules.append(try makeMolecule(name: moleculeName,
                                                  atomLines: atomLines,
                                                  fallbackName: "XYZ Molecule \(molecules.count + 1)"))
                lineIndex = scan
                continue
            }

            let remainder = Array(lines[lineIndex...])
            let atomLines = remainder.filter { !$0.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty }
            guard !atomLines.isEmpty else {
                throw ChemError.emptyInput
            }
            molecules.append(try makeMolecule(name: "",
                                              atomLines: atomLines,
                                              fallbackName: "XYZ Molecule"))
            break
        }

        guard !molecules.isEmpty else { throw ChemError.emptyInput }
        return molecules
    }

    private static func makeMolecule(name: String,
                                     atomLines: [String],
                                     fallbackName: String) throws -> Molecule {
        var atoms: [Atom] = []
        atoms.reserveCapacity(atomLines.count)

        for (index, line) in atomLines.enumerated() {
            guard let parsed = parseAtomLine(line, atomID: index + 1) else {
                throw ChemError.parseFailed("Invalid XYZ atom row: \(line)")
            }
            atoms.append(parsed)
        }

        var molecule = Molecule(name: name.isEmpty ? fallbackName : name,
                                atoms: atoms,
                                bonds: CDKBondPerception.inferSingleBonds(for: atoms))

        if let box = molecule.boundingBox(), box.width <= 0.0001 && box.height <= 0.0001 {
            molecule = Depiction2DGenerator.generate(for: molecule)
        }
        return molecule
    }

    private static func parseAtomLine(_ line: String, atomID: Int) -> Atom? {
        let parts = line.split(whereSeparator: \.isWhitespace).map(String.init)
        guard parts.count >= 4 else { return nil }

        let firstIsNumeric = Double(parts[0]) != nil

        let symbolToken: String
        let xToken: String
        let yToken: String

        if firstIsNumeric {
            guard parts.count >= 4 else { return nil }
            xToken = parts[0]
            yToken = parts[1]
            symbolToken = parts[3]
        } else {
            symbolToken = parts[0]
            xToken = parts[1]
            yToken = parts[2]
        }

        guard let x = Double(xToken), let y = Double(yToken) else {
            return nil
        }

        let symbol = normalizeElementSymbol(from: symbolToken)
        return Atom(id: atomID,
                    element: symbol,
                    position: CGPoint(x: x, y: y))
    }

    private static func normalizeElementSymbol(from raw: String) -> String {
        let cleaned = raw.trimmingCharacters(in: .whitespacesAndNewlines)
        if cleaned.isEmpty { return "C" }
        return CDKDescriptorSupport.canonicalElementSymbol(cleaned)
    }
}
