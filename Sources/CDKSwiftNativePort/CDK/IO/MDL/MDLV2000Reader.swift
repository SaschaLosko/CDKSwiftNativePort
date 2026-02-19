import Foundation
import CoreGraphics

public enum CDKMDLV2000Reader {

    public static func read(text: String) throws -> Molecule {
        let normalized = text
            .replacingOccurrences(of: "\r\n", with: "\n")
            .replacingOccurrences(of: "\r", with: "\n")
        return try read(lines: normalized.components(separatedBy: "\n"))
    }

    public static func read(lines: [String]) throws -> Molecule {
        let trimmed = dropTrailingEmptyLines(lines)
        guard trimmed.count >= 4 else {
            throw ChemError.parseFailed("Molfile too short.")
        }

        let countsLine = trimmed[3]
        if countsLine.uppercased().contains("V3000") {
            return try CDKMDLV3000Reader.read(text: trimmed.joined(separator: "\n"))
        }

        let atomCount = parseCountsValue(in: countsLine, at: 0)
        let bondCount = parseCountsValue(in: countsLine, at: 1)

        guard atomCount > 0 else {
            throw ChemError.parseFailed("Atom count missing/invalid.")
        }
        guard bondCount >= 0 else {
            throw ChemError.parseFailed("Bond count missing/invalid.")
        }

        let firstAtomLine = 4
        let firstBondLine = firstAtomLine + atomCount
        let firstPropertiesLine = firstBondLine + bondCount

        guard trimmed.count >= firstPropertiesLine else {
            throw ChemError.parseFailed("Molfile truncated.")
        }

        var atoms: [Atom] = []
        atoms.reserveCapacity(atomCount)

        for atomIndex in 0..<atomCount {
            let line = trimmed[firstAtomLine + atomIndex]
            let parsedAtom = try parseAtomLine(line, atomID: atomIndex + 1)
            atoms.append(parsedAtom)
        }

        var bonds: [Bond] = []
        bonds.reserveCapacity(bondCount)

        for bondIndex in 0..<bondCount {
            let line = trimmed[firstBondLine + bondIndex]
            let parsedBond = try parseBondLine(line,
                                               bondID: bondIndex + 1,
                                               atomCount: atomCount)
            bonds.append(parsedBond)
        }

        let name = trimmed[0].trimmingCharacters(in: .whitespacesAndNewlines)
        var molecule = Molecule(name: name.isEmpty ? "Molecule" : name, atoms: atoms, bonds: bonds)

        var foundMEnd = false
        var lineIndex = firstPropertiesLine
        while lineIndex < trimmed.count {
            let line = trimmed[lineIndex]
            if line.hasPrefix("M  END") {
                foundMEnd = true
                break
            }

            applyPropertyLine(line, to: &molecule)
            lineIndex += 1
        }

        guard foundMEnd else {
            throw ChemError.parseFailed("Molfile missing M  END record.")
        }

        if let box = molecule.boundingBox(), box.width <= 0.0001 && box.height <= 0.0001 {
            molecule = Depiction2DGenerator.generate(for: molecule)
        }

        return molecule
    }

    private static func parseCountsValue(in countsLine: String, at index: Int) -> Int {
        if let fixed = parseFixedInt(countsLine, start: index * 3, length: 3) {
            return fixed
        }
        let parts = countsLine.split(whereSeparator: \.isWhitespace)
        guard index < parts.count else { return -1 }
        return Int(parts[index]) ?? -1
    }

    private static func parseAtomLine(_ line: String, atomID: Int) throws -> Atom {
        let parts = line.split(whereSeparator: \.isWhitespace).map(String.init)
        guard parts.count >= 4 else {
            throw ChemError.parseFailed("Invalid atom line at atom index \(atomID).")
        }

        let x = parseFixedDouble(line, start: 0, length: 10) ?? Double(parts[0]) ?? 0
        let y = parseFixedDouble(line, start: 10, length: 10) ?? Double(parts[1]) ?? 0

        let symbol = parseElementSymbol(line: line, parts: parts)
        let chargeCode = parseFixedInt(line, start: 36, length: 3) ?? parseInt(parts, index: 5) ?? 0
        let (charge, radicalFromChargeCode) = decodeChargeAndRadical(chargeCode)

        return Atom(
            id: atomID,
            element: symbol,
            position: CGPoint(x: x, y: y),
            charge: charge,
            queryType: queryType(for: symbol),
            radical: radicalFromChargeCode
        )
    }

    private static func parseBondLine(_ line: String,
                                      bondID: Int,
                                      atomCount: Int) throws -> Bond {
        let parts = line.split(whereSeparator: \.isWhitespace).map(String.init)
        guard parts.count >= 3 else {
            throw ChemError.parseFailed("Invalid bond line at bond index \(bondID).")
        }

        let a1 = parseFixedInt(line, start: 0, length: 3) ?? parseInt(parts, index: 0) ?? 0
        let a2 = parseFixedInt(line, start: 3, length: 3) ?? parseInt(parts, index: 1) ?? 0
        let typeCode = parseFixedInt(line, start: 6, length: 3) ?? parseInt(parts, index: 2) ?? 1
        let stereoCode = parseFixedInt(line, start: 9, length: 3) ?? parseInt(parts, index: 3) ?? 0

        guard (1...atomCount).contains(a1), (1...atomCount).contains(a2) else {
            throw ChemError.parseFailed("Bond references unknown atom index.")
        }

        let (order, queryType) = bondFromTypeCode(typeCode)

        return Bond(id: bondID,
                    a1: a1,
                    a2: a2,
                    order: order,
                    stereo: BondStereo.fromMolfile(stereoCode),
                    queryType: queryType)
    }

    private static func applyPropertyLine(_ line: String, to molecule: inout Molecule) {
        guard line.hasPrefix("M  ") else { return }

        let fields = line.split(whereSeparator: \.isWhitespace).map(String.init)
        guard fields.count >= 2 else { return }

        switch fields[1] {
        case "CHG":
            applyChargeLine(fields, to: &molecule)
        case "ISO":
            applyIsotopeLine(fields, to: &molecule)
        case "RAD":
            applyRadicalLine(fields, to: &molecule)
        case "RGP":
            applyRGroupLine(fields, to: &molecule)
        case "ALS":
            applyAtomListLine(fields, to: &molecule)
        case "APO":
            applyAtomIntProperty(fields, to: &molecule) { atom, value in
                atom.attachmentPoint = value
            }
        case "SUB":
            applyAtomIntProperty(fields, to: &molecule) { atom, value in
                atom.substitutionCount = value
            }
        case "UNS":
            applyAtomIntProperty(fields, to: &molecule) { atom, value in
                atom.unsaturated = value
            }
        case "RBC":
            applyAtomIntProperty(fields, to: &molecule) { atom, value in
                atom.ringBondCount = value
            }
        default:
            break
        }
    }

    private static func applyChargeLine(_ fields: [String], to molecule: inout Molecule) {
        for (atomID, charge) in parseAtomValuePairs(fields) {
            guard let atomIndex = molecule.indexOfAtom(id: atomID) else { continue }
            molecule.atoms[atomIndex].charge = charge
        }
    }

    private static func applyIsotopeLine(_ fields: [String], to molecule: inout Molecule) {
        for (atomID, mass) in parseAtomValuePairs(fields) {
            guard let atomIndex = molecule.indexOfAtom(id: atomID) else { continue }
            molecule.atoms[atomIndex].isotopeMassNumber = mass
        }
    }

    private static func applyRadicalLine(_ fields: [String], to molecule: inout Molecule) {
        for (atomID, radical) in parseAtomValuePairs(fields) {
            guard let atomIndex = molecule.indexOfAtom(id: atomID) else { continue }
            molecule.atoms[atomIndex].radical = radical
        }
    }

    private static func applyRGroupLine(_ fields: [String], to molecule: inout Molecule) {
        for (atomID, rLabel) in parseAtomValuePairs(fields) {
            guard let atomIndex = molecule.indexOfAtom(id: atomID) else { continue }
            molecule.atoms[atomIndex].rGroupLabel = rLabel
        }
    }

    private static func applyAtomListLine(_ fields: [String], to molecule: inout Molecule) {
        guard fields.count >= 6,
              let atomID = Int(fields[2]),
              let entryCount = Int(fields[3]),
              entryCount > 0,
              let atomIndex = molecule.indexOfAtom(id: atomID) else {
            return
        }

        let negated = fields[4].uppercased() == "T"
        let rawEntries = fields.dropFirst(5)
        var entries = rawEntries.map { normalizeElementToken($0) }
        if entries.count > entryCount {
            entries = Array(entries.prefix(entryCount))
        }
        guard entries.count == entryCount else { return }

        molecule.atoms[atomIndex].atomList = entries
        molecule.atoms[atomIndex].atomListIsNegated = negated
        molecule.atoms[atomIndex].queryType = .anyAtom
        molecule.atoms[atomIndex].element = "L"
    }

    private static func applyAtomIntProperty(_ fields: [String],
                                             to molecule: inout Molecule,
                                             assign: (inout Atom, Int) -> Void) {
        for (atomID, value) in parseAtomValuePairs(fields) {
            guard let atomIndex = molecule.indexOfAtom(id: atomID) else { continue }
            assign(&molecule.atoms[atomIndex], value)
        }
    }

    private static func parseAtomValuePairs(_ fields: [String]) -> [(Int, Int)] {
        guard fields.count >= 4 else { return [] }

        let pairCount = Int(fields[2]) ?? 0
        guard pairCount > 0 else { return [] }

        var pairs: [(Int, Int)] = []
        pairs.reserveCapacity(pairCount)

        for pair in 0..<pairCount {
            let atomTokenIndex = 3 + pair * 2
            let valueTokenIndex = atomTokenIndex + 1
            guard valueTokenIndex < fields.count,
                  let atomID = Int(fields[atomTokenIndex]),
                  let value = Int(fields[valueTokenIndex]) else {
                continue
            }
            pairs.append((atomID, value))
        }

        return pairs
    }

    private static func parseElementSymbol(line: String, parts: [String]) -> String {
        if let fixed = substring(line, start: 31, length: 3)?.trimmingCharacters(in: .whitespaces),
           !fixed.isEmpty {
            return normalizeElementToken(fixed)
        }

        return normalizeElementToken(parts[3])
    }

    private static func normalizeElementToken(_ token: String) -> String {
        let trimmed = token.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { return "C" }
        if trimmed == "*" { return trimmed }
        if trimmed.allSatisfy(\.isLetter) {
            guard let first = trimmed.first else { return trimmed }
            let head = String(first).uppercased()
            let tail = String(trimmed.dropFirst()).lowercased()
            return head + tail
        }
        return trimmed
    }

    private static func queryType(for symbol: String) -> AtomQueryType? {
        switch symbol.uppercased() {
        case "*":
            return .anyAtom
        case "A":
            return .anyNonHydrogen
        case "Q":
            return .anyHetero
        default:
            return nil
        }
    }

    private static func decodeChargeAndRadical(_ code: Int) -> (charge: Int, radical: Int?) {
        switch code {
        case 1: return (3, nil)
        case 2: return (2, nil)
        case 3: return (1, nil)
        case 4: return (0, 2)
        case 5: return (-1, nil)
        case 6: return (-2, nil)
        case 7: return (-3, nil)
        default: return (0, nil)
        }
    }

    private static func bondFromTypeCode(_ code: Int) -> (BondOrder, BondQueryType?) {
        switch code {
        case 1:
            return (.single, nil)
        case 2:
            return (.double, nil)
        case 3:
            return (.triple, nil)
        case 4:
            return (.aromatic, nil)
        case 5:
            return (.single, .singleOrDouble)
        case 6:
            return (.single, .singleOrAromatic)
        case 7:
            return (.double, .doubleOrAromatic)
        case 8:
            return (.single, .any)
        default:
            return (.single, nil)
        }
    }

    private static func parseInt(_ parts: [String], index: Int) -> Int? {
        guard index < parts.count else { return nil }
        return Int(parts[index])
    }

    private static func parseFixedInt(_ s: String, start: Int, length: Int) -> Int? {
        guard let slice = substring(s, start: start, length: length) else { return nil }
        return Int(slice.trimmingCharacters(in: .whitespaces))
    }

    private static func parseFixedDouble(_ s: String, start: Int, length: Int) -> Double? {
        guard let slice = substring(s, start: start, length: length) else { return nil }
        return Double(slice.trimmingCharacters(in: .whitespaces))
    }

    private static func substring(_ s: String, start: Int, length: Int) -> String? {
        guard start >= 0, length > 0 else { return nil }
        guard s.count >= start + length else { return nil }

        let lower = s.index(s.startIndex, offsetBy: start)
        let upper = s.index(lower, offsetBy: length)
        return String(s[lower..<upper])
    }

    private static func dropTrailingEmptyLines(_ lines: [String]) -> [String] {
        var out = lines
        while let last = out.last, last.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty {
            out.removeLast()
        }
        return out
    }
}
