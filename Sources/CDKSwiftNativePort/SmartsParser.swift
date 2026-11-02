import Foundation

/// Molecular SMARTS parser. Unsupported extensions throw; they are never read as plain SMILES.
public enum CDKSmartsParser {
    public static func parse(_ text: String) throws -> CDKQueryPattern { try parse(text, depth: 0) }
    static func parse(_ text: String, depth: Int) throws -> CDKQueryPattern {
        guard depth <= 16, text.utf8.count <= 32_768 else {
            throw CDKQueryError.limit("SMARTS size or recursion limit.")
        }
        var reader = SmartsGraphReader(chars: Array(text), depth: depth)
        return try reader.read()
    }
}

private struct SmartsGraphReader {
    let chars: [Character]
    let depth: Int
    var position = 0
    var molecule = Molecule()
    var atoms: [Int: CDKAtomExpression] = [:]
    var bonds: [Int: CDKBondExpression] = [:]
    var stereo: [CDKQueryAtomStereo] = []
    var orders: [Int: [Int]] = [:]
    var current: Int?
    var branches: [(Int, Int)] = []
    var rings: [Int: (Int, Int, CDKBondExpression?, BondStereo)] = [:]
    var pending: CDKBondExpression?
    var direction = BondStereo.none
    var peek: Character? { position < chars.count ? chars[position] : nil }
    let defaultBond = CDKBondExpression.or(.order(.single), .order(.aromatic))

    mutating func read() throws -> CDKQueryPattern {
        guard !chars.isEmpty else { throw CDKQueryError.invalid("Empty SMARTS.") }
        while let token = peek {
            if molecule.atoms.count > 1_000 { throw CDKQueryError.limit("More than 1000 query atoms.") }
            switch token {
            case "(":
                guard let current, pending == nil else {
                    throw CDKQueryError.unsupported("Component grouping or invalid branch.")
                }
                guard branches.count < 128 else { throw CDKQueryError.limit("Branch nesting limit.") }
                branches.append((current, molecule.atoms.count)); position += 1
            case ")":
                guard let branch = branches.popLast(), molecule.atoms.count > branch.1, current != nil, pending == nil
                else {
                    throw CDKQueryError.invalid("Unmatched or empty branch.")
                }
                current = branch.0; position += 1
            case ".":
                guard current != nil, pending == nil else {
                    throw CDKQueryError.invalid("Invalid component separator.")
                }
                current = nil; position += 1
            case "0"..."9", "%": try ring()
            case "-", "=", "#", ":", "~", "@", "!", "/", "\\": try bond()
            default: try atom()
            }
        }
        guard current != nil, pending == nil, branches.isEmpty, rings.isEmpty else {
            throw CDKQueryError.invalid("Incomplete SMARTS graph.")
        }
        for index in stereo.indices {
            stereo[index].carriers = orders[stereo[index].focus, default: []]
            guard stereo[index].carriers.count == stereo[index].geometry.carrierCount else {
                throw CDKQueryError.invalid("Stereo centre has the wrong number of ordered ligands.")
            }
        }
        return CDKQueryPattern(molecule: molecule, atoms: atoms, bonds: bonds, stereo: stereo)
    }
    mutating func atom() throws {
        let content: [Character]
        if peek == "[" {
            position += 1; let start = position; var nesting = 1
            while position < chars.count && nesting > 0 {
                if chars[position] == "[" { nesting += 1 }
                if chars[position] == "]" { nesting -= 1 }
                if nesting > 0 { position += 1 }
            }
            guard nesting == 0 else { throw CDKQueryError.invalid("Unclosed atom bracket.") }
            content = Array(chars[start..<position]); position += 1
        } else {
            guard let token = peek, "BCNOPSFIbcnopsaA*".contains(token) else {
                throw CDKQueryError.unsupported("Unsupported SMARTS token at character \(position).")
            }
            position += 1
            if (token == "C" && peek == "l") || (token == "B" && peek == "r") {
                content = [token, chars[position]]; position += 1
            } else {
                content = [token]
            }
        }
        var reader = SmartsAtomReader(chars: content, depth: depth)
        let expression = try reader.read()
        let id = molecule.atoms.count
        molecule.atoms.append(Atom(id: id, element: "*", position: .zero))
        atoms[id] = expression
        if let previous = current { try connect(previous, id, pending ?? defaultBond, direction) }
        if let (geometry, configuration, optional) = reader.stereo {
            if reader.hydrogenLigand {
                guard geometry == .tetrahedral else {
                    throw CDKQueryError.unsupported("Implicit coordination ligand needs explicit metadata.")
                }
                orders[id, default: []].append(id)
            }
            stereo.append(
                .init(
                    focus: id, carriers: [], geometry: geometry, configuration: configuration,
                    allowsUnspecified: optional))
        }
        current = id; pending = nil; direction = .none
    }
    mutating func bond() throws {
        guard current != nil, pending == nil else { throw CDKQueryError.invalid("Bond without an atom.") }
        let start = position
        while let token = peek, "-=#:~@!&,;/\\?".contains(token) { position += 1 }
        let text = Array(chars[start..<position])
        if text.contains("/") || text.contains("\\") {
            guard text.count == 1 else { throw CDKQueryError.unsupported("Compound directional bond expression.") }
            direction = text[0] == "/" ? .up : .down
        }
        var reader = SmartsBondReader(chars: text); pending = try reader.read()
    }
    mutating func ring() throws {
        guard let current else { throw CDKQueryError.invalid("Ring label without an atom.") }
        let label: Int
        if peek == "%" {
            position += 1
            guard position + 1 < chars.count, let value = Int(String(chars[position...position + 1])) else {
                throw CDKQueryError.invalid("Invalid two-digit ring label.")
            }
            label = value; position += 2
        } else {
            label = chars[position].wholeNumberValue!; position += 1
        }
        if let (other, orderIndex, stored, storedDirection) = rings.removeValue(forKey: label) {
            guard pending == nil || stored == nil || pending == stored else {
                throw CDKQueryError.invalid("Conflicting ring bonds.")
            }
            guard direction == .none && storedDirection == .none else {
                throw CDKQueryError.unsupported("Directional ring closure.")
            }
            try connect(other, current, pending ?? stored ?? defaultBond, .none)
            orders[other]!.removeLast(); orders[other]![orderIndex] = current
        } else {
            let index = orders[current, default: []].count
            orders[current, default: []].append(-1)
            rings[label] = (current, index, pending, direction)
        }
        pending = nil; direction = .none
    }
    mutating func connect(_ a: Int, _ b: Int, _ expression: CDKBondExpression, _ stereo: BondStereo) throws {
        guard a != b, molecule.bond(between: a, and: b) == nil else {
            throw CDKQueryError.invalid("Self bond or duplicate edge.")
        }
        let id = molecule.bonds.count
        let order: BondOrder
        if case .order(let specified) = expression { order = specified } else { order = .single }
        molecule.bonds.append(Bond(id: id, a1: a, a2: b, order: order, stereo: stereo))
        bonds[id] = expression; orders[a, default: []].append(b); orders[b, default: []].append(a)
    }
}
