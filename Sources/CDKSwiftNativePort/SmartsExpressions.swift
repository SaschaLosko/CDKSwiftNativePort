import Foundation

// Boolean precedence follows Daylight SMARTS: !, &/adjacency, comma, semicolon.
struct SmartsAtomReader {
    let chars: [Character]
    let depth: Int
    var position = 0
    var stereo: (CDKQueryAtomStereo.Geometry, Int, Bool)?
    var hydrogenLigand = false
    var sawElement = false
    var peek: Character? { position < chars.count ? chars[position] : nil }

    mutating func read() throws -> CDKAtomExpression {
        guard !chars.isEmpty else { throw CDKQueryError.invalid("Empty atom expression.") }
        guard chars.count <= 4096 else { throw CDKQueryError.limit("Atom expression size limit.") }
        let expression = try low()
        guard position == chars.count else { throw CDKQueryError.unsupported("Unsupported atom expression.") }
        return expression
    }
    mutating func low() throws -> CDKAtomExpression {
        var result = try disjunction()
        while peek == ";" { position += 1; result = .and(result, try disjunction()) }
        return result
    }
    mutating func disjunction() throws -> CDKAtomExpression {
        var result = try conjunction()
        while peek == "," { position += 1; result = .or(result, try conjunction()) }
        return result
    }
    mutating func conjunction() throws -> CDKAtomExpression {
        var result = try unary()
        while let token = peek, token != ",", token != ";" {
            if token == "&" { position += 1 }
            result = .and(result, try unary())
        }
        return result
    }
    mutating func unary(_ depth: Int = 0) throws -> CDKAtomExpression {
        guard depth < 64 else { throw CDKQueryError.limit("Predicate nesting limit.") }
        if peek == "!" {
            position += 1
            guard peek != "@" else {
                throw CDKQueryError.unsupported("Negated stereo requires an explicit stereo expression.")
            }
            return .not(try unary(depth + 1))
        }
        guard let token = peek else { throw CDKQueryError.invalid("Missing atom predicate.") }
        position += 1
        if token == "$" {
            guard peek == "(" else { throw CDKQueryError.invalid("Recursive SMARTS requires parentheses.") }
            position += 1; let start = position; var nesting = 1
            while position < chars.count && nesting > 0 {
                if chars[position] == "(" { nesting += 1 }
                if chars[position] == ")" { nesting -= 1 }
                if nesting > 0 { position += 1 }
            }
            guard nesting == 0 else { throw CDKQueryError.invalid("Unclosed recursive query.") }
            let text = String(chars[start..<position]); position += 1
            return .recursive(try CDKSmartsParser.parse(text, depth: depth + 1))
        }
        if token == "@" {
            guard stereo == nil, !chars.contains(","), !chars.contains("!") else {
                throw CDKQueryError.unsupported("Boolean stereo alternatives require separate patterns.")
            }
            var geometry = CDKQueryAtomStereo.Geometry.tetrahedral
            var configuration = 1
            if peek == "@" {
                configuration = 2; position += 1
            } else if let next = peek, ["T", "S", "O", "A"].contains(next) {
                let end = min(position + 2, chars.count)
                let code = String(chars[position..<end])
                let geometries: [String: CDKQueryAtomStereo.Geometry] = [
                    "TH": .tetrahedral, "SP": .squarePlanar, "TB": .trigonalBipyramidal, "OH": .octahedral,
                ]
                guard let selected = geometries[code] else {
                    throw CDKQueryError.unsupported("Unsupported stereo class.")
                }
                geometry = selected; position = end; configuration = try number(required: true)!
            }
            guard (1...geometry.configurationCount).contains(configuration) else {
                throw CDKQueryError.invalid("Stereo configuration is outside its class.")
            }
            let optional = peek == "?"; if optional { position += 1 }
            stereo = (geometry, configuration, optional); return .any
        }
        if token == "*" { sawElement = true; return .any }
        if token == "#" {
            let value = try number(required: true)!
            guard (1..<QueryElements.symbols.count).contains(value) else {
                throw CDKQueryError.invalid("Invalid atomic number.")
            }
            sawElement = true; return .element(value)
        }
        if token.isNumber {
            position -= 1; return .isotope(try number(required: true)!)
        }
        if token == "+" || token == "-" {
            var value = try number(required: false) ?? 1
            while peek == token { value += 1; position += 1 }
            return .charge(token == "+" ? value : -value)
        }
        if token == ":" { _ = try number(required: true); return .any }  // Molecule atom maps are annotations.
        // Recognize two-letter elements before single-letter predicate symbols (He, Rn, Xe, Dy...).
        if let next = peek, next.isLowercase {
            let symbol = String(token).uppercased() + String(next)
            if let atomicNumber = QueryElements.number(symbol) {
                position += 1; sawElement = true
                guard token.isUppercase || ["as", "se"].contains(String([token, next])) else {
                    throw CDKQueryError.invalid("Invalid aromatic element.")
                }
                return .and(.element(atomicNumber), .aromatic(token.isLowercase))
            }
        }
        let value: Int?
        switch token {
        case "D", "X", "v", "H", "h", "R", "r", "x": value = try number(required: false)
        default: value = nil
        }
        switch token {
        case "D": return .degree(value ?? 1)
        case "X": return .totalDegree(value ?? 1)
        case "v": return .valence(value ?? 1)
        case "H":
            if !sawElement && value == nil && (peek == nil || peek == "+" || peek == "-" || peek == ":") {
                sawElement = true; return .element(1)
            }
            hydrogenLigand = (value ?? 1) == 1; return .hydrogens(value ?? 1)
        case "h": return value.map(CDKAtomExpression.implicitHydrogens) ?? .hasImplicitHydrogen
        case "R": return value.map(CDKAtomExpression.ringCount) ?? .inRing
        case "r": return value.map(CDKAtomExpression.smallestRing) ?? .inRing
        case "x": return value.map(CDKAtomExpression.ringBonds) ?? .inRing
        case "a": sawElement = true; return .aromatic(true)
        case "A": sawElement = true; return .aromatic(false)
        default:
            guard let atomicNumber = QueryElements.number(String(token).uppercased()),
                token.isUppercase || "bcnops".contains(token)
            else { throw CDKQueryError.unsupported("Unsupported SMARTS atom primitive: \(token).") }
            sawElement = true
            return .and(.element(atomicNumber), .aromatic(token.isLowercase))
        }
    }
    mutating func number(required: Bool) throws -> Int? {
        let start = position
        while peek?.isNumber == true { position += 1 }
        if start == position {
            if required { throw CDKQueryError.invalid("Expected an integer.") }
            return nil
        }
        guard let number = Int(String(chars[start..<position])) else {
            throw CDKQueryError.invalid("Integer overflow.")
        }
        return number
    }
}

struct SmartsBondReader {
    let chars: [Character]
    var position = 0
    var peek: Character? { position < chars.count ? chars[position] : nil }
    mutating func read() throws -> CDKBondExpression {
        let result = try low()
        guard position == chars.count else { throw CDKQueryError.invalid("Invalid bond expression.") }
        return result
    }
    mutating func low() throws -> CDKBondExpression {
        var result = try disjunction()
        while peek == ";" { position += 1; result = .and(result, try disjunction()) }
        return result
    }
    mutating func disjunction() throws -> CDKBondExpression {
        var result = try conjunction()
        while peek == "," { position += 1; result = .or(result, try conjunction()) }
        return result
    }
    mutating func conjunction() throws -> CDKBondExpression {
        var result = try unary()
        while let token = peek, token != ",", token != ";" {
            if token == "&" { position += 1 }
            result = .and(result, try unary())
        }
        return result
    }
    mutating func unary(_ depth: Int = 0) throws -> CDKBondExpression {
        guard depth < 64 else { throw CDKQueryError.limit("Predicate nesting limit.") }
        guard let token = peek else { throw CDKQueryError.invalid("Missing bond predicate.") }
        position += 1
        switch token {
        case "!": return .not(try unary(depth + 1))
        case "~": return .any
        case "-", "/", "\\": return .order(.single)
        case "=": return .order(.double)
        case "#": return .order(.triple)
        case ":": return .order(.aromatic)
        case "@": return .inRing
        default: throw CDKQueryError.unsupported("Unsupported bond predicate.")
        }
    }
}
