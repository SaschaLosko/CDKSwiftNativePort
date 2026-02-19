import Foundation

/// Swift counterpart of CDK's `org.openscience.cdk.smiles.SmilesParser`.
public final class CDKSmilesParser {
    public let flavor: CDKSmiFlavor

    public init(flavor: CDKSmiFlavor = .cdkDefault) {
        self.flavor = flavor
    }

    public func parseSmiles(_ smiles: String) throws -> Molecule {
        let split = try CDKCxSmilesParser.split(smiles, enabled: flavor.contains(.cxsmiles))
        var molecule = try parseCoreSmiles(split.coreSmiles)
        CDKCxSmilesParser.applyAtomLabels(to: &molecule, state: split.state)
        if let title = split.title, !title.isEmpty {
            molecule.name = title
        }
        return molecule
    }

    public func parseCoreSmiles(_ smiles: String) throws -> Molecule {
        let source = smiles.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !source.isEmpty else { throw ChemError.emptyInput }

        struct BracketAtomDescriptor {
            var element: String
            var aromatic: Bool
            var charge: Int
            var isotopeMass: Int?
            var chirality: AtomChirality
            var explicitHydrogenCount: Int
            var hasExplicitHydrogenSpecification: Bool
            var exactDegree: Int?
            var totalConnectivity: Int?
            var valence: Int?
            var ringMembership: Int?
            var ringSize: Int?
            var unsaturated: Int?
            var atomClass: Int?
            var atomMapNumber: Int?
        }

        struct RingClosureState {
            let atomID: Int
            let explicitOrder: BondOrder?
            let directionalStereo: BondStereo?
        }

        var molecule = Molecule(name: "SMILES")
        var nextAtomID = 0
        var nextBondID = 0

        var currentAtom: Int? = nil
        var branchStack: [Int] = []
        var ringClosures: [Int: RingClosureState] = [:]

        var pendingOrder: BondOrder? = nil
        var pendingDirectionalStereo: BondStereo? = nil

        var i = source.startIndex

        func peek(_ offset: Int = 0) -> Character? {
            var cursor = i
            for _ in 0..<offset {
                guard cursor < source.endIndex else { return nil }
                cursor = source.index(after: cursor)
            }
            return cursor < source.endIndex ? source[cursor] : nil
        }

        func advance() {
            i = source.index(after: i)
        }

        func readDigits(maxCount: Int? = nil) -> String {
            var digits = ""
            while let ch = peek(), ch.isNumber {
                if let maxCount, digits.count >= maxCount { break }
                digits.append(ch)
                advance()
            }
            return digits
        }

        func normalizeElement(symbol: String,
                              aromatic: Bool,
                              isotopeMass: Int?) -> (symbol: String, aromatic: Bool, isotopeMass: Int?) {
            switch symbol.uppercased() {
            case "D":
                return ("H", aromatic, isotopeMass ?? 2)
            case "T":
                return ("H", aromatic, isotopeMass ?? 3)
            default:
                let normalized = symbol.prefix(1).uppercased() + symbol.dropFirst().lowercased()
                return (normalized, aromatic, isotopeMass)
            }
        }

        func appendAtom(element: String,
                        aromatic: Bool = false,
                        charge: Int = 0,
                        isotopeMass: Int? = nil,
                        chirality: AtomChirality = .none,
                        explicitHydrogenCount: Int? = nil,
                        exactDegree: Int? = nil,
                        totalConnectivity: Int? = nil,
                        valence: Int? = nil,
                        ringMembership: Int? = nil,
                        ringSize: Int? = nil,
                        unsaturated: Int? = nil,
                        atomClass: Int? = nil,
                        atomMapNumber: Int? = nil) -> Int {
            nextAtomID += 1
            var atom = Atom(id: nextAtomID,
                            element: element,
                            position: .zero,
                            charge: charge,
                            isotopeMassNumber: isotopeMass,
                            aromatic: aromatic,
                            chirality: chirality,
                            explicitHydrogenCount: explicitHydrogenCount,
                            substitutionCount: exactDegree ?? totalConnectivity,
                            unsaturated: unsaturated ?? valence,
                            ringBondCount: ringMembership ?? ringSize,
                            atomClass: atomClass,
                            atomMapNumber: atomMapNumber)
            atom.charge = charge
            molecule.atoms.append(atom)
            return atom.id
        }

        func appendBond(_ a1: Int, _ a2: Int, order: BondOrder, stereo: BondStereo = .none) {
            nextBondID += 1
            molecule.bonds.append(Bond(id: nextBondID, a1: a1, a2: a2, order: order, stereo: stereo))
        }

        func atomIsAromatic(_ atomID: Int) -> Bool {
            molecule.atom(id: atomID)?.aromatic == true
        }

        func resolveBondOrder(a1: Int, a2: Int, explicit: BondOrder?) -> BondOrder {
            if let explicit { return explicit }
            if flavor.contains(.useAromaticSymbols) && atomIsAromatic(a1) && atomIsAromatic(a2) {
                return .aromatic
            }
            return .single
        }

        func parseRingIndex() throws -> Int? {
            guard let ch = peek() else { return nil }
            if ch.isNumber {
                advance()
                return Int(String(ch))
            }
            if ch == "%" {
                advance()
                let digits = readDigits(maxCount: 2)
                guard digits.count == 2, let idx = Int(digits) else {
                    throw ChemError.parseFailed("Invalid ring index after '%' in SMILES.")
                }
                return idx
            }
            return nil
        }

        func parsePlainAtom() -> (element: String, aromatic: Bool, isotopeMass: Int?)? {
            guard let c0 = peek() else { return nil }
            if c0 == "*" {
                advance()
                return ("*", false, nil)
            }

            if c0.isLowercase {
                if let c1 = peek(1) {
                    let pair = "\(c0)\(c1)"
                    if pair == "se" || pair == "as" {
                        advance()
                        advance()
                        let normalized = normalizeElement(symbol: pair, aromatic: true, isotopeMass: nil)
                        return (normalized.symbol, normalized.aromatic, normalized.isotopeMass)
                    }
                }
                if ["b", "c", "n", "o", "p", "s"].contains(c0) {
                    advance()
                    let normalized = normalizeElement(symbol: String(c0), aromatic: true, isotopeMass: nil)
                    return (normalized.symbol, normalized.aromatic, normalized.isotopeMass)
                }
                return nil
            }

            guard c0.isUppercase else { return nil }
            if let c1 = peek(1) {
                let pair = "\(c0)\(c1)"
                if pair == "Cl" || pair == "Br" {
                    advance()
                    advance()
                    let normalized = normalizeElement(symbol: pair, aromatic: false, isotopeMass: nil)
                    return (normalized.symbol, normalized.aromatic, normalized.isotopeMass)
                }
            }

            advance()
            let normalized = normalizeElement(symbol: String(c0), aromatic: false, isotopeMass: nil)
            return (normalized.symbol, normalized.aromatic, normalized.isotopeMass)
        }

        func parseCharge() -> Int {
            guard let signChar = peek(), signChar == "+" || signChar == "-" else { return 0 }
            let sign = (signChar == "+") ? 1 : -1
            advance()

            if let first = peek(), first.isNumber {
                let digits = readDigits()
                return sign * (Int(digits) ?? 1)
            }

            var magnitude = 1
            while let repeated = peek(), repeated == signChar {
                magnitude += 1
                advance()
            }
            return sign * magnitude
        }

        func parseBracketAtom() throws -> BracketAtomDescriptor {
            guard peek() == "[" else {
                throw ChemError.parseFailed("Invalid bracket atom start.")
            }
            advance() // '['

            let isotopeDigits = readDigits()
            let isotopeMass = isotopeDigits.isEmpty ? nil : Int(isotopeDigits)

            guard let head = peek() else {
                throw ChemError.parseFailed("Unterminated bracket atom.")
            }

            let rawSymbol: String
            let aromatic: Bool

            if head == "*" {
                rawSymbol = "*"
                aromatic = false
                advance()
            } else if head.isUppercase {
                var symbol = String(head)
                advance()
                if let tail = peek(), tail.isLowercase {
                    symbol.append(tail)
                    advance()
                }
                rawSymbol = symbol
                aromatic = false
            } else if head.isLowercase {
                if let c1 = peek(1), "\(head)\(c1)" == "se" || "\(head)\(c1)" == "as" {
                    rawSymbol = "\(head)\(c1)"
                    advance()
                    advance()
                } else {
                    rawSymbol = String(head)
                    advance()
                }
                aromatic = true
            } else {
                throw ChemError.parseFailed("Invalid bracket atom element token.")
            }

            let normalized = normalizeElement(symbol: rawSymbol, aromatic: aromatic, isotopeMass: isotopeMass)
            var descriptor = BracketAtomDescriptor(element: normalized.symbol,
                                                   aromatic: normalized.aromatic,
                                                   charge: 0,
                                                   isotopeMass: normalized.isotopeMass,
                                                   chirality: .none,
                                                   explicitHydrogenCount: 0,
                                                   hasExplicitHydrogenSpecification: false,
                                                   exactDegree: nil,
                                                   totalConnectivity: nil,
                                                   valence: nil,
                                                   ringMembership: nil,
                                                   ringSize: nil,
                                                   unsaturated: nil,
                                                   atomClass: nil,
                                                   atomMapNumber: nil)

            while let ch = peek(), ch != "]" {
                if ch == "@" {
                    advance()
                    var atCount = 1
                    while peek() == "@" {
                        atCount += 1
                        advance()
                    }
                    descriptor.chirality = atCount >= 2 ? .anticlockwise : .clockwise

                    // Consume optional chiral class suffix (TH1, AL2, SP1, TBx, OHx).
                    let classStart = i
                    var classToken = ""
                    while let letter = peek(), letter.isLetter {
                        classToken.append(letter)
                        advance()
                    }
                    let rankToken = readDigits()
                    let known = ["TH", "AL", "SP", "TB", "OH"]
                    if !classToken.isEmpty && !known.contains(classToken.uppercased()) {
                        i = classStart
                    } else {
                        _ = rankToken
                    }
                    continue
                }

                if ch == "D" {
                    advance()
                    let digits = readDigits()
                    descriptor.exactDegree = digits.isEmpty ? 1 : Int(digits)
                    continue
                }

                if ch == "X" {
                    advance()
                    let digits = readDigits()
                    descriptor.totalConnectivity = digits.isEmpty ? 1 : Int(digits)
                    continue
                }

                if ch == "v" {
                    advance()
                    let digits = readDigits()
                    descriptor.valence = digits.isEmpty ? 0 : Int(digits)
                    continue
                }

                if ch == "R" {
                    advance()
                    let digits = readDigits()
                    descriptor.ringMembership = digits.isEmpty ? 1 : Int(digits)
                    continue
                }

                if ch == "r" {
                    advance()
                    let digits = readDigits()
                    descriptor.ringSize = digits.isEmpty ? 0 : Int(digits)
                    continue
                }

                if ch == "u" {
                    advance()
                    let digits = readDigits()
                    descriptor.unsaturated = digits.isEmpty ? 1 : Int(digits)
                    continue
                }

                if ch == "H" {
                    advance()
                    let countDigits = readDigits()
                    let hCount = countDigits.isEmpty ? 1 : (Int(countDigits) ?? 1)
                    descriptor.hasExplicitHydrogenSpecification = true
                    descriptor.explicitHydrogenCount += max(0, hCount)
                    continue
                }

                if ch == "+" || ch == "-" {
                    descriptor.charge += parseCharge()
                    continue
                }

                if ch == ":" {
                    advance()
                    let digits = readDigits()
                    if let value = Int(digits), !digits.isEmpty {
                        descriptor.atomClass = value
                        descriptor.atomMapNumber = value
                    } else if flavor.contains(.strict) {
                        throw ChemError.parseFailed("Missing atom-class/map index after ':' in bracket atom.")
                    }
                    continue
                }

                if ch == ";" || ch == "&" || ch == "," {
                    advance()
                    continue
                }

                if flavor.contains(.strict) {
                    throw ChemError.parseFailed("Unsupported bracket decorator '\(ch)' in strict SMILES mode.")
                }
                advance()
            }

            guard peek() == "]" else {
                throw ChemError.parseFailed("Unterminated bracket atom.")
            }
            advance() // ']'
            return descriptor
        }

        func connectCurrentTo(_ newAtomID: Int, aromatic: Bool) {
            guard let cur = currentAtom else {
                currentAtom = newAtomID
                pendingOrder = nil
                pendingDirectionalStereo = nil
                return
            }

            let order = resolveBondOrder(a1: cur, a2: newAtomID, explicit: pendingOrder)
            let stereo = (order == .single) ? (pendingDirectionalStereo ?? .none) : .none
            appendBond(cur, newAtomID, order: order, stereo: stereo)

            if aromatic && order == .single && flavor.contains(.useAromaticSymbols) {
                // Keep aromatic bracket atoms connected consistently.
                if let idx = molecule.bonds.indices.last {
                    if atomIsAromatic(cur) && atomIsAromatic(newAtomID) {
                        molecule.bonds[idx].order = .aromatic
                    }
                }
            }

            currentAtom = newAtomID
            pendingOrder = nil
            pendingDirectionalStereo = nil
        }

        while i < source.endIndex {
            guard let ch = peek() else { break }

            switch ch {
            case "(":
                guard let cur = currentAtom else {
                    throw ChemError.parseFailed("Branch start '(' has no parent atom.")
                }
                branchStack.append(cur)
                advance()

            case ")":
                guard let parent = branchStack.popLast() else {
                    throw ChemError.parseFailed("Branch close ')' without matching '('.")
                }
                currentAtom = parent
                advance()

            case "-", "=", "#", ":", "/", "\\":
                switch ch {
                case "-":
                    pendingOrder = .single
                    pendingDirectionalStereo = nil
                case "=":
                    pendingOrder = .double
                    pendingDirectionalStereo = nil
                case "#":
                    pendingOrder = .triple
                    pendingDirectionalStereo = nil
                case ":":
                    pendingOrder = .aromatic
                    pendingDirectionalStereo = nil
                case "/":
                    pendingOrder = pendingOrder ?? .single
                    pendingDirectionalStereo = .up
                case "\\":
                    pendingOrder = pendingOrder ?? .single
                    pendingDirectionalStereo = .down
                default:
                    break
                }
                advance()

            case ".":
                currentAtom = nil
                pendingOrder = nil
                pendingDirectionalStereo = nil
                advance()

            case "0"..."9", "%":
                guard let ringIndex = try parseRingIndex() else {
                    throw ChemError.parseFailed("Invalid ring closure token in SMILES.")
                }
                guard let cur = currentAtom else {
                    throw ChemError.parseFailed("Ring closure '\(ringIndex)' has no current atom.")
                }

                if let open = ringClosures[ringIndex] {
                    ringClosures[ringIndex] = nil

                    if let left = open.explicitOrder, let right = pendingOrder, left != right {
                        throw ChemError.parseFailed("Conflicting ring bond order for closure \(ringIndex).")
                    }

                    let explicitOrder = pendingOrder ?? open.explicitOrder
                    let order = resolveBondOrder(a1: open.atomID, a2: cur, explicit: explicitOrder)
                    let stereo = (order == .single) ? (pendingDirectionalStereo ?? open.directionalStereo ?? .none) : .none
                    appendBond(open.atomID, cur, order: order, stereo: stereo)
                } else {
                    ringClosures[ringIndex] = RingClosureState(atomID: cur,
                                                               explicitOrder: pendingOrder,
                                                               directionalStereo: pendingDirectionalStereo)
                }

                pendingOrder = nil
                pendingDirectionalStereo = nil

            case "[":
                let bracket = try parseBracketAtom()
                let newID = appendAtom(element: bracket.element,
                                       aromatic: bracket.aromatic,
                                       charge: bracket.charge,
                                       isotopeMass: bracket.isotopeMass,
                                       chirality: bracket.chirality,
                                       explicitHydrogenCount: bracket.hasExplicitHydrogenSpecification ? bracket.explicitHydrogenCount : nil,
                                       exactDegree: bracket.exactDegree,
                                       totalConnectivity: bracket.totalConnectivity,
                                       valence: bracket.valence,
                                       ringMembership: bracket.ringMembership,
                                       ringSize: bracket.ringSize,
                                       unsaturated: bracket.unsaturated,
                                       atomClass: bracket.atomClass,
                                       atomMapNumber: bracket.atomMapNumber)
                connectCurrentTo(newID, aromatic: bracket.aromatic)

            default:
                if let atom = parsePlainAtom() {
                    let newID = appendAtom(element: atom.element,
                                           aromatic: atom.aromatic,
                                           isotopeMass: atom.isotopeMass)
                    connectCurrentTo(newID, aromatic: atom.aromatic)
                } else {
                    throw ChemError.parseFailed("Unexpected token '\(ch)' in SMILES.")
                }
            }
        }

        guard !molecule.atoms.isEmpty else {
            throw ChemError.parseFailed("No atoms found.")
        }
        guard branchStack.isEmpty else {
            throw ChemError.parseFailed("Unterminated branch in SMILES.")
        }
        guard ringClosures.isEmpty else {
            let missing = ringClosures.keys.sorted().map(String.init).joined(separator: ",")
            throw ChemError.parseFailed("Unterminated ring closure(s): \(missing).")
        }
        if flavor.contains(.strict), pendingOrder != nil || pendingDirectionalStereo != nil {
            throw ChemError.parseFailed("Dangling bond token at end of SMILES.")
        }

        try validateAromaticConstraints(in: molecule)
        annotateDirectionalDoubleBonds(in: &molecule)
        molecule = Depiction2DGenerator.generate(for: molecule)
        molecule.assignWedgeHashFromChiralCenters()
        return molecule
    }

    private func validateAromaticConstraints(in molecule: Molecule) throws {
        for atom in molecule.atoms where atom.aromatic {
            let degree = molecule.neighbors(of: atom.id).count
            let element = atom.element.uppercased()

            let hasTripleBond = molecule.bonds(forAtom: atom.id).contains { $0.order == .triple }
            if hasTripleBond {
                throw ChemError.parseFailed("Invalid aromatic atom '\(atom.element)' with triple bond.")
            }

            switch element {
            case "B", "C":
                if degree > 3 {
                    throw ChemError.parseFailed("Aromatic atom '\(atom.element)' exceeds valence-like degree.")
                }
            case "N", "O":
                if element == "O", degree > 2 && atom.charge <= 0 {
                    throw ChemError.parseFailed("Aromatic atom '\(atom.element)' has invalid degree without positive charge.")
                }
            case "P", "S", "SE", "AS":
                if degree > 3 && atom.charge <= 0 {
                    throw ChemError.parseFailed("Aromatic atom '\(atom.element)' has unsupported high degree.")
                }
            case "*":
                continue
            default:
                throw ChemError.parseFailed("Unsupported aromatic atom '\(atom.element)'.")
            }
        }

        let aromaticFiveRings = molecule.simpleCycles(maxSize: 5).filter { ring in
            ring.count == 5 && ring.allSatisfy { molecule.atom(id: $0)?.aromatic == true }
        }
        for ring in aromaticFiveRings {
            let ringAtoms = ring.compactMap { molecule.atom(id: $0) }
            guard ringAtoms.count == ring.count else { continue }

            let hasAromaticNitrogen = ringAtoms.contains { $0.element.uppercased() == "N" }
            guard hasAromaticNitrogen else { continue }

            let hasPyrrolicDonor = ringAtoms.contains { atom in
                let element = atom.element.uppercased()
                if let count = atom.explicitHydrogenCount, count > 0 {
                    return true
                }
                if atom.charge < 0 && ["N", "O", "S", "SE", "AS", "P"].contains(element) {
                    return true
                }
                if ["O", "S", "SE", "AS"].contains(element) {
                    return true
                }
                if element == "N" {
                    let degree = molecule.neighbors(of: atom.id).count
                    if degree > 2 {
                        return true
                    }
                }
                return false
            }

            if !hasPyrrolicDonor {
                throw ChemError.parseFailed("Aromatic five-member N-heterocycle lacks [nH]-like donor.")
            }
        }
    }

    private func annotateDirectionalDoubleBonds(in molecule: inout Molecule) {
        func isDirectional(_ stereo: BondStereo) -> Bool {
            switch stereo {
            case .up, .down, .upReversed, .downReversed:
                return true
            case .none, .either:
                return false
            }
        }

        for idx in molecule.bonds.indices where molecule.bonds[idx].order == .double {
            let db = molecule.bonds[idx]

            let sideAHasDirection = molecule.bonds.contains { bond in
                guard bond.id != db.id, bond.order == .single, isDirectional(bond.stereo) else { return false }
                if bond.a1 == db.a1 && bond.a2 != db.a2 { return true }
                if bond.a2 == db.a1 && bond.a1 != db.a2 { return true }
                return false
            }

            let sideBHasDirection = molecule.bonds.contains { bond in
                guard bond.id != db.id, bond.order == .single, isDirectional(bond.stereo) else { return false }
                if bond.a1 == db.a2 && bond.a2 != db.a1 { return true }
                if bond.a2 == db.a2 && bond.a1 != db.a1 { return true }
                return false
            }

            if sideAHasDirection && sideBHasDirection {
                molecule.bonds[idx].stereo = .either
            }
        }
    }
}

private extension Character {
    var isUppercase: Bool { String(self).uppercased() == String(self) && isLetter }
    var isLowercase: Bool { String(self).lowercased() == String(self) && isLetter }
}
