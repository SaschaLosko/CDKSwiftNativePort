import Foundation
import CoreGraphics

/// Swift counterpart of CDK's InChI-to-structure adapter.
public final class CDKInChIToStructure {
    private let input: String
    private var molecule: Molecule?
    private(set) var status: CDKInChIStatus = .success
    private(set) var message: String = ""

    public init(inchi: String) {
        self.input = inchi

        do {
            let result = try CDKInChIParser.parse(inchi: inchi)
            molecule = result.molecule
            if !result.ignoredLayers.isEmpty {
                status = .warning
                let list = result.ignoredLayers.map(String.init).sorted().joined(separator: ", ")
                message = "Parsed InChI core layers; ignored unsupported layers: \(list)."
            }
            if !result.ignoredTokens.isEmpty {
                status = .warning
                let tokenSummary = result.ignoredTokens.prefix(5).joined(separator: ", ")
                if message.isEmpty {
                    message = "Parsed with partial fidelity; ignored tokens: \(tokenSummary)."
                } else {
                    message += " Ignored tokens: \(tokenSummary)."
                }
            }
        } catch {
            status = .error
            message = (error as? LocalizedError)?.errorDescription ?? error.localizedDescription
        }
    }

    public func getStatus() -> CDKInChIStatus {
        status
    }

    public func getMessage() -> String {
        message
    }

    public func getAtomContainer() throws -> Molecule {
        guard let molecule else {
            throw ChemError.parseFailed(message.isEmpty ? "Could not parse InChI: \(input)" : message)
        }
        return molecule
    }
}

private struct CDKInChIParseResult {
    let molecule: Molecule
    let ignoredLayers: Set<Character>
    let ignoredTokens: [String]
}

private struct InChIEdgeKey: Hashable {
    let a: Int
    let b: Int

    init(_ u: Int, _ v: Int) {
        a = min(u, v)
        b = max(u, v)
    }
}

private struct InChIMobileHydrogenGroup {
    let count: Int
    let candidateAtoms: [Int]
}

private struct InChIHydrogenLayerResult {
    var fixedByAtom: [Int: Int]
    var mobileGroups: [InChIMobileHydrogenGroup]
    var ignoredTokens: [String]
}

private struct InChIIsotopeLayerResult {
    var massNumberByAtom: [Int: Int]
    var isotopicHydrogenByAtom: [Int: Int]
    var ignoredTokens: [String]
}

private enum CDKInChIParser {
    static func parse(inchi: String) throws -> CDKInChIParseResult {
        let trimmed = inchi.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { throw ChemError.emptyInput }
        guard trimmed.hasPrefix("InChI=") else {
            throw ChemError.parseFailed("Input is not an InChI string (missing 'InChI=' prefix).")
        }

        let payload = String(trimmed.dropFirst("InChI=".count))
        let parts = payload.split(separator: "/", omittingEmptySubsequences: false).map(String.init)
        guard parts.count >= 2 else {
            throw ChemError.parseFailed("InChI is missing required version/formula layers.")
        }
        guard parts[0].hasPrefix("1") else {
            throw ChemError.unsupported("Only InChI version 1.x is supported in this Swift CDK port.")
        }

        let formulaLayer = parts[1]
        let heavyElements = parseHeavyAtomSequence(fromFormula: formulaLayer)
        guard !heavyElements.isEmpty else {
            throw ChemError.parseFailed("InChI formula does not contain heavy atoms supported by this parser.")
        }

        var layers: [Character: String] = [:]
        for segment in parts.dropFirst(2) {
            guard let key = segment.first else { continue }
            let content = String(segment.dropFirst())
            if let existing = layers[key], !existing.isEmpty {
                layers[key] = existing + ";" + content
            } else {
                layers[key] = content
            }
        }

        var atoms: [Atom] = []
        atoms.reserveCapacity(heavyElements.count)
        for (idx, element) in heavyElements.enumerated() {
            atoms.append(Atom(id: idx + 1, element: element, position: .zero))
        }

        let atomCount = atoms.count
        let connectivityLayer = layers["c"] ?? ""
        let edges = try parseConnectivityLayer(connectivityLayer, atomCount: atomCount)

        var bonds: [Bond] = []
        bonds.reserveCapacity(edges.count)
        let sortedEdges = edges.sorted { lhs, rhs in
            if lhs.a != rhs.a { return lhs.a < rhs.a }
            return lhs.b < rhs.b
        }
        for (idx, edge) in sortedEdges.enumerated() {
            bonds.append(Bond(id: idx + 1, a1: edge.a, a2: edge.b, order: .single))
        }

        var molecule = Molecule(name: "InChI", atoms: atoms, bonds: bonds)

        var ignoredTokens: [String] = []

        let hydrogenInfo = parseHydrogenLayer(layers["h"], atomCount: atomCount)
        ignoredTokens.append(contentsOf: hydrogenInfo.ignoredTokens)

        let isotopeInfo = parseIsotopeLayer(layers["i"], atomCount: atomCount, atoms: molecule.atoms)
        ignoredTokens.append(contentsOf: isotopeInfo.ignoredTokens)

        applyIsotopeAssignments(isotopeInfo.massNumberByAtom, to: &molecule)

        let chargeInfo = parseChargeLayer(layers["q"])
        if !chargeInfo.ignoredTokens.isEmpty {
            ignoredTokens.append(contentsOf: chargeInfo.ignoredTokens)
        }
        applyChargeDistribution(chargeInfo.componentCharges, to: &molecule)

        let protonInfo = parseProtonLayer(layers["p"])
        if !protonInfo.ignoredTokens.isEmpty {
            ignoredTokens.append(contentsOf: protonInfo.ignoredTokens)
        }

        var hydrogenByAtom = hydrogenInfo.fixedByAtom
        for (atomID, count) in isotopeInfo.isotopicHydrogenByAtom where count > 0 {
            hydrogenByAtom[atomID, default: 0] += count
        }
        applyProtonLayerAdjustment(protonInfo.delta, to: &molecule, hydrogenByAtom: &hydrogenByAtom)
        distributeMobileHydrogens(hydrogenInfo.mobileGroups, in: molecule, hydrogenByAtom: &hydrogenByAtom)

        assignInferredBondOrders(to: &molecule, hydrogenByAtom: hydrogenByAtom)
        let doubleBondStereoInfo = parseDoubleBondStereoLayer(layers["b"], in: &molecule)
        ignoredTokens.append(contentsOf: doubleBondStereoInfo.ignoredTokens)
        applyTetrahedralLayer(layers: layers, to: &molecule)
        molecule = Depiction2DGenerator.generate(for: molecule)
        molecule.assignWedgeHashFromChiralCenters()

        let supportedKeys: Set<Character> = ["c", "h", "t", "m", "s", "q", "i", "p", "b"]
        var ignoredLayers: Set<Character> = []
        for key in layers.keys where !supportedKeys.contains(key) {
            ignoredLayers.insert(key)
        }

        return CDKInChIParseResult(
            molecule: molecule,
            ignoredLayers: ignoredLayers,
            ignoredTokens: Array(Set(ignoredTokens)).sorted()
        )
    }

    private static func parseHeavyAtomSequence(fromFormula formula: String) -> [String] {
        let normalized = formula
            .replacingOccurrences(of: ".", with: "")
            .replacingOccurrences(of: ";", with: "")
            .replacingOccurrences(of: "+", with: "")
            .replacingOccurrences(of: "-", with: "")

        let pattern = "([A-Z][a-z]?)(\\d*)"
        let regex = try? NSRegularExpression(pattern: pattern)
        let nsRange = NSRange(normalized.startIndex..<normalized.endIndex, in: normalized)

        var out: [String] = []
        regex?.enumerateMatches(in: normalized, options: [], range: nsRange) { match, _, _ in
            guard let match,
                  let symbolRange = Range(match.range(at: 1), in: normalized),
                  let countRange = Range(match.range(at: 2), in: normalized) else {
                return
            }

            let symbol = String(normalized[symbolRange])
            let upper = symbol.uppercased()
            if upper == "H" || upper == "D" || upper == "T" {
                return
            }

            let countToken = String(normalized[countRange])
            let count = max(1, Int(countToken) ?? 1)
            for _ in 0..<count {
                out.append(symbol)
            }
        }

        return out
    }

    private static func parseConnectivityLayer(_ layer: String, atomCount: Int) throws -> Set<InChIEdgeKey> {
        guard !layer.isEmpty else { return [] }

        var edges: Set<InChIEdgeKey> = []
        var idx = layer.startIndex
        var current: Int?
        var stack: [Int] = []

        while idx < layer.endIndex {
            let ch = layer[idx]

            if ch.isWhitespace {
                idx = layer.index(after: idx)
                continue
            }

            if ch.isNumber {
                let atomIndex = parseInteger(layer, from: &idx)
                guard (1...atomCount).contains(atomIndex) else {
                    throw ChemError.parseFailed("InChI connectivity references atom \(atomIndex) outside formula range 1...\(atomCount).")
                }
                if let cur = current, cur != atomIndex {
                    edges.insert(InChIEdgeKey(cur, atomIndex))
                }
                current = atomIndex
                continue
            }

            switch ch {
            case "-", ".":
                break
            case "(":
                guard let cur = current else {
                    throw ChemError.parseFailed("Malformed InChI connectivity branch near '('.")
                }
                stack.append(cur)
            case ")":
                guard let branchRoot = stack.popLast() else {
                    throw ChemError.parseFailed("Unbalanced InChI connectivity parentheses.")
                }
                current = branchRoot
            case ",":
                current = stack.last
            case ";":
                stack.removeAll()
                current = nil
            case "*":
                // InChI repeated-connectivity marker. We parse the local path and ignore multiplicity
                // if present to keep structure extraction robust.
                break
            default:
                if !ch.isLetter {
                    throw ChemError.parseFailed("Unsupported token '\(ch)' in InChI connectivity layer.")
                }
            }

            idx = layer.index(after: idx)
        }

        return edges
    }

    private static func parseHydrogenLayer(_ layer: String?, atomCount: Int) -> InChIHydrogenLayerResult {
        guard let layer, !layer.isEmpty else {
            return InChIHydrogenLayerResult(fixedByAtom: [:], mobileGroups: [], ignoredTokens: [])
        }

        var fixedByAtom: [Int: Int] = [:]
        var mobileGroups: [InChIMobileHydrogenGroup] = []
        var ignoredTokens: [String] = []

        let tokens = splitTopLevel(layer, separators: [",", ";"])
        for rawToken in tokens {
            let token = rawToken.trimmingCharacters(in: .whitespaces)
            guard !token.isEmpty else { continue }

            if token.hasPrefix("(") && token.hasSuffix(")") {
                if let group = parseMobileHydrogenGroup(token, atomCount: atomCount) {
                    mobileGroups.append(group)
                } else {
                    ignoredTokens.append(token)
                }
                continue
            }

            if !parseFixedHydrogenToken(token, atomCount: atomCount, into: &fixedByAtom) {
                ignoredTokens.append(token)
            }
        }

        return InChIHydrogenLayerResult(
            fixedByAtom: fixedByAtom,
            mobileGroups: mobileGroups,
            ignoredTokens: ignoredTokens
        )
    }

    private static func parseFixedHydrogenToken(_ token: String,
                                                atomCount: Int,
                                                into hydrogenByAtom: inout [Int: Int]) -> Bool {
        guard let hIndex = token.firstIndex(where: { $0 == "H" || $0 == "D" || $0 == "T" }) else {
            // InChI allows "h1,2H2" where "1" implies one hydrogen.
            let atoms = parseAtomSpec(token, atomCount: atomCount)
            guard !atoms.isEmpty else { return false }
            for atomID in atoms {
                hydrogenByAtom[atomID, default: 0] += 1
            }
            return true
        }

        let atomSpec = String(token[..<hIndex])
        guard !atomSpec.isEmpty else { return false }

        let symbol = token[hIndex]
        let remainder = String(token[token.index(after: hIndex)...])
        let count = parseLeadingInt(remainder) ?? 1
        guard count > 0 else { return false }

        let atoms = parseAtomSpec(atomSpec, atomCount: atomCount)
        guard !atoms.isEmpty else { return false }

        // D/T are still hydrogens for valence and depiction.
        _ = symbol
        for atomID in atoms {
            hydrogenByAtom[atomID, default: 0] += count
        }
        return true
    }

    private static func parseMobileHydrogenGroup(_ token: String, atomCount: Int) -> InChIMobileHydrogenGroup? {
        guard token.count >= 2 else { return nil }
        let inner = String(token.dropFirst().dropLast())
        let pieces = inner.split(separator: ",").map {
            $0.trimmingCharacters(in: .whitespaces)
        }
        guard pieces.count >= 2 else { return nil }

        let first = pieces[0]
        guard let hChar = first.first, hChar == "H" || hChar == "D" || hChar == "T" else {
            return nil
        }
        let count = parseLeadingInt(String(first.dropFirst())) ?? 1
        guard count > 0 else { return nil }

        let atomSpec = pieces.dropFirst().joined(separator: ",")
        let candidates = parseAtomSpec(atomSpec, atomCount: atomCount)
        guard !candidates.isEmpty else { return nil }

        return InChIMobileHydrogenGroup(count: count, candidateAtoms: candidates)
    }

    private static func parseIsotopeLayer(_ layer: String?,
                                          atomCount: Int,
                                          atoms: [Atom]) -> InChIIsotopeLayerResult {
        guard let layer, !layer.isEmpty else {
            return InChIIsotopeLayerResult(massNumberByAtom: [:], isotopicHydrogenByAtom: [:], ignoredTokens: [])
        }

        let atomByID = Dictionary(uniqueKeysWithValues: atoms.map { ($0.id, $0) })

        var massNumberByAtom: [Int: Int] = [:]
        var isotopicHydrogenByAtom: [Int: Int] = [:]
        var ignoredTokens: [String] = []

        let tokens = splitTopLevel(layer, separators: [",", ";"])
        for rawToken in tokens {
            let token = rawToken.trimmingCharacters(in: .whitespaces)
            guard !token.isEmpty else { continue }

            if let parsed = parseMassShiftToken(token, atomCount: atomCount) {
                for (atomID, shift) in parsed {
                    guard let atom = atomByID[atomID], let base = baseMassNumber(for: atom.element) else { continue }
                    massNumberByAtom[atomID] = max(1, base + shift)
                }
                continue
            }

            if let parsed = parseIsotopicHydrogenToken(token, atomCount: atomCount) {
                for atomID in parsed.atomIDs {
                    isotopicHydrogenByAtom[atomID, default: 0] += parsed.count
                }
                continue
            }

            ignoredTokens.append(token)
        }

        return InChIIsotopeLayerResult(
            massNumberByAtom: massNumberByAtom,
            isotopicHydrogenByAtom: isotopicHydrogenByAtom,
            ignoredTokens: ignoredTokens
        )
    }

    private static func parseMassShiftToken(_ token: String, atomCount: Int) -> [(Int, Int)]? {
        guard let signPos = token.dropFirst().firstIndex(where: { $0 == "+" || $0 == "-" }) else {
            return nil
        }

        let atomSpec = String(token[..<signPos])
        let shiftPart = String(token[signPos...])
        guard let shift = Int(shiftPart) else { return nil }

        let atoms = parseAtomSpec(atomSpec, atomCount: atomCount)
        guard !atoms.isEmpty else { return nil }

        return atoms.map { ($0, shift) }
    }

    private static func parseIsotopicHydrogenToken(_ token: String,
                                                   atomCount: Int) -> (atomIDs: [Int], count: Int)? {
        guard let marker = token.firstIndex(where: { $0 == "D" || $0 == "T" }) else {
            return nil
        }

        let atomSpec = String(token[..<marker])
        guard !atomSpec.isEmpty else { return nil }

        let remainder = String(token[token.index(after: marker)...])
        let count = parseLeadingInt(remainder) ?? 1
        guard count > 0 else { return nil }

        let atoms = parseAtomSpec(atomSpec, atomCount: atomCount)
        guard !atoms.isEmpty else { return nil }
        return (atoms, count)
    }

    private static func baseMassNumber(for element: String) -> Int? {
        switch element.uppercased() {
        case "H": return 1
        case "B": return 11
        case "C": return 12
        case "N": return 14
        case "O": return 16
        case "F": return 19
        case "P": return 31
        case "S": return 32
        case "CL": return 35
        case "BR": return 79
        case "I": return 127
        default: return nil
        }
    }

    private static func applyIsotopeAssignments(_ massByAtom: [Int: Int], to molecule: inout Molecule) {
        guard !massByAtom.isEmpty else { return }
        for idx in molecule.atoms.indices {
            let atomID = molecule.atoms[idx].id
            molecule.atoms[idx].isotopeMassNumber = massByAtom[atomID]
        }
    }

    private struct InChIChargeLayerResult {
        var componentCharges: [Int]
        var ignoredTokens: [String]
    }

    private struct InChIProtonLayerResult {
        var delta: Int
        var ignoredTokens: [String]
    }

    private struct InChIDoubleBondLayerResult {
        var ignoredTokens: [String]
    }

    private static func parseChargeLayer(_ layer: String?) -> InChIChargeLayerResult {
        guard let layer, !layer.isEmpty else {
            return InChIChargeLayerResult(componentCharges: [], ignoredTokens: [])
        }

        var charges: [Int] = []
        var ignored: [String] = []

        let components = splitTopLevel(layer, separators: [";", ","])
        for rawToken in components {
            let token = rawToken.trimmingCharacters(in: .whitespaces)
            guard !token.isEmpty else { continue }

            if let values = parseComponentIntegerToken(token), !values.isEmpty {
                charges.append(contentsOf: values)
                continue
            }
            ignored.append(token)
        }

        return InChIChargeLayerResult(componentCharges: charges, ignoredTokens: ignored)
    }

    private static func parseProtonLayer(_ layer: String?) -> InChIProtonLayerResult {
        guard let layer, !layer.isEmpty else {
            return InChIProtonLayerResult(delta: 0, ignoredTokens: [])
        }

        var delta = 0
        var ignored: [String] = []
        let tokens = splitTopLevel(layer, separators: [";", ","])
        for raw in tokens {
            let token = raw.trimmingCharacters(in: .whitespaces)
            guard !token.isEmpty else { continue }

            if let values = parseComponentIntegerToken(token), !values.isEmpty {
                delta += values.reduce(0, +)
            } else if let scalar = parseSignedInteger(token) {
                delta += scalar
            } else {
                ignored.append(token)
            }
        }

        return InChIProtonLayerResult(delta: delta, ignoredTokens: ignored)
    }

    private static func parseDoubleBondStereoLayer(_ layer: String?,
                                                   in molecule: inout Molecule) -> InChIDoubleBondLayerResult {
        guard let layer, !layer.isEmpty else {
            return InChIDoubleBondLayerResult(ignoredTokens: [])
        }

        var ignored: [String] = []
        let tokens = splitTopLevel(layer, separators: [";", ","])
        for rawToken in tokens {
            let token = rawToken.trimmingCharacters(in: .whitespaces)
            guard !token.isEmpty else { continue }

            guard let pair = parseLeadingAtomPair(in: token) else {
                ignored.append(token)
                continue
            }
            guard let idx = molecule.bonds.firstIndex(where: { bond in
                ((bond.a1 == pair.0 && bond.a2 == pair.1) || (bond.a1 == pair.1 && bond.a2 == pair.0))
                    && bond.order == .double
            }) else {
                ignored.append(token)
                continue
            }
            molecule.bonds[idx].stereo = .either
        }

        return InChIDoubleBondLayerResult(ignoredTokens: ignored)
    }

    private static func parseLeadingAtomPair(in token: String) -> (Int, Int)? {
        var values: [Int] = []
        var current = ""
        for ch in token {
            if ch.isNumber {
                current.append(ch)
                continue
            }
            if !current.isEmpty {
                if let value = Int(current) {
                    values.append(value)
                    if values.count >= 2 { break }
                }
                current = ""
            }
        }
        if values.count < 2, !current.isEmpty, let value = Int(current) {
            values.append(value)
        }
        guard values.count >= 2 else { return nil }
        return (values[0], values[1])
    }

    private static func parseComponentIntegerToken(_ token: String) -> [Int]? {
        var clean = token.trimmingCharacters(in: .whitespacesAndNewlines)
        if clean.hasPrefix("("), clean.hasSuffix(")"), clean.count >= 2 {
            clean.removeFirst()
            clean.removeLast()
            clean = clean.trimmingCharacters(in: .whitespacesAndNewlines)
        }
        guard !clean.isEmpty else { return nil }

        if let value = parseSignedInteger(clean) {
            return [value]
        }

        if let star = clean.firstIndex(of: "*") {
            let lhs = String(clean[..<star]).trimmingCharacters(in: .whitespaces)
            let rhs = String(clean[clean.index(after: star)...]).trimmingCharacters(in: .whitespaces)
            if let count = Int(lhs), count > 0, let value = parseSignedInteger(rhs) {
                return Array(repeating: value, count: count)
            }
        }

        return nil
    }

    private static func applyProtonLayerAdjustment(_ protonDelta: Int,
                                                   to molecule: inout Molecule,
                                                   hydrogenByAtom: inout [Int: Int]) {
        guard protonDelta != 0, !molecule.atoms.isEmpty else { return }

        var atomIndexByID: [Int: Int] = [:]
        for idx in molecule.atoms.indices {
            atomIndexByID[molecule.atoms[idx].id] = idx
        }

        func protonationScore(atomID: Int) -> Int {
            guard let atom = molecule.atom(id: atomID) else { return Int.min }
            let element = atom.element.uppercased()
            let currentH = hydrogenByAtom[atomID] ?? 0
            let bondOrderSum = molecule.bonds(forAtom: atomID).reduce(0) { partial, bond in
                partial + bondOrderValue(bond.order)
            }
            let capacity = max(0, preferredValence(for: atom) - (bondOrderSum + currentH))

            let base: Int
            switch element {
            case "O": base = 108
            case "N": base = 96
            case "S": base = 88
            case "P": base = 82
            case "C": base = 52
            default: base = 36
            }
            let chargeBonus = atom.charge < 0 ? 24 : (atom.charge > 0 ? -10 : 0)
            return base + chargeBonus + capacity * 16 - currentH * 6
        }

        func deprotonationScore(atomID: Int) -> Int {
            guard let atom = molecule.atom(id: atomID) else { return Int.min }
            let currentH = hydrogenByAtom[atomID] ?? 0
            guard currentH > 0 else { return Int.min }

            let element = atom.element.uppercased()
            let base: Int
            switch element {
            case "O": base = 112
            case "S": base = 96
            case "N": base = 84
            case "C": base = 48
            default: base = 30
            }
            let chargeBonus = atom.charge > 0 ? 20 : 0
            return base + chargeBonus + currentH * 8
        }

        let atomIDs = molecule.atoms.map(\.id)
        if protonDelta > 0 {
            for _ in 0..<protonDelta {
                guard let atomID = atomIDs.max(by: { protonationScore(atomID: $0) < protonationScore(atomID: $1) }),
                      let atomIndex = atomIndexByID[atomID] else { break }
                hydrogenByAtom[atomID, default: 0] += 1
                molecule.atoms[atomIndex].charge += 1
            }
        } else {
            for _ in 0..<(-protonDelta) {
                guard let atomID = atomIDs.max(by: { deprotonationScore(atomID: $0) < deprotonationScore(atomID: $1) }),
                      let atomIndex = atomIndexByID[atomID],
                      (hydrogenByAtom[atomID] ?? 0) > 0 else { break }
                hydrogenByAtom[atomID, default: 0] = max(0, (hydrogenByAtom[atomID] ?? 0) - 1)
                molecule.atoms[atomIndex].charge -= 1
            }
        }
    }

    private static func applyChargeDistribution(_ componentCharges: [Int], to molecule: inout Molecule) {
        guard !componentCharges.isEmpty else { return }
        guard !molecule.atoms.isEmpty else { return }

        let components = connectedComponents(in: molecule)
        if componentCharges.count == components.count {
            for (idx, component) in components.enumerated() {
                distributeCharge(componentCharges[idx], in: component, molecule: &molecule)
            }
            return
        }

        let total = componentCharges.reduce(0, +)
        let allAtoms = Set(molecule.atoms.map(\.id))
        distributeCharge(total, in: allAtoms, molecule: &molecule)
    }

    private static func connectedComponents(in molecule: Molecule) -> [Set<Int>] {
        let ids = molecule.atoms.map(\.id).sorted()
        var seen: Set<Int> = []
        var components: [Set<Int>] = []

        for seed in ids where !seen.contains(seed) {
            var stack: [Int] = [seed]
            var comp: Set<Int> = [seed]
            seen.insert(seed)

            while let current = stack.popLast() {
                for nxt in molecule.neighbors(of: current) where !seen.contains(nxt) {
                    seen.insert(nxt)
                    comp.insert(nxt)
                    stack.append(nxt)
                }
            }

            components.append(comp)
        }

        return components.sorted { (lhs, rhs) in
            (lhs.min() ?? 0) < (rhs.min() ?? 0)
        }
    }

    private static func distributeCharge(_ totalCharge: Int,
                                         in component: Set<Int>,
                                         molecule: inout Molecule) {
        guard totalCharge != 0 else { return }
        guard !component.isEmpty else { return }

        var atomIndexByID: [Int: Int] = [:]
        for idx in molecule.atoms.indices {
            atomIndexByID[molecule.atoms[idx].id] = idx
        }

        let degreeByAtom = Dictionary(uniqueKeysWithValues: component.map { ($0, molecule.neighbors(of: $0).count) })

        func chooseCandidate(positive: Bool) -> Int? {
            let candidates = component.compactMap { atomID -> (Int, Int)? in
                guard let idx = atomIndexByID[atomID] else { return nil }
                let atom = molecule.atoms[idx]
                let degree = degreeByAtom[atomID] ?? 0
                let score = chargePlacementScore(atom: atom, degree: degree, positive: positive)
                return (atomID, score)
            }
            return candidates.max { lhs, rhs in
                if lhs.1 != rhs.1 { return lhs.1 < rhs.1 }
                return lhs.0 > rhs.0
            }?.0
        }

        if totalCharge > 0 {
            for _ in 0..<totalCharge {
                guard let atomID = chooseCandidate(positive: true),
                      let idx = atomIndexByID[atomID] else { break }
                molecule.atoms[idx].charge += 1
            }
        } else {
            for _ in 0..<(-totalCharge) {
                guard let atomID = chooseCandidate(positive: false),
                      let idx = atomIndexByID[atomID] else { break }
                molecule.atoms[idx].charge -= 1
            }
        }
    }

    private static func chargePlacementScore(atom: Atom, degree: Int, positive: Bool) -> Int {
        let element = atom.element.uppercased()
        let terminalBonus = degree <= 1 ? 8 : 0
        let chargePenalty = abs(atom.charge) * 28

        let base: Int
        if positive {
            switch element {
            case "N": base = 96
            case "P": base = 92
            case "S": base = 84
            case "O": base = 72
            case "B": base = 60
            case "C": base = 46
            case "F", "CL", "BR", "I": base = 30
            default: base = 24
            }
            let neutralizeBonus = atom.charge < 0 ? 18 : 0
            return base + terminalBonus + neutralizeBonus - chargePenalty
        }

        switch element {
        case "O": base = 100
        case "S": base = 88
        case "N": base = 78
        case "F", "CL", "BR", "I": base = 72
        case "B": base = 58
        case "C": base = 44
        case "P": base = 36
        default: base = 20
        }
        let neutralizeBonus = atom.charge > 0 ? 18 : 0
        return base + terminalBonus + neutralizeBonus - chargePenalty
    }

    private static func distributeMobileHydrogens(_ groups: [InChIMobileHydrogenGroup],
                                                  in molecule: Molecule,
                                                  hydrogenByAtom: inout [Int: Int]) {
        guard !groups.isEmpty else { return }

        for group in groups {
            let candidates = group.candidateAtoms.sorted()
            guard !candidates.isEmpty else { continue }

            for _ in 0..<group.count {
                let chosen = candidates.max { lhs, rhs in
                    let l = mobileHydrogenScore(atomID: lhs, molecule: molecule, hydrogenByAtom: hydrogenByAtom)
                    let r = mobileHydrogenScore(atomID: rhs, molecule: molecule, hydrogenByAtom: hydrogenByAtom)
                    if l != r { return l < r }
                    return lhs > rhs
                }
                guard let atomID = chosen else { continue }
                hydrogenByAtom[atomID, default: 0] += 1
            }
        }
    }

    private static func mobileHydrogenScore(atomID: Int,
                                            molecule: Molecule,
                                            hydrogenByAtom: [Int: Int]) -> Int {
        guard let atom = molecule.atom(id: atomID) else { return Int.min }

        let element = atom.element.uppercased()
        let currentHydrogenCount = hydrogenByAtom[atomID] ?? 0
        let bondOrderSum = molecule.bonds(forAtom: atomID)
            .reduce(0) { partial, bond in
                partial + bondOrderValue(bond.order)
            }

        let target = preferredValence(for: atom)
        let need = max(0, target - (bondOrderSum + currentHydrogenCount))

        let elementScore: Int = {
            switch element {
            case "O": return 96
            case "N": return 90
            case "S": return 82
            case "P": return 70
            case "C": return 44
            default: return 32
            }
        }()

        let chargeBonus = atom.charge < 0 ? 20 : (atom.charge > 0 ? -20 : 0)
        return elementScore + need * 18 + chargeBonus - currentHydrogenCount * 6
    }

    private static func applyTetrahedralLayer(layers: [Character: String], to molecule: inout Molecule) {
        guard let tLayer = layers["t"], !tLayer.isEmpty else { return }

        var chiralityByAtom: [Int: AtomChirality] = [:]
        for token in splitTopLevel(tLayer, separators: [",", ";"]) {
            let piece = token.trimmingCharacters(in: .whitespaces)
            guard !piece.isEmpty else { continue }

            let atomDigits = piece.prefix { $0.isNumber }
            guard !atomDigits.isEmpty, let atomID = Int(atomDigits) else { continue }

            let suffix = piece.dropFirst(atomDigits.count)
            if suffix.contains("+") {
                chiralityByAtom[atomID] = .clockwise
            } else if suffix.contains("-") {
                chiralityByAtom[atomID] = .anticlockwise
            }
        }

        let shouldInvert = shouldInvertTetrahedralParity(mLayer: layers["m"])
        for i in molecule.atoms.indices {
            let id = molecule.atoms[i].id
            guard var chirality = chiralityByAtom[id] else { continue }

            if shouldInvert {
                switch chirality {
                case .clockwise:
                    chirality = .anticlockwise
                case .anticlockwise:
                    chirality = .clockwise
                case .none:
                    break
                }
            }

            molecule.atoms[i].chirality = chirality
        }
    }

    private static func shouldInvertTetrahedralParity(mLayer: String?) -> Bool {
        guard let mLayer, !mLayer.isEmpty else { return false }
        let firstToken = splitTopLevel(mLayer, separators: [",", ";"]).first ?? mLayer
        let value = firstToken.trimmingCharacters(in: .whitespaces)
        return value == "0"
    }

    private static func assignInferredBondOrders(to molecule: inout Molecule, hydrogenByAtom: [Int: Int]) {
        guard !molecule.bonds.isEmpty else { return }

        let degreeByAtom = atomDegrees(in: molecule)
        let ringEdges = ringEdges(in: molecule)

        var deficits = valenceDeficits(in: molecule, hydrogenByAtom: hydrogenByAtom)

        while let idx = bestTripleCandidate(in: molecule,
                                            deficits: deficits,
                                            degreeByAtom: degreeByAtom,
                                            ringEdges: ringEdges) {
            let a = molecule.bonds[idx].a1
            let b = molecule.bonds[idx].a2
            molecule.bonds[idx].order = .triple
            deficits[a, default: 0] = max(0, deficits[a, default: 0] - 2)
            deficits[b, default: 0] = max(0, deficits[b, default: 0] - 2)
        }

        while let idx = bestDoubleCandidate(in: molecule,
                                            deficits: deficits,
                                            degreeByAtom: degreeByAtom,
                                            ringEdges: ringEdges) {
            let a = molecule.bonds[idx].a1
            let b = molecule.bonds[idx].a2
            molecule.bonds[idx].order = .double
            deficits[a, default: 0] = max(0, deficits[a, default: 0] - 1)
            deficits[b, default: 0] = max(0, deficits[b, default: 0] - 1)
        }
    }

    private static func atomDegrees(in molecule: Molecule) -> [Int: Int] {
        var degree: [Int: Int] = [:]
        for bond in molecule.bonds {
            degree[bond.a1, default: 0] += 1
            degree[bond.a2, default: 0] += 1
        }
        return degree
    }

    private static func ringEdges(in molecule: Molecule) -> Set<InChIEdgeKey> {
        let rings = molecule.simpleCycles(maxSize: 12)
        var edges: Set<InChIEdgeKey> = []
        for ring in rings where ring.count >= 3 {
            for i in 0..<ring.count {
                edges.insert(InChIEdgeKey(ring[i], ring[(i + 1) % ring.count]))
            }
        }
        return edges
    }

    private static func valenceDeficits(in molecule: Molecule, hydrogenByAtom: [Int: Int]) -> [Int: Int] {
        var out: [Int: Int] = [:]
        for atom in molecule.atoms {
            let target = preferredValence(for: atom)
            guard target > 0 else {
                out[atom.id] = 0
                continue
            }

            let usedFromBonds = molecule.bonds(forAtom: atom.id).reduce(0) { partial, bond in
                partial + bondOrderValue(bond.order)
            }
            let used = usedFromBonds + (hydrogenByAtom[atom.id] ?? 0)
            out[atom.id] = max(0, target - used)
        }
        return out
    }

    private static func preferredValence(for atom: Atom) -> Int {
        switch atom.element.uppercased() {
        case "C":
            return atom.charge == 0 ? 4 : 3
        case "N":
            return atom.charge > 0 ? 4 : (atom.charge < 0 ? 2 : 3)
        case "O":
            return atom.charge > 0 ? 3 : (atom.charge < 0 ? 1 : 2)
        case "S":
            return atom.charge > 0 ? 3 : (atom.charge < 0 ? 1 : 2)
        case "P":
            return atom.charge > 0 ? 4 : 3
        case "B":
            return 3
        case "F", "CL", "BR", "I":
            return 1
        default:
            return 0
        }
    }

    private static func bondOrderValue(_ order: BondOrder) -> Int {
        switch order {
        case .single:
            return 1
        case .double:
            return 2
        case .triple:
            return 3
        case .aromatic:
            return 2
        }
    }

    private static func bestTripleCandidate(in molecule: Molecule,
                                            deficits: [Int: Int],
                                            degreeByAtom: [Int: Int],
                                            ringEdges: Set<InChIEdgeKey>) -> Int? {
        let elementByAtom = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0.element.uppercased()) })

        var best: (idx: Int, score: Int)?
        for idx in molecule.bonds.indices {
            let bond = molecule.bonds[idx]
            guard bond.order == .single else { continue }

            let a = bond.a1
            let b = bond.a2
            guard deficits[a, default: 0] >= 2, deficits[b, default: 0] >= 2 else { continue }

            let degA = degreeByAtom[a, default: 0]
            let degB = degreeByAtom[b, default: 0]
            let eA = elementByAtom[a] ?? ""
            let eB = elementByAtom[b] ?? ""

            var score = 0
            if ringEdges.contains(InChIEdgeKey(a, b)) {
                score -= 8
            }
            if degA == 1 || degB == 1 {
                score += 8
            }
            if (eA == "C" && eB == "N") || (eA == "N" && eB == "C") {
                score += 4
            }
            if eA == "C" && eB == "C" {
                score += 2
            }

            if score <= 0 { continue }
            if best == nil || score > best!.score {
                best = (idx, score)
            }
        }

        return best?.idx
    }

    private static func bestDoubleCandidate(in molecule: Molecule,
                                            deficits: [Int: Int],
                                            degreeByAtom: [Int: Int],
                                            ringEdges: Set<InChIEdgeKey>) -> Int? {
        let elementByAtom = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0.element.uppercased()) })

        var best: (idx: Int, score: Int)?
        for idx in molecule.bonds.indices {
            let bond = molecule.bonds[idx]
            guard bond.order == .single else { continue }

            let a = bond.a1
            let b = bond.a2
            guard deficits[a, default: 0] >= 1, deficits[b, default: 0] >= 1 else { continue }

            let key = InChIEdgeKey(a, b)
            let degA = degreeByAtom[a, default: 0]
            let degB = degreeByAtom[b, default: 0]
            let eA = elementByAtom[a] ?? ""
            let eB = elementByAtom[b] ?? ""

            var score = 0
            if ringEdges.contains(key) {
                score += 7
            }

            if (eA == "C" && (eB == "O" || eB == "S" || eB == "N") && degB == 1)
                || (eB == "C" && (eA == "O" || eA == "S" || eA == "N") && degA == 1) {
                score += 10
            }

            if eA == "C" && eB == "C" && degA == 2 && degB == 2 {
                score += 5
            }

            if degA == 1 && degB == 1 {
                score -= 10
            }

            if (eA == "O" && degA > 1) || (eB == "O" && degB > 1) {
                score -= 3
            }

            score += min(2, deficits[a, default: 0]) + min(2, deficits[b, default: 0])
            if score <= 0 { continue }

            if best == nil || score > best!.score {
                best = (idx, score)
            }
        }

        return best?.idx
    }

    private static func splitTopLevel(_ text: String, separators: Set<Character>) -> [String] {
        var out: [String] = []
        var level = 0
        var chunk = ""

        for ch in text {
            if ch == "(" {
                level += 1
                chunk.append(ch)
                continue
            }
            if ch == ")" {
                level = max(0, level - 1)
                chunk.append(ch)
                continue
            }
            if level == 0 && separators.contains(ch) {
                out.append(chunk)
                chunk = ""
                continue
            }
            chunk.append(ch)
        }

        if !chunk.isEmpty {
            out.append(chunk)
        }
        return out
    }

    private static func parseAtomSpec(_ atomSpec: String, atomCount: Int) -> [Int] {
        var clean = atomSpec.trimmingCharacters(in: .whitespaces)
        while clean.hasPrefix("("), clean.hasSuffix(")"), clean.count >= 2 {
            clean.removeFirst()
            clean.removeLast()
            clean = clean.trimmingCharacters(in: .whitespaces)
        }
        guard !clean.isEmpty else { return [] }

        var out: [Int] = []
        for part in clean.split(separator: ",") {
            var token = part.trimmingCharacters(in: .whitespaces)
            if token.hasSuffix("?") || token.hasSuffix("+") || token.hasSuffix("-") {
                token = String(token.dropLast())
            }
            token = token.trimmingCharacters(in: .whitespaces)
            guard !token.isEmpty else { continue }

            if let dash = token.firstIndex(of: "-") {
                let left = token[..<dash]
                let right = token[token.index(after: dash)...]
                if let start = Int(left), let end = Int(right) {
                    if start <= end {
                        for id in start...end where (1...atomCount).contains(id) {
                            out.append(id)
                        }
                    } else {
                        for id in stride(from: start, through: end, by: -1) where (1...atomCount).contains(id) {
                            out.append(id)
                        }
                    }
                    continue
                }
            }

            if let id = Int(token), (1...atomCount).contains(id) {
                out.append(id)
            }
        }

        return Array(Set(out)).sorted()
    }

    private static func parseLeadingInt(_ token: String) -> Int? {
        let digits = token.prefix { $0.isNumber }
        guard !digits.isEmpty else { return nil }
        return Int(digits)
    }

    private static func parseSignedInteger(_ token: String) -> Int? {
        let clean = token.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !clean.isEmpty else { return nil }

        if let value = Int(clean) {
            return value
        }

        guard let first = clean.first, (first == "+" || first == "-") else { return nil }
        let magnitudeToken = String(clean.dropFirst())
        guard let magnitude = Int(magnitudeToken) else { return nil }
        return first == "-" ? -magnitude : magnitude
    }

    private static func parseInteger(_ text: String, from index: inout String.Index) -> Int {
        var number = 0
        var idx = index
        while idx < text.endIndex, text[idx].isNumber {
            number = number * 10 + Int(String(text[idx]))!
            idx = text.index(after: idx)
        }
        index = idx
        return number
    }
}
