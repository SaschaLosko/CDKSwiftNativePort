import Foundation
#if canImport(CryptoKit)
import CryptoKit
#endif

struct CDKInChINativeGenerationResult {
    let inchi: String
    let inchiKey: String
    let status: CDKInChIStatus
    let message: String
}

private struct InChIEdgeKey: Hashable {
    let a: Int
    let b: Int

    init(_ u: Int, _ v: Int) {
        a = min(u, v)
        b = max(u, v)
    }
}

/// Pure Swift InChI/InChIKey generator used by the CDKSwiftNativePort package.
///
/// This follows CDK-style layering conventions (`/c`, `/h`, `/q`, `/i`, `/t`, `/m`, `/s`, `/b`)
/// but remains an approximation, not the official IUPAC reference implementation.
enum CDKInChINativeGenerator {
    static func generate(for molecule: Molecule) throws -> CDKInChINativeGenerationResult {
        let normalized = normalizeInput(molecule)
        guard !normalized.atoms.isEmpty else {
            throw ChemError.emptyInput
        }

        if normalized.atoms.contains(where: { $0.queryType != nil || $0.atomList != nil }) {
            throw ChemError.unsupported("InChI generation does not support query atoms in this CDK Swift port.")
        }

        let heavyAtomIDs = normalized.atoms
            .filter { !isHydrogenSymbol($0.element) }
            .map(\.id)
            .sorted()

        guard !heavyAtomIDs.isEmpty else {
            throw ChemError.unsupported("InChI generation requires at least one non-hydrogen atom.")
        }

        let hydrogenByHeavyAtom = hydrogenCountByHeavyAtom(in: normalized, heavyAtomIDs: heavyAtomIDs)
        let canonicalization = canonicalizeHeavyAtoms(in: normalized,
                                                      heavyAtomIDs: heavyAtomIDs,
                                                      hydrogenByHeavyAtom: hydrogenByHeavyAtom)

        let formula = try buildFormula(in: normalized,
                                       heavyAtomIDs: heavyAtomIDs,
                                       hydrogenByHeavyAtom: hydrogenByHeavyAtom)
        let connectivityLayer = buildConnectivityLayer(in: normalized, canonicalization: canonicalization)
        let hydrogenLayer = buildHydrogenLayer(canonicalization: canonicalization,
                                               hydrogenByHeavyAtom: hydrogenByHeavyAtom)
        let chargeLayer = buildChargeLayer(in: normalized)
        let isotopeLayer = buildIsotopeLayer(in: normalized, canonicalization: canonicalization)
        let doubleBondLayer = buildDoubleBondLayer(in: normalized, canonicalization: canonicalization)
        let tetrahedralLayer = buildTetrahedralLayer(in: normalized, canonicalization: canonicalization)

        var segments: [String] = [formula]
        if !connectivityLayer.isEmpty { segments.append("c\(connectivityLayer)") }
        if !hydrogenLayer.isEmpty { segments.append("h\(hydrogenLayer)") }
        if !chargeLayer.isEmpty { segments.append("q\(chargeLayer)") }
        if !isotopeLayer.isEmpty { segments.append("i\(isotopeLayer)") }
        if !doubleBondLayer.isEmpty { segments.append("b\(doubleBondLayer)") }
        if !tetrahedralLayer.isEmpty {
            segments.append("t\(tetrahedralLayer)")
            segments.append("m1")
            segments.append("s1")
        }

        let inchi = "InChI=1S/" + segments.joined(separator: "/")
        let inchiKey = makePseudoInchiKey(from: inchi)
        return CDKInChINativeGenerationResult(inchi: inchi,
                                              inchiKey: inchiKey,
                                              status: .success,
                                              message: "")
    }

    private struct CanonicalizationResult {
        let heavyAtomOrder: [Int]
        let newIDByOldID: [Int: Int]
    }

    private static func canonicalizeHeavyAtoms(in molecule: Molecule,
                                               heavyAtomIDs: [Int],
                                               hydrogenByHeavyAtom: [Int: Int]) -> CanonicalizationResult {
        let atomByID = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0) })
        let bondByEdge = Dictionary(uniqueKeysWithValues: molecule.bonds.map {
            (InChIEdgeKey($0.a1, $0.a2), $0)
        })

        var invariants: [Int: UInt64] = [:]
        for atomID in heavyAtomIDs {
            guard let atom = atomByID[atomID] else { continue }
            let neighbors = heavyNeighbors(of: atomID, in: molecule)
            let heavyDegree = neighbors.count
            let valence = molecule.bonds(forAtom: atomID).reduce(0) { partial, bond in
                partial + bondRank(bond.order)
            }

            var state: UInt64 = 0xcbf29ce484222325
            state = mix(state, UInt64(atomicNumber(for: atom.element)))
            state = mix(state, UInt64(bitPattern: Int64(atom.charge + 128)))
            state = mix(state, UInt64(atom.isotopeMassNumber ?? 0))
            state = mix(state, UInt64(heavyDegree))
            state = mix(state, UInt64(valence))
            state = mix(state, UInt64(hydrogenByHeavyAtom[atomID] ?? 0))
            state = mix(state, atom.aromatic ? 1 : 0)
            invariants[atomID] = state
        }

        for _ in 0..<8 {
            var next = invariants
            for atomID in heavyAtomIDs {
                var state = mix(0x9e3779b97f4a7c15, invariants[atomID] ?? 0)
                let neighbors = heavyNeighbors(of: atomID, in: molecule)
                    .sorted { lhs, rhs in
                        let li = invariants[lhs] ?? 0
                        let ri = invariants[rhs] ?? 0
                        if li != ri { return li > ri }
                        return lhs < rhs
                    }

                for neighbor in neighbors {
                    let edge = InChIEdgeKey(atomID, neighbor)
                    let bond = bondByEdge[edge]
                    state = mix(state, invariants[neighbor] ?? 0)
                    state = mix(state, UInt64(bondRank(bond?.order ?? .single)))
                }

                next[atomID] = state
            }
            invariants = next
        }

        let sorted = heavyAtomIDs.sorted { lhs, rhs in
            let li = invariants[lhs] ?? 0
            let ri = invariants[rhs] ?? 0
            if li != ri { return li > ri }

            let lAtom = atomByID[lhs]
            let rAtom = atomByID[rhs]
            let lSymbol = normalizedElementSymbol(lAtom?.element ?? "")
            let rSymbol = normalizedElementSymbol(rAtom?.element ?? "")
            if lSymbol != rSymbol { return lSymbol < rSymbol }

            let lCharge = lAtom?.charge ?? 0
            let rCharge = rAtom?.charge ?? 0
            if lCharge != rCharge { return lCharge > rCharge }

            let lHydrogen = hydrogenByHeavyAtom[lhs] ?? 0
            let rHydrogen = hydrogenByHeavyAtom[rhs] ?? 0
            if lHydrogen != rHydrogen { return lHydrogen > rHydrogen }

            return lhs < rhs
        }

        var map: [Int: Int] = [:]
        for (index, atomID) in sorted.enumerated() {
            map[atomID] = index + 1
        }

        return CanonicalizationResult(heavyAtomOrder: sorted, newIDByOldID: map)
    }

    private static func buildFormula(in molecule: Molecule,
                                     heavyAtomIDs: [Int],
                                     hydrogenByHeavyAtom: [Int: Int]) throws -> String {
        var composition: [String: Int] = [:]
        let atomByID = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0) })

        for atomID in heavyAtomIDs {
            guard let atom = atomByID[atomID] else { continue }
            let symbol = normalizedElementSymbol(atom.element)
            guard !symbol.isEmpty, symbol != "*" else {
                throw ChemError.unsupported("Unsupported atom symbol '\(atom.element)' for InChI generation.")
            }
            composition[symbol, default: 0] += 1
        }

        let attachedHydrogenCount = hydrogenByHeavyAtom.values.reduce(0, +)
        let attachedExplicitHydrogenIDs = Set(heavyAtomIDs.flatMap { heavyID in
            molecule.neighbors(of: heavyID).filter { neighborID in
                guard let atom = molecule.atom(id: neighborID) else { return false }
                return isHydrogenSymbol(atom.element)
            }
        })

        let detachedHydrogenCount = molecule.atoms.filter { atom in
            isHydrogenSymbol(atom.element) && !attachedExplicitHydrogenIDs.contains(atom.id)
        }.count

        let hydrogenTotal = attachedHydrogenCount + detachedHydrogenCount
        if hydrogenTotal > 0 {
            composition["H", default: 0] += hydrogenTotal
        }

        guard !composition.isEmpty else {
            throw ChemError.parseFailed("Failed to build molecular formula for InChI generation.")
        }

        var tokens: [String] = []
        if let carbon = composition.removeValue(forKey: "C"), carbon > 0 {
            tokens.append("C" + (carbon == 1 ? "" : "\(carbon)"))
        }
        if let hydrogen = composition.removeValue(forKey: "H"), hydrogen > 0 {
            tokens.append("H" + (hydrogen == 1 ? "" : "\(hydrogen)"))
        }
        for element in composition.keys.sorted() {
            guard let count = composition[element], count > 0 else { continue }
            tokens.append(element + (count == 1 ? "" : "\(count)"))
        }

        return tokens.joined()
    }

    private static func buildConnectivityLayer(in molecule: Molecule,
                                               canonicalization: CanonicalizationResult) -> String {
        var edgeSet = Set<InChIEdgeKey>()
        for bond in molecule.bonds {
            guard let a = canonicalization.newIDByOldID[bond.a1],
                  let b = canonicalization.newIDByOldID[bond.a2] else {
                continue
            }
            edgeSet.insert(InChIEdgeKey(a, b))
        }

        let tokens = edgeSet.sorted { lhs, rhs in
            if lhs.a != rhs.a { return lhs.a < rhs.a }
            return lhs.b < rhs.b
        }.map { "\($0.a)-\($0.b)" }

        return tokens.joined(separator: ";")
    }

    private static func buildHydrogenLayer(canonicalization: CanonicalizationResult,
                                           hydrogenByHeavyAtom: [Int: Int]) -> String {
        let tokens = canonicalization.heavyAtomOrder.compactMap { oldAtomID -> String? in
            let count = hydrogenByHeavyAtom[oldAtomID] ?? 0
            guard count > 0, let atomID = canonicalization.newIDByOldID[oldAtomID] else { return nil }
            if count == 1 {
                return "\(atomID)H"
            }
            return "\(atomID)H\(count)"
        }
        return tokens.joined(separator: ",")
    }

    private static func buildChargeLayer(in molecule: Molecule) -> String {
        let totalCharge = molecule.atoms.reduce(0) { $0 + $1.charge }
        guard totalCharge != 0 else { return "" }
        return signedInteger(totalCharge)
    }

    private static func buildIsotopeLayer(in molecule: Molecule,
                                          canonicalization: CanonicalizationResult) -> String {
        let atomByID = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0) })
        var tokens: [String] = []

        for oldAtomID in canonicalization.heavyAtomOrder {
            guard let mappedID = canonicalization.newIDByOldID[oldAtomID],
                  let atom = atomByID[oldAtomID],
                  let isotope = atom.isotopeMassNumber,
                  let base = baseMassNumber(for: atom.element) else {
                continue
            }
            let shift = isotope - base
            guard shift != 0 else { continue }
            tokens.append("\(mappedID)\(signedInteger(shift))")
        }

        for oldAtomID in canonicalization.heavyAtomOrder {
            guard let mappedID = canonicalization.newIDByOldID[oldAtomID] else { continue }
            var dCount = 0
            var tCount = 0
            for neighborID in molecule.neighbors(of: oldAtomID) {
                guard let neighbor = atomByID[neighborID], isHydrogenSymbol(neighbor.element) else { continue }
                let isotope = neighbor.isotopeMassNumber ?? 1
                if isotope == 2 {
                    dCount += 1
                } else if isotope == 3 {
                    tCount += 1
                }
            }

            if dCount > 0 {
                tokens.append("\(mappedID)D" + (dCount == 1 ? "" : "\(dCount)"))
            }
            if tCount > 0 {
                tokens.append("\(mappedID)T" + (tCount == 1 ? "" : "\(tCount)"))
            }
        }

        return tokens.joined(separator: ",")
    }

    private static func buildDoubleBondLayer(in molecule: Molecule,
                                             canonicalization: CanonicalizationResult) -> String {
        var pairs: [(Int, Int)] = []
        for bond in molecule.bonds where bond.order == .double && bond.stereo != .none {
            guard let a = canonicalization.newIDByOldID[bond.a1],
                  let b = canonicalization.newIDByOldID[bond.a2] else {
                continue
            }
            pairs.append((min(a, b), max(a, b)))
        }

        pairs.sort { lhs, rhs in
            if lhs.0 != rhs.0 { return lhs.0 < rhs.0 }
            return lhs.1 < rhs.1
        }

        return pairs.map { "\($0.0)-\($0.1)" }.joined(separator: ",")
    }

    private static func buildTetrahedralLayer(in molecule: Molecule,
                                              canonicalization: CanonicalizationResult) -> String {
        let atomByID = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0) })
        let tokens = canonicalization.heavyAtomOrder.compactMap { oldAtomID -> String? in
            guard let mappedID = canonicalization.newIDByOldID[oldAtomID],
                  let atom = atomByID[oldAtomID] else {
                return nil
            }
            switch atom.chirality {
            case .clockwise:
                return "\(mappedID)+"
            case .anticlockwise:
                return "\(mappedID)-"
            case .none:
                return nil
            }
        }
        return tokens.joined(separator: ",")
    }

    private static func hydrogenCountByHeavyAtom(in molecule: Molecule,
                                                 heavyAtomIDs: [Int]) -> [Int: Int] {
        var counts: [Int: Int] = [:]
        for atomID in heavyAtomIDs {
            guard let atom = molecule.atom(id: atomID) else { continue }
            let base: Int
            if let explicit = atom.explicitHydrogenCount {
                base = max(0, explicit)
            } else {
                base = max(0, CDKDescriptorSupport.implicitHydrogenCount(on: atomID, in: molecule))
            }
            let explicitNeighbors = CDKDescriptorSupport.explicitHydrogenNeighborCount(on: atomID, in: molecule)
            counts[atomID] = max(0, base + explicitNeighbors)
        }
        return counts
    }

    private static func heavyNeighbors(of atomID: Int, in molecule: Molecule) -> [Int] {
        molecule.neighbors(of: atomID).filter { neighborID in
            guard let neighbor = molecule.atom(id: neighborID) else { return false }
            return !isHydrogenSymbol(neighbor.element)
        }
    }

    private static func normalizeInput(_ molecule: Molecule) -> Molecule {
        var copy = molecule
        for index in copy.atoms.indices {
            var atom = copy.atoms[index]
            var normalized = normalizedElementSymbol(atom.element)
            if atom.charge == 0 {
                atom.charge += inferredCharge(fromRawElement: atom.element)
            }

            switch normalized.uppercased() {
            case "D":
                normalized = "H"
                if atom.isotopeMassNumber == nil {
                    atom.isotopeMassNumber = 2
                }
            case "T":
                normalized = "H"
                if atom.isotopeMassNumber == nil {
                    atom.isotopeMassNumber = 3
                }
            default:
                break
            }

            atom.element = normalized
            copy.atoms[index] = atom
        }
        return copy
    }

    private static func normalizedElementSymbol(_ raw: String) -> String {
        let trimmed = raw.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { return "" }

        let letters = trimmed.prefix { $0.isLetter }
        guard !letters.isEmpty else { return trimmed }

        let source = String(letters)
        let first = String(source.prefix(1)).uppercased()
        if source.count == 1 {
            return first
        }

        let second = String(source.dropFirst().prefix(1)).lowercased()
        let candidate = first + second
        if knownElementSymbols.contains(candidate.uppercased()) {
            return candidate
        }
        return first
    }

    private static func inferredCharge(fromRawElement raw: String) -> Int {
        let trimmed = raw.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { return 0 }

        let letterPrefixLength = trimmed.prefix { $0.isLetter }.count
        guard letterPrefixLength < trimmed.count else { return 0 }
        let suffix = String(trimmed.dropFirst(letterPrefixLength))
        guard !suffix.isEmpty else { return 0 }

        var idx = suffix.startIndex
        var total = 0
        while idx < suffix.endIndex {
            let signChar = suffix[idx]
            guard signChar == "+" || signChar == "-" else {
                idx = suffix.index(after: idx)
                continue
            }

            let sign = signChar == "+" ? 1 : -1
            idx = suffix.index(after: idx)

            let digitStart = idx
            while idx < suffix.endIndex, suffix[idx].isNumber {
                idx = suffix.index(after: idx)
            }
            if digitStart != idx {
                let digits = String(suffix[digitStart..<idx])
                total += sign * (Int(digits) ?? 1)
                continue
            }

            var repeated = 1
            while idx < suffix.endIndex, suffix[idx] == signChar {
                repeated += 1
                idx = suffix.index(after: idx)
            }
            total += sign * repeated
        }
        return total
    }

    private static func isHydrogenSymbol(_ symbol: String) -> Bool {
        let upper = normalizedElementSymbol(symbol).uppercased()
        return upper == "H" || upper == "D" || upper == "T"
    }

    private static func baseMassNumber(for element: String) -> Int? {
        switch normalizedElementSymbol(element).uppercased() {
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

    private static func atomicNumber(for element: String) -> Int {
        atomicNumberByElement[normalizedElementSymbol(element).uppercased()] ?? 0
    }

    private static func bondRank(_ order: BondOrder) -> Int {
        switch order {
        case .single: return 1
        case .double: return 2
        case .triple: return 3
        case .aromatic: return 4
        }
    }

    private static func mix(_ state: UInt64, _ value: UInt64) -> UInt64 {
        var x = state ^ (value &+ 0x9e3779b97f4a7c15 &+ (state << 6) &+ (state >> 2))
        x ^= x >> 33
        x &*= 0xff51afd7ed558ccd
        x ^= x >> 33
        x &*= 0xc4ceb9fe1a85ec53
        x ^= x >> 33
        return x
    }

    private static func signedInteger(_ value: Int) -> String {
        value >= 0 ? "+\(value)" : "\(value)"
    }

    private static func makePseudoInchiKey(from inchi: String) -> String {
        let digest = digestBytes(for: inchi)
        let firstBlock = letterBlock(from: digest, seedOffset: 0, length: 14)
        let secondCore = letterBlock(from: digest, seedOffset: 14, length: 8)
        let secondBlock = secondCore + "SA"
        let protonation = inchi.contains("/p-") ? "M" : (inchi.contains("/p+") ? "O" : "N")
        return "\(firstBlock)-\(secondBlock)-\(protonation)"
    }

    private static func letterBlock(from digest: [UInt8], seedOffset: Int, length: Int) -> String {
        guard !digest.isEmpty, length > 0 else { return "" }
        var state: UInt64 = 0x6a09e667f3bcc909 ^ UInt64(seedOffset &* 1315423911)
        var chars: [Character] = []
        chars.reserveCapacity(length)

        for index in 0..<length {
            let byte = UInt64(digest[(seedOffset + index) % digest.count])
            state = mix(state, byte &+ UInt64(index + 1))
            let letter = Character(UnicodeScalar(65 + Int(state % 26))!)
            chars.append(letter)
        }

        return String(chars)
    }

    private static func digestBytes(for text: String) -> [UInt8] {
        let bytes = [UInt8](text.utf8)
        #if canImport(CryptoKit)
        let hash = SHA256.hash(data: Data(bytes))
        return Array(hash)
        #else
        // Deterministic fallback for platforms without CryptoKit.
        var state: UInt64 = 0xcbf29ce484222325
        for byte in bytes {
            state ^= UInt64(byte)
            state &*= 0x100000001b3
        }
        var out: [UInt8] = []
        out.reserveCapacity(32)
        var current = state
        for i in 0..<32 {
            current = mix(current, UInt64(i + 1))
            out.append(UInt8(truncatingIfNeeded: current & 0xff))
        }
        return out
        #endif
    }

    private static let knownElementSymbols: Set<String> = [
        "H", "HE", "LI", "BE", "B", "C", "N", "O", "F", "NE",
        "NA", "MG", "AL", "SI", "P", "S", "CL", "AR", "K", "CA",
        "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN",
        "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y", "ZR",
        "NB", "MO", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB",
        "TE", "I", "XE", "CS", "BA", "PT", "AU", "HG", "PB", "BI"
    ]

    private static let atomicNumberByElement: [String: Int] = [
        "H": 1, "HE": 2, "LI": 3, "BE": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "NE": 10,
        "NA": 11, "MG": 12, "AL": 13, "SI": 14, "P": 15, "S": 16, "CL": 17, "AR": 18, "K": 19, "CA": 20,
        "SC": 21, "TI": 22, "V": 23, "CR": 24, "MN": 25, "FE": 26, "CO": 27, "NI": 28, "CU": 29, "ZN": 30,
        "GA": 31, "GE": 32, "AS": 33, "SE": 34, "BR": 35, "KR": 36, "RB": 37, "SR": 38, "Y": 39, "ZR": 40,
        "NB": 41, "MO": 42, "RU": 44, "RH": 45, "PD": 46, "AG": 47, "CD": 48, "IN": 49, "SN": 50, "SB": 51,
        "TE": 52, "I": 53, "XE": 54, "CS": 55, "BA": 56, "PT": 78, "AU": 79, "HG": 80, "PB": 82, "BI": 83
    ]
}
