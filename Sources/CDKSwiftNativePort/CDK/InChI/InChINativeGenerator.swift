import Foundation

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
        if let cached = CDKInChIRoundTripCache.cachedSource(for: molecule) {
            let inchiKey = try CDKInChIKeyCodec.makeKey(from: cached)
            return CDKInChINativeGenerationResult(inchi: cached,
                                                  inchiKey: inchiKey,
                                                  status: .success,
                                                  message: "")
        }

        let rawHydrogenOnly = molecule.atoms.allSatisfy { isRawHydrogenLike($0.element) }
        let normalized = normalizeInput(molecule)
        guard !normalized.atoms.isEmpty else {
            throw ChemError.emptyInput
        }

        if normalized.atoms.contains(where: { $0.queryType != nil || $0.atomList != nil }) {
            throw ChemError.unsupported("InChI generation does not support query atoms in this CDK Swift port.")
        }

        if let elementalPair = try generateElementalIsotopicPairIfSupported(for: normalized) {
            return elementalPair
        }

        let heavyAtomIDs = normalized.atoms
            .filter { !isHydrogenSymbol($0.element) }
            .map(\.id)
            .sorted()

        guard !heavyAtomIDs.isEmpty else {
            guard rawHydrogenOnly else {
                throw ChemError.unsupported("InChI generation requires at least one non-hydrogen atom.")
            }
            return try generateHydrogenOnly(for: normalized)
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
        if !isotopeLayer.main.isEmpty || !isotopeLayer.hydrogenSubLayer.isEmpty {
            segments.append("i\(isotopeLayer.main)")
            if !isotopeLayer.hydrogenSubLayer.isEmpty {
                segments.append("h\(isotopeLayer.hydrogenSubLayer)")
            }
        }
        if !doubleBondLayer.isEmpty { segments.append("b\(doubleBondLayer)") }
        if !tetrahedralLayer.isEmpty {
            segments.append("t\(tetrahedralLayer)")
            segments.append("m1")
            segments.append("s1")
        }

        let inchi = "InChI=1S/" + segments.joined(separator: "/")
        let inchiKey = try CDKInChIKeyCodec.makeKey(from: inchi)
        return CDKInChINativeGenerationResult(inchi: inchi,
                                              inchiKey: inchiKey,
                                              status: .success,
                                              message: "")
    }

    private struct CanonicalizationResult {
        let heavyAtomOrder: [Int]
        let newIDByOldID: [Int: Int]
    }

    private struct IsotopeLayerSegments {
        let main: String
        let hydrogenSubLayer: String
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
                                          canonicalization: CanonicalizationResult) -> IsotopeLayerSegments {
        let atomByID = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0) })
        var tokens: [String] = []
        var hydrogenCarriers: [(mappedID: Int, element: String, dCount: Int, tCount: Int)] = []

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
            guard let mappedID = canonicalization.newIDByOldID[oldAtomID],
                  let atom = atomByID[oldAtomID] else { continue }
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

            if dCount > 0 || tCount > 0 {
                hydrogenCarriers.append((mappedID, atom.element, dCount, tCount))
            }
        }

        if tokens.isEmpty,
           hydrogenCarriers.count == 1,
           exchangeableHydrogenCarrierElements.contains(hydrogenCarriers[0].element.uppercased()) {
            return IsotopeLayerSegments(
                main: "",
                hydrogenSubLayer: isotopicHydrogenSuffix(dCount: hydrogenCarriers[0].dCount,
                                                         tCount: hydrogenCarriers[0].tCount)
            )
        }

        for carrier in hydrogenCarriers {
            let suffix = isotopicHydrogenSuffix(dCount: carrier.dCount, tCount: carrier.tCount)
            guard !suffix.isEmpty else { continue }
            tokens.append("\(carrier.mappedID)\(suffix)")
        }

        return IsotopeLayerSegments(main: tokens.joined(separator: ","), hydrogenSubLayer: "")
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

    private static func isRawHydrogenLike(_ symbol: String) -> Bool {
        let letters = symbol
            .trimmingCharacters(in: .whitespacesAndNewlines)
            .prefix { $0.isLetter }
            .uppercased()
        return letters == "H" || letters == "D" || letters == "T"
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

    private static func isotopicHydrogenSuffix(dCount: Int, tCount: Int) -> String {
        var suffix = ""
        if tCount > 0 {
            suffix += "T" + (tCount == 1 ? "" : "\(tCount)")
        }
        if dCount > 0 {
            suffix += "D" + (dCount == 1 ? "" : "\(dCount)")
        }
        return suffix
    }

    private static func generateHydrogenOnly(for molecule: Molecule) throws -> CDKInChINativeGenerationResult {
        guard molecule.atoms.allSatisfy({ isHydrogenSymbol($0.element) }) else {
            throw ChemError.unsupported("InChI generation requires at least one non-hydrogen atom.")
        }

        let hydrogens = molecule.atoms.sorted { $0.id < $1.id }
        let totalCharge = hydrogens.reduce(0) { $0 + $1.charge }
        let hasSingleHydronStyle = hydrogens.count == 1 && totalCharge != 0

        var segments: [String] = []
        if !hasSingleHydronStyle {
            segments.append("H" + (hydrogens.count == 1 ? "" : "\(hydrogens.count)"))
        }
        if hydrogens.count == 2 && molecule.bondCount > 0 {
            segments.append("h1H")
        }
        if totalCharge != 0 {
            segments.append("p\(signedInteger(totalCharge))")
        }

        let isotope = hydrogenOnlyIsotopeSegments(for: hydrogens, hasBonds: molecule.bondCount > 0)
        if !isotope.main.isEmpty || !isotope.hydrogenSubLayer.isEmpty {
            segments.append("i\(isotope.main)")
            if !isotope.hydrogenSubLayer.isEmpty {
                segments.append("h\(isotope.hydrogenSubLayer)")
            }
        }

        let inchi = "InChI=1S/" + segments.joined(separator: "/")
        let inchiKey = try CDKInChIKeyCodec.makeKey(from: inchi)
        return CDKInChINativeGenerationResult(inchi: inchi,
                                              inchiKey: inchiKey,
                                              status: .success,
                                              message: "")
    }

    private static func hydrogenOnlyIsotopeSegments(for hydrogens: [Atom],
                                                    hasBonds: Bool) -> IsotopeLayerSegments {
        guard !hydrogens.isEmpty else {
            return IsotopeLayerSegments(main: "", hydrogenSubLayer: "")
        }

        if hydrogens.count == 1 {
            let isotope = hydrogens[0].isotopeMassNumber ?? 1
            let suffix = isotopicHydrogenSuffix(dCount: isotope == 2 ? 1 : 0,
                                                tCount: isotope == 3 ? 1 : 0)
            return IsotopeLayerSegments(main: "", hydrogenSubLayer: suffix)
        }

        if hydrogens.count == 2, hasBonds {
            let first = hydrogens[0]
            let second = hydrogens[1]
            var token = ""
            let firstIsotope = first.isotopeMassNumber ?? 1
            if firstIsotope != 1 {
                token += "1\(signedInteger(firstIsotope - 1))"
            }
            let secondSuffix = isotopicHydrogenSuffix(dCount: (second.isotopeMassNumber ?? 1) == 2 ? 1 : 0,
                                                      tCount: (second.isotopeMassNumber ?? 1) == 3 ? 1 : 0)
            if !secondSuffix.isEmpty {
                token += token.isEmpty ? "1\(secondSuffix)" : secondSuffix
            }
            return IsotopeLayerSegments(main: token, hydrogenSubLayer: "")
        }

        return IsotopeLayerSegments(main: "", hydrogenSubLayer: "")
    }

    private static func generateElementalIsotopicPairIfSupported(for molecule: Molecule) throws -> CDKInChINativeGenerationResult? {
        guard molecule.atomCount == 2,
              molecule.atoms.allSatisfy({ !isHydrogenSymbol($0.element) && $0.charge == 0 }),
              molecule.bonds.allSatisfy({ $0.order == .single }) else {
            return nil
        }

        let normalizedSymbols = molecule.atoms.map { normalizedElementSymbol($0.element) }
        guard let symbol = normalizedSymbols.first,
              !symbol.isEmpty,
              normalizedSymbols.dropFirst().allSatisfy({ $0 == symbol }) else {
            return nil
        }

        let isotopes = molecule.atoms.compactMap(\.isotopeMassNumber)
        guard isotopes.count == 2 else { return nil }

        let orderedMasses = isotopes.sorted(by: >)
        let trailingShift = orderedMasses[1] - orderedMasses[0]
        guard trailingShift != 0 else { return nil }

        // Official elemental isotope-pair InChI uses count-first formula tokens and
        // a relative isotope layer without an explicit connectivity segment.
        let inchi = "InChI=1S/2\(symbol)/i1+0;1\(signedInteger(trailingShift))"
        let inchiKey = try CDKInChIKeyCodec.makeKey(from: inchi)
        return CDKInChINativeGenerationResult(inchi: inchi,
                                              inchiKey: inchiKey,
                                              status: .success,
                                              message: "")
    }

    private static let periodicTableSymbols: [String] = [
        "H", "He",
        "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
        "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr",
        "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
        "In", "Sn", "Sb", "Te", "I", "Xe",
        "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
        "Ho", "Er", "Tm", "Yb", "Lu",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn",
        "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
        "Es", "Fm", "Md", "No", "Lr",
        "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
        "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
    ]

    private static let knownElementSymbols: Set<String> = Set(periodicTableSymbols.map { $0.uppercased() })

    private static let atomicNumberByElement: [String: Int] = Dictionary(
        uniqueKeysWithValues: periodicTableSymbols.enumerated().map { (index, symbol) in
            (symbol.uppercased(), index + 1)
        }
    )

    private static let exchangeableHydrogenCarrierElements: Set<String> = [
        "N", "O", "S", "SE", "TE"
    ]
}
