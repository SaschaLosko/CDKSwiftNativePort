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

        if let multicomponent = try generateMultiComponentIfNeeded(for: normalized) {
            return multicomponent
        }

        return try generateSingleComponent(for: normalized)
    }

    private struct CanonicalizationResult {
        let heavyAtomOrder: [Int]
        let newIDByOldID: [Int: Int]
    }

    private struct IsotopeLayerSegments {
        let main: String
        let hydrogenSubLayer: String
    }

    private struct ParsedComponentInChI {
        let inchi: String
        let formula: String
        let layers: [Character: String]
        let isotopeHydrogenLayer: String
    }

    private static func generateSingleComponent(for molecule: Molecule) throws -> CDKInChINativeGenerationResult {
        let rawHydrogenOnly = molecule.atoms.allSatisfy { isHydrogenSymbol($0.element) }

        let heavyAtomIDs = molecule.atoms
            .filter { !isHydrogenSymbol($0.element) }
            .map(\.id)
            .sorted()

        guard !heavyAtomIDs.isEmpty else {
            guard rawHydrogenOnly else {
                throw ChemError.unsupported("InChI generation requires at least one non-hydrogen atom.")
            }
            return try generateHydrogenOnly(for: molecule)
        }

        let hydrogenByHeavyAtom = hydrogenCountByHeavyAtom(in: molecule, heavyAtomIDs: heavyAtomIDs)
        let canonicalization = canonicalizeHeavyAtoms(in: molecule,
                                                      heavyAtomIDs: heavyAtomIDs,
                                                      hydrogenByHeavyAtom: hydrogenByHeavyAtom)

        let formula = try buildFormula(in: molecule,
                                       heavyAtomIDs: heavyAtomIDs,
                                       hydrogenByHeavyAtom: hydrogenByHeavyAtom)
        let connectivityLayer = buildConnectivityLayer(in: molecule, canonicalization: canonicalization)
        let hydrogenLayer = buildHydrogenLayer(in: molecule,
                                               canonicalization: canonicalization,
                                               hydrogenByHeavyAtom: hydrogenByHeavyAtom)
        let chargeLayer = buildChargeLayer(in: molecule)
        let isotopeLayer = buildIsotopeLayer(in: molecule, canonicalization: canonicalization)
        let doubleBondLayer = buildDoubleBondLayer(in: molecule, canonicalization: canonicalization)
        let tetrahedralLayer = buildTetrahedralLayer(in: molecule, canonicalization: canonicalization)

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

    private static func generateMultiComponentIfNeeded(for molecule: Molecule) throws -> CDKInChINativeGenerationResult? {
        let components = disconnectedReferenceComponents(in: molecule)
        guard components.count > 1,
              !molecule.atoms.allSatisfy({ isHydrogenSymbol($0.element) }) else {
            return nil
        }

        let heavyAtomIDs = molecule.atoms
            .filter { !isHydrogenSymbol($0.element) }
            .map(\.id)
            .sorted()
        let originalHydrogenByHeavyAtom = hydrogenCountByHeavyAtom(in: molecule, heavyAtomIDs: heavyAtomIDs)
        let preservedImplicitHydrogenCountByAtom = Dictionary(uniqueKeysWithValues: heavyAtomIDs.map { atomID in
            let explicitNeighbors = CDKDescriptorSupport.explicitHydrogenNeighborCount(on: atomID, in: molecule)
            let total = originalHydrogenByHeavyAtom[atomID] ?? 0
            return (atomID, max(0, total - explicitNeighbors))
        })

        let parsedComponents = try components
            .map { submolecule(from: molecule, atomIDs: $0, preservedImplicitHydrogenCountByAtom: preservedImplicitHydrogenCountByAtom) }
            .map(generateSingleComponent(for:))
            .map { try parseGeneratedComponentInChI($0.inchi) }
            .sorted(by: componentSortOrder)

        var segments: [String] = [compressFormulaComponents(parsedComponents.map(\.formula))]

        if let connectivityLayer = joinComponentLayer(parsedComponents.map { $0.layers["c"] ?? "" }) {
            segments.append("c\(connectivityLayer)")
        }
        if let hydrogenLayer = joinComponentLayer(parsedComponents.map { $0.layers["h"] ?? "" }) {
            segments.append("h\(hydrogenLayer)")
        }
        if let chargeLayer = joinComponentLayer(parsedComponents.map { $0.layers["q"] ?? "" }) {
            segments.append("q\(chargeLayer)")
        }

        let isotopeMain = joinComponentLayer(parsedComponents.map { $0.layers["i"] ?? "" })
        let isotopeHydrogen = joinComponentLayer(parsedComponents.map(\.isotopeHydrogenLayer))
        if let isotopeMain {
            segments.append("i\(isotopeMain)")
            if let isotopeHydrogen {
                segments.append("h\(isotopeHydrogen)")
            }
        }

        if let doubleBondLayer = joinComponentLayer(parsedComponents.map { $0.layers["b"] ?? "" }) {
            segments.append("b\(doubleBondLayer)")
        }
        if let tetrahedralLayer = joinComponentLayer(parsedComponents.map { $0.layers["t"] ?? "" }) {
            segments.append("t\(tetrahedralLayer)")
        }
        if let parityLayer = joinComponentLayer(parsedComponents.map { $0.layers["m"] ?? "" }) {
            segments.append("m\(parityLayer)")
        }
        if let stereoFlagLayer = joinComponentLayer(parsedComponents.map { $0.layers["s"] ?? "" }) {
            segments.append("s\(stereoFlagLayer)")
        }

        let inchi = "InChI=1S/" + segments.joined(separator: "/")
        let inchiKey = try CDKInChIKeyCodec.makeKey(from: inchi)
        return CDKInChINativeGenerationResult(inchi: inchi,
                                              inchiKey: inchiKey,
                                              status: .success,
                                              message: "")
    }

    private static func disconnectedReferenceComponents(in molecule: Molecule) -> [Set<Int>] {
        let atomIDs = molecule.atoms.map(\.id).sorted()
        guard !atomIDs.isEmpty else { return [] }

        var adjacency: [Int: Set<Int>] = [:]
        for atomID in atomIDs {
            adjacency[atomID] = []
        }

        for bond in molecule.bonds where !shouldDisconnectBondForReferenceInChI(bond, in: molecule) {
            adjacency[bond.a1, default: []].insert(bond.a2)
            adjacency[bond.a2, default: []].insert(bond.a1)
        }

        var seen: Set<Int> = []
        var components: [Set<Int>] = []
        for seed in atomIDs where !seen.contains(seed) {
            var stack = [seed]
            var component: Set<Int> = [seed]
            seen.insert(seed)

            while let current = stack.popLast() {
                for neighbor in adjacency[current, default: []] where !seen.contains(neighbor) {
                    seen.insert(neighbor)
                    component.insert(neighbor)
                    stack.append(neighbor)
                }
            }

            components.append(component)
        }

        return components
    }

    private static func shouldDisconnectBondForReferenceInChI(_ bond: Bond, in molecule: Molecule) -> Bool {
        guard let left = molecule.atom(id: bond.a1),
              let right = molecule.atom(id: bond.a2) else {
            return false
        }
        let leftIsMetal = isMetalElement(left.element)
        let rightIsMetal = isMetalElement(right.element)
        return leftIsMetal != rightIsMetal
    }

    private static func isMetalElement(_ element: String) -> Bool {
        let upper = normalizedElementSymbol(element).uppercased()
        guard !upper.isEmpty else { return false }
        return knownElementSymbols.contains(upper) && !nonMetalOrMetalloidSymbols.contains(upper)
    }

    private static func submolecule(from molecule: Molecule,
                                    atomIDs: Set<Int>,
                                    preservedImplicitHydrogenCountByAtom: [Int: Int]) -> Molecule {
        let atoms = molecule.atoms
            .filter { atomIDs.contains($0.id) }
            .map { atom -> Atom in
                guard !isHydrogenSymbol(atom.element) else { return atom }
                var updated = atom
                updated.explicitHydrogenCount = preservedImplicitHydrogenCountByAtom[atom.id] ?? atom.explicitHydrogenCount
                return updated
            }
        let bonds = molecule.bonds.filter { atomIDs.contains($0.a1) && atomIDs.contains($0.a2) }
        return Molecule(name: molecule.name,
                        externalID: molecule.externalID,
                        atoms: atoms,
                        bonds: bonds)
    }

    private static func parseGeneratedComponentInChI(_ inchi: String) throws -> ParsedComponentInChI {
        guard inchi.hasPrefix("InChI=1S/") else {
            throw ChemError.parseFailed("Expected generated component InChI, got '\(inchi)'.")
        }

        let parts = String(inchi.dropFirst("InChI=1S/".count))
            .split(separator: "/", omittingEmptySubsequences: false)
            .map(String.init)

        guard !parts.isEmpty else {
            throw ChemError.parseFailed("Generated component InChI is missing formula/layer segments.")
        }

        let formula: String
        let layerSegments: ArraySlice<String>
        if componentLayerSegment(parts[0]) {
            formula = ""
            layerSegments = parts[0...]
        } else {
            formula = parts[0]
            layerSegments = parts.dropFirst()
        }

        var layers: [Character: String] = [:]
        var isotopeHydrogenLayer = ""
        var previousKey: Character?
        for segment in layerSegments {
            guard let key = segment.first else { continue }
            let content = String(segment.dropFirst())
            if key == "h", previousKey == "i" {
                isotopeHydrogenLayer = content
                previousKey = key
                continue
            }
            layers[key] = content
            previousKey = key
        }

        return ParsedComponentInChI(inchi: inchi,
                                    formula: formula,
                                    layers: layers,
                                    isotopeHydrogenLayer: isotopeHydrogenLayer)
    }

    private static func componentLayerSegment(_ segment: String) -> Bool {
        guard let key = segment.first else { return false }
        return key.isLowercase && ["c", "h", "q", "i", "b", "t", "m", "s"].contains(key)
    }

    private static func componentSortOrder(_ lhs: ParsedComponentInChI, _ rhs: ParsedComponentInChI) -> Bool {
        let lhsHasCarbon = formulaContainsCarbon(lhs.formula)
        let rhsHasCarbon = formulaContainsCarbon(rhs.formula)
        if lhsHasCarbon != rhsHasCarbon { return lhsHasCarbon && !rhsHasCarbon }
        if lhs.formula != rhs.formula { return lhs.formula < rhs.formula }
        return lhs.inchi < rhs.inchi
    }

    private static func formulaContainsCarbon(_ formula: String) -> Bool {
        var index = formula.startIndex
        while index < formula.endIndex {
            let char = formula[index]
            if char.isNumber {
                index = formula.index(after: index)
                continue
            }
            if char == "C" {
                let next = formula.index(after: index)
                if next == formula.endIndex || !formula[next].isLowercase {
                    return true
                }
            }
            index = formula.index(after: index)
        }
        return false
    }

    private static func compressFormulaComponents(_ formulas: [String]) -> String {
        guard !formulas.isEmpty else { return "" }

        var tokens: [String] = []
        var index = 0
        while index < formulas.count {
            let formula = formulas[index]
            var end = index
            while end + 1 < formulas.count, formulas[end + 1] == formula {
                end += 1
            }
            let count = end - index + 1
            tokens.append(count == 1 ? formula : "\(count)\(formula)")
            index = end + 1
        }
        return tokens.joined(separator: ".")
    }

    private static func joinComponentLayer(_ tokens: [String]) -> String? {
        guard tokens.contains(where: { !$0.isEmpty }) else { return nil }

        var out: [String] = []
        var index = 0
        while index < tokens.count {
            let token = tokens[index]
            if token.isEmpty {
                out.append("")
                index += 1
                continue
            }

            var end = index
            while end + 1 < tokens.count, tokens[end + 1] == token, !tokens[end + 1].isEmpty {
                end += 1
            }
            let count = end - index + 1
            out.append(count == 1 ? token : "\(count)*\(token)")
            index = end + 1
        }

        return out.joined(separator: ";")
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
        if composition["C"] != nil {
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
        } else {
            for element in composition.keys.sorted() {
                guard let count = composition[element], count > 0 else { continue }
                tokens.append(element + (count == 1 ? "" : "\(count)"))
            }
        }

        return tokens.joined()
    }

    private static func buildConnectivityLayer(in molecule: Molecule,
                                               canonicalization: CanonicalizationResult) -> String {
        var edgeSet = Set<InChIEdgeKey>()
        var adjacency: [Int: Set<Int>] = [:]
        for bond in molecule.bonds {
            guard let a = canonicalization.newIDByOldID[bond.a1],
                  let b = canonicalization.newIDByOldID[bond.a2] else {
                continue
            }
            edgeSet.insert(InChIEdgeKey(a, b))
            adjacency[a, default: []].insert(b)
            adjacency[b, default: []].insert(a)
        }

        let atomIDs = Set(canonicalization.newIDByOldID.values)
        let components = connectivityComponents(atomIDs: atomIDs, adjacency: adjacency)

        let tokens = components.compactMap { component -> String? in
            guard component.count > 1 else { return nil }

            let componentEdgeCount = edgeSet.reduce(into: 0) { count, edge in
                if component.contains(edge.a), component.contains(edge.b) {
                    count += 1
                }
            }

            if componentEdgeCount == component.count - 1 {
                return renderTreeConnectivity(component: component, adjacency: adjacency)
            }

            let componentEdges = edgeSet.filter { component.contains($0.a) && component.contains($0.b) }
            return componentEdges.sorted { lhs, rhs in
                if lhs.a != rhs.a { return lhs.a < rhs.a }
                return lhs.b < rhs.b
            }.map { "\($0.a)-\($0.b)" }.joined(separator: ";")
        }

        return tokens.joined(separator: ";")
    }

    private static func buildHydrogenLayer(in molecule: Molecule,
                                           canonicalization: CanonicalizationResult,
                                           hydrogenByHeavyAtom: [Int: Int]) -> String {
        let atomByID = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0) })
        let tokens = canonicalization.heavyAtomOrder.compactMap { oldAtomID -> (atomID: Int, count: Int, exchangeable: Bool)? in
            let count = hydrogenByHeavyAtom[oldAtomID] ?? 0
            guard count > 0,
                  let atomID = canonicalization.newIDByOldID[oldAtomID],
                  let atom = atomByID[oldAtomID] else {
                return nil
            }
            let exchangeable = exchangeableHydrogenCarrierElements.contains(normalizedElementSymbol(atom.element).uppercased())
            return (atomID, count, exchangeable)
        }

        let sorted = tokens.sorted { lhs, rhs in
            if lhs.exchangeable != rhs.exchangeable {
                return lhs.exchangeable && !rhs.exchangeable
            }
            return lhs.atomID < rhs.atomID
        }

        var grouped: [String] = []
        var index = 0
        while index < sorted.count {
            let current = sorted[index]
            var end = index
            while end + 1 < sorted.count,
                  sorted[end + 1].exchangeable == current.exchangeable,
                  sorted[end + 1].count == current.count,
                  sorted[end + 1].atomID == sorted[end].atomID + 1 {
                end += 1
            }

            let atomToken = end == index
                ? "\(current.atomID)"
                : "\(current.atomID)-\(sorted[end].atomID)"
            grouped.append(current.count == 1 ? "\(atomToken)H" : "\(atomToken)H\(current.count)")
            index = end + 1
        }

        return grouped.joined(separator: ",")
    }

    private static func connectivityComponents(atomIDs: Set<Int>,
                                               adjacency: [Int: Set<Int>]) -> [[Int]] {
        var visited: Set<Int> = []
        var components: [[Int]] = []

        for atomID in atomIDs.sorted() {
            if visited.contains(atomID) { continue }
            var stack = [atomID]
            var component: [Int] = []
            visited.insert(atomID)

            while let current = stack.popLast() {
                component.append(current)
                for neighbor in adjacency[current, default: []].sorted() where !visited.contains(neighbor) {
                    visited.insert(neighbor)
                    stack.append(neighbor)
                }
            }

            components.append(component.sorted())
        }

        return components
    }

    private static func renderTreeConnectivity(component: [Int],
                                               adjacency: [Int: Set<Int>]) -> String {
        let componentSet = Set(component)
        let leaves = component.filter { degree(of: $0, in: componentSet, adjacency: adjacency) <= 1 }.sorted()
        let start = leaves.first ?? component.first ?? 0
        let end = leaves.last ?? start
        let mainPath = path(between: start, and: end, in: componentSet, adjacency: adjacency) ?? [start]

        var nextOnMainPath: [Int: Int] = [:]
        for index in 0..<(mainPath.count - 1) {
            nextOnMainPath[mainPath[index]] = mainPath[index + 1]
        }

        return renderTreeConnectivityNode(start,
                                          parent: nil,
                                          forcedNext: nextOnMainPath,
                                          component: componentSet,
                                          adjacency: adjacency)
    }

    private static func renderTreeConnectivityNode(_ atomID: Int,
                                                   parent: Int?,
                                                   forcedNext: [Int: Int],
                                                   component: Set<Int>,
                                                   adjacency: [Int: Set<Int>]) -> String {
        let children = adjacency[atomID, default: []]
            .filter { $0 != parent && component.contains($0) }
            .sorted()
        let continuation = forcedNext[atomID] ?? chooseContinuationChild(from: children,
                                                                         parent: atomID,
                                                                         component: component,
                                                                         adjacency: adjacency)
        let branches = children.filter { $0 != continuation }

        var out = "\(atomID)"
        if !branches.isEmpty {
            let branchText = branches
                .sorted { lhs, rhs in
                    let left = subtreeLeafExtremes(root: lhs,
                                                   parent: atomID,
                                                   component: component,
                                                   adjacency: adjacency)
                    let right = subtreeLeafExtremes(root: rhs,
                                                    parent: atomID,
                                                    component: component,
                                                    adjacency: adjacency)
                    if left.min != right.min { return left.min < right.min }
                    return left.max < right.max
                }
                .map { renderTreeConnectivityNode($0,
                                                  parent: atomID,
                                                  forcedNext: forcedNext,
                                                  component: component,
                                                  adjacency: adjacency) }
                .joined(separator: ",")
            out += "(\(branchText))"
        }

        if let continuation {
            let continuationText = renderTreeConnectivityNode(continuation,
                                                              parent: atomID,
                                                              forcedNext: forcedNext,
                                                              component: component,
                                                              adjacency: adjacency)
            if branches.isEmpty {
                out += "-\(continuationText)"
            } else {
                out += continuationText
            }
        }

        return out
    }

    private static func chooseContinuationChild(from children: [Int],
                                                parent: Int,
                                                component: Set<Int>,
                                                adjacency: [Int: Set<Int>]) -> Int? {
        guard !children.isEmpty else { return nil }
        return children.max { lhs, rhs in
            let left = subtreeLeafExtremes(root: lhs, parent: parent, component: component, adjacency: adjacency)
            let right = subtreeLeafExtremes(root: rhs, parent: parent, component: component, adjacency: adjacency)
            if left.max != right.max { return left.max < right.max }
            if left.min != right.min { return left.min < right.min }
            return lhs < rhs
        }
    }

    private static func subtreeLeafExtremes(root: Int,
                                            parent: Int,
                                            component: Set<Int>,
                                            adjacency: [Int: Set<Int>]) -> (min: Int, max: Int) {
        let children = adjacency[root, default: []].filter { $0 != parent && component.contains($0) }
        if children.isEmpty {
            return (root, root)
        }

        var minLeaf = Int.max
        var maxLeaf = Int.min
        for child in children {
            let childExtremes = subtreeLeafExtremes(root: child,
                                                    parent: root,
                                                    component: component,
                                                    adjacency: adjacency)
            minLeaf = min(minLeaf, childExtremes.min)
            maxLeaf = max(maxLeaf, childExtremes.max)
        }
        return (minLeaf, maxLeaf)
    }

    private static func degree(of atomID: Int,
                               in component: Set<Int>,
                               adjacency: [Int: Set<Int>]) -> Int {
        adjacency[atomID, default: []].filter { component.contains($0) }.count
    }

    private static func path(between start: Int,
                             and end: Int,
                             in component: Set<Int>,
                             adjacency: [Int: Set<Int>]) -> [Int]? {
        guard start != end else { return [start] }

        var queue = [start]
        var visited: Set<Int> = [start]
        var parentByNode: [Int: Int] = [:]

        while !queue.isEmpty {
            let current = queue.removeFirst()
            for neighbor in adjacency[current, default: []].sorted() where component.contains(neighbor) && !visited.contains(neighbor) {
                visited.insert(neighbor)
                parentByNode[neighbor] = current

                if neighbor == end {
                    var out = [end]
                    var cursor = end
                    while let parent = parentByNode[cursor] {
                        out.append(parent)
                        if parent == start { break }
                        cursor = parent
                    }
                    return out.reversed()
                }

                queue.append(neighbor)
            }
        }

        return nil
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

    private static let nonMetalOrMetalloidSymbols: Set<String> = [
        "H", "HE", "B", "C", "N", "O", "F", "NE",
        "SI", "P", "S", "CL", "AR",
        "GE", "AS", "SE", "BR", "KR",
        "SB", "TE", "I", "XE",
        "AT", "RN", "TS", "OG"
    ]

    private static let exchangeableHydrogenCarrierElements: Set<String> = [
        "N", "O", "S", "SE", "TE"
    ]
}
