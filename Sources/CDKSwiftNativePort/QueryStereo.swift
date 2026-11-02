import Foundation

/// A specified graph plus ordered stereo elements, including coordination geometries.
/// Use this reader when coordination stereo must be retained for query matching.
public struct CDKQueryStructure: Hashable, Codable, Sendable {
    public var molecule: Molecule
    public var stereo: [CDKQueryAtomStereo]
    public init(molecule: Molecule, stereo: [CDKQueryAtomStereo] = []) {
        self.molecule = molecule; self.stereo = stereo
    }
    public static func parseSmiles(_ text: String) throws -> Self {
        // Keep target chemistry distinct from SMARTS. The legacy SMILES reader accepts
        // some query decorators, so validate concrete bracket syntax before using it.
        let brackets = try NSRegularExpression(pattern: "\\[[^\\]]*\\]")
        let concreteAtom = try NSRegularExpression(
            pattern:
                "^\\[(?:[0-9]+)?(?:[A-Z][a-z]?|[bcnops]|se|as)(?:@(?:TH[12]|SP[1-3]|TB(?:20|1[0-9]|[1-9])|OH(?:30|[12][0-9]|[1-9])|@)?)?(?:H[0-9]*)?(?:[+-](?:[0-9]+|[+-]*))?(?::[0-9]+)?\\]$"
        )
        for match in brackets.matches(in: text, range: NSRange(text.startIndex..., in: text)) {
            let atom = (text as NSString).substring(with: match.range)
            guard concreteAtom.firstMatch(in: atom, range: NSRange(atom.startIndex..., in: atom)) != nil else {
                throw CDKQueryError.invalid("Target SMILES contains a query or unsupported atom decorator.")
            }
        }
        let query = try CDKSmartsParser.parse(text)
        // The ordinary molecule reader only models tetrahedral atom parity. Remove
        // atom stereo tokens before reading and retain all classes in the sidecar.
        let regex = try NSRegularExpression(
            pattern: "@(?:TH[12]|SP[1-3]|TB(?:[12][0-9]|[1-9])|OH(?:[123][0-9]|[1-9])|@)?")
        let stripped = regex.stringByReplacingMatches(
            in: text, range: NSRange(text.startIndex..., in: text), withTemplate: "")
        let molecule = try CDKSmilesParser(generateCoordinates: false).parseSmiles(stripped)
        guard molecule.atoms.allSatisfy({ (QueryElements.number($0.element) ?? 0) > 0 }) else {
            throw CDKQueryError.invalid("Target SMILES requires specified elements.")
        }
        guard molecule.atoms.count == query.molecule.atoms.count else {
            throw CDKQueryError.invalid("Stereo graph differs from molecule graph.")
        }
        let map = Dictionary(uniqueKeysWithValues: zip(query.molecule.atoms, molecule.atoms).map { ($0.id, $1.id) })
        var stereo = query.stereo
        for index in stereo.indices {
            stereo[index].focus = map[stereo[index].focus]!
            stereo[index].carriers = stereo[index].carriers.map { map[$0]! }
            guard !stereo[index].allowsUnspecified else {
                throw CDKQueryError.invalid("A target SMILES cannot contain optional stereo.")
            }
        }
        return .init(molecule: molecule, stereo: stereo)
    }
}

extension CDKQueryMatcher {
    public static func match(
        _ query: CDKQueryPattern, in target: CDKQueryStructure,
        mode: Mode = .substructure, options: CDKQueryMatchOptions = .init()
    ) -> CDKQueryMatchResult {
        match(query, in: target.molecule, mode: mode, targetStereo: target.stereo, options: options)
    }
}

enum QueryStereoSupport {
    static func elements(_ molecule: Molecule, explicit: [CDKQueryAtomStereo]) throws -> [CDKQueryAtomStereo] {
        guard molecule.cxState?.racemic != true, molecule.cxState?.racemicFragments.isEmpty != false else {
            throw CDKQueryError.unsupported("Racemic fragment mixtures require explicit grouped stereo metadata.")
        }
        for (index, group) in molecule.cxState?.stereoGroups ?? [:] {
            guard molecule.atoms.indices.contains(index), molecule.atoms[index].cxStereoGroup == group else {
                throw CDKQueryError.invalid("Inconsistent enhanced stereo metadata.")
            }
        }
        var result = explicit
        let explicitIDs = Set(explicit.map(\.focus))
        guard explicitIDs.count == explicit.count else { throw CDKQueryError.invalid("Duplicate stereo focus.") }
        for atom in molecule.atoms where atom.chirality != .none && !explicitIDs.contains(atom.id) {
            guard let carriers = atom.ligandOrderingAtomIDs else {
                throw CDKQueryError.unsupported("Stereo requires ordered ligands.")
            }
            var group: Int?
            if let encoded = atom.cxStereoGroup, let decoded = CDKCxSmilesParser.decodeStereoGroup(encoded),
                decoded.kind != "abs"
            {
                group = encoded
            }
            result.append(
                .init(
                    focus: atom.id, carriers: carriers, geometry: .tetrahedral,
                    configuration: atom.chirality == .clockwise ? 1 : 2, relativeGroup: group))
        }
        for element in result {
            let neighbors = Set(molecule.neighbors(of: element.focus))
            let actual = element.carriers.filter { $0 != element.focus }
            guard molecule.atom(id: element.focus) != nil,
                element.carriers.count == element.geometry.carrierCount,
                (1...element.geometry.configurationCount).contains(element.configuration),
                Set(actual).count == actual.count, Set(actual) == neighbors,
                element.relativeGroup == nil || element.geometry == .tetrahedral
            else {
                throw CDKQueryError.invalid("Invalid ordered stereo element.")
            }
            if element.geometry == .tetrahedral && element.carriers.filter({ $0 == element.focus }).count > 1 {
                throw CDKQueryError.invalid("Tetrahedral stereo has multiple implicit ligands.")
            }
        }
        return result
    }

    static func representatives(_ geometry: CDKQueryAtomStereo.Geometry) -> [[Int]] {
        switch geometry {
        case .tetrahedral: [[0, 1, 2, 3], [0, 1, 3, 2]]
        case .squarePlanar: QueryStereoPermutations.squarePlanar
        case .trigonalBipyramidal: QueryStereoPermutations.trigonalBipyramidal
        case .octahedral: QueryStereoPermutations.octahedral
        }
    }
    static func rotations(_ geometry: CDKQueryAtomStereo.Geometry) -> [[Int]] {
        switch geometry {
        case .tetrahedral:
            return [
                [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2], [1, 0, 3, 2],
                [1, 2, 0, 3], [1, 3, 2, 0], [2, 0, 1, 3], [2, 1, 3, 0],
                [2, 3, 0, 1], [3, 0, 2, 1], [3, 1, 0, 2], [3, 2, 1, 0],
            ]
        case .squarePlanar: return QueryStereoPermutations.squarePlanarRotations
        case .trigonalBipyramidal: return QueryStereoPermutations.trigonalBipyramidalRotations
        case .octahedral: return QueryStereoPermutations.octahedralRotations
        }
    }
    static func equivalent(_ a: CDKQueryAtomStereo, _ b: CDKQueryAtomStereo, mapping: [Int: Int]) -> Bool {
        guard a.geometry == b.geometry else { return false }
        let orders = representatives(a.geometry)
        let left = orders[a.configuration - 1].compactMap { mapping[a.carriers[$0]] }
        let right = orders[b.configuration - 1].map { b.carriers[$0] }
        guard left.count == right.count else { return false }
        return rotations(a.geometry).contains { rotation in rotation.map { left[$0] } == right }
    }

    static func compare(
        query: [CDKQueryAtomStereo], target: [CDKQueryAtomStereo], mapping: [Int: Int],
        targetMolecule: Molecule
    ) throws -> QueryTruth {
        var missing = false
        // Equations qGroup XOR tGroup = observed inversion. Node 0 is absolute.
        var equations: [(String, String, Bool)] = []
        for original in query {
            var element = original
            guard let focus = mapping[element.focus] else { throw CDKQueryError.invalid("Unmapped stereo focus.") }
            guard var other = target.first(where: { $0.focus == focus }) else {
                if !element.allowsUnspecified { missing = true }
                continue
            }
            if element.geometry != .tetrahedral {
                if !coordinationEquivalent(element, other, mapping: mapping) { return .no }
                continue
            }
            guard other.geometry == .tetrahedral else { return .no }
            var stereoMap = mapping
            if element.geometry == .tetrahedral {
                func normalize(_ id: Int) -> Int {
                    guard let atom = targetMolecule.atom(id: id), atom.element == "H",
                        atom.isotopeMassNumber == nil, atom.charge == 0,
                        targetMolecule.bonds(forAtom: id).count == 1
                    else { return id }
                    return focus
                }
                element.focus = focus
                element.carriers = try element.carriers.map {
                    guard let mapped = mapping[$0] else { throw CDKQueryError.invalid("Unmapped stereo carrier.") }
                    return normalize(mapped)
                }
                other.carriers = other.carriers.map(normalize)
                stereoMap = Dictionary(uniqueKeysWithValues: Set(element.carriers + [focus]).map { ($0, $0) })
            }
            let same = equivalent(element, other, mapping: stereoMap)
            var inverted = element; inverted.configuration = 3 - element.configuration
            guard same || equivalent(inverted, other, mapping: stereoMap) else { return .no }
            let left = element.relativeGroup.map { "q\($0)" } ?? "absolute"
            let right = other.relativeGroup.map { "t\($0)" } ?? "absolute"
            equations.append((left, right, !same))
        }
        var values: [String: Bool] = ["absolute": false]
        for seed in ["absolute"] + equations.flatMap({ [$0.0, $0.1] }) {
            if values[seed] == nil { values[seed] = false }
            var queue = [seed]; var cursor = 0
            while cursor < queue.count {
                let current = queue[cursor]; cursor += 1
                for (a, b, flip) in equations where a == current || b == current {
                    let next = a == current ? b : a
                    let value = values[current]! != flip
                    if let existing = values[next] {
                        if existing != value { return .no }
                    } else {
                        values[next] = value; queue.append(next)
                    }
                }
            }
        }
        return missing ? .unknown : .yes
    }

    struct Alkene {
        let a: Int; let b: Int; let left: Int; let right: Int; let cis: Bool
    }
    static func alkenes(_ molecule: Molecule) throws -> [Alkene] {
        var result: [Alkene] = []
        var consumed: Set<Int> = []
        func direction(_ bond: Bond, from center: Int) -> Int? {
            let sign: Int
            switch bond.stereo {
            case .up, .downReversed: sign = 1
            case .down, .upReversed: sign = -1
            default: return nil
            }
            return bond.a1 == center ? sign : -sign
        }
        for bond in molecule.bonds where bond.order == .double {
            if let stereo = bond.doubleBondStereo {
                guard let refs = bond.stereoReferenceAtomIDs, refs.count == 4,
                    Set(refs[1...2]) == Set([bond.a1, bond.a2]),
                    molecule.bond(between: refs[0], and: refs[1]) != nil,
                    molecule.bond(between: refs[2], and: refs[3]) != nil
                else {
                    throw CDKQueryError.invalid("Invalid geometric stereo references.")
                }
                result.append(.init(a: refs[1], b: refs[2], left: refs[0], right: refs[3], cis: stereo == .cis))
                // Directional bonds may also be retained by an importer.
                consumed.formUnion(
                    molecule.bonds(forAtom: bond.a1).filter { direction($0, from: bond.a1) != nil }.map(\.id))
                consumed.formUnion(
                    molecule.bonds(forAtom: bond.a2).filter { direction($0, from: bond.a2) != nil }.map(\.id))
                continue
            }
            func sides(_ center: Int) -> [Bond] {
                molecule.bonds(forAtom: center).filter {
                    $0.order == .single && direction($0, from: center) != nil
                        && (bond.stereo == .either
                            || (molecule.atom(id: $0.a1)?.chirality == AtomChirality.none
                                && molecule.atom(id: $0.a2)?.chirality == AtomChirality.none))
                }
            }
            let left = sides(bond.a1); let right = sides(bond.a2)
            if left.isEmpty && right.isEmpty { continue }
            guard let a = left.first, let b = right.first else {
                throw CDKQueryError.unsupported("Incomplete geometric stereo directions.")
            }
            for (sides, center) in [(left, bond.a1), (right, bond.a2)] {
                guard sides.count <= 2,
                    sides.count != 2 || direction(sides[0], from: center) != direction(sides[1], from: center)
                else {
                    throw CDKQueryError.invalid("Conflicting geometric stereo directions.")
                }
                consumed.formUnion(sides.map(\.id))
            }
            result.append(
                .init(
                    a: bond.a1, b: bond.a2, left: a.a1 == bond.a1 ? a.a2 : a.a1,
                    right: b.a1 == bond.a2 ? b.a2 : b.a1,
                    cis: direction(a, from: bond.a1) == direction(b, from: bond.a2)))
        }
        for bond in molecule.bonds where bond.stereo != .none && bond.order != .double && !consumed.contains(bond.id) {
            guard
                molecule.atom(id: bond.a1)?.chirality != AtomChirality.none
                    || molecule.atom(id: bond.a2)?.chirality != AtomChirality.none
            else {
                throw CDKQueryError.unsupported("Uninterpreted directional stereo bond.")
            }
        }
        return result
    }
    static func doubleBonds(_ query: Molecule, target: Molecule, mapping: [Int: Int]) throws -> QueryTruth {
        let expected = try alkenes(query)
        if expected.isEmpty { return .yes }
        let actual = try alkenes(target)
        var missing = false
        for element in expected {
            guard let a = mapping[element.a], let b = mapping[element.b],
                let left = mapping[element.left], let right = mapping[element.right]
            else {
                throw CDKQueryError.invalid("Unmapped geometric stereo reference.")
            }
            guard let other = actual.first(where: { Set([$0.a, $0.b]) == Set([a, b]) }) else {
                missing = true; continue
            }
            let targetLeft = other.a == a ? other.left : other.right
            let targetRight = other.a == a ? other.right : other.left
            if (element.cis == other.cis) != ((left == targetLeft) == (right == targetRight)) { return .no }
        }
        return missing ? .unknown : .yes
    }
}
