import Foundation

/// CDK-style Kekulization support for assigning localized single/double bond
/// orders to aromatic systems.
public enum CDKKekulization {
    public static func kekulized(_ molecule: Molecule) throws -> Molecule {
        var copy = molecule
        try kekulize(&copy)
        return copy
    }

    public static func kekulize(_ molecule: inout Molecule) throws {
        let context = Context(molecule: molecule)
        let availableAtomIDs = availableAtomIDs(in: molecule, context: context)
        let matching = try perfectMatching(availableAtomIDs: availableAtomIDs, context: context)

        for index in molecule.bonds.indices where molecule.bonds[index].order == .aromatic {
            molecule.bonds[index].order = .single
        }

        for edge in matching {
            guard let bondIndex = context.bondIndexByEdge[edge] else { continue }
            guard molecule.bonds[bondIndex].order == .single else {
                throw ChemError.unsupported("Cannot assign Kekulé structure; non-sigma bond order is already assigned.")
            }
            molecule.bonds[bondIndex].order = .double
            molecule.bonds[bondIndex].stereo = .none
        }
    }

    private static func availableAtomIDs(in molecule: Molecule, context: Context) -> Set<Int> {
        var available = Set<Int>()
        for atom in molecule.atoms where context.aromaticAtomIDs.contains(atom.id) {
            guard let atomicNumber = atomicNumber(for: atom.element), atomicNumber > 0 else { continue }

            var piBondCount = 0
            var skipAtom = false
            for bond in molecule.bonds(forAtom: atom.id) {
                switch bond.order {
                case .double:
                    piBondCount += 1
                case .triple:
                    skipAtom = true
                case .single, .aromatic:
                    break
                }
            }
            guard !skipAtom else { continue }

            let valence = molecule.neighbors(of: atom.id).count
                + molecule.implicitHydrogenCount(for: atom.id)
                + piBondCount
            if isAvailableForPiBond(atomicNumber: atomicNumber, charge: atom.charge, valence: valence) {
                available.insert(atom.id)
            }
        }
        return available
    }

    private static func perfectMatching(availableAtomIDs: Set<Int>, context: Context) throws -> Set<EdgeKey> {
        guard !availableAtomIDs.isEmpty else { return [] }
        guard availableAtomIDs.count.isMultiple(of: 2) else {
            throw ChemError.unsupported("Cannot assign Kekulé structure without randomly creating radicals.")
        }

        let adjacency = context.neighborIDsByAtomID.mapValues { neighbors in
            neighbors.filter { availableAtomIDs.contains($0) }
        }
        guard let matching = backtrackingPerfectMatching(unmatched: availableAtomIDs, adjacency: adjacency) else {
            throw ChemError.unsupported("Cannot assign Kekulé structure without randomly creating radicals.")
        }
        return matching
    }

    private static func backtrackingPerfectMatching(
        unmatched: Set<Int>,
        adjacency: [Int: [Int]]
    ) -> Set<EdgeKey>? {
        guard !unmatched.isEmpty else { return [] }
        guard let atomID = unmatched.min(by: { lhs, rhs in
            let leftDegree = adjacency[lhs, default: []].filter(unmatched.contains).count
            let rightDegree = adjacency[rhs, default: []].filter(unmatched.contains).count
            if leftDegree != rightDegree { return leftDegree < rightDegree }
            return lhs < rhs
        }) else {
            return []
        }

        let candidates = adjacency[atomID, default: []]
            .filter { unmatched.contains($0) && $0 != atomID }
            .sorted()
        guard !candidates.isEmpty else { return nil }

        for neighborID in candidates {
            var nextUnmatched = unmatched
            nextUnmatched.remove(atomID)
            nextUnmatched.remove(neighborID)
            guard var downstream = backtrackingPerfectMatching(unmatched: nextUnmatched, adjacency: adjacency) else {
                continue
            }
            downstream.insert(EdgeKey(atomID, neighborID))
            return downstream
        }
        return nil
    }

    private static func isAvailableForPiBond(atomicNumber: Int, charge: Int, valence: Int) -> Bool {
        switch atomicNumber {
        case 5:
            if charge == 0 && valence <= 2 { return true }
            if charge == -1 && valence <= 3 { return true }
        case 6, 14, 32, 50:
            if charge == 0 && valence <= 3 { return true }
        case 7, 15, 33, 51:
            if charge == 0 { return valence <= 2 || valence == 4 }
            if charge == 1 { return valence <= 3 }
        case 8, 16, 34, 52:
            if charge == 0 { return valence <= 1 || valence == 3 || valence == 5 }
            if charge == 1 { return valence <= 2 || valence == 4 }
        default:
            break
        }
        return false
    }

    private static func atomicNumber(for symbol: String) -> Int? {
        atomicNumberBySymbol[canonicalElementSymbol(symbol).uppercased()]
    }

    private static func canonicalElementSymbol(_ raw: String) -> String {
        let trimmed = raw.trimmingCharacters(in: .whitespacesAndNewlines)
        guard let first = trimmed.first, first.isLetter else { return trimmed }
        return String(first).uppercased() + String(trimmed.dropFirst()).lowercased()
    }
}

private struct Context {
    let aromaticAtomIDs: Set<Int>
    let neighborIDsByAtomID: [Int: [Int]]
    let bondIndexByEdge: [EdgeKey: Int]

    init(molecule: Molecule) {
        let explicitAromaticAtomIDs = Set(molecule.atoms.compactMap { $0.aromatic ? $0.id : nil })
        let aromaticBondAtomIDs = molecule.bonds.flatMap { bond -> [Int] in
            bond.order == .aromatic ? [bond.a1, bond.a2] : []
        }
        aromaticAtomIDs = explicitAromaticAtomIDs.union(aromaticBondAtomIDs)

        var neighbors: [Int: Set<Int>] = [:]
        var bondIndices: [EdgeKey: Int] = [:]
        for atom in molecule.atoms {
            neighbors[atom.id] = []
        }
        for index in molecule.bonds.indices {
            let bond = molecule.bonds[index]
            neighbors[bond.a1, default: []].insert(bond.a2)
            neighbors[bond.a2, default: []].insert(bond.a1)
            bondIndices[EdgeKey(bond.a1, bond.a2)] = index
        }
        neighborIDsByAtomID = neighbors.mapValues { $0.sorted() }
        bondIndexByEdge = bondIndices
    }
}

private struct EdgeKey: Hashable {
    let a: Int
    let b: Int

    init(_ lhs: Int, _ rhs: Int) {
        a = min(lhs, rhs)
        b = max(lhs, rhs)
    }
}

private let atomicNumberBySymbol: [String: Int] = [
    "H": 1, "He": 2,
    "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
    "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26,
    "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34,
    "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42,
    "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58,
    "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66,
    "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72, "Ta": 73, "W": 74,
    "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82,
    "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90,
    "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98,
    "Es": 99, "Fm": 100, "Md": 101, "No": 102, "Lr": 103, "Rf": 104, "Db": 105,
    "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111, "Cn": 112,
    "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118
]
