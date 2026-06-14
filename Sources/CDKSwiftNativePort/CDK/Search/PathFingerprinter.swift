import Foundation

public enum CDKPathFingerprinterError: Error, Equatable, Sendable {
    case pathLimitExceeded(limit: Int)
}

/// Path-based fingerprint compatible with CDK's standard `Fingerprinter`.
///
/// CDK uses this fingerprint as a fast substructure screen: if a query bit is
/// absent in the target, the exact substructure match can be skipped safely.
/// The final decision still belongs to graph matching because hashed path
/// fingerprints can produce false positives.
public final class CDKPathFingerprinter: @unchecked Sendable {
    public static let defaultSize = 1024
    public static let defaultSearchDepth = 7
    public static let defaultPathLimit = 42_000

    public let size: Int
    public let searchDepth: Int
    public let pathLimit: Int
    public let hashPseudoAtoms: Bool
    public let hashExplicitHydrogens: Bool

    public init(
        size: Int = CDKPathFingerprinter.defaultSize,
        searchDepth: Int = CDKPathFingerprinter.defaultSearchDepth,
        pathLimit: Int = CDKPathFingerprinter.defaultPathLimit,
        hashPseudoAtoms: Bool = false,
        hashExplicitHydrogens: Bool = false
    ) {
        precondition(size > 0, "Fingerprint size must be positive.")
        precondition(searchDepth >= 0, "Search depth cannot be negative.")
        precondition(pathLimit > 0, "Path limit must be positive.")
        self.size = size
        self.searchDepth = searchDepth
        self.pathLimit = pathLimit
        self.hashPseudoAtoms = hashPseudoAtoms
        self.hashExplicitHydrogens = hashExplicitHydrogens
    }

    public func bitFingerprint(for molecule: Molecule) throws -> CDKBitFingerprint {
        var calculator = CDKPathFingerprintCalculator(
            molecule: molecule,
            size: size,
            searchDepth: searchDepth,
            pathLimit: pathLimit,
            hashPseudoAtoms: hashPseudoAtoms,
            hashExplicitHydrogens: hashExplicitHydrogens)
        return try calculator.bitFingerprint()
    }

    public func screensIn(query: Molecule, target: Molecule) throws -> Bool {
        let queryFingerprint = try bitFingerprint(for: query)
        let targetFingerprint = try bitFingerprint(for: target)
        return targetFingerprint.isSuperset(of: queryFingerprint)
    }
}

public typealias CDKSubstructureScreenFingerprinter = CDKPathFingerprinter

private struct CDKPathFingerprintCalculator {
    private let size: Int
    private let maxDepth: Int
    private let pathLimit: Int
    private let hashPseudoAtoms: Bool
    private let hashExplicitHydrogens: Bool
    private let atoms: [Atom]
    private let bonds: [Bond]
    private let atomIDToIndex: [Int: Int]
    private var adjacency: [[PathEdge]]
    private var visitedAtomIndices: Set<Int> = []
    private var atomPath: [Int] = []
    private var bondPath: [Int] = []
    private var pathCount = 0
    private var bitIndices: Set<Int> = []

    init(
        molecule: Molecule,
        size: Int,
        searchDepth: Int,
        pathLimit: Int,
        hashPseudoAtoms: Bool,
        hashExplicitHydrogens: Bool
    ) {
        self.size = size
        self.maxDepth = searchDepth + 1
        self.pathLimit = pathLimit
        self.hashPseudoAtoms = hashPseudoAtoms
        self.hashExplicitHydrogens = hashExplicitHydrogens
        self.atoms = molecule.atoms
        self.bonds = molecule.bonds
        self.atomIDToIndex = Dictionary(
            uniqueKeysWithValues: molecule.atoms.enumerated().map { ($0.element.id, $0.offset) })
        self.adjacency = Array(repeating: [], count: molecule.atoms.count)

        for (bondIndex, bond) in molecule.bonds.enumerated() {
            guard let firstIndex = atomIDToIndex[bond.a1],
                let secondIndex = atomIDToIndex[bond.a2]
            else {
                continue
            }
            adjacency[firstIndex].append(PathEdge(atomIndex: secondIndex, bondIndex: bondIndex))
            adjacency[secondIndex].append(PathEdge(atomIndex: firstIndex, bondIndex: bondIndex))
        }
    }

    mutating func bitFingerprint() throws -> CDKBitFingerprint {
        bitIndices.removeAll(keepingCapacity: true)

        for atomIndex in atoms.indices {
            pathCount = 0
            visitedAtomIndices = [atomIndex]
            try traverse(from: atomIndex, through: nil)
            visitedAtomIndices.removeAll(keepingCapacity: true)
        }

        return CDKBitFingerprint(size: size, bits: bitIndices)
    }

    private mutating func traverse(
        from atomIndex: Int,
        through bondIndex: Int?
    ) throws {
        guard !skipAtom(atoms[atomIndex]) else {
            return
        }

        atomPath.append(atomIndex)
        if let bondIndex {
            bondPath.append(bondIndex)
        }
        defer {
            _ = atomPath.popLast()
            if bondIndex != nil {
                _ = bondPath.popLast()
            }
        }

        try storeCurrentPath()

        guard atomPath.count < maxDepth else {
            return
        }

        for edge in adjacency[atomIndex] {
            guard edge.bondIndex != bondIndex,
                !visitedAtomIndices.contains(edge.atomIndex)
            else {
                continue
            }

            visitedAtomIndices.insert(edge.atomIndex)
            try traverse(from: edge.atomIndex, through: edge.bondIndex)
            visitedAtomIndices.remove(edge.atomIndex)
        }
    }

    private mutating func storeCurrentPath() throws {
        pathCount += 1
        guard pathCount <= pathLimit else {
            throw CDKPathFingerprinterError.pathLimitExceeded(limit: pathLimit)
        }

        let hash: Int32
        if bondPath.isEmpty {
            hash = javaStringHash(atomSymbol(for: atoms[atomPath[0]]))
        } else if comparePath() >= 0 {
            hash = pathHash(reversed: false)
        } else {
            hash = pathHash(reversed: true)
        }

        var random = CDKJavaRandom(seed: hash)
        bitIndices.insert(random.nextInt(bound: size))
    }

    private func comparePath() -> Int {
        var forwardIndex = 0
        var reverseIndex = atomPath.count - 1

        var comparison = compareAtoms(atomPath[forwardIndex], atomPath[reverseIndex])
        if comparison != 0 {
            return comparison
        }

        forwardIndex += 1
        reverseIndex -= 1
        while reverseIndex != 0 {
            comparison = compareBonds(bondPath[forwardIndex - 1], bondPath[reverseIndex])
            if comparison != 0 {
                return comparison
            }
            comparison = compareAtoms(atomPath[forwardIndex], atomPath[reverseIndex])
            if comparison != 0 {
                return comparison
            }
            forwardIndex += 1
            reverseIndex -= 1
        }
        return 0
    }

    private func pathHash(reversed: Bool) -> Int32 {
        var hash: Int32 = 0

        if reversed {
            var atomIndex = atomPath.count - 1
            hash = appendCDKPathHash(hash, atomSymbol(for: atoms[atomPath[atomIndex]]))
            while atomIndex > 0 {
                let bondIndex = atomIndex - 1
                hash = appendCDKPathHash(hash, bondSymbol(for: bonds[bondPath[bondIndex]]))
                hash = appendCDKPathHash(hash, atomSymbol(for: atoms[atomPath[bondIndex]]))
                atomIndex -= 1
            }
        } else {
            hash = appendCDKPathHash(hash, atomSymbol(for: atoms[atomPath[0]]))
            for index in 1..<atomPath.count {
                hash = appendCDKPathHash(hash, bondSymbol(for: bonds[bondPath[index - 1]]))
                hash = appendCDKPathHash(hash, atomSymbol(for: atoms[atomPath[index]]))
            }
        }

        return hash
    }

    private func compareAtoms(_ lhs: Int, _ rhs: Int) -> Int {
        let lhsElement = atomicNumber(for: atoms[lhs])
        let rhsElement = atomicNumber(for: atoms[rhs])
        guard lhsElement != rhsElement else {
            return 0
        }
        return compareUTF16(atomSymbol(for: atoms[lhs]), atomSymbol(for: atoms[rhs]))
    }

    private func compareBonds(_ lhs: Int, _ rhs: Int) -> Int {
        compareUTF16(bondSymbol(for: bonds[lhs]), bondSymbol(for: bonds[rhs]))
    }

    private func skipAtom(_ atom: Atom) -> Bool {
        let element = atomicNumber(for: atom)
        return (!hashPseudoAtoms && element == 0) || (!hashExplicitHydrogens && element == 1)
    }

    private func atomSymbol(for atom: Atom) -> String {
        switch atomicNumber(for: atom) {
        case 0:
            return "*"
        case 6:
            return "C"
        case 7:
            return "N"
        case 8:
            return "O"
        case 17:
            return "X"
        case 35:
            return "Z"
        case 14:
            return "Y"
        case 33:
            return "D"
        case 3:
            return "L"
        case 34:
            return "E"
        case 11:
            return "G"
        case 20:
            return "J"
        case 13:
            return "A"
        default:
            return CDKDescriptorSupport.canonicalElementSymbol(atom.element)
        }
    }

    private func bondSymbol(for bond: Bond) -> String {
        if bond.order == .aromatic {
            return ":"
        }
        switch bond.order {
        case .single:
            return "-"
        case .double:
            return "="
        case .triple:
            return "#"
        case .aromatic:
            return ":"
        }
    }

    private func atomicNumber(for atom: Atom) -> Int {
        Self.atomicNumberBySymbol[CDKDescriptorSupport.canonicalElementSymbol(atom.element).uppercased()] ?? 0
    }

    private static let atomicSymbols: [String] = [
        "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
        "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
        "Rg", "Cn",
    ]

    private static let atomicNumberBySymbol: [String: Int] = {
        var values: [String: Int] = [:]
        for (index, symbol) in atomicSymbols.enumerated() where index > 0 {
            values[symbol.uppercased()] = index
        }
        return values
    }()

    private struct PathEdge: Hashable {
        let atomIndex: Int
        let bondIndex: Int
    }
}

private func appendCDKPathHash(_ currentHash: Int32, _ symbol: String) -> Int32 {
    var hash = currentHash
    let utf16 = Array(symbol.utf16)
    guard let first = utf16.first else {
        return hash
    }
    let firstValue = Int32(first)
    for _ in utf16 {
        hash = hash &* 31 &+ firstValue
    }
    return hash
}

private func javaStringHash(_ value: String) -> Int32 {
    var hash: Int32 = 0
    for codeUnit in value.utf16 {
        hash = hash &* 31 &+ Int32(codeUnit)
    }
    return hash
}

private func compareUTF16(_ lhs: String, _ rhs: String) -> Int {
    let left = Array(lhs.utf16)
    let right = Array(rhs.utf16)
    let count = min(left.count, right.count)
    for index in 0..<count {
        if left[index] < right[index] {
            return -1
        }
        if left[index] > right[index] {
            return 1
        }
    }
    if left.count < right.count {
        return -1
    }
    if left.count > right.count {
        return 1
    }
    return 0
}

private struct CDKJavaRandom {
    private static let multiplier: Int64 = 0x5DEEC_E66D
    private static let addend: Int64 = 0xB
    private static let mask: Int64 = (1 << 48) - 1
    private var seed: Int64

    init(seed: Int32) {
        let signedSeed = Int64(seed)
        self.seed = (signedSeed ^ Self.multiplier) & Self.mask
    }

    mutating func nextInt(bound: Int) -> Int {
        precondition(bound > 0, "Bound must be positive.")

        if (bound & -bound) == bound {
            return Int((Int64(bound) * Int64(next(bits: 31))) >> 31)
        }

        var bits: Int32
        var value: Int32
        repeat {
            bits = next(bits: 31)
            value = bits % Int32(bound)
        } while bits - value + (Int32(bound) - 1) < 0
        return Int(value)
    }

    private mutating func next(bits: Int) -> Int32 {
        seed = (seed &* Self.multiplier &+ Self.addend) & Self.mask
        return Int32(seed >> (48 - bits))
    }
}
