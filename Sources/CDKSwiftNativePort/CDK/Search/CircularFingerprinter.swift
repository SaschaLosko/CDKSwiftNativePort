import Foundation

public enum CDKCircularFingerprintClass: String, CaseIterable, Codable, Sendable {
    case ecfp0 = "ECFP0"
    case ecfp2 = "ECFP2"
    case ecfp4 = "ECFP4"
    case ecfp6 = "ECFP6"
    case fcfp0 = "FCFP0"
    case fcfp2 = "FCFP2"
    case fcfp4 = "FCFP4"
    case fcfp6 = "FCFP6"

    var iterationCount: Int {
        switch self {
        case .ecfp0, .fcfp0:
            return 0
        case .ecfp2, .fcfp2:
            return 1
        case .ecfp4, .fcfp4:
            return 2
        case .ecfp6, .fcfp6:
            return 3
        }
    }

    var usesFunctionalAtomClasses: Bool {
        switch self {
        case .fcfp0, .fcfp2, .fcfp4, .fcfp6:
            return true
        case .ecfp0, .ecfp2, .ecfp4, .ecfp6:
            return false
        }
    }
}

public struct CDKCircularFingerprintFeature: Hashable, Codable, Sendable {
    public let hashCode: Int32
    public let iteration: Int
    public let atomIDs: [Int]

    public init(hashCode: Int32, iteration: Int, atomIDs: [Int]) {
        self.hashCode = hashCode
        self.iteration = iteration
        self.atomIDs = atomIDs
    }
}

public struct CDKBitFingerprint: Hashable, Codable, Sendable {
    public let size: Int
    public let bits: Set<Int>

    public var cardinality: Int { bits.count }
    public var sortedBitIndices: [Int] { bits.sorted() }

    public init(size: Int, bits: Set<Int>) {
        precondition(size > 0, "Fingerprint size must be positive.")
        self.size = size
        self.bits = bits
    }

    public func contains(_ bit: Int) -> Bool {
        bits.contains(bit)
    }

    public func isSuperset(of other: CDKBitFingerprint) -> Bool {
        precondition(size == other.size, "Fingerprint sizes must match.")
        return bits.isSuperset(of: other.bits)
    }
}

public struct CDKCountFingerprint: Hashable, Codable, Sendable {
    public let countsByHash: [Int32: Int]

    public var populatedBinCount: Int { countsByHash.count }
    public var sortedBins: [(hash: Int32, count: Int)] {
        countsByHash.keys.sorted().map { hash in
            (hash, countsByHash[hash] ?? 0)
        }
    }

    public init(countsByHash: [Int32: Int]) {
        self.countsByHash = countsByHash.filter { $0.value > 0 }
    }

    public func count(forHash hash: Int32) -> Int {
        countsByHash[hash] ?? 0
    }
}

/// Swift counterpart of CDK's circular ECFP/FCFP fingerprinter.
///
/// The hash generation follows CDK's CRC32-based circular iteration so folded
/// bit positions and count fingerprints are stable across platforms for the
/// supported Swift `Molecule` model.
public final class CDKCircularFingerprinter {
    public let fingerprintClass: CDKCircularFingerprintClass
    public let foldedSize: Int

    public init(
        fingerprintClass: CDKCircularFingerprintClass = .ecfp6,
        foldedSize: Int = 1024
    ) {
        precondition(foldedSize > 0, "Fingerprint size must be positive.")
        self.fingerprintClass = fingerprintClass
        self.foldedSize = foldedSize
    }

    public func calculate(_ molecule: Molecule) -> [CDKCircularFingerprintFeature] {
        var calculator = CircularFingerprintCalculator(
            molecule: molecule,
            fingerprintClass: fingerprintClass)
        return calculator.calculate()
    }

    public func bitFingerprint(for molecule: Molecule) -> CDKBitFingerprint {
        let features = calculate(molecule)
        var bits = Set<Int>()
        bits.reserveCapacity(features.count)
        for feature in features {
            let unsigned = UInt32(bitPattern: feature.hashCode)
            bits.insert(Int(unsigned % UInt32(foldedSize)))
        }
        return CDKBitFingerprint(size: foldedSize, bits: bits)
    }

    public func countFingerprint(for molecule: Molecule) -> CDKCountFingerprint {
        let features = calculate(molecule)
        var counts: [Int32: Int] = [:]
        counts.reserveCapacity(features.count)
        for feature in features {
            counts[feature.hashCode, default: 0] += 1
        }
        return CDKCountFingerprint(countsByHash: counts)
    }
}

private struct CircularFingerprintCalculator {
    private let molecule: Molecule
    private let fingerprintClass: CDKCircularFingerprintClass
    private let atoms: [Atom]
    private let bonds: [Bond]
    private let atomIDToIndex: [Int: Int]

    private var heavyAtomMask: [Bool] = []
    private var hydrogenCounts: [Int] = []
    private var atomAdjacency: [[Int]] = []
    private var bondAdjacency: [[Int]] = []
    private var ringBlock: [Int] = []
    private var smallRings: [[Int]] = []
    private var bondOrder: [Int] = []
    private var atomAromatic: [Bool] = []
    private var bondAromatic: [Bool] = []
    private var tetrahedralRubrics: [[Int]?] = []

    private var donorMask: [Bool] = []
    private var acceptorMask: [Bool] = []
    private var positiveMask: [Bool] = []
    private var negativeMask: [Bool] = []
    private var biotypeAromaticMask: [Bool] = []
    private var halogenMask: [Bool] = []
    private var bondSum: [Int] = []
    private var hasDoubleBond: [Bool] = []
    private var aliphatic: [Bool] = []
    private var isOxide: [Bool] = []
    private var lonePair: [Bool] = []
    private var tetrazole: [Bool] = []

    private var identity: [Int32] = []
    private var resolvedChiral: [Bool] = []
    private var atomGroups: [[Int]] = []
    private var features: [InternalFeature] = []

    init(molecule: Molecule, fingerprintClass: CDKCircularFingerprintClass) {
        self.molecule = molecule
        self.fingerprintClass = fingerprintClass
        self.atoms = molecule.atoms
        self.bonds = molecule.bonds
        self.atomIDToIndex = Dictionary(
            uniqueKeysWithValues: molecule.atoms.enumerated().map { ($0.element.id, $0.offset) })
    }

    mutating func calculate() -> [CDKCircularFingerprintFeature] {
        features.removeAll(keepingCapacity: true)
        excavateMolecule()
        if fingerprintClass.usesFunctionalAtomClasses {
            calculateBioTypes()
        }

        let atomCount = atoms.count
        identity = Array(repeating: 0, count: atomCount)
        resolvedChiral = Array(repeating: false, count: atomCount)
        atomGroups = Array(repeating: [], count: atomCount)

        for atomIndex in 0..<atomCount where heavyAtomMask[atomIndex] {
            identity[atomIndex] =
                fingerprintClass.usesFunctionalAtomClasses
                ? initialIdentityFCFP(atomIndex)
                : initialIdentityECFP(atomIndex)
            atomGroups[atomIndex] = [atomIndex]
            features.append(
                InternalFeature(
                    hashCode: identity[atomIndex],
                    iteration: 0,
                    atomIndices: [atomIndex]))
        }

        if fingerprintClass.iterationCount > 0 {
            for iteration in 1...fingerprintClass.iterationCount {
                var newIdentity = Array(repeating: Int32(0), count: atomCount)
                for atomIndex in 0..<atomCount where heavyAtomMask[atomIndex] {
                    newIdentity[atomIndex] = circularIterate(iteration: iteration, atomIndex: atomIndex)
                }
                identity = newIdentity

                for atomIndex in 0..<atomCount where heavyAtomMask[atomIndex] {
                    atomGroups[atomIndex] = growAtoms(atomGroups[atomIndex])
                    considerNewFeature(
                        InternalFeature(
                            hashCode: identity[atomIndex],
                            iteration: iteration,
                            atomIndices: atomGroups[atomIndex]))
                }
            }
        }

        return features.map { feature in
            CDKCircularFingerprintFeature(
                hashCode: feature.hashCode,
                iteration: feature.iteration,
                atomIDs: feature.atomIndices.map { atoms[$0].id })
        }
    }

    private mutating func initialIdentityECFP(_ atomIndex: Int) -> Int32 {
        let atomicNumber = atomicNumber(for: atoms[atomIndex])
        let hydrogenCount = hydrogenCounts[atomIndex]
        let elementBonding = Self.elementBondingValue(forAtomicNumber: atomicNumber)
        let degree = elementBonding - hydrogenCount
        let inRing = ringBlock[atomIndex] > 0 ? 1 : 0

        var crc = CRC32()
        crc.updateByte((atomAdjacency[atomIndex].count << 4) | degree)
        crc.updateByte(atomicNumber)
        crc.updateByte(atoms[atomIndex].charge + 0x80)
        crc.updateByte((hydrogenCount << 4) | inRing)
        return Int32(bitPattern: crc.checksum)
    }

    private func initialIdentityFCFP(_ atomIndex: Int) -> Int32 {
        var identity = 0
        if donorMask[atomIndex] { identity |= 0x01 }
        if acceptorMask[atomIndex] { identity |= 0x02 }
        if positiveMask[atomIndex] { identity |= 0x04 }
        if negativeMask[atomIndex] { identity |= 0x08 }
        if atomAromatic[atomIndex] { identity |= 0x10 }
        if halogenMask[atomIndex] { identity |= 0x20 }
        return Int32(identity)
    }

    private mutating func circularIterate(iteration: Int, atomIndex: Int) -> Int32 {
        let adjacentAtoms = atomAdjacency[atomIndex]
        let adjacentBonds = bondAdjacency[atomIndex]
        var pairs: [(bondOrder: Int, identity: Int32)] = []
        pairs.reserveCapacity(adjacentAtoms.count)

        for offset in adjacentAtoms.indices {
            let bondIndex = adjacentBonds[offset]
            let order = bondAromatic[bondIndex] ? 0xF : bondOrder[bondIndex]
            pairs.append((order, identity[adjacentAtoms[offset]]))
        }

        pairs.sort {
            if $0.bondOrder != $1.bondOrder {
                return $0.bondOrder < $1.bondOrder
            }
            return $0.identity < $1.identity
        }

        var crc = CRC32()
        crc.updateByte(iteration)
        crc.updateInt32(identity[atomIndex])
        for pair in pairs {
            crc.updateByte(pair.bondOrder)
            crc.updateInt32(pair.identity)
        }

        if !resolvedChiral[atomIndex], let rubric = tetrahedralRubrics[atomIndex] {
            let priorities = rubric.map { index -> Int32 in
                index < 0 ? 0 : identity[index]
            }
            if Set(priorities).count == priorities.count {
                var parity = 0
                if priorities[0] < priorities[1] { parity += 1 }
                if priorities[0] < priorities[2] { parity += 1 }
                if priorities[0] < priorities[3] { parity += 1 }
                if priorities[1] < priorities[2] { parity += 1 }
                if priorities[1] < priorities[3] { parity += 1 }
                if priorities[2] < priorities[3] { parity += 1 }
                crc.updateByte((parity & 1) + 1)
                resolvedChiral[atomIndex] = true
            }
        }

        return Int32(bitPattern: crc.checksum)
    }

    private func growAtoms(_ atomIndices: [Int]) -> [Int] {
        var mask = Set(atomIndices)
        for index in atomIndices {
            for adjacent in atomAdjacency[index] {
                mask.insert(adjacent)
            }
        }
        return mask.sorted()
    }

    private mutating func considerNewFeature(_ newFeature: InternalFeature) {
        guard let hitIndex = features.firstIndex(where: { $0.atomIndices == newFeature.atomIndices }) else {
            features.append(newFeature)
            return
        }

        let existing = features[hitIndex]
        if existing.iteration < newFeature.iteration || existing.hashCode < newFeature.hashCode {
            return
        }
        features[hitIndex] = newFeature
    }

    private mutating func excavateMolecule() {
        let atomCount = atoms.count
        let bondCount = bonds.count

        heavyAtomMask = atoms.map { atomicNumber(for: $0) > 1 }
        atomAdjacency = Array(repeating: [], count: atomCount)
        bondAdjacency = Array(repeating: [], count: atomCount)
        bondOrder = Array(repeating: 0, count: bondCount)
        hydrogenCounts = Array(repeating: 0, count: atomCount)

        for bondIndex in bonds.indices {
            let bond = bonds[bondIndex]
            guard let first = atomIDToIndex[bond.a1], let second = atomIDToIndex[bond.a2] else { continue }
            bondOrder[bondIndex] = numericBondOrder(bond.order)

            if heavyAtomMask[first] && heavyAtomMask[second] {
                atomAdjacency[first].append(second)
                bondAdjacency[first].append(bondIndex)
                atomAdjacency[second].append(first)
                bondAdjacency[second].append(bondIndex)
            } else {
                if !heavyAtomMask[first], heavyAtomMask[second] { hydrogenCounts[second] += 1 }
                if !heavyAtomMask[second], heavyAtomMask[first] { hydrogenCounts[first] += 1 }
            }
        }

        for atomIndex in 0..<atomCount where heavyAtomMask[atomIndex] {
            let atom = atoms[atomIndex]
            let element = normalizedElement(atom).uppercased()
            var valence = 0
            switch element {
            case "C": valence = 4
            case "N": valence = 3
            case "O": valence = 2
            case "S": valence = 2
            case "P": valence = 3
            default: break
            }
            guard valence > 0 else { continue }
            var charge = atom.charge
            if element == "C" {
                charge = -abs(charge)
            }
            var implicit = Double(valence + charge)
            for bond in molecule.bonds(forAtom: atom.id) {
                implicit -= bond.order.valenceContribution
            }
            hydrogenCounts[atomIndex] += max(0, Int(implicit.rounded()))
        }

        detectRings()
        detectStrictAromaticity()
        detectTetrahedralRubrics()
    }

    private mutating func detectRings() {
        ringBlock = Array(repeating: 0, count: atoms.count)
        smallRings = []

        for ring in molecule.simpleCycles(maxSize: 7) {
            let indices = ring.compactMap { atomIDToIndex[$0] }
            guard indices.count == ring.count, indices.allSatisfy({ heavyAtomMask[$0] }) else { continue }
            smallRings.append(indices)
            for index in indices {
                ringBlock[index] = 1
            }
        }
    }

    private mutating func detectStrictAromaticity() {
        atomAromatic = atoms.map(\.aromatic)
        bondAromatic = bonds.map { $0.order == .aromatic }

        guard !smallRings.isEmpty else { return }

        var piAtom = Array(repeating: false, count: atoms.count)
        for bondIndex in bonds.indices where bondOrder[bondIndex] == 2 {
            guard let first = atomIDToIndex[bonds[bondIndex].a1],
                let second = atomIDToIndex[bonds[bondIndex].a2]
            else { continue }
            piAtom[first] = true
            piAtom[second] = true
        }

        var maybe = smallRings.filter { ring in
            guard ring.count == 6 else { return false }
            for offset in ring.indices {
                let first = ring[offset]
                let second = ring[(offset + 1) % ring.count]
                guard piAtom[first], let bondIndex = findBond(first, second) else {
                    return false
                }
                if bondOrder[bondIndex] != 1 && bondOrder[bondIndex] != 2 && !bondAromatic[bondIndex] {
                    return false
                }
            }
            return true
        }

        while true {
            var changed = false
            for candidateIndex in maybe.indices.reversed() {
                let ring = maybe[candidateIndex]
                var phaseOne = true
                var phaseTwo = true
                for offset in ring.indices {
                    guard let bondIndex = findBond(ring[offset], ring[(offset + 1) % ring.count]) else {
                        phaseOne = false
                        phaseTwo = false
                        break
                    }
                    if bondAromatic[bondIndex] { continue }
                    phaseOne = phaseOne && bondOrder[bondIndex] == (2 - (offset & 1))
                    phaseTwo = phaseTwo && bondOrder[bondIndex] == (1 + (offset & 1))
                }
                guard phaseOne || phaseTwo else { continue }

                for offset in ring.indices {
                    atomAromatic[ring[offset]] = true
                    if let bondIndex = findBond(ring[offset], ring[(offset + 1) % ring.count]) {
                        bondAromatic[bondIndex] = true
                    }
                }
                maybe.remove(at: candidateIndex)
                changed = true
            }
            if !changed { break }
        }
    }

    private mutating func detectTetrahedralRubrics() {
        tetrahedralRubrics = Array(repeating: nil, count: atoms.count)
        for atomIndex in atoms.indices {
            let atom = atoms[atomIndex]
            guard atom.chirality != .none else { continue }
            var rubric = atomAdjacency[atomIndex]
            guard rubric.count == 3 || rubric.count == 4 else { continue }
            if rubric.count == 3, hydrogenCounts[atomIndex] == 1 {
                rubric.append(-1)
            }
            guard rubric.count == 4 else { continue }
            if atom.chirality == .anticlockwise {
                rubric.swapAt(2, 3)
            }
            tetrahedralRubrics[atomIndex] = rubric
        }
    }

    private mutating func calculateBioTypes() {
        let atomCount = atoms.count
        donorMask = Array(repeating: false, count: atomCount)
        acceptorMask = Array(repeating: false, count: atomCount)
        positiveMask = Array(repeating: false, count: atomCount)
        negativeMask = Array(repeating: false, count: atomCount)
        biotypeAromaticMask = Array(repeating: false, count: atomCount)
        halogenMask = Array(repeating: false, count: atomCount)
        aliphatic = Array(repeating: false, count: atomCount)
        bondSum = Array(repeating: 0, count: atomCount)

        for atomIndex in atoms.indices where heavyAtomMask[atomIndex] {
            aliphatic[atomIndex] = atomicNumber(for: atoms[atomIndex]) == 6
            bondSum[atomIndex] = hydrogenCounts[atomIndex]
        }

        hasDoubleBond = Array(repeating: false, count: atomCount)
        isOxide = Array(repeating: false, count: atomCount)
        for bondIndex in bonds.indices {
            let bond = bonds[bondIndex]
            guard let first = atomIDToIndex[bond.a1], let second = atomIDToIndex[bond.a2],
                heavyAtomMask[first], heavyAtomMask[second]
            else { continue }

            let order = bondOrder[bondIndex]
            bondSum[first] += order
            bondSum[second] += order
            if order == 2 {
                hasDoubleBond[first] = true
                hasDoubleBond[second] = true
                if atomicNumber(for: atoms[first]) == 8 { isOxide[second] = true }
                if atomicNumber(for: atoms[second]) == 8 { isOxide[first] = true }
            }
            if order != 1 {
                aliphatic[first] = false
                aliphatic[second] = false
            }
        }

        lonePair = Array(repeating: false, count: atomCount)
        for atomIndex in atoms.indices {
            let element = normalizedElement(atoms[atomIndex]).uppercased()
            let valence: Int
            switch element {
            case "N": valence = 3
            case "O", "S": valence = 2
            default: valence = 0
            }
            if valence > 0 && bondSum[atomIndex] + atoms[atomIndex].charge <= valence {
                lonePair[atomIndex] = true
            }
        }

        tetrazole = Array(repeating: false, count: atomCount)
        for ring in smallRings where ring.count >= 5 && ring.count <= 7 {
            considerBioTypeAromaticity(ring)
            if ring.count == 5 {
                considerBioTypeTetrazole(ring)
            }
        }

        for atomIndex in atoms.indices where heavyAtomMask[atomIndex] {
            donorMask[atomIndex] = determineDonor(atomIndex)
            acceptorMask[atomIndex] = determineAcceptor(atomIndex)
            positiveMask[atomIndex] = determinePositive(atomIndex)
            negativeMask[atomIndex] = determineNegative(atomIndex)
            halogenMask[atomIndex] = determineHalide(atomIndex)
        }
    }

    private mutating func considerBioTypeAromaticity(_ ring: [Int]) {
        var doubleCount = 0
        for atomIndex in ring {
            if hasDoubleBond[atomIndex] {
                doubleCount += 1
            } else if !lonePair[atomIndex] {
                return
            }
        }
        guard doubleCount >= ring.count - 2 else { return }
        for atomIndex in ring {
            biotypeAromaticMask[atomIndex] = true
        }
    }

    private mutating func considerBioTypeTetrazole(_ ring: [Int]) {
        var carbonCount = 0
        var nitrogenCount = 0
        var doubleCount = 0
        for offset in ring.indices {
            let atomIndex = ring[offset]
            guard atoms[atomIndex].charge == 0 else { return }
            switch normalizedElement(atoms[atomIndex]).uppercased() {
            case "C": carbonCount += 1
            case "N": nitrogenCount += 1
            default: break
            }
            if let bondIndex = findBond(atomIndex, ring[(offset + 1) % ring.count]), bondOrder[bondIndex] == 2 {
                doubleCount += 1
            }
        }
        guard carbonCount == 1, nitrogenCount == 4, doubleCount == 2 else { return }
        for atomIndex in ring where atomicNumber(for: atoms[atomIndex]) == 7 {
            tetrazole[atomIndex] = true
        }
    }

    private func determineDonor(_ atomIndex: Int) -> Bool {
        guard hydrogenCounts[atomIndex] > 0 else { return false }
        let element = normalizedElement(atoms[atomIndex]).uppercased()
        switch element {
        case "N", "O":
            if tetrazole[atomIndex] { return false }
            for neighbor in atomAdjacency[atomIndex] where isOxide[neighbor] {
                if atomicNumber(for: atoms[neighbor]) != 6 || element != "N" {
                    return false
                }
            }
            return true
        case "S":
            for neighbor in atomAdjacency[atomIndex] where hasDoubleBond[neighbor] {
                return false
            }
            return true
        case "C":
            return bondAdjacency[atomIndex].contains { bondOrderBioType($0) == 3 }
        default:
            return false
        }
    }

    private func determineAcceptor(_ atomIndex: Int) -> Bool {
        guard lonePair[atomIndex], atoms[atomIndex].charge <= 0 else { return false }
        if atomicNumber(for: atoms[atomIndex]) == 7 {
            let basic = atomAdjacency[atomIndex].allSatisfy { aliphatic[$0] }
            if basic { return false }
        }
        return true
    }

    private func determinePositive(_ atomIndex: Int) -> Bool {
        let atom = atoms[atomIndex]
        if atom.charge < 0 { return false }
        if atom.charge > 0 {
            return !atomAdjacency[atomIndex].contains { atoms[$0].charge < 0 }
        }

        let element = normalizedElement(atom).uppercased()
        if element == "N" {
            if atomAdjacency[atomIndex].allSatisfy({ aliphatic[$0] }) {
                return true
            }

            if hasDoubleBond[atomIndex] && hydrogenCounts[atomIndex] == 0,
                let other = bondAdjacency[atomIndex].first(where: { bondOrderBioType($0) == 2 })
                    .flatMap({ otherAtomIndex(for: bonds[$0], from: atomIndex) })
            {
                var amines = 0
                for neighbor in atomAdjacency[other] where neighbor != atomIndex {
                    if hasDoubleBond[neighbor] { return false }
                    let neighborElement = normalizedElement(atoms[neighbor]).uppercased()
                    if neighborElement == "N" {
                        if hydrogenCounts[neighbor] > 0 { return false }
                        amines += 1
                    } else if neighborElement != "C" {
                        return false
                    }
                }
                if amines > 0 { return true }
            }
        } else if element == "C" {
            var imine = false
            var amine = false
            for offset in atomAdjacency[atomIndex].indices {
                let neighbor = atomAdjacency[atomIndex][offset]
                if tetrazole[neighbor] { return false }
                guard atomicNumber(for: atoms[neighbor]) == 7 else { continue }
                if bondOrderBioType(bondAdjacency[atomIndex][offset]) == 2 {
                    imine = true
                } else if hydrogenCounts[neighbor] == 1 {
                    amine = true
                }
            }
            return imine && amine
        }

        return false
    }

    private func determineNegative(_ atomIndex: Int) -> Bool {
        let atom = atoms[atomIndex]
        if atom.charge > 0 { return false }
        if atom.charge < 0 {
            return !atomAdjacency[atomIndex].contains { atoms[$0].charge > 0 }
        }

        let element = normalizedElement(atom).uppercased()
        if tetrazole[atomIndex], element == "N" {
            return true
        }

        if isOxide[atomIndex], element == "C" || element == "S" || element == "P" {
            for offset in atomAdjacency[atomIndex].indices where bondOrderBioType(bondAdjacency[atomIndex][offset]) == 1
            {
                let neighbor = atomAdjacency[atomIndex][offset]
                if atomicNumber(for: atoms[neighbor]) == 8, hydrogenCounts[neighbor] > 0 {
                    return true
                }
            }
        }

        return false
    }

    private func determineHalide(_ atomIndex: Int) -> Bool {
        switch normalizedElement(atoms[atomIndex]).uppercased() {
        case "F", "CL", "BR", "I":
            return true
        default:
            return false
        }
    }

    private func bondOrderBioType(_ bondIndex: Int) -> Int {
        guard let first = atomIDToIndex[bonds[bondIndex].a1],
            let second = atomIDToIndex[bonds[bondIndex].a2]
        else { return 0 }
        if biotypeAromaticMask[first] && biotypeAromaticMask[second] {
            return -1
        }
        return bondOrder[bondIndex]
    }

    private func findBond(_ first: Int, _ second: Int) -> Int? {
        for offset in atomAdjacency[first].indices where atomAdjacency[first][offset] == second {
            return bondAdjacency[first][offset]
        }
        return nil
    }

    private func otherAtomIndex(for bond: Bond, from atomIndex: Int) -> Int? {
        let atomID = atoms[atomIndex].id
        if bond.a1 == atomID { return atomIDToIndex[bond.a2] }
        if bond.a2 == atomID { return atomIDToIndex[bond.a1] }
        return nil
    }

    private func normalizedElement(_ atom: Atom) -> String {
        CDKDescriptorSupport.canonicalElementSymbol(atom.element)
    }

    private func atomicNumber(for atom: Atom) -> Int {
        Self.atomicNumberBySymbol[normalizedElement(atom).uppercased()] ?? 0
    }

    private func numericBondOrder(_ order: BondOrder) -> Int {
        switch order {
        case .single: return 1
        case .double: return 2
        case .triple: return 3
        case .aromatic: return 1
        }
    }

    private static func elementBondingValue(forAtomicNumber atomicNumber: Int) -> Int {
        guard atomicNumber > 0, atomicNumber < elementBonding.count else { return 0 }
        return elementBonding[atomicNumber]
    }

    private struct InternalFeature: Hashable {
        let hashCode: Int32
        let iteration: Int
        let atomIndices: [Int]
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

    private static let elementBonding: [Int] = [
        0, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2,
        3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 3, 4, 3, 2, 1, 0, 1, 2, 3, 4,
        5, 6, 7, 8, 9, 10, 11, 12, 3, 4, 3, 2, 1, 0, 1, 2, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
        3, 4, 5, 6, 7, 8, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
    ]
}

private struct CRC32 {
    private var state: UInt32 = 0xFFFF_FFFF

    var checksum: UInt32 {
        state ^ 0xFFFF_FFFF
    }

    mutating func updateByte(_ value: Int) {
        let byte = UInt8(truncatingIfNeeded: value)
        let index = Int((state ^ UInt32(byte)) & 0xFF)
        state = (state >> 8) ^ Self.table[index]
    }

    mutating func updateInt32(_ value: Int32) {
        let unsigned = UInt32(bitPattern: value)
        updateByte(Int((unsigned >> 24) & 0xFF))
        updateByte(Int((unsigned >> 16) & 0xFF))
        updateByte(Int((unsigned >> 8) & 0xFF))
        updateByte(Int(unsigned & 0xFF))
    }

    private static let table: [UInt32] = {
        (0..<256).map { index in
            var value = UInt32(index)
            for _ in 0..<8 {
                if value & 1 == 1 {
                    value = 0xEDB8_8320 ^ (value >> 1)
                } else {
                    value >>= 1
                }
            }
            return value
        }
    }()
}
