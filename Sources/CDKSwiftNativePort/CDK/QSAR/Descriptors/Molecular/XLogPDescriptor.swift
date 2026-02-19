import Foundation

/// CDK port of `org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor`.
///
/// Notes:
/// - This path mirrors CDK 2.11 atom-fragment and correction-factor logic.
/// - The descriptor builds a local explicit-hydrogen working graph by expanding
///   implicit hydrogens into virtual nodes, matching CDK's "explicit H required"
///   execution model.
public enum CDKXLogPDescriptor {
    public static func calculate(for molecule: Molecule,
                                 checkAromaticity: Bool = true,
                                 salicylFlag: Bool = false) -> Double? {
        guard !molecule.atoms.isEmpty else { return nil }
        let context = XLogPContext(molecule: molecule, checkAromaticity: checkAromaticity)
        return context.calculate(salicylFlag: salicylFlag)
    }
}

private struct XLogPContext {
    private struct Node {
        let id: Int
        let symbol: String
        let charge: Int
        let aromatic: Bool
        let isVirtualHydrogen: Bool
    }

    private enum MaxBondOrder {
        case single
        case double
        case triple
    }

    private let sourceMolecule: Molecule
    private let nodes: [Int: Node]
    private let nodeOrder: [Int]
    private let neighborsByID: [Int: [Int]]
    private let bondOrderByKey: [Int64: BondOrder]
    private let ringSizesByAtomID: [Int: Int]
    private let aromaticRingMembership: Set<Int>
    private let aromaticAtomIDs: Set<Int>

    init(molecule: Molecule, checkAromaticity: Bool) {
        sourceMolecule = molecule

        var aromaticIDs = Set<Int>()
        aromaticIDs.formUnion(molecule.atoms.filter(\.aromatic).map(\.id))
        if checkAromaticity {
            aromaticIDs.formUnion(molecule.aromaticDisplayRings().flatMap { $0 })
        }
        aromaticAtomIDs = aromaticIDs

        let cycles = molecule.simpleCycles(maxSize: 12)
        var ringSizes: [Int: Int] = [:]
        var aromaticRingSet = Set<Int>()
        for cycle in cycles {
            let size = cycle.count
            for atomID in cycle {
                ringSizes[atomID] = min(ringSizes[atomID] ?? Int.max, size)
                if size >= 6, aromaticAtomIDs.contains(atomID) {
                    aromaticRingSet.insert(atomID)
                }
            }
        }
        ringSizesByAtomID = ringSizes
        aromaticRingMembership = aromaticRingSet

        var nextVirtualID = -1
        var builtNodes: [Int: Node] = [:]
        var order: [Int] = []
        var adjacency: [Int: [Int]] = [:]
        var bondOrders: [Int64: BondOrder] = [:]

        for atom in molecule.atoms {
            let canonical = CDKDescriptorSupport.canonicalElementSymbol(atom.element)
            let aromatic = aromaticAtomIDs.contains(atom.id)
            builtNodes[atom.id] = Node(id: atom.id,
                                       symbol: canonical,
                                       charge: atom.charge,
                                       aromatic: aromatic,
                                       isVirtualHydrogen: false)
            order.append(atom.id)
            adjacency[atom.id] = adjacency[atom.id] ?? []
        }

        for bond in molecule.bonds {
            let key = Self.bondKey(bond.a1, bond.a2)
            bondOrders[key] = bond.order
            adjacency[bond.a1, default: []].append(bond.a2)
            adjacency[bond.a2, default: []].append(bond.a1)
        }

        for atom in molecule.atoms {
            let canonical = CDKDescriptorSupport.canonicalElementSymbol(atom.element)
            guard canonical.uppercased() != "H" else { continue }
            let implicitHydrogens = max(0, CDKDescriptorSupport.implicitHydrogenCount(on: atom.id, in: molecule))
            guard implicitHydrogens > 0 else { continue }
            for _ in 0..<implicitHydrogens {
                let virtualID = nextVirtualID
                nextVirtualID -= 1

                builtNodes[virtualID] = Node(id: virtualID,
                                             symbol: "H",
                                             charge: 0,
                                             aromatic: false,
                                             isVirtualHydrogen: true)
                order.append(virtualID)

                adjacency[virtualID] = [atom.id]
                adjacency[atom.id, default: []].append(virtualID)
                bondOrders[Self.bondKey(virtualID, atom.id)] = .single
            }
        }

        nodes = builtNodes
        nodeOrder = order
        neighborsByID = adjacency
        bondOrderByKey = bondOrders
    }

    func calculate(salicylFlag: Bool) -> Double {
        var xlogP = 0.0
        var previous = 0.0
        var previousSymbol = "H"
        var hBondAcceptors: [Int] = []
        var hBondDonors: [Int] = []

        for atomID in nodeOrder {
            let symbol = canonicalSymbol(atomID)
            if previous == xlogP, !previousSymbol.isEmpty, previousSymbol != "H" {
                // Keep CDK port behavior where assignment mismatches are tolerated.
            }
            previous = xlogP
            previousSymbol = symbol

            let bondCount = connectedBondsCount(atomID)
            let hsCount = hydrogenCount(atomID)
            let maxBondOrder = maximumBondOrder(atomID)

            if symbol == "C" {
                if bondCount == 2 {
                    if hsCount >= 1 {
                        xlogP += 0.209
                    } else if maxBondOrder == .double {
                        xlogP += 2.073
                    } else if maxBondOrder == .triple {
                        xlogP += 0.33
                    }
                }

                if bondCount == 3 {
                    if isInAromaticRing(atomID) {
                        if aromaticCarbonsCount(atomID) >= 2, aromaticNitrogensCount(atomID) == 0 {
                            if hsCount == 0 {
                                if atomTypeXCount(atomID) == 0 {
                                    xlogP += 0.296
                                } else {
                                    xlogP -= 0.151
                                }
                            } else {
                                xlogP += 0.337
                            }
                        } else if aromaticNitrogensCount(atomID) >= 1 {
                            if hsCount == 0 {
                                if atomTypeXCount(atomID) == 0 {
                                    xlogP += 0.174
                                } else {
                                    xlogP += 0.366
                                }
                            } else if hsCount == 1 {
                                xlogP += 0.126
                            }
                        }
                    } else {
                        if hsCount == 0 {
                            if atomTypeXCount(atomID) == 0 {
                                if piSystemsCount(atomID) <= 1 {
                                    xlogP += 0.05
                                } else {
                                    xlogP += 0.013
                                }
                            } else if atomTypeXCount(atomID) == 1 {
                                if piSystemsCount(atomID) == 0 {
                                    xlogP -= 0.03
                                } else {
                                    xlogP -= 0.027
                                }
                            } else if atomTypeXCount(atomID) == 2 {
                                if piSystemsCount(atomID) == 0 {
                                    xlogP += 0.005
                                } else {
                                    xlogP -= 0.315
                                }
                            }
                        }

                        if hsCount == 1 {
                            if atomTypeXCount(atomID) == 0 {
                                if piSystemsCount(atomID) == 0 {
                                    xlogP += 0.466
                                }
                                if piSystemsCount(atomID) == 1 {
                                    xlogP += 0.136
                                }
                            } else {
                                if piSystemsCount(atomID) == 0 {
                                    xlogP += 0.001
                                }
                                if piSystemsCount(atomID) == 1 {
                                    xlogP -= 0.31
                                }
                            }
                        }

                        if hsCount == 2 {
                            xlogP += 0.42
                        }
                        if ifCarbonIsHydrophobic(atomID) {
                            xlogP += 0.211
                        }
                    }
                }

                if bondCount == 4 {
                    if hsCount == 0 {
                        if atomTypeXCount(atomID) == 0 {
                            if piSystemsCount(atomID) == 0 {
                                xlogP -= 0.006
                            }
                            if piSystemsCount(atomID) == 1 {
                                xlogP -= 0.57
                            }
                            if piSystemsCount(atomID) >= 2 {
                                xlogP -= 0.317
                            }
                        } else {
                            if piSystemsCount(atomID) == 0 {
                                xlogP -= 0.316
                            } else {
                                xlogP -= 0.723
                            }
                        }
                    }

                    if hsCount == 1 {
                        if atomTypeXCount(atomID) == 0 {
                            if piSystemsCount(atomID) == 0 {
                                xlogP += 0.127
                            }
                            if piSystemsCount(atomID) == 1 {
                                xlogP -= 0.243
                            }
                            if piSystemsCount(atomID) >= 2 {
                                xlogP -= 0.499
                            }
                        } else {
                            if piSystemsCount(atomID) == 0 {
                                xlogP -= 0.205
                            }
                            if piSystemsCount(atomID) == 1 {
                                xlogP -= 0.305
                            }
                            if piSystemsCount(atomID) >= 2 {
                                xlogP -= 0.709
                            }
                        }
                    }

                    if hsCount == 2 {
                        if atomTypeXCount(atomID) == 0 {
                            if piSystemsCount(atomID) == 0 {
                                xlogP += 0.358
                            }
                            if piSystemsCount(atomID) == 1 {
                                xlogP -= 0.008
                            }
                            if piSystemsCount(atomID) == 2 {
                                xlogP -= 0.185
                            }
                        } else {
                            if piSystemsCount(atomID) == 0 {
                                xlogP -= 0.137
                            }
                            if piSystemsCount(atomID) == 1 {
                                xlogP -= 0.303
                            }
                            if piSystemsCount(atomID) == 2 {
                                xlogP -= 0.815
                            }
                        }
                    }

                    if hsCount > 2 {
                        if atomTypeXCount(atomID) == 0 {
                            if piSystemsCount(atomID) == 0 {
                                xlogP += 0.528
                            }
                            if piSystemsCount(atomID) == 1 {
                                xlogP += 0.267
                            }
                        } else {
                            xlogP -= 0.032
                        }
                    }
                    if ifCarbonIsHydrophobic(atomID) {
                        xlogP += 0.211
                    }
                }
            }

            if symbol == "N" {
                if bondOrderSum(atomID) >= 3.0, oxygenCount(atomID) >= 2, maxBondOrder == .double {
                    xlogP += 1.178
                } else if presenceOfCarbonyl(atomID) >= 1 {
                    if hsCount == 0 {
                        if atomTypeXCount(atomID) == 0 {
                            xlogP += 0.078
                        }
                        if atomTypeXCount(atomID) == 1 {
                            xlogP -= 0.118
                        }
                    }
                    if hsCount == 1 {
                        if atomTypeXCount(atomID) == 0 {
                            xlogP -= 0.096
                            hBondDonors.append(atomID)
                        } else {
                            xlogP -= 0.044
                            hBondDonors.append(atomID)
                        }
                    }
                    if hsCount == 2 {
                        xlogP -= 0.646
                        hBondDonors.append(atomID)
                    }
                } else {
                    if bondCount == 1 {
                        if carbonsCount(atomID) == 1 {
                            xlogP -= 0.566
                        }
                    } else if bondCount == 2 {
                        if isInAromaticRing(atomID) {
                            xlogP -= 0.493
                        } else if doubleBondedCarbonsCount(atomID) == 0 {
                            if doubleBondedNitrogensCount(atomID) == 0 {
                                if doubleBondedOxygensCount(atomID) == 1 {
                                    xlogP += 0.427
                                }
                            }
                            if doubleBondedNitrogensCount(atomID) == 1 {
                                if atomTypeXCount(atomID) == 0 {
                                    xlogP += 0.536
                                }
                                if atomTypeXCount(atomID) == 1 {
                                    xlogP -= 0.597
                                }
                            }
                        } else if doubleBondedCarbonsCount(atomID) == 1 {
                            if atomTypeXCount(atomID) == 0 {
                                if piSystemsCount(atomID) == 0 {
                                    xlogP += 0.007
                                }
                                if piSystemsCount(atomID) == 1 {
                                    xlogP -= 0.275
                                }
                            } else if atomTypeXCount(atomID) == 1 {
                                if piSystemsCount(atomID) == 0 {
                                    xlogP += 0.366
                                }
                                if piSystemsCount(atomID) == 1 {
                                    xlogP += 0.251
                                }
                            }
                        }
                    } else if bondCount == 3 {
                        if hsCount == 0 {
                            if isAtomAromatic(atomID) || (ringContains(atomID) && smallestRingSize(atomID) > 3 && piSystemsCount(atomID) >= 1) {
                                if atomTypeXCount(atomID) == 0 {
                                    xlogP += 0.881
                                } else {
                                    xlogP -= 0.01
                                }
                            } else if atomTypeXCount(atomID) == 0 {
                                if piSystemsCount(atomID) == 0 {
                                    xlogP += 0.159
                                }
                                if piSystemsCount(atomID) > 0 {
                                    xlogP += 0.761
                                }
                            } else {
                                xlogP -= 0.239
                            }
                        } else if hsCount == 1 {
                            if atomTypeXCount(atomID) == 0 {
                                if isAtomAromatic(atomID) || (ringContains(atomID) && smallestRingSize(atomID) > 3 && piSystemsCount(atomID) >= 2) {
                                    xlogP += 0.545
                                    hBondDonors.append(atomID)
                                } else {
                                    if piSystemsCount(atomID) == 0 {
                                        xlogP -= 0.112
                                        hBondDonors.append(atomID)
                                    }
                                    if piSystemsCount(atomID) > 0 {
                                        xlogP += 0.166
                                        hBondDonors.append(atomID)
                                    }
                                }
                            } else if ringContains(atomID) {
                                xlogP += 0.153
                                hBondDonors.append(atomID)
                            } else {
                                xlogP += 0.324
                                hBondDonors.append(atomID)
                            }
                        } else if hsCount == 2 {
                            if atomTypeXCount(atomID) == 0 {
                                if piSystemsCount(atomID) == 0 {
                                    xlogP -= 0.534
                                    hBondDonors.append(atomID)
                                }
                                if piSystemsCount(atomID) == 1 {
                                    xlogP -= 0.329
                                    hBondDonors.append(atomID)
                                }
                            } else {
                                xlogP -= 1.082
                                hBondDonors.append(atomID)
                            }
                        }
                    }
                }
            }

            if symbol == "O" {
                if bondCount == 1, maxBondOrder == .double {
                    xlogP -= 0.399
                    if !presenceOfHydroxy(atomID) {
                        hBondAcceptors.append(atomID)
                    }
                } else if (bondCount == 1 && hsCount == 0 && (presenceOfNitro(atomID) || presenceOfCarbonyl(atomID) == 1)) || presenceOfSulfat(atomID) {
                    xlogP -= 0.399
                    if !presenceOfHydroxy(atomID) {
                        hBondAcceptors.append(atomID)
                    }
                } else if bondCount >= 1 {
                    if hsCount == 0, bondCount == 2 {
                        if atomTypeXCount(atomID) == 0 {
                            if piSystemsCount(atomID) == 0 {
                                xlogP += 0.084
                            }
                            if piSystemsCount(atomID) > 0 {
                                xlogP += 0.435
                            }
                        } else if atomTypeXCount(atomID) == 1 {
                            xlogP += 0.105
                        }
                    } else if atomTypeXCount(atomID) == 0 {
                        if piSystemsCount(atomID) == 0 {
                            xlogP -= 0.467
                            hBondDonors.append(atomID)
                            hBondAcceptors.append(atomID)
                        }
                        if piSystemsCount(atomID) == 1 {
                            xlogP += 0.082
                            hBondDonors.append(atomID)
                            hBondAcceptors.append(atomID)
                        }
                    } else if atomTypeXCount(atomID) == 1 {
                        xlogP -= 0.522
                        hBondDonors.append(atomID)
                        hBondAcceptors.append(atomID)
                    }
                }
            }

            if symbol == "S" {
                if (bondCount == 1 && maxBondOrder == .double) || (bondCount == 1 && charge(atomID) == -1) {
                    xlogP -= 0.148
                } else if bondCount == 2 {
                    if hsCount == 0 {
                        xlogP += 0.255
                    } else {
                        xlogP += 0.419
                    }
                } else if bondCount == 3 {
                    if oxygenCount(atomID) >= 1 {
                        xlogP -= 1.375
                    }
                } else if bondCount == 4 {
                    if doubleBondedOxygensCount(atomID) >= 2 {
                        xlogP -= 0.168
                    }
                }
            }

            if symbol == "P" {
                if doubleBondedSulfurCount(atomID) >= 1, bondCount >= 4 {
                    xlogP += 1.253
                } else if (oxygenCount(atomID) >= 1 || doubleBondedOxygensCount(atomID) == 1), bondCount >= 4 {
                    xlogP -= 0.447
                }
            }

            if symbol == "F" {
                if piSystemsCount(atomID) == 0 {
                    xlogP += 0.375
                } else if piSystemsCount(atomID) == 1 {
                    xlogP += 0.202
                }
            }
            if symbol == "Cl" {
                if piSystemsCount(atomID) == 0 {
                    xlogP += 0.512
                } else if piSystemsCount(atomID) >= 1 {
                    xlogP += 0.663
                }
            }
            if symbol == "Br" {
                if piSystemsCount(atomID) == 0 {
                    xlogP += 0.85
                } else if piSystemsCount(atomID) == 1 {
                    xlogP += 0.839
                }
            }
            if symbol == "I" {
                if piSystemsCount(atomID) == 0 {
                    xlogP += 1.05
                } else if piSystemsCount(atomID) == 1 {
                    xlogP += 1.109
                }
            }

            let halCount = halogenCount(atomID)
            if halCount == 2 {
                xlogP += 0.137
            } else if halCount == 3 {
                xlogP += 3 * 0.137
            } else if halCount == 4 {
                xlogP += 6 * 0.137
            }

            if presenceOfCarbonyl(atomID) == 2, !ringContains(atomID) {
                xlogP += 0.580
            }
        }

        if !hBondAcceptors.isEmpty, !hBondDonors.isEmpty {
            var seenPairs = Set<Int64>()
            for acceptor in hBondAcceptors {
                for donor in hBondDonors {
                    guard checkRingLink(acceptor) || checkRingLink(donor) else { continue }
                    guard let distance = shortestPathDistance(from: acceptor, to: donor) else { continue }
                    let pairKey = Self.bondKey(acceptor, donor)
                    guard !seenPairs.contains(pairKey) else { continue }

                    let bothRingLinked = checkRingLink(acceptor) && checkRingLink(donor)
                    if bothRingLinked && distance == 3 {
                        xlogP += 0.429
                        seenPairs.insert(pairKey)
                    } else if !bothRingLinked && distance == 4 {
                        xlogP += 0.429
                        seenPairs.insert(pairKey)
                    }
                }
            }
        }

        if matchesTerminalAminoAcid() {
            xlogP -= 2.166
        }
        if matchesPAminoSulphonicAcid() {
            xlogP -= 0.501
        }
        if salicylFlag, matchesSalicylicAcid() {
            xlogP += 0.554
        }
        if hasOrthoOxygenPair() {
            xlogP -= 0.268
        }

        return xlogP
    }

    private func connectedBondsCount(_ atomID: Int) -> Int {
        neighborsByID[atomID]?.count ?? 0
    }

    private func neighbors(of atomID: Int) -> [Int] {
        neighborsByID[atomID] ?? []
    }

    private func canonicalSymbol(_ atomID: Int) -> String {
        nodes[atomID]?.symbol ?? ""
    }

    private func isAtomAromatic(_ atomID: Int) -> Bool {
        nodes[atomID]?.aromatic ?? false
    }

    private func charge(_ atomID: Int) -> Int {
        nodes[atomID]?.charge ?? 0
    }

    private func bondOrder(_ a: Int, _ b: Int) -> BondOrder? {
        bondOrderByKey[Self.bondKey(a, b)]
    }

    private func bondOrderSum(_ atomID: Int) -> Double {
        neighbors(of: atomID).reduce(0.0) { partial, neighbor in
            guard let order = bondOrder(atomID, neighbor) else { return partial }
            switch order {
            case .single:
                return partial + 1.0
            case .double:
                return partial + 2.0
            case .triple:
                return partial + 3.0
            case .aromatic:
                return partial + 1.5
            }
        }
    }

    private func maximumBondOrder(_ atomID: Int) -> MaxBondOrder {
        let rank = neighbors(of: atomID).reduce(1) { current, neighbor in
            guard let order = bondOrder(atomID, neighbor) else { return current }
            let nextRank: Int
            switch order {
            case .single:
                nextRank = 1
            case .double, .aromatic:
                nextRank = 2
            case .triple:
                nextRank = 3
            }
            return max(current, nextRank)
        }
        switch rank {
        case 3:
            return .triple
        case 2:
            return .double
        default:
            return .single
        }
    }

    private func hydrogenCount(_ atomID: Int) -> Int {
        neighbors(of: atomID).reduce(0) { partial, neighbor in
            partial + (canonicalSymbol(neighbor) == "H" ? 1 : 0)
        }
    }

    private func halogenCount(_ atomID: Int) -> Int {
        neighbors(of: atomID).reduce(0) { partial, neighbor in
            switch canonicalSymbol(neighbor).uppercased() {
            case "F", "CL", "BR", "I":
                return partial + 1
            default:
                return partial
            }
        }
    }

    private func atomTypeXCount(_ atomID: Int) -> Int {
        neighbors(of: atomID).reduce(0) { partial, neighbor in
            let symbol = canonicalSymbol(neighbor).uppercased()
            guard (symbol == "N" || symbol == "O"),
                  !isInAromaticRing(neighbor),
                  bondOrder(atomID, neighbor) != .double else {
                return partial
            }
            return partial + 1
        }
    }

    private func aromaticCarbonsCount(_ atomID: Int) -> Int {
        neighbors(of: atomID).reduce(0) { partial, neighbor in
            let isAromaticCarbon = canonicalSymbol(neighbor).uppercased() == "C" && isAtomAromatic(neighbor)
            return partial + (isAromaticCarbon ? 1 : 0)
        }
    }

    private func carbonsCount(_ atomID: Int) -> Int {
        neighbors(of: atomID).reduce(0) { partial, neighbor in
            guard canonicalSymbol(neighbor).uppercased() == "C", !isAtomAromatic(neighbor) else {
                return partial
            }
            return partial + 1
        }
    }

    private func oxygenCount(_ atomID: Int) -> Int {
        neighbors(of: atomID).reduce(0) { partial, neighbor in
            guard canonicalSymbol(neighbor).uppercased() == "O", !isAtomAromatic(neighbor) else {
                return partial
            }
            return partial + 1
        }
    }

    private func doubleBondedCarbonsCount(_ atomID: Int) -> Int {
        neighbors(of: atomID).reduce(0) { partial, neighbor in
            guard canonicalSymbol(neighbor).uppercased() == "C",
                  bondOrder(atomID, neighbor) == .double else {
                return partial
            }
            return partial + 1
        }
    }

    private func doubleBondedOxygensCount(_ atomID: Int) -> Int {
        let positiveCenter = charge(atomID) >= 1
        return neighbors(of: atomID).reduce(0) { partial, neighbor in
            guard canonicalSymbol(neighbor).uppercased() == "O" else {
                return partial
            }
            let order = bondOrder(atomID, neighbor)
            let neighborCharge = charge(neighbor)

            if positiveCenter, neighborCharge == -1, order == .single {
                return partial + 1
            }
            if !isAtomAromatic(neighbor), order == .double {
                return partial + 1
            }
            return partial
        }
    }

    private func doubleBondedSulfurCount(_ atomID: Int) -> Int {
        neighbors(of: atomID).reduce(0) { partial, neighbor in
            guard canonicalSymbol(neighbor).uppercased() == "S" else {
                return partial
            }
            if charge(atomID) == 1, charge(neighbor) == -1 {
                return partial + 1
            }
            if !isAtomAromatic(neighbor), bondOrder(atomID, neighbor) == .double {
                return partial + 1
            }
            return partial
        }
    }

    private func doubleBondedNitrogensCount(_ atomID: Int) -> Int {
        neighbors(of: atomID).reduce(0) { partial, neighbor in
            guard canonicalSymbol(neighbor).uppercased() == "N",
                  !isAtomAromatic(neighbor),
                  bondOrder(atomID, neighbor) == .double else {
                return partial
            }
            return partial + 1
        }
    }

    private func aromaticNitrogensCount(_ atomID: Int) -> Int {
        neighbors(of: atomID).reduce(0) { partial, neighbor in
            let isAromaticNitrogen = canonicalSymbol(neighbor).uppercased() == "N" && isInAromaticRing(neighbor)
            return partial + (isAromaticNitrogen ? 1 : 0)
        }
    }

    private func piSystemsCount(_ atomID: Int) -> Int {
        neighbors(of: atomID).reduce(0) { partial, neighbor in
            let symbol = canonicalSymbol(neighbor).uppercased()
            guard symbol != "P", symbol != "S" else {
                return partial
            }
            var piCount = 0
            var aromaticCounted = false
            for secondNeighbor in neighbors(of: neighbor) where secondNeighbor != atomID {
                guard let order = bondOrder(neighbor, secondNeighbor) else { continue }
                switch order {
                case .double, .triple:
                    piCount += 1
                case .aromatic:
                    // CDK executes on atom-typed/kekulized forms where aromatic neighbors
                    // typically contribute one effective pi relationship.
                    if !aromaticCounted {
                        piCount += 1
                        aromaticCounted = true
                    }
                case .single:
                    break
                }
            }
            return partial + piCount
        }
    }

    private func presenceOfHydroxy(_ atomID: Int) -> Bool {
        let connected = neighbors(of: atomID)
        guard let firstNeighbor = connected.first else { return false }
        guard canonicalSymbol(firstNeighbor).uppercased() == "C" else { return false }

        for candidate in neighbors(of: firstNeighbor) where canonicalSymbol(candidate).uppercased() == "O" {
            guard bondOrder(firstNeighbor, candidate) == .single else { continue }
            if connectedBondsCount(candidate) > 1, hydrogenCount(candidate) == 0 {
                return false
            }
            return true
        }
        return false
    }

    private func presenceOfNitro(_ atomID: Int) -> Bool {
        for neighbor in neighbors(of: atomID) where canonicalSymbol(neighbor).uppercased() == "N" {
            for oxygen in neighbors(of: neighbor) where canonicalSymbol(oxygen).uppercased() == "O" {
                if bondOrder(neighbor, oxygen) == .double {
                    return true
                }
            }
        }
        return false
    }

    private func presenceOfSulfat(_ atomID: Int) -> Bool {
        for neighbor in neighbors(of: atomID) where canonicalSymbol(neighbor).uppercased() == "S" {
            if oxygenCount(neighbor) >= 2, connectedBondsCount(neighbor) == 4 {
                return true
            }
        }
        return false
    }

    private func presenceOfCarbonyl(_ atomID: Int) -> Int {
        var count = 0
        for carbon in neighbors(of: atomID) where canonicalSymbol(carbon).uppercased() == "C" {
            for oxygen in neighbors(of: carbon) where canonicalSymbol(oxygen).uppercased() == "O" {
                if bondOrder(carbon, oxygen) == .double {
                    count += 1
                }
            }
        }
        return count
    }

    private func ifCarbonIsHydrophobic(_ atomID: Int) -> Bool {
        let first = neighbors(of: atomID)
        guard !first.isEmpty else { return false }

        for firstAtom in first {
            let firstSymbol = canonicalSymbol(firstAtom).uppercased()
            guard firstSymbol == "C" || firstSymbol == "H" else { return false }

            let second = neighbors(of: firstAtom)
            guard !second.isEmpty else { return false }
            for secondAtom in second {
                let secondSymbol = canonicalSymbol(secondAtom).uppercased()
                guard secondSymbol == "C" || secondSymbol == "H" else { return false }

                let third = neighbors(of: secondAtom)
                guard !third.isEmpty else { return false }
                for thirdAtom in third {
                    let thirdSymbol = canonicalSymbol(thirdAtom).uppercased()
                    guard thirdSymbol == "C" || thirdSymbol == "H" else { return false }
                }
            }
        }
        return true
    }

    private func ringContains(_ atomID: Int) -> Bool {
        ringSizesByAtomID[atomID] != nil
    }

    private func smallestRingSize(_ atomID: Int) -> Int {
        ringSizesByAtomID[atomID] ?? 0
    }

    private func isInAromaticRing(_ atomID: Int) -> Bool {
        aromaticRingMembership.contains(atomID)
    }

    private func checkRingLink(_ atomID: Int) -> Bool {
        if ringContains(atomID) {
            return true
        }
        return neighbors(of: atomID).contains(where: { ringContains($0) })
    }

    private func shortestPathDistance(from source: Int, to target: Int) -> Int? {
        if source == target { return 0 }
        var queue: [(Int, Int)] = [(source, 0)]
        var index = 0
        var seen: Set<Int> = [source]

        while index < queue.count {
            let (current, distance) = queue[index]
            index += 1
            for neighbor in neighbors(of: current) {
                if neighbor == target {
                    return distance + 1
                }
                if seen.insert(neighbor).inserted {
                    queue.append((neighbor, distance + 1))
                }
            }
        }
        return nil
    }

    private func matchesTerminalAminoAcid() -> Bool {
        for nitrogen in nodeOrder where canonicalSymbol(nitrogen).uppercased() == "N" {
            for alphaCarbon in neighbors(of: nitrogen) where canonicalSymbol(alphaCarbon).uppercased() == "C" {
                guard !isBondInRing(nitrogen, alphaCarbon) else { continue }
                for carbonylCarbon in neighbors(of: alphaCarbon) where carbonylCarbon != nitrogen {
                    guard canonicalSymbol(carbonylCarbon).uppercased() == "C" else { continue }
                    guard hasCarbonylOxygen(carbonylCarbon) else { continue }
                    if hasTerminalAcidOxygen(carbonylCarbon) {
                        return true
                    }
                }
            }
        }
        return false
    }

    private func hasCarbonylOxygen(_ atomID: Int) -> Bool {
        neighbors(of: atomID).contains { neighbor in
            canonicalSymbol(neighbor).uppercased() == "O" && bondOrder(atomID, neighbor) == .double
        }
    }

    private func hasTerminalAcidOxygen(_ atomID: Int) -> Bool {
        for oxygen in neighbors(of: atomID) where canonicalSymbol(oxygen).uppercased() == "O" && bondOrder(atomID, oxygen) == .single {
            let bonds = connectedBondsCount(oxygen)
            let hs = hydrogenCount(oxygen)
            let q = charge(oxygen)
            if (bonds == 2 && hs == 1 && q == 0) || (bonds == 1 && hs == 0 && q == -1) {
                return true
            }
        }
        return false
    }

    private func matchesPAminoSulphonicAcid() -> Bool {
        for sulfur in nodeOrder where canonicalSymbol(sulfur).uppercased() == "S" {
            guard connectedBondsCount(sulfur) == 4 else { continue }
            guard doubleBondedOxygensCount(sulfur) >= 2 else { continue }

            let neighbors = self.neighbors(of: sulfur)
            let carbonNeighbors = neighbors.filter { canonicalSymbol($0).uppercased() == "C" }
            guard carbonNeighbors.count >= 2 else { continue }

            let aromaticNeighbor = carbonNeighbors.first(where: { isAtomAromatic($0) })
            guard aromaticNeighbor != nil else { continue }

            if aromaticRingHasExocyclicAmino(aromaticSeed: aromaticNeighbor!) {
                return true
            }
        }
        return false
    }

    private func aromaticRingHasExocyclicAmino(aromaticSeed: Int) -> Bool {
        let cycles = sourceMolecule.aromaticDisplayRings().filter { $0.contains(aromaticSeed) }
        for cycle in cycles {
            for atomID in cycle {
                for neighbor in neighbors(of: atomID) where !cycle.contains(neighbor) {
                    guard canonicalSymbol(neighbor).uppercased() == "N" else { continue }
                    if hydrogenCount(neighbor) >= 1 {
                        return true
                    }
                }
            }
        }
        return false
    }

    private func matchesSalicylicAcid() -> Bool {
        let aromaticCycles = sourceMolecule.aromaticDisplayRings().filter { $0.count == 6 }
        guard !aromaticCycles.isEmpty else { return false }

        for carbonylCarbon in nodeOrder where canonicalSymbol(carbonylCarbon).uppercased() == "C" {
            guard hasCarbonylOxygen(carbonylCarbon) else { continue }
            guard hasTerminalAcidOxygen(carbonylCarbon) else { continue }

            let aromaticAttachment = neighbors(of: carbonylCarbon).first { canonicalSymbol($0).uppercased() == "C" && isAtomAromatic($0) }
            guard let aromaticAttachment else { continue }

            for cycle in aromaticCycles where cycle.contains(aromaticAttachment) {
                for ringAtom in cycle {
                    if ringAtom == aromaticAttachment { continue }
                    guard areAdjacentInCycle(a: ringAtom, b: aromaticAttachment, cycle: cycle) else { continue }
                    if hasHydroxySubstituent(on: ringAtom, excluding: cycle) {
                        return true
                    }
                }
            }
        }
        return false
    }

    private func hasOrthoOxygenPair() -> Bool {
        for carbonA in nodeOrder
        where canonicalSymbol(carbonA).uppercased() == "C" &&
            isInAromaticRing(carbonA) {
            guard hasOxygenSubstituent(on: carbonA) else { continue }
            for carbonB in neighbors(of: carbonA)
            where canonicalSymbol(carbonB).uppercased() == "C" &&
                isInAromaticRing(carbonB) {
                if hasOxygenSubstituent(on: carbonB) {
                    return true
                }
            }
        }
        return false
    }

    private func hasOxygenSubstituent(on carbonID: Int) -> Bool {
        for oxygen in neighbors(of: carbonID) where canonicalSymbol(oxygen).uppercased() == "O" {
            if !hasCarbonylOxygen(oxygen), bondOrder(carbonID, oxygen) == .single {
                return true
            }
        }
        return false
    }

    private func hasHydroxySubstituent(on carbonID: Int, excluding cycle: [Int]) -> Bool {
        for oxygen in neighbors(of: carbonID) where !cycle.contains(oxygen) && canonicalSymbol(oxygen).uppercased() == "O" {
            if bondOrder(carbonID, oxygen) == .single, hydrogenCount(oxygen) >= 1 {
                return true
            }
        }
        return false
    }

    private func areAdjacentInCycle(a: Int, b: Int, cycle: [Int]) -> Bool {
        guard let idxA = cycle.firstIndex(of: a),
              let idxB = cycle.firstIndex(of: b) else {
            return false
        }
        let n = cycle.count
        return abs(idxA - idxB) == 1 || abs(idxA - idxB) == n - 1
    }

    private func isBondInRing(_ a: Int, _ b: Int) -> Bool {
        guard let nodeA = nodes[a], let nodeB = nodes[b],
              !nodeA.isVirtualHydrogen, !nodeB.isVirtualHydrogen,
              let sourceBond = sourceMolecule.bond(between: a, and: b) else {
            return false
        }
        return CDKDescriptorSupport.isBondInRing(sourceBond, in: sourceMolecule)
    }

    private static func bondKey(_ a: Int, _ b: Int) -> Int64 {
        let lo = Int64(min(a, b))
        let hi = Int64(max(a, b))
        return (lo << 32) ^ (hi & 0xFFFF_FFFF)
    }
}
