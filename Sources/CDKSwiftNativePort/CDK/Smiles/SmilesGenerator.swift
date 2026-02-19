import Foundation

/// Swift counterpart of CDK's `org.openscience.cdk.smiles.SmilesGenerator`.
public final class CDKSmilesGenerator {
    public let flavor: CDKSmiFlavor

    public init(flavor: CDKSmiFlavor = .cdkDefault) {
        self.flavor = flavor
    }

    public func create(_ molecule: Molecule) -> String {
        guard !molecule.atoms.isEmpty else { return "" }

        let atomByID = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0) })
        let adjacency = buildAdjacency(molecule)
        let bondByEdge = buildBondLookup(molecule)
        let aromaticBondIDs = flavor.contains(.useAromaticSymbols) ? molecule.aromaticDisplayBondIDs() : Set<Int>()

        let components = connectedComponents(molecule, adjacency: adjacency)

        var componentStrings: [String] = []
        componentStrings.reserveCapacity(components.count)

        for component in components {
            let componentSet = Set(component)
            guard let root = component.min() else { continue }

            let (treeAdjacency, treeEdges) = buildTree(component: componentSet,
                                                       root: root,
                                                       adjacency: adjacency)

            let componentEdges = Set(
                bondByEdge.keys.filter { componentSet.contains($0.a) && componentSet.contains($0.b) }
            )
            let ringEdges = componentEdges.subtracting(treeEdges)
            let ringNumbers = assignRingNumbers(to: ringEdges)

            var ringClosuresByAtom: [Int: [RingClosure]] = [:]
            for edge in ringEdges {
                guard let bond = bondByEdge[edge], let ringNumber = ringNumbers[edge] else { continue }
                let symbol = bondToken(bond, atomByID: atomByID, aromaticBondIDs: aromaticBondIDs)
                ringClosuresByAtom[edge.a, default: []].append(
                    RingClosure(number: ringNumber, bondSymbol: symbol, partner: edge.b)
                )
                ringClosuresByAtom[edge.b, default: []].append(
                    RingClosure(number: ringNumber, bondSymbol: symbol, partner: edge.a)
                )
            }

            for atomID in ringClosuresByAtom.keys {
                ringClosuresByAtom[atomID]?.sort { lhs, rhs in
                    if lhs.number != rhs.number { return lhs.number < rhs.number }
                    return lhs.partner < rhs.partner
                }
            }

            var visitedTreeAtoms: Set<Int> = []
            let rendered = renderTreeAtom(root,
                                          parent: nil,
                                          treeAdjacency: treeAdjacency,
                                          ringClosuresByAtom: ringClosuresByAtom,
                                          bondByEdge: bondByEdge,
                                          atomByID: atomByID,
                                          aromaticBondIDs: aromaticBondIDs,
                                          visited: &visitedTreeAtoms)
            componentStrings.append(rendered)
        }

        return componentStrings.joined(separator: ".")
    }

    private struct EdgeKey: Hashable, Comparable {
        let a: Int
        let b: Int

        init(_ u: Int, _ v: Int) {
            a = min(u, v)
            b = max(u, v)
        }

        static func < (lhs: EdgeKey, rhs: EdgeKey) -> Bool {
            if lhs.a != rhs.a { return lhs.a < rhs.a }
            return lhs.b < rhs.b
        }
    }

    private struct RingClosure {
        let number: Int
        let bondSymbol: String
        let partner: Int
    }

    private func buildAdjacency(_ molecule: Molecule) -> [Int: [Int]] {
        var map: [Int: Set<Int>] = [:]
        for atom in molecule.atoms {
            map[atom.id] = []
        }
        for bond in molecule.bonds {
            map[bond.a1, default: []].insert(bond.a2)
            map[bond.a2, default: []].insert(bond.a1)
        }
        return map.mapValues { Array($0).sorted() }
    }

    private func buildBondLookup(_ molecule: Molecule) -> [EdgeKey: Bond] {
        var lookup: [EdgeKey: Bond] = [:]
        for bond in molecule.bonds {
            lookup[EdgeKey(bond.a1, bond.a2)] = bond
        }
        return lookup
    }

    private func connectedComponents(_ molecule: Molecule,
                                     adjacency: [Int: [Int]]) -> [[Int]] {
        var visited: Set<Int> = []
        var components: [[Int]] = []

        for atomID in molecule.atoms.map(\.id).sorted() {
            if visited.contains(atomID) { continue }

            var stack = [atomID]
            var component: [Int] = []
            visited.insert(atomID)

            while let node = stack.popLast() {
                component.append(node)
                for neighbor in adjacency[node, default: []] where !visited.contains(neighbor) {
                    visited.insert(neighbor)
                    stack.append(neighbor)
                }
            }

            component.sort()
            components.append(component)
        }

        components.sort { lhs, rhs in
            (lhs.first ?? Int.max) < (rhs.first ?? Int.max)
        }

        return components
    }

    private func buildTree(component: Set<Int>,
                           root: Int,
                           adjacency: [Int: [Int]]) -> (adjacency: [Int: [Int]], edges: Set<EdgeKey>) {
        var treeAdjacency: [Int: [Int]] = [:]
        var treeEdges: Set<EdgeKey> = []
        var visited: Set<Int> = []

        func dfs(_ atomID: Int) {
            visited.insert(atomID)
            let neighbors = adjacency[atomID, default: []]
                .filter { component.contains($0) }
                .sorted()

            for neighbor in neighbors where !visited.contains(neighbor) {
                treeAdjacency[atomID, default: []].append(neighbor)
                treeAdjacency[neighbor, default: []].append(atomID)
                treeEdges.insert(EdgeKey(atomID, neighbor))
                dfs(neighbor)
            }
        }

        dfs(root)

        for atomID in component {
            treeAdjacency[atomID, default: []].sort()
        }

        return (treeAdjacency, treeEdges)
    }

    private func assignRingNumbers(to ringEdges: Set<EdgeKey>) -> [EdgeKey: Int] {
        var assignments: [EdgeKey: Int] = [:]
        var next = 1

        for edge in ringEdges.sorted() {
            assignments[edge] = next
            next += 1
            if next > 99 { next = 1 }
        }

        return assignments
    }

    private func renderTreeAtom(_ atomID: Int,
                                parent: Int?,
                                treeAdjacency: [Int: [Int]],
                                ringClosuresByAtom: [Int: [RingClosure]],
                                bondByEdge: [EdgeKey: Bond],
                                atomByID: [Int: Atom],
                                aromaticBondIDs: Set<Int>,
                                visited: inout Set<Int>) -> String {
        visited.insert(atomID)

        guard let atom = atomByID[atomID] else { return "" }

        var out = atomToken(atom)

        for closure in ringClosuresByAtom[atomID] ?? [] {
            out += closure.bondSymbol + ringToken(closure.number)
        }

        let children = treeAdjacency[atomID, default: []]
            .filter { $0 != parent }
            .sorted()

        for (index, childID) in children.enumerated() {
            let edge = EdgeKey(atomID, childID)
            guard let bond = bondByEdge[edge] else { continue }

            let branchBond = bondToken(bond, atomByID: atomByID, aromaticBondIDs: aromaticBondIDs)
            let branchText = branchBond + renderTreeAtom(childID,
                                                         parent: atomID,
                                                         treeAdjacency: treeAdjacency,
                                                         ringClosuresByAtom: ringClosuresByAtom,
                                                         bondByEdge: bondByEdge,
                                                         atomByID: atomByID,
                                                         aromaticBondIDs: aromaticBondIDs,
                                                         visited: &visited)

            if index == 0 {
                out += branchText
            } else {
                out += "(\(branchText))"
            }
        }

        return out
    }

    private func ringToken(_ number: Int) -> String {
        if number < 10 {
            return "\(number)"
        }
        return "%\(number)"
    }

    private func atomToken(_ atom: Atom) -> String {
        if atom.queryType != nil || atom.atomList != nil {
            return "*"
        }

        let element = atom.element.trimmingCharacters(in: .whitespacesAndNewlines)
        if element == "*" {
            return "*"
        }

        let plain = plainSymbol(for: atom)
        if !requiresBracket(atom), let plain {
            return plain
        }

        var token = "["

        if flavor.contains(.isomeric), let isotope = atom.isotopeMassNumber {
            token += "\(isotope)"
        }

        token += bracketElementSymbol(for: atom)

        if flavor.contains(.isomeric) {
            switch atom.chirality {
            case .clockwise:
                token += "@"
            case .anticlockwise:
                token += "@@"
            case .none:
                break
            }
        }

        if let h = atom.explicitHydrogenCount {
            token += hydrogenToken(h)
        }

        if atom.charge != 0 {
            token += chargeToken(atom.charge)
        }

        token += "]"
        return token
    }

    private func requiresBracket(_ atom: Atom) -> Bool {
        if atom.queryType != nil || atom.atomList != nil { return true }
        if atom.charge != 0 { return true }
        if atom.explicitHydrogenCount != nil { return true }
        if flavor.contains(.isomeric) {
            if atom.isotopeMassNumber != nil { return true }
            if atom.chirality != .none { return true }
        }
        return plainSymbol(for: atom) == nil
    }

    private func plainSymbol(for atom: Atom) -> String? {
        let upper = atom.element.uppercased()
        if upper == "*" {
            return "*"
        }

        let aromaticEnabled = flavor.contains(.useAromaticSymbols)
        if aromaticEnabled && atom.aromatic {
            switch upper {
            case "B", "C", "N", "O", "P", "S":
                return upper.lowercased()
            case "SE":
                return "se"
            case "AS":
                return "as"
            default:
                return nil
            }
        }

        switch upper {
        case "B", "C", "N", "O", "P", "S", "F", "CL", "BR", "I":
            return normalizedElementSymbol(upper)
        default:
            return nil
        }
    }

    private func bracketElementSymbol(for atom: Atom) -> String {
        let upper = atom.element.uppercased()
        if upper == "*" { return "*" }

        if flavor.contains(.useAromaticSymbols), atom.aromatic {
            switch upper {
            case "B", "C", "N", "O", "P", "S":
                return upper.lowercased()
            case "SE":
                return "se"
            case "AS":
                return "as"
            default:
                break
            }
        }

        return normalizedElementSymbol(upper)
    }

    private func normalizedElementSymbol(_ upper: String) -> String {
        guard let first = upper.first else { return upper }
        let head = String(first)
        let tail = String(upper.dropFirst()).lowercased()
        return head + tail
    }

    private func hydrogenToken(_ count: Int) -> String {
        if count <= 0 {
            return "H0"
        }
        if count == 1 {
            return "H"
        }
        return "H\(count)"
    }

    private func chargeToken(_ charge: Int) -> String {
        if charge == 1 { return "+" }
        if charge == -1 { return "-" }

        let sign = charge > 0 ? "+" : "-"
        let magnitude = abs(charge)
        if magnitude <= 3 {
            return String(repeating: sign, count: magnitude)
        }
        return "\(sign)\(magnitude)"
    }

    private func bondToken(_ bond: Bond,
                           atomByID: [Int: Atom],
                           aromaticBondIDs: Set<Int>) -> String {
        if bond.queryType != nil {
            return "-"
        }

        switch bond.order {
        case .single:
            if flavor.contains(.isomeric) {
                switch bond.stereo {
                case .up, .upReversed:
                    return "/"
                case .down, .downReversed:
                    return "\\"
                case .either, .none:
                    break
                }
            }

            if flavor.contains(.useAromaticSymbols), aromaticBondIDs.contains(bond.id) {
                let leftAromatic = atomByID[bond.a1]?.aromatic == true
                let rightAromatic = atomByID[bond.a2]?.aromatic == true
                if leftAromatic && rightAromatic {
                    return ""
                }
            }

            return ""

        case .double:
            return "="

        case .triple:
            return "#"

        case .aromatic:
            if flavor.contains(.useAromaticSymbols) {
                let leftAromatic = atomByID[bond.a1]?.aromatic == true
                let rightAromatic = atomByID[bond.a2]?.aromatic == true
                if leftAromatic && rightAromatic {
                    return ""
                }
            }
            return ":"
        }
    }
}

/// Factory facade mirroring CDK-style generator construction patterns.
public final class CDKSmilesGeneratorFactory {
    public static let shared = CDKSmilesGeneratorFactory()

    private init() {}

    public func newSmilesGenerator(flavor: CDKSmiFlavor = .cdkDefault) -> CDKSmilesGenerator {
        CDKSmilesGenerator(flavor: flavor)
    }
}
