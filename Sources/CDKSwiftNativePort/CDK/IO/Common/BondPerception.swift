import Foundation
import CoreGraphics

enum CDKBondPerception {
    private struct Candidate {
        let a1: Int
        let a2: Int
        let distance: Double
        let threshold: Double
    }

    static func inferSingleBonds(for atoms: [Atom], existingBonds: [Bond] = []) -> [Bond] {
        guard atoms.count >= 2 else { return existingBonds }

        let atomByID = Dictionary(uniqueKeysWithValues: atoms.map { ($0.id, $0) })
        var bonds = existingBonds
        var occupiedPairs = Set(existingBonds.map(canonicalPair(of:)))

        var valenceLoad: [Int: Int] = [:]
        for atom in atoms {
            valenceLoad[atom.id] = 0
        }
        for bond in existingBonds {
            let load = bondOrderLoad(bond.order)
            valenceLoad[bond.a1, default: 0] += load
            valenceLoad[bond.a2, default: 0] += load
        }

        var candidates: [Candidate] = []
        for i in 0..<(atoms.count - 1) {
            for j in (i + 1)..<atoms.count {
                let a = atoms[i]
                let b = atoms[j]
                let pair = canonicalPair(a.id, b.id)
                if occupiedPairs.contains(pair) { continue }

                let symbolA = CDKDescriptorSupport.canonicalElementSymbol(a.element).uppercased()
                let symbolB = CDKDescriptorSupport.canonicalElementSymbol(b.element).uppercased()
                if symbolA == "H" && symbolB == "H" { continue }

                let distance = hypot(a.position.x - b.position.x, a.position.y - b.position.y)
                if distance < 0.20 || distance > 2.60 { continue }

                let threshold = bondThreshold(symbolA: symbolA, symbolB: symbolB)
                if distance <= threshold {
                    candidates.append(Candidate(a1: a.id,
                                                a2: b.id,
                                                distance: Double(distance),
                                                threshold: threshold))
                }
            }
        }

        candidates.sort {
            if abs($0.distance - $1.distance) > 0.0001 {
                return $0.distance < $1.distance
            }
            let slackA = $0.threshold - $0.distance
            let slackB = $1.threshold - $1.distance
            if abs(slackA - slackB) > 0.0001 {
                return slackA > slackB
            }
            if min($0.a1, $0.a2) != min($1.a1, $1.a2) {
                return min($0.a1, $0.a2) < min($1.a1, $1.a2)
            }
            return max($0.a1, $0.a2) < max($1.a1, $1.a2)
        }

        var nextBondID = (bonds.map(\.id).max() ?? 0) + 1
        for candidate in candidates {
            let pair = canonicalPair(candidate.a1, candidate.a2)
            if occupiedPairs.contains(pair) { continue }

            guard let atomA = atomByID[candidate.a1], let atomB = atomByID[candidate.a2] else {
                continue
            }

            let maxValenceA = maxValence(for: atomA)
            let maxValenceB = maxValence(for: atomB)
            let currentA = valenceLoad[candidate.a1, default: 0]
            let currentB = valenceLoad[candidate.a2, default: 0]

            if currentA + 1 > maxValenceA || currentB + 1 > maxValenceB {
                continue
            }

            bonds.append(Bond(id: nextBondID,
                              a1: candidate.a1,
                              a2: candidate.a2,
                              order: .single,
                              stereo: .none))
            nextBondID += 1
            occupiedPairs.insert(pair)
            valenceLoad[candidate.a1, default: 0] += 1
            valenceLoad[candidate.a2, default: 0] += 1
        }

        return bonds
    }

    private static func bondOrderLoad(_ order: BondOrder) -> Int {
        switch order {
        case .single, .aromatic:
            return 1
        case .double:
            return 2
        case .triple:
            return 3
        }
    }

    private static func bondThreshold(symbolA: String, symbolB: String) -> Double {
        let rA = covalentRadius(for: symbolA)
        let rB = covalentRadius(for: symbolB)
        let scaled = (rA + rB) * 1.22
        return min(2.30, max(0.65, scaled + 0.06))
    }

    private static func maxValence(for atom: Atom) -> Int {
        let symbol = CDKDescriptorSupport.canonicalElementSymbol(atom.element).uppercased()
        switch symbol {
        case "H":
            return 1
        case "B":
            return 3
        case "C":
            return 4
        case "N":
            return atom.charge > 0 ? 4 : 3
        case "O":
            return atom.charge > 0 ? 3 : 2
        case "P":
            return 5
        case "S":
            return 6
        case "F", "CL", "BR", "I":
            return 1
        case "SI":
            return 4
        default:
            return 8
        }
    }

    private static func covalentRadius(for symbol: String) -> Double {
        covalentRadiusByElement[symbol] ?? 0.77
    }

    private static func canonicalPair(of bond: Bond) -> String {
        canonicalPair(bond.a1, bond.a2)
    }

    private static func canonicalPair(_ a: Int, _ b: Int) -> String {
        if a < b { return "\(a)-\(b)" }
        return "\(b)-\(a)"
    }

    private static let covalentRadiusByElement: [String: Double] = [
        "H": 0.31,
        "B": 0.85,
        "C": 0.76,
        "N": 0.71,
        "O": 0.66,
        "F": 0.57,
        "SI": 1.11,
        "P": 1.07,
        "S": 1.05,
        "CL": 1.02,
        "BR": 1.20,
        "I": 1.39,
        "NA": 1.66,
        "MG": 1.41,
        "K": 2.03,
        "CA": 1.76,
        "FE": 1.24,
        "ZN": 1.22,
        "CU": 1.32
    ]
}
