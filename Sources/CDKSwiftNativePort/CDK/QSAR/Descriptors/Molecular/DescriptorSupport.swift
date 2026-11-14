import Foundation

enum CDKDescriptorSupport {
    static func canonicalElementSymbol(_ raw: String) -> String {
        let trimmed = raw.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { return "" }

        guard let first = trimmed.first else { return "" }
        if !first.isLetter {
            return trimmed
        }

        let head = String(first).uppercased()
        let tail = String(trimmed.dropFirst()).lowercased()
        return head + tail
    }

    static func explicitHydrogenNeighborCount(on atomID: Int, in molecule: Molecule) -> Int {
        molecule.bonds(forAtom: atomID).reduce(0) { partial, bond in
            let neighborID = bond.a1 == atomID ? bond.a2 : bond.a1
            guard let neighbor = molecule.atom(id: neighborID) else { return partial }
            return partial + (canonicalElementSymbol(neighbor.element).uppercased() == "H" ? 1 : 0)
        }
    }

    static func implicitHydrogenCount(on atomID: Int, in molecule: Molecule) -> Int {
        guard let atom = molecule.atom(id: atomID) else { return 0 }
        let element = canonicalElementSymbol(atom.element).uppercased()
        guard element != "H" else { return 0 }
        if let explicit = atom.explicitHydrogenCount {
            return max(0, explicit)
        }

        let targetValence = preferredValence(for: atom)
        guard targetValence > 0 else { return 0 }

        let bondOrderSum = molecule.bonds(forAtom: atomID)
            .reduce(0.0) { partial, bond in
                partial + bondOrderContribution(for: bond.order)
            }
        return max(0, Int(round(targetValence - bondOrderSum)))
    }

    static func totalHydrogenCount(on atomID: Int, in molecule: Molecule) -> Int {
        explicitHydrogenNeighborCount(on: atomID, in: molecule) + implicitHydrogenCount(on: atomID, in: molecule)
    }

    static func heavyDegree(of atomID: Int, in molecule: Molecule) -> Int {
        molecule.bonds(forAtom: atomID).reduce(0) { partial, bond in
            let neighborID = bond.a1 == atomID ? bond.a2 : bond.a1
            guard let neighbor = molecule.atom(id: neighborID) else { return partial }
            return partial + (canonicalElementSymbol(neighbor.element).uppercased() == "H" ? 0 : 1)
        }
    }

    static func heavyAtomIDs(in molecule: Molecule) -> Set<Int> {
        Set(
            molecule.atoms.compactMap { atom in
                canonicalElementSymbol(atom.element).uppercased() == "H" ? nil : atom.id
            })
    }

    static func heavyAtomConnectedComponentCount(in molecule: Molecule) -> Int {
        let heavyIDs = heavyAtomIDs(in: molecule)
        guard !heavyIDs.isEmpty else { return 0 }

        var seen: Set<Int> = []
        var components = 0
        for seed in heavyIDs where !seen.contains(seed) {
            components += 1
            var stack = [seed]
            seen.insert(seed)
            while let current = stack.popLast() {
                for neighbor in molecule.neighbors(of: current)
                where heavyIDs.contains(neighbor) && !seen.contains(neighbor) {
                    seen.insert(neighbor)
                    stack.append(neighbor)
                }
            }
        }
        return components
    }

    static func isBondInRing(_ bond: Bond, in molecule: Molecule) -> Bool {
        guard bond.a1 != bond.a2 else { return true }
        return hasPath(from: bond.a1, to: bond.a2, in: molecule, excludingBondID: bond.id)
    }

    private static func hasPath(
        from source: Int,
        to target: Int,
        in molecule: Molecule,
        excludingBondID: Int
    ) -> Bool {
        var seen: Set<Int> = [source]
        var stack: [Int] = [source]

        while let current = stack.popLast() {
            for bond in molecule.bonds(forAtom: current) where bond.id != excludingBondID {
                let next = bond.a1 == current ? bond.a2 : bond.a1
                if next == target { return true }
                if seen.insert(next).inserted {
                    stack.append(next)
                }
            }
        }
        return false
    }

    static func isAmideLikeCNBond(_ bond: Bond, in molecule: Molecule) -> Bool {
        guard bond.order == .single else { return false }
        guard let atomA = molecule.atom(id: bond.a1), let atomB = molecule.atom(id: bond.a2) else { return false }

        let symbolA = canonicalElementSymbol(atomA.element).uppercased()
        let symbolB = canonicalElementSymbol(atomB.element).uppercased()

        let carbonID: Int
        let nitrogenID: Int
        if symbolA == "C" && symbolB == "N" {
            carbonID = atomA.id
            nitrogenID = atomB.id
        } else if symbolA == "N" && symbolB == "C" {
            carbonID = atomB.id
            nitrogenID = atomA.id
        } else {
            return false
        }

        guard molecule.bond(between: carbonID, and: nitrogenID)?.id == bond.id else { return false }

        for carbonBond in molecule.bonds(forAtom: carbonID)
        where carbonBond.id != bond.id && carbonBond.order == .double {
            let neighborID = carbonBond.a1 == carbonID ? carbonBond.a2 : carbonBond.a1
            guard let neighbor = molecule.atom(id: neighborID) else { continue }
            let neighborSymbol = canonicalElementSymbol(neighbor.element).uppercased()
            if neighborSymbol == "O" || neighborSymbol == "S" || neighborSymbol == "N" {
                return true
            }
        }
        return false
    }

    private static func bondOrderContribution(for order: BondOrder) -> Double {
        switch order {
        case .single: return 1.0
        case .double: return 2.0
        case .triple: return 3.0
        case .aromatic: return 1.0
        }
    }

    private static func preferredValence(for atom: Atom) -> Double {
        let element = canonicalElementSymbol(atom.element).uppercased()
        if atom.aromatic {
            switch element {
            case "C": return 3.0
            case "N", "O", "S", "P": return 2.0
            default: break
            }
        }

        switch element {
        case "C":
            return 4.0
        case "N":
            return atom.charge > 0 ? 4.0 : 3.0
        case "O":
            if atom.charge > 0 { return 3.0 }
            if atom.charge < 0 { return 1.0 }
            return 2.0
        case "S":
            return atom.charge > 0 ? 3.0 : 2.0
        case "P":
            return atom.charge > 0 ? 4.0 : 3.0
        case "B":
            return atom.charge < 0 ? 4.0 : 3.0
        case "F", "CL", "BR", "I":
            return 1.0
        default:
            return 0.0
        }
    }

}
