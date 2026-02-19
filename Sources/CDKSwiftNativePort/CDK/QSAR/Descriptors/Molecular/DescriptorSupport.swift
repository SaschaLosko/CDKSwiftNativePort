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

    static func averageAtomicMass(for atom: Atom) -> Double {
        if let isotope = atom.isotopeMassNumber {
            return Double(isotope)
        }
        return averageAtomicMass(forElementSymbol: atom.element)
    }

    static func averageAtomicMass(forElementSymbol symbol: String) -> Double {
        let key = canonicalElementSymbol(symbol).uppercased()
        return averageAtomicMassByElement[key] ?? 0.0
    }

    static func monoisotopicAtomicMass(for atom: Atom) -> Double {
        if let isotope = atom.isotopeMassNumber {
            return Double(isotope)
        }
        return monoisotopicAtomicMass(forElementSymbol: atom.element)
    }

    static func monoisotopicAtomicMass(forElementSymbol symbol: String) -> Double {
        let key = canonicalElementSymbol(symbol).uppercased()
        return monoisotopicMassByElement[key] ?? averageAtomicMassByElement[key] ?? 0.0
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
        Set(molecule.atoms.compactMap { atom in
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
                for neighbor in molecule.neighbors(of: current) where heavyIDs.contains(neighbor) && !seen.contains(neighbor) {
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

    private static func hasPath(from source: Int,
                                to target: Int,
                                in molecule: Molecule,
                                excludingBondID: Int) -> Bool {
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

        for carbonBond in molecule.bonds(forAtom: carbonID) where carbonBond.id != bond.id && carbonBond.order == .double {
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
            return 3.0
        case "F", "CL", "BR", "I":
            return 1.0
        default:
            return 0.0
        }
    }

    private static let averageAtomicMassByElement: [String: Double] = [
        "H": 1.00794,
        "HE": 4.002602,
        "LI": 6.941,
        "BE": 9.012182,
        "B": 10.811,
        "C": 12.0107,
        "N": 14.0067,
        "O": 15.9994,
        "F": 18.9984032,
        "NE": 20.1797,
        "NA": 22.98976928,
        "MG": 24.305,
        "AL": 26.9815386,
        "SI": 28.0855,
        "P": 30.973762,
        "S": 32.065,
        "CL": 35.453,
        "AR": 39.948,
        "K": 39.0983,
        "CA": 40.078,
        "SC": 44.955912,
        "TI": 47.867,
        "V": 50.9415,
        "CR": 51.9961,
        "MN": 54.938045,
        "FE": 55.845,
        "CO": 58.933195,
        "NI": 58.6934,
        "CU": 63.546,
        "ZN": 65.38,
        "GA": 69.723,
        "GE": 72.64,
        "AS": 74.9216,
        "SE": 78.96,
        "BR": 79.904,
        "KR": 83.798,
        "RB": 85.4678,
        "SR": 87.62,
        "Y": 88.90585,
        "ZR": 91.224,
        "NB": 92.90638,
        "MO": 95.96,
        "RU": 101.07,
        "RH": 102.9055,
        "PD": 106.42,
        "AG": 107.8682,
        "CD": 112.411,
        "IN": 114.818,
        "SN": 118.71,
        "SB": 121.76,
        "TE": 127.6,
        "I": 126.90447,
        "XE": 131.293,
        "CS": 132.9054519,
        "BA": 137.327,
        "PT": 195.084,
        "AU": 196.966569,
        "HG": 200.59,
        "PB": 207.2,
        "BI": 208.9804
    ]

    private static let monoisotopicMassByElement: [String: Double] = [
        "H": 1.00782503223,
        "B": 11.00930536,
        "C": 12.0,
        "N": 14.00307400443,
        "O": 15.99491461957,
        "F": 18.99840316273,
        "NA": 22.989769282,
        "MG": 23.985041697,
        "AL": 26.98153853,
        "SI": 27.97692653465,
        "P": 30.97376199842,
        "S": 31.9720711744,
        "CL": 34.968852682,
        "K": 38.9637064864,
        "CA": 39.962590863,
        "FE": 55.93493633,
        "CU": 62.92959772,
        "ZN": 63.92914201,
        "SE": 79.9165218,
        "BR": 78.9183376,
        "I": 126.9044719,
        "PT": 194.9647917,
        "AU": 196.96656879,
        "HG": 201.9706434,
        "PB": 207.9766525
    ]
}
