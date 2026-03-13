import Foundation

/// CDK port of `org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor`.
///
/// The implementation mirrors the original fragment-profile approach from CDK
/// and accepts molecules with either implicit or explicit hydrogens.
public enum CDKTPSADescriptor {
    public static func calculate(for molecule: Molecule, checkAromaticity: Bool = false) -> Double {
        let aromaticBondIDs = CDKDescriptorRingSupport.aromaticBondIDs(in: molecule,
                                                                       perceiveAromaticity: checkAromaticity)

        var tpsa = 0.0
        for atom in molecule.atoms {
            let symbol = CDKDescriptorSupport.canonicalElementSymbol(atom.element)
            let uppercasedSymbol = symbol.uppercased()
            guard uppercasedSymbol == "N" || uppercasedSymbol == "O" || uppercasedSymbol == "S" || uppercasedSymbol == "P" else {
                continue
            }

            var singleBondCount = 0
            var doubleBondCount = 0
            var tripleBondCount = 0
            var aromaticBondCount = 0

            for bond in molecule.bonds(forAtom: atom.id) {
                if aromaticBondIDs.contains(bond.id) {
                    aromaticBondCount += 1
                    continue
                }

                switch bond.order {
                case .single:
                    singleBondCount += 1
                case .double:
                    doubleBondCount += 1
                case .triple:
                    tripleBondCount += 1
                case .aromatic:
                    aromaticBondCount += 1
                }
            }

            let implicitHydrogenCount = CDKDescriptorSupport.implicitHydrogenCount(on: atom.id, in: molecule)
            let explicitHydrogenCount = CDKDescriptorSupport.explicitHydrogenNeighborCount(on: atom.id, in: molecule)

            let hydrogenCount = explicitHydrogenCount + implicitHydrogenCount
            let numberOfNeighbors = molecule.neighbors(of: atom.id).count + implicitHydrogenCount
            singleBondCount += implicitHydrogenCount

            let bondOrderSum =
                Double(singleBondCount) +
                (Double(doubleBondCount) * 2.0) +
                (Double(tripleBondCount) * 3.0) +
                (Double(aromaticBondCount) * 1.5)

            var maxBondOrder = 0.0
            if singleBondCount > 0 { maxBondOrder = 1.0 }
            if aromaticBondCount > 0 { maxBondOrder = 1.5 }
            if doubleBondCount > 0 { maxBondOrder = 2.0 }
            if tripleBondCount > 0 { maxBondOrder = 3.0 }

            let isInThreeMemberedRing = CDKDescriptorRingSupport.hasThreeMemberedRing(containing: atom.id, in: molecule) ? 1 : 0
            let profile = makeProfile(symbol: symbol,
                                      maxBondOrder: maxBondOrder,
                                      bondOrderSum: bondOrderSum,
                                      numberOfNeighbors: numberOfNeighbors,
                                      hydrogenCount: hydrogenCount,
                                      formalCharge: atom.charge,
                                      aromaticBondCount: aromaticBondCount,
                                      isInThreeMemberedRing: isInThreeMemberedRing,
                                      singleBondCount: singleBondCount,
                                      doubleBondCount: doubleBondCount,
                                      tripleBondCount: tripleBondCount)
            tpsa += contributionsByProfile[profile] ?? 0.0
        }

        return tpsa
    }

    private static func makeProfile(symbol: String,
                                    maxBondOrder: Double,
                                    bondOrderSum: Double,
                                    numberOfNeighbors: Int,
                                    hydrogenCount: Int,
                                    formalCharge: Int,
                                    aromaticBondCount: Int,
                                    isInThreeMemberedRing: Int,
                                    singleBondCount: Int,
                                    doubleBondCount: Int,
                                    tripleBondCount: Int) -> String {
        "\(symbol)+\(formatted(maxBondOrder))+\(formatted(bondOrderSum))+\(numberOfNeighbors)+\(hydrogenCount)+\(formalCharge)+\(aromaticBondCount)+\(isInThreeMemberedRing)+\(singleBondCount)+\(doubleBondCount)+\(tripleBondCount)"
    }

    private static func formatted(_ value: Double) -> String {
        String(format: "%.1f", locale: Locale(identifier: "en_US_POSIX"), value)
    }

    // Fragment profile contributions copied from CDK's TPSADescriptor.
    private static let contributionsByProfile: [String: Double] = [
        "N+1.0+3.0+3+0+0+0+0+3+0+0": 3.24,
        "N+2.0+3.0+2+0+0+0+0+1+1+0": 12.36,
        "N+3.0+3.0+1+0+0+0+0+0+0+1": 23.79,
        "N+2.0+5.0+3+0+0+0+0+1+2+0": 11.68,
        "N+3.0+5.0+2+0+0+0+0+0+1+1": 13.6,
        "N+1.0+3.0+3+0+0+0+1+3+0+0": 3.01,
        "N+1.0+3.0+3+1+0+0+0+3+0+0": 12.03,
        "N+1.0+3.0+3+1+0+0+1+3+0+0": 21.94,
        "N+2.0+3.0+2+1+0+0+0+1+1+0": 23.85,
        "N+1.0+3.0+3+2+0+0+0+3+0+0": 26.02,
        "N+1.0+4.0+4+0+1+0+0+4+0+0": 0.0,
        "N+2.0+4.0+3+0+1+0+0+2+1+0": 3.01,
        "N+3.0+4.0+2+0+1+0+0+1+0+1": 4.36,
        "N+1.0+4.0+4+1+1+0+0+4+0+0": 4.44,
        "N+2.0+4.0+3+1+1+0+0+2+1+0": 13.97,
        "N+1.0+4.0+4+2+1+0+0+4+0+0": 16.61,
        "N+2.0+4.0+3+2+1+0+0+2+1+0": 25.59,
        "N+1.0+4.0+4+3+1+0+0+4+0+0": 27.64,
        "N+1.5+3.0+2+0+0+2+0+0+0+0": 12.89,
        "N+1.5+4.5+3+0+0+3+0+0+0+0": 4.41,
        "N+1.5+4.0+3+0+0+2+0+1+0+0": 4.93,
        "N+2.0+5.0+3+0+0+2+0+0+1+0": 8.39,
        "N+1.5+4.0+3+1+0+2+0+1+0+0": 15.79,
        "N+1.5+4.5+3+0+1+3+0+0+0+0": 4.1,
        "N+1.5+4.0+3+0+1+2+0+1+0+0": 3.88,
        "N+1.5+4.0+3+1+1+2+0+1+0+0": 14.14,
        "O+1.0+2.0+2+0+0+0+0+2+0+0": 9.23,
        "O+1.0+2.0+2+0+0+0+1+2+0+0": 12.53,
        "O+2.0+2.0+1+0+0+0+0+0+1+0": 17.07,
        "O+1.0+1.0+1+0+-1+0+0+1+0+0": 23.06,
        "O+1.0+2.0+2+1+0+0+0+2+0+0": 20.23,
        "O+1.5+3.0+2+0+0+2+0+0+0+0": 13.14,
        "S+1.0+2.0+2+0+0+0+0+2+0+0": 25.3,
        "S+2.0+2.0+1+0+0+0+0+0+1+0": 32.09,
        "S+2.0+4.0+3+0+0+0+0+2+1+0": 19.21,
        "S+2.0+6.0+4+0+0+0+0+2+2+0": 8.38,
        "S+1.0+2.0+2+1+0+0+0+2+0+0": 38.8,
        "S+1.5+3.0+2+0+0+2+0+0+0+0": 28.24,
        "S+2.0+5.0+3+0+0+2+0+0+1+0": 21.7,
        "P+1.0+3.0+3+0+0+0+0+3+0+0": 13.59,
        "P+2.0+3.0+3+0+0+0+0+1+1+0": 34.14,
        "P+2.0+5.0+4+0+0+0+0+3+1+0": 9.81,
        "P+2.0+5.0+4+1+0+0+0+3+1+0": 23.47
    ]
}
