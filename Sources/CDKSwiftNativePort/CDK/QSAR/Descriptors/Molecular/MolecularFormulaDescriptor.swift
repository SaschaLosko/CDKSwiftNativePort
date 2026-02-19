import Foundation

public enum CDKMolecularFormulaDescriptor {
    public static func calculate(for molecule: Molecule,
                                 includeImplicitHydrogens: Bool = true) -> String {
        var composition: [String: Int] = [:]

        for atom in molecule.atoms {
            let symbol = CDKDescriptorSupport.canonicalElementSymbol(atom.element)
            guard !symbol.isEmpty else { continue }
            composition[symbol, default: 0] += 1

            if includeImplicitHydrogens && symbol.uppercased() != "H" {
                let implicitHydrogenCount = CDKDescriptorSupport.implicitHydrogenCount(on: atom.id, in: molecule)
                if implicitHydrogenCount > 0 {
                    composition["H", default: 0] += implicitHydrogenCount
                }
            }
        }

        guard !composition.isEmpty else { return "Unknown" }

        var tokens: [String] = []
        if let carbon = composition.removeValue(forKey: "C"), carbon > 0 {
            tokens.append("C" + (carbon == 1 ? "" : "\(carbon)"))
        }
        if let hydrogen = composition.removeValue(forKey: "H"), hydrogen > 0 {
            tokens.append("H" + (hydrogen == 1 ? "" : "\(hydrogen)"))
        }
        for element in composition.keys.sorted() {
            guard let count = composition[element], count > 0 else { continue }
            tokens.append(element + (count == 1 ? "" : "\(count)"))
        }

        return tokens.joined()
    }
}
