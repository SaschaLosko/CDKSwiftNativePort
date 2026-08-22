import Foundation

/// Extracts InChI canonical atom labels for deterministic Universal SMILES
/// traversal. Fixed-H (`/F:`) labels take precedence over standard (`/N:`)
/// labels, matching CDK 2.13.
public enum CDKInChIAuxiliaryNumbering {
    public static func universalSmilesNumbers(from auxInfo: String,
                                              molecule: Molecule) throws -> [Int] {
        guard !molecule.atoms.isEmpty else { return [] }
        let atomCount = molecule.atomCount
        let numberingStart = numberingLayerStart(in: auxInfo)
        guard let numberingStart,
              let baseLayer = layerValue(in: auxInfo, marker: "/N:", from: numberingStart) else {
            throw ChemError.parseFailed("InChI AuxInfo does not contain canonical /N: numbering.")
        }

        let baseComponents = baseLayer.split(separator: ";", omittingEmptySubsequences: false).map(String.init)
        var numbers = Array(repeating: 0, count: atomCount)
        var firstAtomIndices = Array(repeating: -1, count: baseComponents.count)
        var nextLabel = 1

        func assign(_ component: String, baseComponentIndex: Int?) throws {
            let atomNumbers = try parseAtomNumbers(component, atomCount: atomCount)
            if let baseComponentIndex,
               baseComponentIndex < firstAtomIndices.count,
               let first = atomNumbers.first {
                firstAtomIndices[baseComponentIndex] = first - 1
            }
            for atomNumber in atomNumbers {
                numbers[atomNumber - 1] = nextLabel
                nextLabel += 1
            }
        }

        if let fixedLayer = layerValue(in: auxInfo, marker: "/F:", from: numberingStart) {
            let fixedComponents = fixedLayer.split(separator: ";", omittingEmptySubsequences: false).map(String.init)
            var baseComponentIndex = 0
            for component in fixedComponents where !component.isEmpty {
                if component.hasSuffix("m") {
                    let multiplierText = component.dropLast()
                    let multiplier = multiplierText.isEmpty ? 1 : Int(multiplierText) ?? 0
                    guard multiplier > 0 else {
                        throw ChemError.parseFailed("Invalid InChI fixed-H component '\(component)'.")
                    }
                    for offset in 0..<multiplier where baseComponentIndex + offset < baseComponents.count {
                        try assign(baseComponents[baseComponentIndex + offset],
                                   baseComponentIndex: baseComponentIndex + offset)
                    }
                    baseComponentIndex += multiplier
                } else {
                    try assign(component, baseComponentIndex: nil)
                    baseComponentIndex += 1
                }
            }
        } else {
            for (index, component) in baseComponents.enumerated() where !component.isEmpty {
                try assign(component, baseComponentIndex: index)
            }
        }

        applyAnionicOxygenStartCorrection(numbers: &numbers,
                                          firstAtomIndices: firstAtomIndices,
                                          molecule: molecule)
        for index in numbers.indices where numbers[index] == 0 {
            numbers[index] = nextLabel
            nextLabel += 1
        }
        return numbers
    }

    private static func numberingLayerStart(in auxInfo: String) -> String.Index? {
        guard let reconnectRange = auxInfo.range(of: "/R:") else {
            return auxInfo.startIndex
        }
        return auxInfo.range(of: "/N:", range: reconnectRange.upperBound..<auxInfo.endIndex)?.lowerBound
            ?? auxInfo.startIndex
    }

    private static func layerValue(in auxInfo: String,
                                   marker: String,
                                   from start: String.Index) -> String? {
        guard let markerRange = auxInfo.range(of: marker, range: start..<auxInfo.endIndex) else {
            return nil
        }
        let valueStart = markerRange.upperBound
        let valueEnd = auxInfo[valueStart...].firstIndex(of: "/") ?? auxInfo.endIndex
        return String(auxInfo[valueStart..<valueEnd])
    }

    private static func parseAtomNumbers(_ component: String,
                                         atomCount: Int) throws -> [Int] {
        try component.split(separator: ",").map { token in
            guard let atomNumber = Int(token), (1...atomCount).contains(atomNumber) else {
                throw ChemError.parseFailed("Invalid atom number '\(token)' in InChI AuxInfo.")
            }
            return atomNumber
        }
    }

    private static func applyAnionicOxygenStartCorrection(numbers: inout [Int],
                                                           firstAtomIndices: [Int],
                                                           molecule: Molecule) {
        for atomIndex in firstAtomIndices where molecule.atoms.indices.contains(atomIndex) {
            let atom = molecule.atoms[atomIndex]
            guard atom.element.uppercased() == "O",
                  atom.charge == -1,
                  molecule.neighbors(of: atom.id).count == 1,
                  let neighborID = molecule.neighbors(of: atom.id).first else {
                continue
            }
            guard let neutralOxygen = molecule.bonds(forAtom: neighborID).compactMap({ bond -> Atom? in
                guard bond.order == .double else { return nil }
                let candidateID = bond.a1 == neighborID ? bond.a2 : bond.a1
                guard let candidate = molecule.atom(id: candidateID),
                      candidate.element.uppercased() == "O",
                      candidate.charge == 0 else {
                    return nil
                }
                return candidate
            }).first,
                let neutralIndex = molecule.indexOfAtom(id: neutralOxygen.id) else {
                continue
            }
            numbers.swapAt(atomIndex, neutralIndex)
        }
    }
}
