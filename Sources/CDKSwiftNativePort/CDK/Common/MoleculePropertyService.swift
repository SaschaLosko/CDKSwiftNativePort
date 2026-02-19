import Foundation

public struct CDKRuleOfFiveResult: Equatable {
    public let molecularWeight: Double
    public let hBondDonorCount: Int
    public let hBondAcceptorCount: Int
    public let xlogP: Double?
    public let violations: Int
    public let evaluatedCriteriaCount: Int
    public let isCompliant: Bool?

    public init(molecularWeight: Double,
                hBondDonorCount: Int,
                hBondAcceptorCount: Int,
                xlogP: Double?,
                violations: Int,
                evaluatedCriteriaCount: Int,
                isCompliant: Bool?) {
        self.molecularWeight = molecularWeight
        self.hBondDonorCount = hBondDonorCount
        self.hBondAcceptorCount = hBondAcceptorCount
        self.xlogP = xlogP
        self.violations = violations
        self.evaluatedCriteriaCount = evaluatedCriteriaCount
        self.isCompliant = isCompliant
    }

    public var statusText: String {
        let suffix = violations == 1 ? "" : "s"
        if let isCompliant {
            let status = isCompliant ? "Pass" : "Fail"
            return "\(status) (\(violations) violation\(suffix))"
        }
        return "Partial (\(violations) violation\(suffix), XLogP unavailable)"
    }
}

public struct CDKMolecularProperties: Equatable {
    public let formula: String
    public let molecularWeight: Double
    public let molarMass: Double
    public let monoisotopicMass: Double
    public let heavyAtomCount: Int
    public let hBondDonorCount: Int
    public let hBondAcceptorCount: Int
    public let rotatableBondCount: Int
    public let ringCount: Int
    public let ruleOfFive: CDKRuleOfFiveResult

    public init(formula: String,
                molecularWeight: Double,
                molarMass: Double,
                monoisotopicMass: Double,
                heavyAtomCount: Int,
                hBondDonorCount: Int,
                hBondAcceptorCount: Int,
                rotatableBondCount: Int,
                ringCount: Int,
                ruleOfFive: CDKRuleOfFiveResult) {
        self.formula = formula
        self.molecularWeight = molecularWeight
        self.molarMass = molarMass
        self.monoisotopicMass = monoisotopicMass
        self.heavyAtomCount = heavyAtomCount
        self.hBondDonorCount = hBondDonorCount
        self.hBondAcceptorCount = hBondAcceptorCount
        self.rotatableBondCount = rotatableBondCount
        self.ringCount = ringCount
        self.ruleOfFive = ruleOfFive
    }
}

public enum CDKRuleOfFiveDescriptor {
    public static func evaluate(for molecule: Molecule, xlogP: Double? = nil) -> CDKRuleOfFiveResult {
        let molecularWeight = CDKMolecularWeightDescriptor.calculate(for: molecule)
        let donors = CDKHBondDonorCountDescriptor.calculate(for: molecule)
        let acceptors = CDKHBondAcceptorCountDescriptor.calculate(for: molecule)
        let resolvedXlogP = xlogP ?? CDKXLogPDescriptor.calculate(for: molecule)
        return evaluate(molecularWeight: molecularWeight,
                        hBondDonorCount: donors,
                        hBondAcceptorCount: acceptors,
                        xlogP: resolvedXlogP)
    }

    public static func evaluate(molecularWeight: Double,
                                hBondDonorCount: Int,
                                hBondAcceptorCount: Int,
                                xlogP: Double? = nil) -> CDKRuleOfFiveResult {
        var violations = 0
        if molecularWeight > 500.0 { violations += 1 }
        if hBondDonorCount > 5 { violations += 1 }
        if hBondAcceptorCount > 10 { violations += 1 }
        if let xlogP, xlogP > 5.0 { violations += 1 }

        let criteriaCount = xlogP == nil ? 3 : 4
        let isCompliant = xlogP == nil ? nil : (violations <= 1)

        return CDKRuleOfFiveResult(molecularWeight: molecularWeight,
                                   hBondDonorCount: hBondDonorCount,
                                   hBondAcceptorCount: hBondAcceptorCount,
                                   xlogP: xlogP,
                                   violations: violations,
                                   evaluatedCriteriaCount: criteriaCount,
                                   isCompliant: isCompliant)
    }
}

public enum CDKMoleculePropertyService {
    public static func compute(for molecule: Molecule, xlogP: Double? = nil) -> CDKMolecularProperties {
        let formula = CDKMolecularFormulaDescriptor.calculate(for: molecule)
        let molecularWeight = CDKMolecularWeightDescriptor.calculate(for: molecule)
        let monoisotopicMass = CDKExactMassDescriptor.calculate(for: molecule)
        let heavyAtomCount = CDKHeavyAtomCountDescriptor.calculate(for: molecule)
        let hBondDonorCount = CDKHBondDonorCountDescriptor.calculate(for: molecule)
        let hBondAcceptorCount = CDKHBondAcceptorCountDescriptor.calculate(for: molecule)
        let rotatableBondCount = CDKRotatableBondsCountDescriptor.calculate(for: molecule)
        let ringCount = CDKRingCountDescriptor.calculate(for: molecule)
        let resolvedXlogP = xlogP ?? CDKXLogPDescriptor.calculate(for: molecule)
        let ruleOfFive = CDKRuleOfFiveDescriptor.evaluate(molecularWeight: molecularWeight,
                                                          hBondDonorCount: hBondDonorCount,
                                                          hBondAcceptorCount: hBondAcceptorCount,
                                                          xlogP: resolvedXlogP)

        return CDKMolecularProperties(formula: formula,
                                      molecularWeight: molecularWeight,
                                      molarMass: molecularWeight,
                                      monoisotopicMass: monoisotopicMass,
                                      heavyAtomCount: heavyAtomCount,
                                      hBondDonorCount: hBondDonorCount,
                                      hBondAcceptorCount: hBondAcceptorCount,
                                      rotatableBondCount: rotatableBondCount,
                                      ringCount: ringCount,
                                      ruleOfFive: ruleOfFive)
    }
}
