// Isotope values derived from CDK 2.14.1 / BODR 10.
// Copyright (C) 2001-2007 Christoph Steinbeck; 2013,2016 Egon Willighagen;
// 2026 John Mayfield. CDK adaptation licensed under LGPL-2.1-or-later.
// Dataset provenance and precision are documented in Documentation/ISOTOPES.md.
import Foundation

/// Isotopic mass in Da and natural abundance in percent (0...100).
/// Values retain the precision of CDK 2.14.1's BODR 10 dataset.
public struct CDKIsotope: Codable, Equatable, Sendable {
    public let atomicNumber: Int
    public let massNumber: Int
    public let exactMass: Double
    public let naturalAbundance: Double
}

public enum CDKIsotopeError: Error, Equatable, Sendable, LocalizedError {
    case dataUnavailable
    case unsupportedMass(element: String, massNumber: Int?)

    public var errorDescription: String? {
        switch self {
        case .dataUnavailable: "The bundled isotope dataset is unavailable or invalid."
        case .unsupportedMass(let element, let number):
            "No supported mass data for \(number.map(String.init) ?? "natural") \(element)."
        }
    }
}

/// Immutable, shared isotope data. Missing physical data returns nil; a missing
/// or invalid packaged dataset throws, rather than masquerading as an unknown isotope.
/// CTfile mass-difference reference numbers are a separate format-specific table.
public enum CDKIsotopes {
    public static let datasetVersion = "CDK 2.14.1 / BODR 10"

    public static func isotope(atomicNumber: Int, massNumber: Int) throws -> CDKIsotope? {
        try catalog.get().byElement[atomicNumber]?[massNumber]
    }

    public static func isotopes(atomicNumber: Int) throws -> [CDKIsotope] {
        try (catalog.get().byElement[atomicNumber]?.values.map { $0 } ?? []).sorted {
            $0.massNumber < $1.massNumber
        }
    }

    /// Most abundant naturally occurring isotope; nil where no natural abundance is supplied.
    public static func majorIsotope(atomicNumber: Int) throws -> CDKIsotope? {
        try catalog.get().major[atomicNumber]
    }

    public static func averageMass(atomicNumber: Int) throws -> Double? {
        try catalog.get().averages[atomicNumber]
    }

    public static func atomicNumber(symbol: String) -> Int? {
        symbolIndices[CDKDescriptorSupport.canonicalElementSymbol(symbol)]
    }

    /// Explicit labels use their isotope's exact mass in both mass conventions.
    public static func mass(element: String, massNumber: Int? = nil, average: Bool = false) throws -> Double {
        guard let number = atomicNumber(symbol: element) else {
            throw CDKIsotopeError.unsupportedMass(element: element, massNumber: massNumber)
        }
        let value: Double?
        if let massNumber {
            value = try isotope(atomicNumber: number, massNumber: massNumber)?.exactMass
        } else if average {
            value = try averageMass(atomicNumber: number)
        } else {
            value = try majorIsotope(atomicNumber: number)?.exactMass
        }
        guard let value else { throw CDKIsotopeError.unsupportedMass(element: element, massNumber: massNumber) }
        return value
    }

    private struct Source: Decodable {
        struct Average: Decodable { let symbol: String; let average: Double }
        let source: String
        let averageMasses: [Average]
        let isotopes: [CDKIsotope]
    }

    private struct Catalog: Sendable {
        var byElement: [Int: [Int: CDKIsotope]] = [:]
        var major: [Int: CDKIsotope] = [:]
        var averages: [Int: Double] = [:]
    }

    // One initialization shared safely by all threads. No runtime network or mutable cache.
    private static let catalog: Result<Catalog, CDKIsotopeError> = {
        do {
            guard let url = Bundle.module.url(forResource: "Isotopes", withExtension: "json") else {
                return .failure(.dataUnavailable)
            }
            let source = try JSONDecoder().decode(Source.self, from: Data(contentsOf: url))
            guard source.source == datasetVersion, source.averageMasses.count == 92, !source.isotopes.isEmpty else {
                return .failure(.dataUnavailable)
            }
            var result = Catalog()
            for (index, entry) in source.averageMasses.enumerated() {
                guard atomicNumber(symbol: entry.symbol) == index + 1, entry.average.isFinite, entry.average >= 0 else {
                    return .failure(.dataUnavailable)
                }
                if entry.average > 0 { result.averages[index + 1] = entry.average }
            }
            for isotope in source.isotopes {
                guard (1...118).contains(isotope.atomicNumber), isotope.massNumber > 0,
                    isotope.exactMass.isFinite, isotope.exactMass > 0,
                    (0...100).contains(isotope.naturalAbundance),
                    result.byElement[isotope.atomicNumber]?[isotope.massNumber] == nil
                else { return .failure(.dataUnavailable) }
                result.byElement[isotope.atomicNumber, default: [:]][isotope.massNumber] = isotope
                if isotope.naturalAbundance > (result.major[isotope.atomicNumber]?.naturalAbundance ?? 0) {
                    result.major[isotope.atomicNumber] = isotope
                }
            }
            return .success(result)
        } catch { return .failure(.dataUnavailable) }
    }()

    private static let symbolIndices: [String: Int] = Dictionary(
        uniqueKeysWithValues: ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
            + "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd "
            + "Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac "
            + "Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og").split(
                separator: " "
            ).enumerated().map { (String($0.element), $0.offset + 1) })
}
