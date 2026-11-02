import Foundation

/// Explicit, unsubstituted, acyclic saturated classes. No source prose is interpreted.
/// All constitutional attachment isomers are included; stereochemistry is unconstrained.
public struct CDKFiniteChemicalClass: Codable, Equatable, Sendable {
    public enum Kind: String, Codable, CaseIterable, Sendable { case alkyl, alkoxy }
    public var kind: Kind
    public var minimumCarbons: Int
    public var maximumCarbons: Int

    public init(kind: Kind, minimumCarbons: Int, maximumCarbons: Int) {
        self.kind = kind; self.minimumCarbons = minimumCarbons; self.maximumCarbons = maximumCarbons
    }

    /// Rooted fragment SMILES, in carbon-count then canonical-tree order. The first atom is the attachment.
    /// Carbon limits are deliberately finite; broader ranges are rejected, never truncated.
    public func fragments() throws -> [String] {
        guard minimumCarbons >= 1, maximumCarbons <= 6, minimumCarbons <= maximumCarbons else {
            throw CDKRGroupEnumerationError.unsupportedFeature("Finite alkyl and alkoxy classes require 1–6 carbons.")
        }
        // A carbon root has room for three carbon children and its external bond.
        // Every child has the same constraint after accounting for its parent bond.
        // Unordered rooted ternary trees enumerate constitutional alkyl groups (OEIS A000598).
        var trees: [Int: [String]] = [1: ["C"]]
        if maximumCarbons >= 2 {
            for size in 2...maximumCarbons {
                try Task.checkCancellation()
                let available = (1..<size).flatMap { count in
                    trees[count, default: []].map { (count, $0) }
                }
                var generated: [String] = []
                func children(_ remaining: Int, _ start: Int, _ chosen: [String]) throws {
                    try Task.checkCancellation()
                    if remaining == 0 {
                        generated.append("C" + chosen.map { "(\($0))" }.joined())
                        return
                    }
                    guard chosen.count < 3 else { return }
                    for index in start..<available.count where available[index].0 <= remaining {
                        try children(remaining - available[index].0, index, chosen + [available[index].1])
                    }
                }
                try children(size - 1, 0, [])
                trees[size] = generated.sorted()
            }
        }
        return (minimumCarbons...maximumCarbons).flatMap { count in
            trees[count, default: []].map { kind == .alkoxy ? "O\($0)" : $0 }
        }
    }
}
