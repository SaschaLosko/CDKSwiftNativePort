import Foundation

public struct CDKCxSmilesState: Equatable {
    public var atomLabels: [Int: String] = [:]
    public var fragmentGroups: [[Int]] = []
    public var racemic: Bool = false
    public var racemicFragments: [Int] = []

    public init(atomLabels: [Int: String] = [:],
                fragmentGroups: [[Int]] = [],
                racemic: Bool = false,
                racemicFragments: [Int] = []) {
        self.atomLabels = atomLabels
        self.fragmentGroups = fragmentGroups
        self.racemic = racemic
        self.racemicFragments = racemicFragments
    }
}

public enum CDKCxSmilesParser {
    public struct SplitResult: Equatable {
        public let coreSmiles: String
        public let title: String?
        public let state: CDKCxSmilesState

        public init(coreSmiles: String, title: String?, state: CDKCxSmilesState) {
            self.coreSmiles = coreSmiles
            self.title = title
            self.state = state
        }
    }

    public static func split(_ input: String, enabled: Bool) throws -> SplitResult {
        let trimmed = input.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { throw ChemError.emptyInput }

        guard enabled else {
            return SplitResult(coreSmiles: trimmed, title: nil, state: CDKCxSmilesState())
        }

        guard let firstPipe = trimmed.firstIndex(of: "|") else {
            return SplitResult(coreSmiles: trimmed, title: nil, state: CDKCxSmilesState())
        }

        let afterFirst = trimmed.index(after: firstPipe)
        guard let secondPipe = trimmed[afterFirst...].firstIndex(of: "|") else {
            throw ChemError.parseFailed("Unterminated CXSMILES layer (missing closing '|').")
        }

        let core = String(trimmed[..<firstPipe]).trimmingCharacters(in: .whitespacesAndNewlines)
        let cxBody = String(trimmed[afterFirst..<secondPipe])
        let trailing = String(trimmed[trimmed.index(after: secondPipe)...]).trimmingCharacters(in: .whitespacesAndNewlines)

        guard !core.isEmpty else {
            throw ChemError.parseFailed("Missing core SMILES before CXSMILES layer.")
        }
        if trailing.contains("|") {
            throw ChemError.parseFailed("Malformed CXSMILES tail.")
        }

        let state = try parseLayers(cxBody)
        let title = trailing.isEmpty ? nil : trailing
        return SplitResult(coreSmiles: core, title: title, state: state)
    }

    public static func applyAtomLabels(to molecule: inout Molecule, state: CDKCxSmilesState) {
        guard !state.atomLabels.isEmpty else { return }
        guard !molecule.atoms.isEmpty else { return }

        for (atomIndex, label) in state.atomLabels {
            guard atomIndex >= 0, atomIndex < molecule.atoms.count else { continue }
            let atomID = molecule.atoms[atomIndex].id
            guard let idx = molecule.atoms.firstIndex(where: { $0.id == atomID }) else { continue }

            let existing = molecule.atoms[idx]
            molecule.atoms[idx] = Atom(id: existing.id,
                                       element: label,
                                       position: existing.position,
                                       charge: existing.charge,
                                       isotopeMassNumber: existing.isotopeMassNumber,
                                       aromatic: false,
                                       chirality: existing.chirality,
                                       explicitHydrogenCount: existing.explicitHydrogenCount)
        }
    }

    private static func parseLayers(_ raw: String) throws -> CDKCxSmilesState {
        var state = CDKCxSmilesState()
        guard !raw.isEmpty else { return state }

        let layers = splitTopLevelLayers(raw)
        for layer in layers {
            let token = layer.trimmingCharacters(in: .whitespacesAndNewlines)
            guard !token.isEmpty else { continue }

            if token.hasPrefix("$") {
                guard token.hasSuffix("$"), token.count >= 2 else {
                    throw ChemError.parseFailed("Malformed CXSMILES atom-label layer.")
                }
                parseAtomLabels(token, into: &state)
                continue
            }

            if token == "r" {
                state.racemic = true
                continue
            }

            if token.hasPrefix("r:") {
                let body = String(token.dropFirst(2))
                let values = body.split(separator: ",").compactMap { Int($0) }
                if values.isEmpty && !body.isEmpty {
                    throw ChemError.parseFailed("Malformed CXSMILES racemic-fragment layer.")
                }
                state.racemicFragments = values
                continue
            }

            if token.hasPrefix("f:") {
                let body = String(token.dropFirst(2))
                if body.isEmpty { continue }
                let groups = body.split(separator: ",")
                for group in groups {
                    let ids = group.split(separator: ".").compactMap { Int($0) }
                    if ids.isEmpty {
                        throw ChemError.parseFailed("Malformed CXSMILES fragment-group layer.")
                    }
                    state.fragmentGroups.append(ids)
                }
                continue
            }
        }

        return state
    }

    private static func parseAtomLabels(_ token: String, into state: inout CDKCxSmilesState) {
        let content = String(token.dropFirst().dropLast())
        let entries = content.split(separator: ";", omittingEmptySubsequences: false)
        for (idx, rawEntry) in entries.enumerated() {
            guard !rawEntry.isEmpty else { continue }
            var label = unescape(String(rawEntry))
            if label.hasPrefix("_") {
                label = String(label.dropFirst())
            }
            state.atomLabels[idx] = label
        }
    }

    private static func splitTopLevelLayers(_ input: String) -> [String] {
        var out: [String] = []
        var current = ""
        var inAtomLabel = false
        var parenDepth = 0

        for ch in input {
            if ch == "$" {
                inAtomLabel.toggle()
                current.append(ch)
                continue
            }
            if !inAtomLabel {
                if ch == "(" { parenDepth += 1 }
                if ch == ")" { parenDepth = max(0, parenDepth - 1) }
                if ch == "," && parenDepth == 0 {
                    out.append(current)
                    current.removeAll(keepingCapacity: true)
                    continue
                }
            }
            current.append(ch)
        }
        out.append(current)
        return out
    }

    // Minimal numeric HTML entity decoding used by CXSMILES escaped labels.
    private static func unescape(_ input: String) -> String {
        var out = input
        while let start = out.range(of: "&#"), let semi = out[start.upperBound...].firstIndex(of: ";") {
            let numeric = out[start.upperBound..<semi]
            if let value = Int(numeric), let scalar = UnicodeScalar(value) {
                out.replaceSubrange(start.lowerBound...semi, with: String(Character(scalar)))
            } else {
                break
            }
        }
        return out
    }
}
