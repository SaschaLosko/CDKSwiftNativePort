import Foundation

/// A display-only result. Never substitute its graph for a strict parse of the source.
/// This type deliberately has no persistence or export conformance.
public struct CDKSmilesPreview: Sendable {
    public enum Status: String, Sendable { case empty, complete, incomplete, invalid, unsupported, limited }
    public enum Issue: String, Sendable { case atom, ringIndex, ringClosure, branch, bond, fragment }

    public let status: Status
    /// The atoms actually read; no missing atom, branch endpoint or ring bond is invented.
    /// Coordinates are not generated. Incomplete graphs have their stereo assignments cleared.
    public let moleculeForDepiction: Molecule?
    public let issues: [Issue]
    public let diagnostic: String?

    public static let maximumInputBytes = 1024
    public static let maximumAtoms = 96

    /// Bounded preview of one ordinary SMILES record (an optional title is ignored).
    /// Reaction, CXSMILES and multiline sources remain on the normal import path.
    /// Throws only cancellation; syntax failures become explicit invalid results.
    public static func evaluate(_ input: String) throws -> Self {
        try Task.checkCancellation()
        guard input.utf8.count <= maximumInputBytes else { return .init(status: .limited) }
        let text = input.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !text.isEmpty else { return .init(status: .empty) }
        guard !text.contains(where: { $0.isNewline }), !text.contains("|"), !text.contains(">"),
            !text.hasPrefix("InChI="), !text.hasPrefix("<"), !text.hasPrefix("{")
        else { return .init(status: .unsupported) }
        let core = String(text.prefix { !$0.isWhitespace })
        let context = CDKSmilesPreviewContext()
        do {
            let molecule = try CDKSmilesParser(generateCoordinates: false).parseCoreSmiles(core, preview: context)
            guard
                molecule.atoms.allSatisfy({ $0.element == "*" || CDKIsotopes.atomicNumber(symbol: $0.element) != nil })
            else { return .init(status: .invalid, diagnostic: "Unknown element symbol.") }
            return .init(
                status: context.issues.isEmpty ? .complete : .incomplete,
                molecule: molecule.atoms.isEmpty ? nil : molecule, issues: context.issues)
        } catch is CancellationError {
            throw CancellationError()
        } catch CDKSmilesPreviewInterruption.limit {
            return .init(status: .limited)
        } catch {
            return .init(status: .invalid, diagnostic: error.localizedDescription)
        }
    }

    private init(status: Status, molecule: Molecule? = nil, issues: [Issue] = [], diagnostic: String? = nil) {
        self.status = status
        self.moleculeForDepiction = molecule
        self.issues = issues
        self.diagnostic = diagnostic
    }
}

// Internal recovery is available only through the separate preview entry point.
final class CDKSmilesPreviewContext {
    var issues: [CDKSmilesPreview.Issue] = []
    func check(atomCount: Int) throws {
        try Task.checkCancellation()
        guard atomCount <= CDKSmilesPreview.maximumAtoms else { throw CDKSmilesPreviewInterruption.limit }
    }
}

enum CDKSmilesPreviewInterruption: Error {
    case incomplete(CDKSmilesPreview.Issue)
    case limit
}
