import Foundation

#if canImport(CoreGraphics)
    import CoreGraphics
#endif

public enum ChemError: LocalizedError {
    case emptyInput
    case unsupported(String)
    case parseFailed(String)

    public var errorDescription: String? {
        switch self {
        case .emptyInput:
            return "No input provided."
        case .unsupported(let s):
            return "Unsupported: \(s)"
        case .parseFailed(let s):
            return "Could not parse: \(s)"
        }
    }
}
