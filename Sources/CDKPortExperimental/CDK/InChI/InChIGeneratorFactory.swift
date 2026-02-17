import Foundation

/// Swift counterpart of CDK's InChIGeneratorFactory facade.
public final class CDKInChIGeneratorFactory {
    public static let shared = CDKInChIGeneratorFactory()

    private init() {}

    public func getInChIToStructure(_ inchi: String) -> CDKInChIToStructure {
        CDKInChIToStructure(inchi: inchi)
    }
}

public enum CDKInChIStatus: Equatable {
    case success
    case warning
    case error
}
