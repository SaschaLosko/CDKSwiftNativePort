import Foundation

/// Swift counterpart of CDK's `SmilesParser` construction/facade entry point.
public final class CDKSmilesParserFactory {
    public static let shared = CDKSmilesParserFactory()

    private init() {}

    public func newSmilesParser(flavor: CDKSmiFlavor = .cdkDefault) -> CDKSmilesParser {
        CDKSmilesParser(flavor: flavor)
    }
}
