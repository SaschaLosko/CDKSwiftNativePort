import Foundation

/// Swift counterpart of CDK's `SmilesParser` construction/facade entry point.
public final class CDKSmilesParserFactory: @unchecked Sendable {
    public static let shared = CDKSmilesParserFactory()

    private init() {}

    public func newSmilesParser(
        flavor: CDKSmiFlavor = .cdkDefault, generateCoordinates: Bool = true
    ) -> CDKSmilesParser {
        CDKSmilesParser(flavor: flavor, generateCoordinates: generateCoordinates)
    }
}
