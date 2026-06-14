import XCTest
@testable import CDKSwiftNativePort

final class PathFingerprinterTests: XCTestCase {
    private let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    func testCarbonAtomBitMatchesCDKRandomFoldedReference() throws {
        let methane = try parser.parseSmiles("C")
        let fingerprint = try CDKPathFingerprinter().bitFingerprint(for: methane)

        XCTAssertEqual(fingerprint.size, CDKPathFingerprinter.defaultSize)
        XCTAssertEqual(fingerprint.bits, [742])
    }

    func testSubstructureCandidateScreensInBeforeExactMatch() throws {
        let query = try parser.parseSmiles("CCO")
        let target = try parser.parseSmiles("CCCO")
        let fingerprinter = CDKSubstructureScreenFingerprinter()

        let queryFingerprint = try fingerprinter.bitFingerprint(for: query)
        let targetFingerprint = try fingerprinter.bitFingerprint(for: target)

        XCTAssertTrue(targetFingerprint.isSuperset(of: queryFingerprint))
        XCTAssertTrue(target.containsSubstructure(query))
    }

    func testMissingPathScreensOutNonCandidate() throws {
        let query = try parser.parseSmiles("C=O")
        let target = try parser.parseSmiles("CCO")
        let fingerprinter = CDKPathFingerprinter()

        let queryFingerprint = try fingerprinter.bitFingerprint(for: query)
        let targetFingerprint = try fingerprinter.bitFingerprint(for: target)

        XCTAssertFalse(targetFingerprint.isSuperset(of: queryFingerprint))
        XCTAssertFalse(target.containsSubstructure(query))
    }

    func testDefaultScreenIgnoresExplicitHydrogens() throws {
        let explicitHydrogenQuery = try parser.parseSmiles("[H]O")
        let target = try parser.parseSmiles("CO")
        let fingerprinter = CDKPathFingerprinter()

        let queryFingerprint = try fingerprinter.bitFingerprint(for: explicitHydrogenQuery)
        let targetFingerprint = try fingerprinter.bitFingerprint(for: target)

        XCTAssertTrue(targetFingerprint.isSuperset(of: queryFingerprint))
    }

    func testPathLimitOverflowThrows() throws {
        let propane = try parser.parseSmiles("CCC")
        let fingerprinter = CDKPathFingerprinter(pathLimit: 1)

        XCTAssertThrowsError(try fingerprinter.bitFingerprint(for: propane)) { error in
            XCTAssertEqual(error as? CDKPathFingerprinterError, .pathLimitExceeded(limit: 1))
        }
    }
}
