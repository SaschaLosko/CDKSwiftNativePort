import XCTest
@testable import CDKSwiftNativePort

final class FingerprintSimilarityTests: XCTestCase {
    func testBitFingerprintTanimotoDiceAndCosine() throws {
        let first = CDKBitFingerprint(size: 16, bits: [1, 2, 3, 8])
        let second = CDKBitFingerprint(size: 16, bits: [2, 3, 4, 9])

        XCTAssertEqual(try CDKFingerprintSimilarity.tanimoto(first, second), 2.0 / 6.0, accuracy: 1e-12)
        XCTAssertEqual(try CDKFingerprintSimilarity.dice(first, second), 4.0 / 8.0, accuracy: 1e-12)
        XCTAssertEqual(try CDKFingerprintSimilarity.cosine(first, second), 2.0 / 4.0, accuracy: 1e-12)
    }

    func testContinuousTanimotoMatchesCDKReferenceCase() throws {
        let features = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
        XCTAssertEqual(try CDKFingerprintSimilarity.tanimoto(features, features), 1.0, accuracy: 0.001)
    }

    func testCountFingerprintTanimotoMethodsMatchCDKReferenceCases() throws {
        let first = CDKCountFingerprint(countsByHash: [65: 3])
        let second = CDKCountFingerprint(countsByHash: [65: 4])

        XCTAssertEqual(try CDKFingerprintSimilarity.countTanimotoMethod1(first, second), 0.923, accuracy: 0.001)
        XCTAssertEqual(try CDKFingerprintSimilarity.countTanimotoMethod2(first, second), 0.75, accuracy: 0.001)
        XCTAssertEqual(try CDKFingerprintSimilarity.tanimoto(first, second), 0.75, accuracy: 0.001)
    }

    func testRejectsMismatchedBitFingerprintSizes() {
        let first = CDKBitFingerprint(size: 16, bits: [1])
        let second = CDKBitFingerprint(size: 32, bits: [1])

        XCTAssertThrowsError(try CDKFingerprintSimilarity.tanimoto(first, second)) { error in
            XCTAssertEqual(error as? CDKFingerprintSimilarity.SimilarityError, .mismatchedBitFingerprintSize)
        }
    }
}
