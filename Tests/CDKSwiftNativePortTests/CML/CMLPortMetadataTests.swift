import XCTest

private struct CDKCMLPortMetadata: Decodable {
    struct PortedCase: Decodable {
        let id: String
        let source: String
    }

    let schemaVersion: Int
    let suite: String
    let cdkReferenceVersion: String
    let cdkReferenceTag: String
    let sourceTests: [String]
    let portedCases: [PortedCase]
}

final class CMLPortMetadataTests: XCTestCase {
    func testPortMetadataExistsAndReferencesCDKTests() throws {
        let url = URL(fileURLWithPath: #filePath)
            .deletingLastPathComponent()
            .appendingPathComponent("port_metadata.json")

        let data = try Data(contentsOf: url)
        let metadata = try JSONDecoder().decode(CDKCMLPortMetadata.self, from: data)

        XCTAssertGreaterThanOrEqual(metadata.schemaVersion, 1)
        XCTAssertFalse(metadata.suite.isEmpty)
        XCTAssertEqual(metadata.cdkReferenceVersion, "2.13")
        XCTAssertEqual(metadata.cdkReferenceTag, "cdk-2.13")

        XCTAssertTrue(metadata.sourceTests.contains(where: { $0.contains("CMLReaderTest.java") }))
        XCTAssertTrue(metadata.sourceTests.contains(where: { $0.contains("CML23FragmentsTest.java") }))
        XCTAssertTrue(metadata.sourceTests.contains(where: { $0.contains("CML2Test.java") }))
        XCTAssertTrue(metadata.sourceTests.contains(where: { $0.contains("CML2WriterTest.java") }))
        XCTAssertTrue(metadata.sourceTests.contains(where: { $0.contains("CMLRoundTripTest.java") }))

        XCTAssertGreaterThanOrEqual(metadata.portedCases.count, 20)
        XCTAssertTrue(metadata.portedCases.allSatisfy { !$0.id.isEmpty && !$0.source.isEmpty })

        let uniqueIDs = Set(metadata.portedCases.map(\.id))
        XCTAssertEqual(uniqueIDs.count, metadata.portedCases.count)
    }
}
