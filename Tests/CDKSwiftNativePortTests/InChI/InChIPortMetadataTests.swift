import XCTest

private struct CDKInChIPortMetadata: Decodable {
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

final class InChIPortMetadataTests: XCTestCase {

    func testPortMetadataExistsAndReferencesCDKTests() throws {
        let url = URL(fileURLWithPath: #filePath)
            .deletingLastPathComponent()
            .appendingPathComponent("port_metadata.json")

        let data = try Data(contentsOf: url)
        let metadata = try JSONDecoder().decode(CDKInChIPortMetadata.self, from: data)

        XCTAssertGreaterThanOrEqual(metadata.schemaVersion, 1)
        XCTAssertFalse(metadata.suite.isEmpty)
        XCTAssertEqual(metadata.cdkReferenceVersion, "2.11")
        XCTAssertEqual(metadata.cdkReferenceTag, "cdk-2.11")

        XCTAssertTrue(metadata.sourceTests.contains(where: { $0.contains("InChIGeneratorTest.java") }))
        XCTAssertTrue(metadata.sourceTests.contains(where: { $0.contains("InChIToStructureTest.java") }))
        XCTAssertGreaterThanOrEqual(metadata.portedCases.count, 8)
        XCTAssertTrue(metadata.portedCases.allSatisfy { !$0.id.isEmpty && !$0.source.isEmpty })

        let uniqueIDs = Set(metadata.portedCases.map(\.id))
        XCTAssertEqual(uniqueIDs.count, metadata.portedCases.count)
    }
}
