import XCTest

private struct CDKReactionPortMetadata: Decodable {
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

final class ReactionPortMetadataTests: XCTestCase {
    func testPortMetadataTracksUpstreamReactionSurface() throws {
        let url = URL(fileURLWithPath: #filePath)
            .deletingLastPathComponent()
            .appendingPathComponent("port_metadata.json")

        let data = try Data(contentsOf: url)
        let metadata = try JSONDecoder().decode(CDKReactionPortMetadata.self, from: data)

        XCTAssertGreaterThanOrEqual(metadata.schemaVersion, 1)
        XCTAssertEqual(metadata.cdkReferenceVersion, "2.13")
        XCTAssertEqual(metadata.cdkReferenceTag, "cdk-2.13")
        XCTAssertTrue(metadata.sourceTests.contains(where: { $0.contains("MDLRXNReaderTest.java") }))
        XCTAssertTrue(
            metadata.sourceTests.contains(where: { $0.contains("MDLRXNV2000ReaderTest.java") }))
        XCTAssertTrue(
            metadata.sourceTests.contains(where: { $0.contains("MDLRXNV3000ReaderTest.java") }))
        XCTAssertTrue(metadata.sourceTests.contains(where: { $0.contains("MDLRXNWriterTest.java") }))
        XCTAssertTrue(
            metadata.sourceTests.contains(where: { $0.contains("ReactionManipulatorTest.java") }))
        XCTAssertTrue(
            metadata.sourceTests.contains(where: { $0.contains("ReactionSetManipulatorTest.java") }))
        XCTAssertTrue(
            metadata.sourceTests.contains(where: { $0.contains("ReactionSchemeManipulatorTest.java") }))

        XCTAssertGreaterThanOrEqual(metadata.portedCases.count, 20)
        XCTAssertEqual(Set(metadata.portedCases.map(\.id)).count, metadata.portedCases.count)
        XCTAssertTrue(metadata.portedCases.allSatisfy { !$0.source.isEmpty })
    }
}
