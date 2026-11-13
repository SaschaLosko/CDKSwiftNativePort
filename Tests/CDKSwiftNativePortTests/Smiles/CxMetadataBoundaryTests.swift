import XCTest
@testable import CDKSwiftNativePort

final class CxMetadataBoundaryTests: XCTestCase {
    func testInvalidLaterLayerDoesNotPublishEarlierMetadata() throws {
        for suffix in ["r:bad", "C:bad", "unknown:value"] {
            let title = "|$;R1$,ha:0,\(suffix)| original title"
            let result = try CDKCxSmilesParser.split("C* \(title)", enabled: true)
            XCTAssertEqual(result.coreSmiles, "C*")
            XCTAssertEqual(result.title, title)
            XCTAssertEqual(result.state, CDKCxSmilesState())
        }
    }

    func testDisabledCXPreservesNestedDefinitionsAsTitle() throws {
        let title = "|$;R1$,RG:_R1={Cl* |$;_AP1$|}| sample"
        let disabled = try CDKCxSmilesParser.split("C* \(title)", enabled: false)
        XCTAssertEqual(disabled.title, title)
        XCTAssertEqual(disabled.state, CDKCxSmilesState())
        let enabled = try CDKCxSmilesParser.split("C* \(title)", enabled: true)
        XCTAssertEqual(enabled.title, "sample")
        XCTAssertEqual(enabled.state.rGroupDefinitions["R1"], ["Cl* |$;_AP1$|"])
        XCTAssertEqual(enabled.state.atomLabels[1], "R1")
    }
}
