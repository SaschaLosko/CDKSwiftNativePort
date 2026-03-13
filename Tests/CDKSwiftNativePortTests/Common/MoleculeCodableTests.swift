import XCTest
@testable import CDKSwiftNativePort

final class MoleculeCodableTests: XCTestCase {

    func testDecodesLegacyPayloadWithoutDataFields() throws {
        let data = Data(#"{"name":"Legacy","atoms":[],"bonds":[]}"#.utf8)

        let molecule = try JSONDecoder().decode(Molecule.self, from: data)

        XCTAssertEqual(molecule.name, "Legacy")
        XCTAssertTrue(molecule.dataFields.isEmpty)
        XCTAssertTrue(molecule.orderedDataFieldNames.isEmpty)
    }
}
