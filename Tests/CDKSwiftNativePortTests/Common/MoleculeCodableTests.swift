import XCTest
@testable import CDKSwiftNativePort

final class MoleculeCodableTests: XCTestCase {

    func testCoordinateProvenanceSurvivesCoding() throws {
        for generated: Bool? in [nil, false, true] {
            var molecule = try CDKSmilesParser().parseSmiles("FC(Cl)Br")
            molecule.coordinatesAreGenerated = generated
            let data = try JSONEncoder().encode(molecule)
            let restored = try JSONDecoder().decode(Molecule.self, from: data)
            XCTAssertEqual(restored.coordinatesAreGenerated, generated)
            XCTAssertEqual(restored, molecule)
        }
    }

    func testDecodesLegacyPayloadWithoutDataFields() throws {
        let data = Data(#"{"name":"Legacy","atoms":[],"bonds":[]}"#.utf8)

        let molecule = try JSONDecoder().decode(Molecule.self, from: data)

        XCTAssertEqual(molecule.name, "Legacy")
        XCTAssertTrue(molecule.dataFields.isEmpty)
        XCTAssertTrue(molecule.orderedDataFieldNames.isEmpty)
        XCTAssertNil(molecule.coordinatesAreGenerated)
    }
}
