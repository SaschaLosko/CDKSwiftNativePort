import XCTest
@testable import CDKSwiftNativePort

final class CMLReaderPortTests: XCTestCase {
    func testParsesBasicCMLMolecule() throws {
        let text = """
        <cml>
          <molecule id="mol1" title="Formaldehyde">
            <atomArray>
              <atom id="a1" elementType="C" x2="0.0" y2="0.0" />
              <atom id="a2" elementType="O" x2="1.2" y2="0.0" />
            </atomArray>
            <bondArray>
              <bond id="b1" atomRefs2="a1 a2" order="2" />
            </bondArray>
          </molecule>
        </cml>
        """

        let molecules = try CDKCMLReader.read(text: text)
        XCTAssertEqual(molecules.count, 1)

        let molecule = molecules[0]
        XCTAssertEqual(molecule.name, "Formaldehyde")
        XCTAssertEqual(molecule.atomCount, 2)
        XCTAssertEqual(molecule.bondCount, 1)
        XCTAssertEqual(molecule.bonds[0].order, .double)
    }

    func testParsesVectorAtomArrayForm() throws {
        let text = """
        <cml>
          <molecule id="vec">
            <atomArray atomID="a1 a2" elementType="C O" x2="0.0 1.3" y2="0.0 0.0" />
            <bondArray atomRef1="a1" atomRef2="a2" order="1" />
          </molecule>
        </cml>
        """

        let molecules = try CDKCMLReader.read(text: text)
        XCTAssertEqual(molecules.count, 1)
        XCTAssertEqual(molecules[0].atomCount, 2)
        XCTAssertEqual(molecules[0].bondCount, 1)
    }

    func testRejectsMalformedCML() {
        XCTAssertThrowsError(try CDKCMLReader.read(text: "<cml><molecule><atomArray></cml>"))
    }
}
