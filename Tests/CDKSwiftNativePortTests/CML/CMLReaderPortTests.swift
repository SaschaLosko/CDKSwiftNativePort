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

    func testParsesListWithMultipleMolecules() throws {
        let text = """
        <list>
          <molecule id="m1">
            <atomArray><atom id="a1" elementType="C" /></atomArray>
          </molecule>
          <molecule id="m2">
            <atomArray><atom id="a1" elementType="O" /></atomArray>
          </molecule>
        </list>
        """

        let molecules = try CDKCMLReader.read(text: text)
        XCTAssertEqual(molecules.count, 2)
        XCTAssertEqual(molecules[0].name, "m1")
        XCTAssertEqual(molecules[1].name, "m2")
    }

    func testParsesVector3DCoordinates() throws {
        let text = """
        <molecule id="m1">
          <atomArray atomID="a1 a2" x3="0.0 0.1" y3="1.2 1.3" z3="2.1 2.5" />
        </molecule>
        """

        let molecules = try CDKCMLReader.read(text: text)
        XCTAssertEqual(molecules[0].atoms[0].position.x, 0.0, accuracy: 0.00001)
        XCTAssertEqual(molecules[0].atoms[0].position.y, 1.2, accuracy: 0.00001)
        XCTAssertEqual(try XCTUnwrap(molecules[0].atoms[0].zPosition), 2.1, accuracy: 0.00001)
        XCTAssertEqual(try XCTUnwrap(molecules[0].atoms[1].zPosition), 2.5, accuracy: 0.00001)
    }

    func testParsesXYZ3OnIndividualAtoms() throws {
        let text = """
        <molecule id="m1">
          <atomArray>
            <atom id="a1" xyz3="0.0 0.1 0.2" />
            <atom id="a2" />
            <atom id="a3" xyz3="0.1 0.0 0.2" />
          </atomArray>
        </molecule>
        """

        let molecules = try CDKCMLReader.read(text: text)
        XCTAssertNotNil(molecules[0].atoms[0].zPosition)
        XCTAssertNil(molecules[0].atoms[1].zPosition)
        XCTAssertNotNil(molecules[0].atoms[2].zPosition)
    }

    func testParsesBondStereoFromDictRefAndCharContent() throws {
        let text = """
        <molecule id="m1">
          <atomArray atomID="a1 a2 a3 a4 a5" />
          <bondArray>
            <bond atomRefs2="a1 a2" order="S"><bondStereo dictRef="cml:H" /></bond>
            <bond atomRefs2="a2 a3" order="S"><bondStereo>W</bondStereo></bond>
            <bond atomRefs2="a3 a4" order="S"><bondStereo dictRef="cml:W">H</bondStereo></bond>
            <bond atomRefs2="a4 a5" order="S"><bondStereo dictRef="cml:">W</bondStereo></bond>
          </bondArray>
        </molecule>
        """

        let molecule = try XCTUnwrap(CDKCMLReader.read(text: text).first)
        XCTAssertEqual(molecule.bonds.map(\.stereo), [.down, .up, .up, .up])
    }

    func testParsesAtomParityAndInfersDisplayStereoAfterLayout() throws {
        let text = """
        <cml>
          <molecule id="trans12dichlorocyclohexane">
            <atomArray>
              <atom id="a1" elementType="C">
                <atomParity atomRefs4="a7 a2 a6 a8">1</atomParity>
              </atom>
              <atom id="a2" elementType="C">
                <atomParity atomRefs4="a8 a3 a1 a7">-1</atomParity>
              </atom>
              <atom id="a3" elementType="C" />
              <atom id="a4" elementType="C" />
              <atom id="a5" elementType="C" />
              <atom id="a6" elementType="C" />
              <atom id="a7" elementType="Cl" />
              <atom id="a8" elementType="Cl" />
            </atomArray>
            <bondArray>
              <bond atomRefs2="a1 a2" order="S" />
              <bond atomRefs2="a2 a3" order="S" />
              <bond atomRefs2="a3 a4" order="S" />
              <bond atomRefs2="a4 a5" order="S" />
              <bond atomRefs2="a5 a6" order="S" />
              <bond atomRefs2="a6 a1" order="S" />
              <bond atomRefs2="a1 a7" order="S" />
              <bond atomRefs2="a2 a8" order="S" />
            </bondArray>
          </molecule>
        </cml>
        """

        let molecule = try XCTUnwrap(CDKCMLReader.read(text: text).first)
        XCTAssertEqual(molecule.atoms.filter { $0.chirality != .none }.count, 2)

        let chlorineStereo = molecule.bonds.filter { bond in
            let a1 = molecule.atom(id: bond.a1)?.element.uppercased()
            let a2 = molecule.atom(id: bond.a2)?.element.uppercased()
            return a1 == "CL" || a2 == "CL"
        }.map(\.stereo)

        XCTAssertEqual(chlorineStereo.count, 2)
        XCTAssertTrue(chlorineStereo.contains(.up) || chlorineStereo.contains(.upReversed))
        XCTAssertTrue(chlorineStereo.contains(.down) || chlorineStereo.contains(.downReversed))
    }

    func testParsesAromaticBondTypeChildAndMarksAtomsAromatic() throws {
        let text = """
        <molecule id="m1">
          <atomArray atomID="a1 a2" />
          <bondArray>
            <bond atomRefs2="a1 a2" order="2">
              <bondType dictRef="cdk:aromaticBond" />
            </bond>
          </bondArray>
        </molecule>
        """

        let molecule = try XCTUnwrap(CDKCMLReader.read(text: text).first)
        XCTAssertEqual(molecule.bonds[0].order, .double)
        XCTAssertTrue(molecule.atoms[0].aromatic)
        XCTAssertTrue(molecule.atoms[1].aromatic)
    }

    func testParsesPseudoAtomFromDummyTitle() throws {
        let text = """
        <molecule id="m1">
          <atomArray>
            <atom id="a1" elementType="Du" title="Glu55" x2="0.0" y2="0.0" />
          </atomArray>
        </molecule>
        """

        let molecule = try XCTUnwrap(CDKCMLReader.read(text: text).first)
        XCTAssertEqual(molecule.atoms[0].element, "R")
        XCTAssertEqual(molecule.atoms[0].aliasLabel, "Glu55")
        XCTAssertEqual(molecule.atoms[0].symbolToDraw, "Glu55")
    }

    func testRejectsMalformedCML() {
        XCTAssertThrowsError(try CDKCMLReader.read(text: "<cml><molecule><atomArray></cml>"))
    }
}
