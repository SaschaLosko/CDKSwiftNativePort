import XCTest
@testable import CDKSwiftNativePort

final class SDFWriterPortTests: XCTestCase {

    func testChoosesV2000AndV3000RecordsAutomatically() throws {
        let molecules = [singleCarbon(name: "small"), longCarbonChain(count: 1000, name: "large")]

        let text = try CDKSDFWriter.write(molecules)

        XCTAssertTrue(text.contains("V2000"))
        XCTAssertTrue(text.contains("V3000"))
        XCTAssertEqual(text.components(separatedBy: "$$$$").count - 1, 2)
    }

    func testAlwaysV3000WritesAllRecordsAsV3000() throws {
        let molecules = [singleCarbon(name: "small"), singleCarbon(name: "second")]

        let text = try CDKSDFWriter.write(molecules,
                                          options: CDKSDFWriterOptions(alwaysV3000: true))

        XCTAssertFalse(text.contains("V2000"))
        XCTAssertEqual(text.components(separatedBy: "V3000").count - 1, 2)
    }

    func testSanitizesInvalidHeaderTags() throws {
        var molecule = singleCarbon(name: "Tagged")
        molecule.appendDataFieldValue("URL", named: "http://not-valid.com")

        let text = try CDKSDFWriter.write(molecule)

        XCTAssertTrue(text.contains("> <http://not_valid_com>"))
    }

    func testCanDisableDataFieldOutput() throws {
        var molecule = singleCarbon(name: "Tagged")
        molecule.appendDataFieldValue("a", named: "one")
        molecule.appendDataFieldValue("b", named: "two")

        let text = try CDKSDFWriter.write(molecule,
                                          options: CDKSDFWriterOptions(writeDataFields: false))

        XCTAssertFalse(text.contains("> <one>"))
        XCTAssertFalse(text.contains("> <two>"))
    }

    func testAcceptedDataFieldNamesFilterOutput() throws {
        var molecule = singleCarbon(name: "Tagged")
        molecule.appendDataFieldValue("a", named: "one")
        molecule.appendDataFieldValue("b", named: "two")

        let text = try CDKSDFWriter.write(molecule,
                                          options: CDKSDFWriterOptions(acceptedDataFieldNames: ["one"]))

        XCTAssertTrue(text.contains("> <one>"))
        XCTAssertFalse(text.contains("> <two>"))
    }

    func testTruncatesLongDataFieldsWhenEnabled() throws {
        var molecule = singleCarbon(name: "Tagged")
        molecule.appendDataFieldValue(String(repeating: "x", count: 250), named: "MyLongField")

        let text = try CDKSDFWriter.write(molecule,
                                          options: CDKSDFWriterOptions(truncateLongDataFields: true))

        let longLine = text
            .components(separatedBy: "\n")
            .first {
                !$0.isEmpty &&
                $0.allSatisfy(\.isASCII) &&
                $0.trimmingCharacters(in: .whitespacesAndNewlines).allSatisfy { $0 == "x" }
            }

        XCTAssertEqual(longLine?.count, 200)
        XCTAssertFalse(text.contains(String(repeating: "x", count: 250)))
    }

    func testUsesProgramNameInHeadersForV2000AndV3000Records() throws {
        let molecules = [singleCarbon(name: "small"), longCarbonChain(count: 1000, name: "large")]

        let text = try CDKSDFWriter.write(molecules,
                                          options: CDKSDFWriterOptions(programName: "Bioclipse"))

        let records = text.components(separatedBy: "$$$$\n").filter { !$0.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty }
        XCTAssertEqual(records.count, 2)
        XCTAssertTrue(records.allSatisfy { $0.contains("  Bioclip") })
    }

    func testKeepsSingleRacemicStereoGroupInV2000() throws {
        let molecule = try CDKMDLReader.read(text: v3000SingleRacemicGroup())

        let text = try CDKSDFWriter.write(molecule)

        XCTAssertTrue(text.contains("V2000"))
        XCTAssertFalse(text.contains("BEGIN COLLECTION"))
    }

    func testPromotesRelativeStereoGroupToV3000() throws {
        let molecule = try CDKMDLReader.read(text: v3000RelativeStereoGroup())

        let text = try CDKSDFWriter.write(molecule)

        XCTAssertTrue(text.contains("V3000"))
        XCTAssertTrue(text.contains("MDLV30/STEREL1"))
    }

    private func singleCarbon(name: String) -> Molecule {
        Molecule(name: name,
                 atoms: [Atom(id: 1, element: "C", position: .zero)],
                 bonds: [])
    }

    private func longCarbonChain(count: Int, name: String) -> Molecule {
        var atoms: [Atom] = []
        atoms.reserveCapacity(count)
        var bonds: [Bond] = []
        bonds.reserveCapacity(max(0, count - 1))

        for index in 0..<count {
            atoms.append(Atom(id: index + 1,
                              element: "C",
                              position: CGPoint(x: Double(index) * 1.4, y: 0)))
            if index > 0 {
                bonds.append(Bond(id: index,
                                  a1: index,
                                  a2: index + 1,
                                  order: .single))
            }
        }

        return Molecule(name: name, atoms: atoms, bonds: bonds)
    }

    private func v3000SingleRacemicGroup() -> String {
        """

          Mrv1810 02052113162D

          0  0  0     0  0            999 V3000
        M  V30 BEGIN CTAB
        M  V30 COUNTS 7 7 0 0 0
        M  V30 BEGIN ATOM
        M  V30 1 C -2.1407 12.3148 0 0 CFG=2
        M  V30 2 C -3.4743 11.5447 0 0
        M  V30 3 C -3.4743 10.0047 0 0
        M  V30 4 C -2.1407 9.2347 0 0
        M  V30 5 C -0.807 10.0047 0 0
        M  V30 6 N -0.807 11.5447 0 0
        M  V30 7 O -2.1407 13.8548 0 0
        M  V30 END ATOM
        M  V30 BEGIN BOND
        M  V30 1 1 1 2
        M  V30 2 1 2 3
        M  V30 3 1 3 4
        M  V30 4 1 4 5
        M  V30 5 1 5 6
        M  V30 6 1 1 6
        M  V30 7 1 1 7 CFG=1
        M  V30 END BOND
        M  V30 BEGIN COLLECTION
        M  V30 MDLV30/STERAC1 ATOMS=(1 1)
        M  V30 END COLLECTION
        M  V30 END CTAB
        M  END
        """
    }

    private func v3000RelativeStereoGroup() -> String {
        """

          Mrv1810 02052113162D

          0  0  0     0  0            999 V3000
        M  V30 BEGIN CTAB
        M  V30 COUNTS 7 7 0 0 0
        M  V30 BEGIN ATOM
        M  V30 1 C -2.1407 12.3148 0 0 CFG=2
        M  V30 2 C -3.4743 11.5447 0 0
        M  V30 3 C -3.4743 10.0047 0 0
        M  V30 4 C -2.1407 9.2347 0 0
        M  V30 5 C -0.807 10.0047 0 0
        M  V30 6 N -0.807 11.5447 0 0
        M  V30 7 O -2.1407 13.8548 0 0
        M  V30 END ATOM
        M  V30 BEGIN BOND
        M  V30 1 1 1 2
        M  V30 2 1 2 3
        M  V30 3 1 3 4
        M  V30 4 1 4 5
        M  V30 5 1 5 6
        M  V30 6 1 1 6
        M  V30 7 1 1 7 CFG=1
        M  V30 END BOND
        M  V30 BEGIN COLLECTION
        M  V30 MDLV30/STEREL5 ATOMS=(1 1)
        M  V30 END COLLECTION
        M  V30 END CTAB
        M  END
        """
    }
}
