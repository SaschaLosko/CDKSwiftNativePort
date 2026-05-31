import XCTest
@testable import CDKSwiftNativePort

private struct InChIOfficialReferenceFixtures: Decodable {
    struct Case: Decodable {
        let id: String
        let expectedInChI: String
        let expectedInChIKey: String
        let molfile: String
    }

    let cases: [Case]
}

final class InChIGeneratorPortTests: XCTestCase {
    private let smilesParser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    func testGeneratesNativeInchiForEthanol() throws {
        let molecule = try smilesParser.parseSmiles("CCO")
        let generator = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)

        XCTAssertEqual(generator.getStatus(), .success)
        XCTAssertEqual(generator.getMessage(), "")

        let inchi = try generator.getInchi()
        let inchiKey = try generator.getInchiKey()
        XCTAssertEqual(inchi, "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
        XCTAssertEqual(inchiKey, "LFQSCWFLJHTTHZ-UHFFFAOYSA-N")
    }

    func testGeneratesNativeInchiForAceticAcid() throws {
        let molecule = try smilesParser.parseSmiles("CC(=O)O")
        let generator = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)

        let inchi = try generator.getInchi()
        XCTAssertEqual(inchi, "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)")
        XCTAssertEqual(try generator.getInchiKey(), "QTBSBXVTEAMEQO-UHFFFAOYSA-N")
    }

    func testGeneratedInchiCanBeParsedBack() throws {
        let molecule = try smilesParser.parseSmiles("CC(=O)O")
        let generator = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)

        let inchi = try generator.getInchi()
        let parser = CDKInChIGeneratorFactory.shared.getInChIToStructure(inchi)
        let reparsed = try parser.getAtomContainer()

        XCTAssertEqual(reparsed.atomCount, molecule.atomCount)
        XCTAssertEqual(reparsed.bondCount, molecule.bondCount)
    }

    func testGeneratesDeterministicInchiKey() throws {
        let molecule = try smilesParser.parseSmiles("c1ccccc1")
        let generatorA = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)
        let generatorB = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)

        let keyA = try generatorA.getInchiKey()
        let keyB = try generatorB.getInchiKey()

        XCTAssertEqual(keyA, keyB)
        assertInchiKeyFormat(keyA)
    }

    func testGeneratesChargeAndIsotopeLayers() throws {
        let molecule = try smilesParser.parseSmiles("[13CH3][NH3+]")
        let generator = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)

        let inchi = try generator.getInchi()
        XCTAssertTrue(inchi.contains("/p+1"))
        XCTAssertTrue(inchi.contains("/i"))
    }

    func testHandlesMdlChargeStyledElementsWithoutFailure() throws {
        var molecule = try smilesParser.parseSmiles("O=[N+]([O-])C")
        for index in molecule.atoms.indices {
            if molecule.atoms[index].element.uppercased() == "N", molecule.atoms[index].charge == 1 {
                molecule.atoms[index].element = "N+"
                molecule.atoms[index].charge = 0
            } else if molecule.atoms[index].element.uppercased() == "O", molecule.atoms[index].charge == -1 {
                molecule.atoms[index].element = "O-"
                molecule.atoms[index].charge = 0
            }
        }
        let generator = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)
        XCTAssertEqual(generator.getStatus(), .success)
        let inchi = try generator.getInchi()
        XCTAssertTrue(inchi.hasPrefix("InChI=1S/"))
        XCTAssertFalse(inchi.contains("N+"))
        XCTAssertFalse(inchi.contains("O-"))
    }

    func testGeneratesInchiFromV2000WithSignedElementTokens() throws {
        let mol = """
Nitro
CDKSwiftNativePort

  3  2  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 N+  0  0  0  0  0  0  0  0  0  0  0  0
    1.2000    0.0000    0.0000 O-  0  0  0  0  0  0  0  0  0  0  0  0
   -1.2000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  2  0  0  0  0
M  END
"""
        let molecule = try CDKMDLV2000Reader.read(text: mol)
        let generator = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)
        XCTAssertEqual(generator.getStatus(), .success)

        let inchi = try generator.getInchi()
        XCTAssertTrue(inchi.hasPrefix("InChI=1S/"))
        XCTAssertTrue(inchi.contains("N"))
        XCTAssertTrue(inchi.contains("O"))
        XCTAssertFalse(inchi.contains("N+"))
        XCTAssertFalse(inchi.contains("O-"))
    }

    func testGeneratesHydrogenOnlyReferenceInchi() throws {
        let mol = """
Deutron
CDKSwiftNativePort

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 H   1  3  0  0  0  0  0  0  0  0  0  0
M  CHG  1   1   1
M  ISO  1   1   2
M  END
"""
        let molecule = try CDKMDLV2000Reader.read(text: mol)
        let generator = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)

        XCTAssertEqual(generator.getStatus(), .success)
        XCTAssertEqual(try generator.getInchi(), "InChI=1S/p+1/i/hD")
        XCTAssertEqual(try generator.getInchiKey(), "GPRLSGONYQIRFK-DYCDLGHISA-N")
    }

    func testGeneratesCompactIsotopicHydrogenTokens() throws {
        let methane = """
Methane Isotopes
CDKSwiftNativePort

  4  3  0  0  0  0  0  0  0  0999 V2000
   10.1667   -6.7458    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.0833   -6.2125    0.0000 H   1  0  0  0  0  0  0  0  0  0  0  0
   10.6958   -7.6625    0.0000 H   2  0  0  0  0  0  0  0  0  0  0  0
    9.1042   -6.7417    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  3  1  0  0  0  0
  1  2  1  0  0  0  0
  1  4  1  1  0  0  0
M  ISO  2   2   2   3   3
M  END
"""
        let methaneMolecule = try CDKMDLV2000Reader.read(text: methane)
        XCTAssertEqual(try CDKInChIGeneratorFactory.shared.getInChIGenerator(methaneMolecule).getInchi(),
                       "InChI=1S/CH4/h1H4/i1TD")

        let water = """
Water Isotopes
CDKSwiftNativePort

  3  2  0  0  0  0  0  0  0  0999 V2000
   18.3720  -19.7020    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   20.2318  -19.7020    0.0000 H   1  0  0  0  0  0  0  0  0  0  0  0
   17.3013  -21.3025    0.0000 H   2  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
  3  1  1  0  0  0  0
M  ISO  2   2   2   3   3
M  END
"""
        let waterMolecule = try CDKMDLV2000Reader.read(text: water)
        XCTAssertEqual(try CDKInChIGeneratorFactory.shared.getInChIGenerator(waterMolecule).getInchi(),
                       "InChI=1S/H2O/h1H2/i/hTD")
    }

    func testGeneratesReferenceInchiForSimpleTreeCases() throws {
        let cases = [
            ("CO", "InChI=1S/CH4O/c1-2/h2H,1H3"),
            ("CC", "InChI=1S/C2H6/c1-2/h1-2H3"),
            ("C=C", "InChI=1S/C2H4/c1-2/h1-2H2"),
            ("C[P](C)C", "InChI=1S/C3H9P/c1-4(2)3/h1-3H3")
        ]

        for (smiles, expectedInChI) in cases {
            let molecule = try smilesParser.parseSmiles(smiles)
            let generator = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)
            XCTAssertEqual(generator.getStatus(), .success)
            XCTAssertEqual(try generator.getInchi(), expectedInChI, "Unexpected InChI for \(smiles)")
        }
    }

    func testGeneratesReferenceInchiForTetrafluoroborate() throws {
        let mol = """
BF4
CDKSwiftNativePort

  5  4  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 B   0  0
    1.3000    0.0000    0.0000 F   0  0
   -1.3000    0.0000    0.0000 F   0  0
    0.0000    1.3000    0.0000 F   0  0
    0.0000   -1.3000    0.0000 F   0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
M  CHG  1   1  -1
M  END
"""
        let molecule = try CDKMDLV2000Reader.read(text: mol)
        let generator = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)

        XCTAssertEqual(generator.getStatus(), .success)
        XCTAssertEqual(try generator.getInchi(), "InChI=1S/BF4/c2-1(3,4)5/q-1")
    }

    func testOfficialLibraryDerivesInchiKeysFromReferenceStrings() throws {
        let references = [
            ("InChI=1S/CH4/h1H4/i1TD", "VNWKTOKETHGBQD-XIGASBNHSA-N"),
            ("InChI=1S/p+1/i/hD", "GPRLSGONYQIRFK-DYCDLGHISA-N"),
            ("InChI=1S/C4H2N2S2/c5-1-3(7)4(8)2-6/h7-8H/p-2/b4-3+", "DMDOIBWPFWJPQJ-ONEGZZNKSA-L")
        ]

        for (inchi, expectedKey) in references {
            XCTAssertEqual(try CDKInChINativeGenerator.inchiKey(for: inchi), expectedKey, "Reference InChIKey mismatch for \(inchi)")
        }
    }

    func testGeneratesReferenceElementalIsotopePairs() throws {
        let dy = """
Dy Pair
CDKSwiftNativePort

  2  1  0  0  0  0  0  0  0  0999 V2000
   -1.2875    0.0000    0.0000 Dy -1  0  0  0  0  0  0  0  0  0  0  0
   -0.2250    0.0000    0.0000 Dy  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  ISO  2   1 162   2 163
M  END
"""
        let ag = """
Ag Pair
CDKSwiftNativePort

  2  1  0  0  0  0  0  0  0  0999 V2000
   -1.2875    0.0000    0.0000 Ag -1  0  0  0  0  0  0  0  0  0  0  0
   -0.2250    0.0000    0.0000 Ag  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  ISO  2   1 107   2 108
M  END
"""
        let cases = [
            (dy, "InChI=1S/2Dy/i1+0;1-1", "WSFQKYAVYHDRER-QCIKTKHTSA-N"),
            (ag, "InChI=1S/2Ag/i1+0;1-1", "OGFYIDCVDSATDC-QCIKTKHTSA-N")
        ]

        for (mol, expectedInChI, expectedKey) in cases {
            let molecule = try CDKMDLV2000Reader.read(text: mol)
            let generator = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)

            XCTAssertEqual(generator.getStatus(), .success)
            XCTAssertEqual(try generator.getInchi(), expectedInChI)
            XCTAssertEqual(try generator.getInchiKey(), expectedKey)
        }
    }

    func testGeneratesReferenceInchiForSelectedOfficialMulticomponentCases() throws {
        let ids = [
            "mcule:4028583412",
            "mcule:4699275687",
            "mcule:2940963212",
            "mcule:8769289721",
            "mcule:9786532355"
        ]

        let fixtures = try loadOfficialReferenceFixturesByID()
        for id in ids {
            guard let fixture = fixtures[id] else {
                XCTFail("Missing official reference fixture for \(id)")
                continue
            }

            let molecule = try CDKMDLReader.read(text: fixture.molfile)
            let generator = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)

            XCTAssertEqual(generator.getStatus(), .success)
            XCTAssertEqual(try generator.getInchi(), fixture.expectedInChI, "Unexpected multicomponent InChI for \(id)")
            XCTAssertEqual(try generator.getInchiKey(), fixture.expectedInChIKey, "Unexpected multicomponent InChIKey for \(id)")
        }
    }

    func testGeneratesReferenceInchiForSelectedOfficialAcyclicTreeCases() throws {
        let ids = [
            "mcule:1292201352",
            "mcule:1803511350",
            "mcule:1837794750",
            "mcule:4053513648",
            "mcule:4525381422",
            "mcule:5227399299",
            "mcule:5269746171",
            "mcule:5369504805",
            "mcule:6601882022",
            "mcule:9646570107",
            "mcule:9872974813"
        ]

        let fixtures = try loadOfficialReferenceFixturesByID()
        for id in ids {
            guard let fixture = fixtures[id] else {
                XCTFail("Missing official reference fixture for \(id)")
                continue
            }

            let molecule = try CDKMDLReader.read(text: fixture.molfile)
            let generator = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)

            XCTAssertEqual(generator.getStatus(), .success)
            XCTAssertEqual(try generator.getInchi(), fixture.expectedInChI, "Unexpected acyclic-tree InChI for \(id)")
            XCTAssertEqual(try generator.getInchiKey(), fixture.expectedInChIKey, "Unexpected acyclic-tree InChIKey for \(id)")
        }
    }

    private func assertInchiKeyFormat(_ key: String,
                                      file: StaticString = #filePath,
                                      line: UInt = #line) {
        XCTAssertEqual(key.count, 27, file: file, line: line)
        let parts = key.split(separator: "-")
        XCTAssertEqual(parts.count, 3, file: file, line: line)
        XCTAssertEqual(parts[0].count, 14, file: file, line: line)
        XCTAssertEqual(parts[1].count, 10, file: file, line: line)
        XCTAssertEqual(parts[2].count, 1, file: file, line: line)
        XCTAssertTrue(key.allSatisfy { $0 == "-" || ("A"..."Z").contains(String($0)) },
                      file: file,
                      line: line)
    }

    private func loadOfficialReferenceFixturesByID() throws -> [String: InChIOfficialReferenceFixtures.Case] {
        let url = URL(fileURLWithPath: #filePath)
            .deletingLastPathComponent()
            .appendingPathComponent("OfficialReference")
            .appendingPathComponent("official_reference_cases.json")
        let data = try Data(contentsOf: url)
        let fixtures = try JSONDecoder().decode(InChIOfficialReferenceFixtures.self, from: data)
        return Dictionary(uniqueKeysWithValues: fixtures.cases.map { ($0.id, $0) })
    }
}
