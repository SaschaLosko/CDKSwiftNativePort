import XCTest
@testable import CDKSwiftNativePort

final class InChIGeneratorPortTests: XCTestCase {
    private let smilesParser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    func testGeneratesNativeInchiForEthanol() throws {
        let molecule = try smilesParser.parseSmiles("CCO")
        let generator = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)

        XCTAssertEqual(generator.getStatus(), .success)
        XCTAssertEqual(generator.getMessage(), "")

        let inchi = try generator.getInchi()
        let inchiKey = try generator.getInchiKey()
        XCTAssertTrue(inchi.hasPrefix("InChI=1S/C2H6O"))
        XCTAssertTrue(inchi.contains("/c"))
        XCTAssertTrue(inchi.contains("/h"))
        assertInchiKeyFormat(inchiKey)
    }

    func testGeneratesNativeInchiForAceticAcid() throws {
        let molecule = try smilesParser.parseSmiles("CC(=O)O")
        let generator = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)

        let inchi = try generator.getInchi()
        XCTAssertTrue(inchi.hasPrefix("InChI=1S/C2H4O2"))
        XCTAssertTrue(inchi.contains("/c"))
        XCTAssertTrue(inchi.contains("/h"))
        assertInchiKeyFormat(try generator.getInchiKey())
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
        XCTAssertTrue(inchi.contains("/q+1"))
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

    func testDerivesOfficialInchiKeysFromReferenceStrings() throws {
        let references = [
            ("InChI=1S/CH4/h1H4/i1TD", "VNWKTOKETHGBQD-XIGASBNHSA-N"),
            ("InChI=1S/p+1/i/hD", "GPRLSGONYQIRFK-DYCDLGHISA-N"),
            ("InChI=1S/C4H2N2S2/c5-1-3(7)4(8)2-6/h7-8H/p-2/b4-3+", "DMDOIBWPFWJPQJ-ONEGZZNKSA-L")
        ]

        for (inchi, expectedKey) in references {
            XCTAssertEqual(try CDKInChIKeyCodec.makeKey(from: inchi), expectedKey, "Reference InChIKey mismatch for \(inchi)")
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
}
