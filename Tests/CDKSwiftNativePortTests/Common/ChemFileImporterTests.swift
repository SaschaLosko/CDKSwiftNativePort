import XCTest
@testable import CDKSwiftNativePort

final class ChemFileImporterTests: XCTestCase {
    func testSupportedExtensionsIncludeNewFormats() {
        let supported = Set(CDKFileImporter.supportedFileExtensions)
        XCTAssertTrue(supported.contains("mol2"))
        XCTAssertTrue(supported.contains("pdb"))
        XCTAssertTrue(supported.contains("xyz"))
        XCTAssertTrue(supported.contains("cml"))
        XCTAssertTrue(supported.contains("rxn"))
        XCTAssertTrue(supported.contains("rdf"))
    }

    func testAutoDetectsInChIFromPlainText() throws {
        let text = "InChI=1S/H2O/h1H2 water"
        let molecules = try CDKFileImporter.readMolecules(text: text, fileExtension: "txt")
        XCTAssertEqual(molecules.count, 1)
        XCTAssertEqual(molecules[0].name, "water")
    }

    func testDispatchesMol2ByExtension() throws {
        let text = """
        @<TRIPOS>MOLECULE
        Test
         1 0 0 0 0
        SMALL
        NO_CHARGES
        @<TRIPOS>ATOM
             1 C1 0.0 0.0 0.0 C.3 1 T 0.0
        """

        let molecules = try CDKFileImporter.readMolecules(text: text, fileExtension: "mol2")
        XCTAssertEqual(molecules.count, 1)
        XCTAssertEqual(molecules[0].atomCount, 1)
    }

    func testReadsFromFileURL() throws {
        let temporaryDirectory = URL(fileURLWithPath: NSTemporaryDirectory(), isDirectory: true)
        let fileURL = temporaryDirectory.appendingPathComponent("cdk_file_importer_test.xyz")
        let text = """
        1
        argon
        Ar 0.0 0.0 0.0
        """
        try text.write(to: fileURL, atomically: true, encoding: .utf8)
        defer { try? FileManager.default.removeItem(at: fileURL) }

        let molecules = try CDKFileImporter.readMolecules(from: fileURL)
        XCTAssertEqual(molecules.count, 1)
        XCTAssertEqual(molecules[0].name, "argon")
    }
}
