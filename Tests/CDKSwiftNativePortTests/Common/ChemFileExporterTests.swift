import XCTest
@testable import CDKSwiftNativePort

final class ChemFileExporterTests: XCTestCase {
    private let smilesParser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    func testSupportedExtensionsIncludeWriterFormats() {
        let supported = Set(CDKFileExporter.supportedFileExtensions)
        XCTAssertTrue(supported.contains("mol"))
        XCTAssertTrue(supported.contains("sdf"))
        XCTAssertTrue(supported.contains("smi"))
        XCTAssertTrue(supported.contains("ism"))
        XCTAssertTrue(supported.contains("inchi"))
        XCTAssertTrue(supported.contains("mol2"))
        XCTAssertTrue(supported.contains("pdb"))
        XCTAssertTrue(supported.contains("xyz"))
        XCTAssertTrue(supported.contains("cml"))
        XCTAssertTrue(supported.contains("rxn"))
        XCTAssertTrue(supported.contains("rdf"))
        XCTAssertTrue(supported.contains("svg"))
    }

    func testWritesMolAndRoundTrips() throws {
        let molecule = try referenceMolecule(name: "Aspirin")
        let text = try CDKFileExporter.write(molecule: molecule, as: .mol)
        let parsed = try CDKMDLReader.read(text: text)
        XCTAssertEqual(parsed.atomCount, molecule.atomCount)
        XCTAssertEqual(parsed.bondCount, molecule.bondCount)
    }

    func testWritingMolWithMultipleMoleculesThrows() throws {
        let molecules = try referenceMolecules()
        XCTAssertThrowsError(try CDKFileExporter.write(molecules: molecules, as: .mol))
    }

    func testWritesSDFAndRoundTrips() throws {
        let molecules = try referenceMolecules()
        let text = try CDKFileExporter.write(molecules: molecules, as: .sdf)
        XCTAssertTrue(text.contains("$$$$"))
        let parsed = try CDKIteratingSDFReader.read(text: text)
        XCTAssertEqual(parsed.count, molecules.count)
    }

    func testWritesSmilesAndRoundTrips() throws {
        let molecules = try referenceMolecules()
        let text = try CDKFileExporter.write(molecules: molecules, as: .smiles)
        let parsed = try CDKSMILESReader.read(text: text)
        XCTAssertEqual(parsed.count, molecules.count)
    }

    func testWritesIsomericSmilesAndRoundTrips() throws {
        let molecule = try referenceMolecule(name: "LacticAcid")
        let text = try CDKFileExporter.write(molecule: molecule, as: .isomericSmiles)
        let parsed = try CDKSMILESReader.read(text: text)
        XCTAssertEqual(parsed.count, 1)
        XCTAssertEqual(parsed[0].atomCount, molecule.atomCount)
    }

    func testWritesInChIAndRoundTrips() throws {
        let molecule = try referenceMolecule(name: "Ethanol")
        let text = try CDKFileExporter.write(molecule: molecule, as: .inchi)
        XCTAssertTrue(text.contains("InChI=1S/"))
        let parsed = try CDKInChIReader.read(text: text)
        XCTAssertEqual(parsed.count, 1)
        XCTAssertEqual(parsed[0].atomCount, molecule.atomCount)
    }

    func testWritesXYZAndRoundTrips() throws {
        let molecule = try referenceMolecule(name: "Aspirin")
        let text = try CDKFileExporter.write(molecule: molecule, as: .xyz)
        let parsed = try CDKXYZReader.read(text: text)
        XCTAssertEqual(parsed.count, 1)
        XCTAssertEqual(parsed[0].atomCount, molecule.atomCount)
    }

    func testWritesPDBAndRoundTrips() throws {
        let molecule = try referenceMolecule(name: "Aspirin")
        let text = try CDKFileExporter.write(molecule: molecule, as: .pdb)
        let parsed = try CDKPDBReader.read(text: text)
        XCTAssertEqual(parsed.count, 1)
        XCTAssertEqual(parsed[0].atomCount, molecule.atomCount)
    }

    func testWritesMol2AndRoundTrips() throws {
        let molecule = try referenceMolecule(name: "Aspirin")
        let text = try CDKFileExporter.write(molecule: molecule, as: .mol2)
        let parsed = try CDKMol2Reader.read(text: text)
        XCTAssertEqual(parsed.count, 1)
        XCTAssertEqual(parsed[0].atomCount, molecule.atomCount)
    }

    func testWritesCMLAndRoundTrips() throws {
        let molecule = try referenceMolecule(name: "Aspirin")
        let text = try CDKFileExporter.write(molecule: molecule, as: .cml)
        let parsed = try CDKCMLReader.read(text: text)
        XCTAssertEqual(parsed.count, 1)
        XCTAssertEqual(parsed[0].atomCount, molecule.atomCount)
    }

    func testWritesRXNAndRoundTrips() throws {
        let molecules = try referenceMolecules()
        let text = try CDKFileExporter.write(molecules: molecules, as: .rxn)
        let parsed = try CDKRXNReader.read(text: text)
        XCTAssertEqual(parsed.count, molecules.count)
    }

    func testWritesRDFAndRoundTrips() throws {
        let molecules = try referenceMolecules()
        let text = try CDKFileExporter.write(molecules: molecules, as: .rdf)
        let parsed = try CDKRDFReader.read(text: text)
        XCTAssertEqual(parsed.count, molecules.count)
    }

    func testWritesSVG() throws {
        let molecule = try referenceMolecule(name: "Aspirin")
        var options = CDKFileExportOptions()
        options.svgIncludeBackground = true
        let text = try CDKFileExporter.write(molecule: molecule, as: .svg, options: options)
        XCTAssertTrue(text.contains("<svg"))
        XCTAssertTrue(text.contains("</svg>"))
    }

    func testWriteToURLInfersFormatFromExtension() throws {
        let molecule = try referenceMolecule(name: "Aspirin")
        let temporaryDirectory = URL(fileURLWithPath: NSTemporaryDirectory(), isDirectory: true)
        let fileURL = temporaryDirectory.appendingPathComponent("cdk_exporter_roundtrip_test.sdf")
        defer { try? FileManager.default.removeItem(at: fileURL) }

        try CDKFileExporter.write(molecule: molecule, to: fileURL)
        let written = try String(contentsOf: fileURL, encoding: .utf8)
        XCTAssertTrue(written.contains("$$$$"))
    }

    private func referenceMolecule(name: String) throws -> Molecule {
        var molecule = try smilesParser.parseSmiles("CC(=O)OC1=CC=CC=C1C(=O)O")
        molecule.name = name
        return molecule
    }

    private func referenceMolecules() throws -> [Molecule] {
        var first = try referenceMolecule(name: "Aspirin")
        var second = try smilesParser.parseSmiles("CCO")
        second.name = "Ethanol"
        first = Depiction2DGenerator.generate(for: first)
        second = Depiction2DGenerator.generate(for: second)
        return [first, second]
    }
}
