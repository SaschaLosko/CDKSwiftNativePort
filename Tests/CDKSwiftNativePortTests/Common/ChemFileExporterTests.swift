import XCTest
@testable import CDKSwiftNativePort

final class ChemFileExporterTests: XCTestCase {
    private let smilesParser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    func testSupportedExtensionsIncludeWriterFormats() {
        let supported = Set(CDKFileExporter.supportedFileExtensions)
        XCTAssertTrue(supported.contains("mol"))
        XCTAssertTrue(supported.contains("rgf"))
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

    func testFormatsIncludeV3000MolfileExport() {
        XCTAssertTrue(CDKFileExporter.formats.contains(where: { $0.format == .molV3000 && $0.displayName.contains("V3000") }))
    }

    func testFormatsIncludeRGfileExport() {
        XCTAssertTrue(CDKFileExporter.formats.contains(where: { $0.format == .rgfile && $0.fileExtensions.contains("rgf") }))
    }

    func testWritesMolAndRoundTrips() throws {
        let molecule = try referenceMolecule(name: "Aspirin")
        let text = try CDKFileExporter.write(molecule: molecule, as: .mol)
        let parsed = try CDKMDLReader.read(text: text)
        XCTAssertEqual(parsed.atomCount, molecule.atomCount)
        XCTAssertEqual(parsed.bondCount, molecule.bondCount)
    }

    func testWritesV3000MolAndRoundTrips() throws {
        let molecule = try referenceMolecule(name: "Aspirin")
        let text = try CDKFileExporter.write(molecule: molecule, as: .molV3000)
        XCTAssertTrue(text.contains("V3000"))
        let parsed = try CDKMDLReader.read(text: text)
        XCTAssertEqual(parsed.atomCount, molecule.atomCount)
        XCTAssertEqual(parsed.bondCount, molecule.bondCount)
    }

    func testWritesRGFileAndRoundTripsAdvancedLogic() throws {
        let query = try CDKRGFileReader.read(text: RGFileFixtures.simpleQuery)
        let flattened = CDKRGroupQueryManipulator.toFlatMolecule(query)

        let text = try CDKFileExporter.write(molecule: flattened, as: .rgfile)
        let reparsed = try CDKRGFileReader.read(text: text)

        XCTAssertTrue(text.contains("$RGP"))
        XCTAssertTrue(text.contains("M  LOG  1   1   0   1   0,1-3"))
        XCTAssertEqual(reparsed.rGroupDefinitions[1]?.occurrence, "0,1-3")
        XCTAssertEqual(reparsed.rGroupDefinitions[1]?.restH, true)
        XCTAssertEqual(reparsed.rGroupDefinitions[1]?.rGroups.count, 3)
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

    func testWritesSDFAsV3000WhenRequested() throws {
        let molecule = try referenceMolecule(name: "Aspirin")
        var exportOptions = CDKFileExportOptions()
        exportOptions.sdfOptions = CDKSDFWriterOptions(alwaysV3000: true)

        let text = try CDKFileExporter.write(molecule: molecule, as: .sdf, options: exportOptions)

        XCTAssertTrue(text.contains("V3000"))
        XCTAssertFalse(text.contains("V2000"))
    }

    func testWritesSDFDataFieldsAndRoundTripsThem() throws {
        var molecules = try referenceMolecules()
        molecules[0].appendDataFieldValue("001", named: "ID")
        molecules[0].appendDataFieldValue("alpha", named: "Tags")
        molecules[0].appendDataFieldValue("beta", named: "Tags")
        molecules[1].ensureDataField(named: "Empty")

        let text = try CDKFileExporter.write(molecules: molecules, as: .sdf)
        let parsed = try CDKIteratingSDFReader.read(text: text)

        XCTAssertEqual(parsed[0].orderedDataFieldNames, ["ID", "Tags"])
        XCTAssertEqual(parsed[0].dataFieldValues(named: "Tags"), ["alpha", "beta"])
        XCTAssertEqual(parsed[1].orderedDataFieldNames, ["Empty"])
        XCTAssertEqual(parsed[1].dataFieldValues(named: "Empty"), [])
    }

    func testSDFExportOptionsPropagateToWriter() throws {
        var molecule = try referenceMolecule(name: "Tagged")
        molecule.appendDataFieldValue("a", named: "one")
        molecule.appendDataFieldValue("b", named: "two")

        var exportOptions = CDKFileExportOptions()
        exportOptions.sdfOptions = CDKSDFWriterOptions(acceptedDataFieldNames: ["one"],
                                                       programName: "Bioclipse")

        let text = try CDKFileExporter.write(molecule: molecule, as: .sdf, options: exportOptions)

        XCTAssertTrue(text.contains("> <one>"))
        XCTAssertFalse(text.contains("> <two>"))
        XCTAssertTrue(text.contains("  Bioclip"))
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

    func testWritesReactionCMLAndRoundTrips() throws {
        let reaction = CDKReaction(reactants: [try referenceMolecule(name: "Reactant")],
                                   agents: [],
                                   products: [try referenceMolecule(name: "Product")],
                                   id: "cml-rxn-1",
                                   properties: ["Ka": "3"])

        let text = try CDKFileExporter.write(reaction: reaction, as: .cml)
        let parsed = try CDKCMLReactionReader.readReaction(text: text)

        XCTAssertEqual(parsed.id, "cml-rxn-1")
        XCTAssertEqual(parsed.reactantCount, 1)
        XCTAssertEqual(parsed.productCount, 1)
        XCTAssertEqual(parsed.properties["Ka"], "3")
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
