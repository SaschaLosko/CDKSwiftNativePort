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
        XCTAssertTrue(supported.contains("rsmi"))
        XCTAssertTrue(supported.contains("cxsmiles"))
    }

    func testReadsCxSmilesByCxSmilesExtension() throws {
        let text = "C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|"

        let molecules = try CDKFileImporter.readMolecules(text: text, fileExtension: "cxsmiles")
        let molecule = try XCTUnwrap(molecules.first)

        XCTAssertEqual(molecules.count, 1)
        XCTAssertEqual(molecule.atomCount, 15)
        XCTAssertTrue(molecule.cxState?.rGroupDefinitions.isEmpty == false)
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

    func testReadsReactionByExtension() throws {
        let text = """
        $RXN
        ImporterReaction
          CDKSwiftNativePort

          1  1
        $MOL
        Reactant
          CDK

          1  0  0  0  0  0            999 V2000
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        M  END
        $MOL
        Product
          CDK

          1  0  0  0  0  0            999 V2000
            0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
        M  END
        """

        let reaction = try CDKFileImporter.readReaction(text: text, fileExtension: "rxn")
        XCTAssertEqual(reaction.reactantCount, 1)
        XCTAssertEqual(reaction.productCount, 1)
        XCTAssertEqual(reaction.agentCount, 0)
    }

    func testReadsReactionSmilesByRSMIExtension() throws {
        let text = "CCO>O>CC=O"

        let molecules = try CDKFileImporter.readMolecules(text: text, fileExtension: "rsmi")
        XCTAssertEqual(molecules.count, 3)

        let reaction = try CDKFileImporter.readReaction(text: text, fileExtension: "rsmi")
        XCTAssertEqual(reaction.reactantCount, 1)
        XCTAssertEqual(reaction.agentCount, 1)
        XCTAssertEqual(reaction.productCount, 1)
    }

    func testReadsReactionCMLByCMLExtension() throws {
        let reaction = try CDKFileImporter.readReaction(text: CMLReactionFixtures.reactionWithProperties,
                                                        fileExtension: "cml")

        XCTAssertEqual(reaction.id, "reaction.2")
        XCTAssertEqual(reaction.reactantCount, 1)
        XCTAssertEqual(reaction.productCount, 1)
        XCTAssertEqual(reaction.properties["Ka"], "3")

        let molecules = try CDKFileImporter.readMolecules(text: CMLReactionFixtures.reactionWithProperties,
                                                          fileExtension: "cml")
        XCTAssertEqual(molecules.count, 2)
        XCTAssertEqual(molecules.map(\.externalID), ["react", "product"])
    }

    func testAutoDetectsReactionCMLFromPlainText() throws {
        let reaction = try CDKFileImporter.readReaction(text: CMLReactionFixtures.fragmentReaction,
                                                        fileExtension: "txt")
        XCTAssertEqual(reaction.reactantCount, 1)
        XCTAssertEqual(reaction.productCount, 1)
        XCTAssertEqual(reaction.agentCount, 1)
    }

    func testReadsMappedReactionSmilesByRSMIExtension() throws {
        let text = "[CH3:1][C:2](=[O:3])[CH2:4][S:5][CH2:6][CH:7]=[CH2:8].[CH3:9][CH:10]([NH2:11])[CH:12]([OH:13])[CH:14]=[CH2:15]>>[CH3:1][C:2](=[O:3])[CH2:4][S:5][CH2:6][CH:7]=[CH:14][CH:12]([CH:10]([CH3:9])[NH2:11])[OH:13]"

        let molecules = try CDKFileImporter.readMolecules(text: text, fileExtension: "rsmi")
        XCTAssertEqual(molecules.count, 3)

        let reaction = try CDKFileImporter.readReaction(text: text, fileExtension: "rsmi")
        XCTAssertEqual(reaction.reactantCount, 2)
        XCTAssertEqual(reaction.agentCount, 0)
        XCTAssertEqual(reaction.productCount, 1)

        let reactant0Maps = Set(reaction.reactants[0].atoms.compactMap(\.atomMapNumber))
        let reactant1Maps = Set(reaction.reactants[1].atoms.compactMap(\.atomMapNumber))
        let productMaps = Set(reaction.products[0].atoms.compactMap(\.atomMapNumber))

        XCTAssertEqual(reactant0Maps, Set(1...8))
        XCTAssertEqual(reactant1Maps, Set(9...15))
        XCTAssertTrue(productMaps.isSuperset(of: [1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14]))
    }

    func testReadsCxSmilesBySmilesExtensionWithoutTreatingCxAsName() throws {
        let text = "C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|"

        let molecules = try CDKFileImporter.readMolecules(text: text, fileExtension: "smi")
        let molecule = try XCTUnwrap(molecules.first)

        XCTAssertEqual(molecules.count, 1)
        XCTAssertEqual(molecule.name, "")
        XCTAssertEqual(molecule.sgroups.count, 1)
        XCTAssertEqual(molecule.sgroups.first?.subscriptText, "1-3")
        XCTAssertEqual(molecule.atoms.filter { $0.rGroupMembership == nil && $0.symbolToDraw == "R1" }.count, 1)
    }

    func testAutoDetectsCxSmilesWithoutTreatingCxAsName() throws {
        let text = "C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|"

        let molecules = try CDKFileImporter.readMolecules(text: text, fileExtension: nil)
        let molecule = try XCTUnwrap(molecules.first)

        XCTAssertEqual(molecules.count, 1)
        XCTAssertEqual(molecule.name, "")
        XCTAssertEqual(molecule.sgroups.count, 1)
        XCTAssertEqual(molecule.sgroups.first?.subscriptText, "1-3")
        XCTAssertEqual(molecule.atoms.filter { $0.rGroupMembership == nil && $0.symbolToDraw == "R1" }.count, 1)
    }

    func testReadsReactionSmilesFromFileURLByExtension() throws {
        let temporaryDirectory = URL(fileURLWithPath: NSTemporaryDirectory(), isDirectory: true)
        let fileURL = temporaryDirectory.appendingPathComponent("cdk_file_importer_test.rsmi")
        let text = "CCO>O>CC=O"
        try text.write(to: fileURL, atomically: true, encoding: .utf8)
        defer { try? FileManager.default.removeItem(at: fileURL) }

        let reaction = try CDKFileImporter.readReaction(from: fileURL)
        XCTAssertEqual(reaction.reactantCount, 1)
        XCTAssertEqual(reaction.agentCount, 1)
        XCTAssertEqual(reaction.productCount, 1)

        let molecules = try CDKFileImporter.readMolecules(from: fileURL)
        XCTAssertEqual(molecules.count, 3)
    }

    func testRejectsNonStandardReactionSmilesStoichiometryPrefixes() {
        XCTAssertThrowsError(try CDKFileImporter.readMolecules(text: "{2}C.O>>CO", fileExtension: nil))
        XCTAssertThrowsError(try CDKFileImporter.readReaction(text: "{2}C.O>>CO", fileExtension: nil))
        XCTAssertThrowsError(try CDKFileImporter.readMolecules(text: "{2}C.O>>CO", fileExtension: "rsmi"))
        XCTAssertThrowsError(try CDKFileImporter.readReaction(text: "{2}C.O>>CO", fileExtension: "rsmi"))
    }

    func testRejectsNonStandardReactionSmilesAgentSlashShorthand() {
        XCTAssertThrowsError(try CDKFileImporter.readReaction(text: "C=CC.O=O>Pd/Cu>CC(C)=O", fileExtension: nil))
        XCTAssertThrowsError(try CDKFileImporter.readReaction(text: "C=CC.O=O>Pd/Cu>CC(C)=O", fileExtension: "rsmi"))
    }

    func testRejectsNonStandardReactionSmilesDiatomicFormulaShorthand() {
        XCTAssertThrowsError(try CDKFileImporter.readReaction(text: "H2.O2>>O", fileExtension: nil))
        XCTAssertThrowsError(try CDKFileImporter.readReaction(text: "H2.O2>>O", fileExtension: "rsmi"))
    }

    func testSupportsImplicitStoichiometryByRepeatedReactionComponents() throws {
        let text = "C.C>>O"

        let molecules = try CDKFileImporter.readMolecules(text: text, fileExtension: "rsmi")
        XCTAssertEqual(molecules.count, 3)

        let reaction = try CDKFileImporter.readReaction(text: text, fileExtension: "rsmi")
        XCTAssertEqual(reaction.reactantCount, 2)
        XCTAssertEqual(reaction.agentCount, 0)
        XCTAssertEqual(reaction.productCount, 1)
        XCTAssertTrue(reaction.reactantParticipants.allSatisfy { $0.stoichiometry == nil })
    }
}
