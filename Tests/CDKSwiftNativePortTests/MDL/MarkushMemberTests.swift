import XCTest

@testable import CDKSwiftNativePort

final class MarkushMemberTests: XCTestCase {
    private func parse(_ text: String) throws -> Molecule {
        try CDKSmilesParser(flavor: .cdkDefault).parseSmiles(text)
    }
    private func identity(_ molecule: Molecule) throws -> String {
        try CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule).getInchi()
    }

    func testAll103M2SBackbonesRemainVariableIncludingElementLikeAliases() throws {
        struct Record: Decodable {
            let id: Int
            let cxsmiles: String
        }
        let url = try XCTUnwrap(
            Bundle.module.url(
                forResource: "m2s-backbones", withExtension: "json",
                subdirectory: "RGroupEnumerationFixtures"))
        let records = try JSONDecoder().decode([Record].self, from: Data(contentsOf: url))
        XCTAssertEqual(records.count, 103)
        for record in records {
            let molecule = try parse(record.cxsmiles)
            XCTAssertTrue(
                CDKMarkushStructure.requiresSpecifiedStructure(molecule), "M2S \(record.id)")
            let encoded = CDKSmilesGenerator(flavor: .cxsmiles).create(molecule)
            XCTAssertTrue(
                CDKMarkushStructure.requiresSpecifiedStructure(try parse(encoded)),
                "Round trip M2S \(record.id)")
        }
        for text in ["B", "[Y]", "[Bi]", "[Pb]", "[Na+].[Cl-]"] {
            let molecule = try parse(text)
            XCTAssertFalse(CDKMarkushStructure.requiresSpecifiedStructure(molecule), text)
            let roundTrip = try parse(
                CDKSmilesGenerator(flavor: [.cxsmiles, .isomeric]).create(molecule))
            XCTAssertEqual(try identity(molecule), try identity(roundTrip), text)
        }
    }

    func testOriginalRingExampleProducesNineDistinctExpectedMembers() throws {
        let source = try parse("C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|")
        let snapshot = source
        let members = try CDKMarkushEnumerator.configurations(of: source)
        XCTAssertEqual(members.count, 9)
        XCTAssertEqual(
            members.map(\.repeatCounts),
            Array(repeating: [1], count: 3) + Array(repeating: [2], count: 3)
                + Array(repeating: [3], count: 3))
        var expected: Set<String> = []
        for prefix in ["C1", "C1C", "C1CC"] {
            for substituent in ["OC", "Cl", "C#N"] {
                expected.insert(try identity(parse(prefix + "CNCC(\(substituent))C1")))
            }
        }
        XCTAssertEqual(try Set(members.map { try identity($0.molecule) }), expected)
        XCTAssertEqual(source, snapshot)
        XCTAssertTrue(
            members.allSatisfy { !CDKMarkushStructure.requiresSpecifiedStructure($0.molecule) })
    }

    func testExplicitOuterAtomIndicesAndBranchReplication() throws {
        // Middle CH has two repeat ports (F and Cl); its methyl branch repeats with it.
        let source = try parse("FC(C)Cl |LN:1:1.3.0.3|")
        XCTAssertEqual(source.sgroups.first?.atomIDs.count, 2)
        let members = try CDKMarkushEnumerator.configurations(of: source)
        let expected = try Set(
            ["FC(C)Cl", "FC(C)C(C)Cl", "FC(C)C(C)C(C)Cl"].map { try identity(parse($0)) })
        XCTAssertEqual(try Set(members.map { try identity($0.molecule) }), expected)
        let roundTrip = try parse(CDKSmilesGenerator(flavor: .cxsmiles).create(source))
        XCTAssertEqual(
            try Set(
                CDKMarkushEnumerator.configurations(of: roundTrip).map { try identity($0.molecule) }
            ),
            expected)
    }

    func testTwoRepeatUnitsComposeAndLimitsAreExplicit() throws {
        let source = try parse("FCCNCCCl |LN:1:1.2,4:1.2|")
        XCTAssertEqual(try CDKMarkushEnumerator.configurations(of: source).count, 4)
        XCTAssertThrowsError(
            try CDKMarkushEnumerator.configurations(of: source, maximumConfigurations: 3)
        ) {
            XCTAssertEqual($0 as? CDKRGroupEnumerationError, .configurationLimitExceeded(3))
        }
        XCTAssertEqual(try CDKMarkushEnumerator.enumerate(source) { _ in false }, 1)
        XCTAssertThrowsError(try CDKMarkushEnumerator.configurations(of: source, maximumAtoms: 7))
    }

    func testUnsupportedAndMalformedFamiliesNeverBecomeConcrete() throws {
        for text in [
            "C* |$;X$|", "CC |LN:9:1.3|", "FC(C)Cl |LN:1:1.3.0.1|", "FCCl |LN:1:0.3|",
            "F[C@H](C)Cl |LN:1:1.3.0.3|", "FCCl |LN:1:1.999999|", "C* |$;R1$|",
        ] {
            XCTAssertThrowsError(try CDKMarkushEnumerator.configurations(of: parse(text)), text)
        }
    }

    func testCXLogicAndRGFileTwoPreserveRulesAndMembers() throws {
        let url = try XCTUnwrap(
            Bundle.module.url(
                forResource: "rgfile.2", withExtension: "mol",
                subdirectory: "RGroupEnumerationFixtures"))
        let original = try CDKRGFileReader.read(text: String(contentsOf: url, encoding: .utf8))
        let flat = CDKRGroupQueryManipulator.toFlatMolecule(original)
        let encoded = try CDKSmilesGenerator(flavor: [.cxsmiles, .isomeric])
            .createPreservingMarkush(
                flat)
        let decoded = try parse(encoded)
        let rebuilt = try CDKRGroupQueryManipulator.toRGroupQuery(decoded)
        XCTAssertEqual(decoded.rGroupLogicDefinitions.count, 3, encoded)
        XCTAssertEqual(rebuilt.rGroupDefinitions.values.reduce(0) { $0 + $1.rGroups.count }, 5)
        for (number, expected) in original.rGroupDefinitions {
            XCTAssertEqual(rebuilt.rGroupDefinitions[number]?.occurrence, expected.occurrence)
            XCTAssertEqual(rebuilt.rGroupDefinitions[number]?.restH, expected.restH)
            XCTAssertEqual(
                rebuilt.rGroupDefinitions[number]?.requiredRGroupNumber,
                expected.requiredRGroupNumber)
        }
        XCTAssertEqual(try rebuilt.getAllConfigurations().count, 12, encoded)
        XCTAssertEqual(
            try Set(rebuilt.getAllConfigurations().map(identity)),
            try Set(original.getAllConfigurations().map(identity)))
        XCTAssertTrue(decoded.rGroupLogicDefinitions.values.allSatisfy { $0.source == .cxsmiles })
    }

    func testPublishedCXLOGSyntaxAndInvalidConditions() throws {
        let source =
            "[*]C1CCCCC1[*] |$_R1;;;;;;;_R2$,RG:_R1={CCC},_R2={N},LOG={_R1:;;>0._R2:_R1;H;0,1}|"
        let molecule = try parse(source)
        XCTAssertFalse(molecule.name.hasPrefix("|"))
        XCTAssertEqual(molecule.rGroupLogicDefinitions[2]?.requiredRGroupNumber, 1)
        XCTAssertEqual(molecule.rGroupLogicDefinitions[2]?.occurrence, "0,1")
        XCTAssertEqual(try CDKMarkushEnumerator.configurations(of: molecule).count, 2)
        let invalid = try parse("C* |$;R1$,RG:_R1={O},LOG={_R1:;banana;>0}|")
        XCTAssertTrue(invalid.name.hasPrefix("|"))
        XCTAssertThrowsError(try CDKMarkushEnumerator.configurations(of: invalid))
    }

    func testCancellation() async throws {
        let molecule = try parse("FCCl |LN:1:1.3|")
        let task = Task {
            withUnsafeCurrentTask { $0?.cancel() }
            return try CDKMarkushEnumerator.configurations(of: molecule)
        }
        do {
            _ = try await task.value
            XCTFail("Expected cancellation")
        } catch is CancellationError {} catch { XCTFail("\(error)") }
    }
}
