// Upstream fixture expectations: CDK 2.13 RGroupQueryReaderTest and
// RGroupQueryManipulatorTest. Copyright 2010 Mark Rijnbeek, 2026 John Mayfield.
// SPDX-License-Identifier: LGPL-2.1-or-later
import XCTest

@testable import CDKSwiftNativePort

final class RGroupEnumerationTests: XCTestCase {
    private func fixture(_ number: Int) throws -> CDKRGroupQuery {
        let url = try XCTUnwrap(
            Bundle.module.url(
                forResource: "rgfile.\(number)", withExtension: "mol",
                subdirectory: "RGroupEnumerationFixtures"))
        return try CDKRGFileReader.read(text: String(contentsOf: url, encoding: .utf8))
    }

    private func parse(_ text: String) throws -> Molecule {
        try CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault).parseSmiles(text)
    }

    private func query(
        _ root: String, _ alternatives: [String], occurrence: String = ">0", restH: Bool = false
    ) throws -> CDKRGroupQuery {
        try CDKRGroupQuery(
            rootStructure: parse(root),
            rGroupDefinitions: [
                1: CDKRGroupList(
                    rGroupNumber: 1, restH: restH, occurrence: occurrence,
                    rGroups: alternatives.map { try CDKRGroup(group: parse($0)) })
            ])
    }

    private func identities(_ molecules: [Molecule]) throws -> Set<String> {
        try Set(molecules.map { try CDKInChIGeneratorFactory.shared.getInChIGenerator($0).getInchi() })
    }

    func testUpstreamRGFileConfigurationCountsAndRestH() throws {
        for (number, expected) in [(1, 4), (2, 12), (3, 8), (4, 2), (5, 17_820), (7, 20)] {
            let input = try fixture(number)
            let original = input
            var count = 0
            var smallestRestH = 0
            try input.enumerateConfigurations { molecule in
                count += 1
                XCTAssertTrue(
                    molecule.atoms.allSatisfy { $0.rGroupLabel == nil && $0.attachmentPoint == nil })
                let ids = Set(molecule.atoms.map(\.id))
                XCTAssertEqual(ids.count, molecule.atomCount)
                XCTAssertTrue(molecule.bonds.allSatisfy { ids.contains($0.a1) && ids.contains($0.a2) })
                if number == 1, molecule.atomCount == 6 {
                    XCTAssertTrue(
                        molecule.atoms.filter { $0.element == "P" && $0.rGroupRestH != nil }.allSatisfy {
                            $0.rGroupRestH == true
                        })
                }
                if number == 2, molecule.atomCount == 13 {
                    smallestRestH += molecule.atoms.filter { $0.rGroupRestH != nil }.count
                    XCTAssertTrue(
                        molecule.atoms.filter { $0.element == "P" && $0.rGroupRestH != nil }.allSatisfy {
                            $0.rGroupRestH == true
                        })
                }
                return true
            }
            XCTAssertEqual(count, expected, "RGfile \(number)")
            if number == 2 { XCTAssertEqual(smallestRestH, 6) }
            XCTAssertEqual(input, original)
        }
        XCTAssertThrowsError(try fixture(6).getAllConfigurations())
    }

    func testFreshJavaCDK213GraphAndRestHSignatures() throws {
        struct Reference: Decodable { let count: Int; let signatureCount: Int; let sha256: String }
        let url = try XCTUnwrap(
            Bundle.module.url(
                forResource: "java-cdk-2.13", withExtension: "json", subdirectory: "RGroupEnumerationFixtures"))
        let references = try JSONDecoder().decode([String: Reference].self, from: Data(contentsOf: url))
        for number in [1, 2, 3, 4, 5, 7] {
            var histogram: [String: Int] = [:]
            try fixture(number).enumerateConfigurations { molecule in
                let atoms = molecule.atoms.map {
                    "\($0.element):\($0.charge):\($0.rGroupRestH.map(String.init) ?? "nil")"
                }.sorted()
                let bonds = molecule.bonds.map { bond in
                    let elements = [molecule.atom(id: bond.a1)!.element, molecule.atom(id: bond.a2)!.element].sorted()
                    return "\(elements[0]):\(bond.order.rawValue):\(elements[1])"
                }.sorted()
                histogram[atoms.joined(separator: ",") + "|" + bonds.joined(separator: ","), default: 0] += 1
                return true
            }
            let serialized = histogram.keys.sorted().map { $0 + "\t" + String(histogram[$0]!) + "\n" }.joined()
            let digest = CDKSHA256.hash(bytes: Array(serialized.utf8)).map { String(format: "%02x", $0) }.joined()
            let expected = try XCTUnwrap(references[String(number)])
            XCTAssertEqual(histogram.values.reduce(0, +), expected.count)
            XCTAssertEqual(histogram.count, expected.signatureCount)
            XCTAssertEqual(digest, expected.sha256, "Java CDK 2.13 RGfile \(number)")
        }
    }

    func testUpstreamExplicitAttachmentAlternativesAndFlatRoundTrip() throws {
        let input = try query(
            "*c1ncccc1 |$R1$|", ["*Cl |$_AP1$|", "*Br |$_AP1$|", "*OC |$_AP1$|", "*C(=O)C |$_AP1$|"])
        let expected = try identities(
            ["Clc1ncccc1", "Brc1ncccc1", "COc1ncccc1", "CC(=O)c1ncccc1"].map(parse))
        XCTAssertEqual(try identities(input.getAllConfigurations()), expected)
        let roundTrip = try CDKRGroupQueryManipulator.toRGroupQuery(
            CDKRGroupQueryManipulator.toFlatMolecule(input))
        XCTAssertEqual(try identities(roundTrip.getAllConfigurations()), expected)
    }

    func testUpstreamTwoPointLigandOrdersAndFlatRoundTrip() throws {
        for (order, expected) in [
            ("0.2", ["NCCOC(=O)O", "NCCNC(=O)O"]), ("2.0", ["NOCCC(=O)O", "NNCCC(=O)O"]),
        ] {
            let input = try query(
                "N*C(=O)O |$;R1$,LO:1:\(order)|", ["*CCO* |$_AP1;;;;_AP2$|", "*CCN* |$_AP1;;;;_AP2$|"])
            let expectedIDs = try identities(expected.map(parse))
            XCTAssertEqual(try identities(input.getAllConfigurations()), expectedIDs)
            let roundTrip = try CDKRGroupQueryManipulator.toRGroupQuery(
                CDKRGroupQueryManipulator.toFlatMolecule(input))
            XCTAssertEqual(try identities(roundTrip.getAllConfigurations()), expectedIDs)
        }
    }

    func testPositionalOccurrencesDependenciesAndSymmetricDuplicates() throws {
        let input = try query("*CC* |$R1;;;R1$|", ["*F |$_AP1$|", "*Cl |$_AP1$|"], occurrence: ">0")
        let configurations = try input.getAllConfigurations()
        XCTAssertEqual(configurations.count, 8)
        XCTAssertEqual(try identities(configurations).count, 5)  // duplicates are intentional positional assignments
        XCTAssertEqual(try input.getAllConfigurations(), configurations)
        var constrained = input
        constrained.rGroupDefinitions[1]?.occurrence = "2"
        XCTAssertEqual(try constrained.getAllConfigurations().count, 4)
        constrained.rGroupDefinitions[1]?.occurrence = "4-6"
        XCTAssertEqual(try constrained.getAllConfigurations().count, 0)
        constrained.rGroupDefinitions[1]?.occurrence = ">\(Int.max)"
        XCTAssertEqual(try constrained.getAllConfigurations().count, 0)
        var dependent = try query("*CC* |$R1;;;R2$|", ["*F |$_AP1$|"], occurrence: "0,1")
        dependent.rGroupDefinitions[1]?.requiredRGroupNumber = 2
        dependent.rGroupDefinitions[2] = try CDKRGroupList(
            rGroupNumber: 2, occurrence: "0,1", rGroups: [CDKRGroup(group: parse("*Cl |$_AP1$|"))])
        XCTAssertEqual(try dependent.getAllConfigurations().count, 3)
        dependent.rGroupDefinitions[2]?.requiredRGroupNumber = 1
        XCTAssertEqual(try dependent.getAllConfigurations().count, 2)
    }

    func testOmittedSiteHydrogensAndRestHTrueFalse() throws {
        for restH in [true, false] {
            let input = try query("[CH2]* |$;R1$|", ["*F |$_AP1$|"], occurrence: "0", restH: restH)
            let product = try XCTUnwrap(input.getAllConfigurations().first)
            XCTAssertEqual(product.atoms[0].explicitHydrogenCount, 3)
            XCTAssertEqual(product.atoms[0].rGroupRestH, restH)
            let decoded = try JSONDecoder().decode(Molecule.self, from: JSONEncoder().encode(product))
            XCTAssertEqual(decoded.atoms[0].rGroupRestH, restH)
        }
    }

    func testScaffoldAndSubstituentStereochemistrySurviveIDRemapping() throws {
        let inputs = [
            ("N[C@H](C)* |$;;;R1$|", "*F |$_AP1$|", "N[C@H](C)F"),
            ("N* |$;R1$|", "*[C@H](F)Cl |$_AP1$|", "N[C@H](F)Cl"),
            ("N* |$;R1$|", "*C/C=C/Cl |$_AP1$|", "NC/C=C/Cl"),
        ]
        for (root, alternative, expected) in inputs {
            let input = try query(root, [alternative])
            XCTAssertEqual(try identities(input.getAllConfigurations()), try identities([parse(expected)]))
            let roundTrip = try CDKRGroupQueryManipulator.toRGroupQuery(CDKRGroupQueryManipulator.toFlatMolecule(input))
            XCTAssertEqual(try identities(roundTrip.getAllConfigurations()), try identities([parse(expected)]))
        }
    }

    func testLimitsStreamingStopAndUnsupportedRepeat() throws {
        let input = try fixture(5)
        XCTAssertThrowsError(try input.getAllConfigurations(maximumConfigurations: 10)) {
            XCTAssertEqual($0 as? CDKRGroupEnumerationError, .configurationLimitExceeded(10))
        }
        var delivered = 0
        XCTAssertEqual(
            try input.enumerateConfigurations(maximumConfigurations: 10) { _ in
                delivered += 1
                return delivered < 3
            }, 3)
        XCTAssertThrowsError(try input.getAllConfigurations(maximumConfigurations: 0))
        let repeated = try CDKRGroupQueryManipulator.toRGroupQuery(
            parse("C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|"))
        XCTAssertThrowsError(try repeated.getAllConfigurations()) {
            guard case .unsupportedFeature = $0 as? CDKRGroupEnumerationError else {
                return XCTFail("Expected unsupported repeat")
            }
        }
    }

    func testMalformedAttachmentsAndMissingDependenciesFail() throws {
        var input = try query("C* |$;R1$|", ["F"])
        XCTAssertEqual(try input.getAllConfigurations().first?.bondCount, 0)
        input.rGroupDefinitions[1]?.rGroups[0].firstAttachmentPointAtomID = 999
        XCTAssertThrowsError(try input.getAllConfigurations())
        input = try query("C* |$;R1$|", ["*F |$_AP1$|"])
        input.rGroupDefinitions[1]?.requiredRGroupNumber = 2
        XCTAssertThrowsError(try input.getAllConfigurations())
        input.rGroupDefinitions[1]?.requiredRGroupNumber = 0
        input.rootStructure.atoms[1].ligandOrderingAtomIDs = [999]
        XCTAssertThrowsError(try input.getAllConfigurations())
    }

    func testCancellation() async throws {
        let input = try fixture(5)
        let task = Task {
            withUnsafeCurrentTask { $0?.cancel() }
            return try input.getAllConfigurations()
        }
        do {
            _ = try await task.value
            XCTFail("Expected cancellation")
        } catch is CancellationError {}
    }
}
