import XCTest
@testable import CDKSwiftNativePort

final class MarkushBoundaryTests: XCTestCase {
    private func parse(_ text: String) throws -> Molecule {
        try CDKSmilesParser(flavor: .cdkDefault).parseSmiles(text)
    }

    func testExpandedAtomLimitIsUnresolvedAndStopsBeforeDelivery() throws {
        let family = try parse("C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|")
        var delivered = 0
        XCTAssertThrowsError(
            try CDKMarkushEnumerator.enumerate(family, maximumAtoms: 2) { _ in
                delivered += 1; return true
            }
        ) { XCTAssertEqual($0 as? CDKRGroupEnumerationError, .atomLimitExceeded(2)) }
        XCTAssertEqual(delivered, 0)
        let query = try CDKRGroupQueryManipulator.toRGroupQuery(parse("*C* |$R1;;R1$,RG:_R1={CCCC}|"))
        XCTAssertThrowsError(
            try query.enumerateConfigurations(maximumAtoms: 4) { _ in
                XCTFail("An oversized configuration was delivered"); return true
            }
        ) { XCTAssertEqual($0 as? CDKRGroupEnumerationError, .atomLimitExceeded(4)) }
    }

    func testNestedCompilationBoundsRetainedAtoms() throws {
        let definition = CDKFiniteMarkushDefinition(
            sourceIdentifier: "synthetic:bounded", backbone: "C* |$;X$|",
            reviewed: true,
            variables: [
                .init(label: "X", sourceText: "O-Y", alternatives: ["O* |$;Y$|"]),
                .init(label: "Y", sourceText: "ethyl or propyl", alternatives: ["CC", "CCC"]),
            ])
        XCTAssertThrowsError(try definition.compile(maximumAtoms: 8))
        XCTAssertThrowsError(try definition.compile(maximumAlternatives: 1))
        XCTAssertEqual(try CDKMarkushEnumerator.configurations(of: definition.compile()).count, 2)
    }

    func testInconsistentMetadataAndDuplicateIDsCannotBeErasedByExpansion() throws {
        var position = try parse("O*.c1ccccc1 |m:1:2.3|")
        position.cxState?.positionalVariations[1] = [999]
        XCTAssertThrowsError(try CDKMarkushEnumerator.configurations(of: position))
        var stereo = try parse("F[C@H](Cl)Br |o1:1|")
        stereo.cxState?.stereoGroups[1] = CDKCxSmilesParser.encodeStereoGroup(kind: "and", number: 2)
        XCTAssertThrowsError(try CDKMarkushEnumerator.configurations(of: stereo))
        stereo.atoms.append(stereo.atoms[0])
        XCTAssertEqual(CDKMarkushRelations.compare(stereo, try parse("C"), relation: .overlaps).verdict, .unresolved)
    }

    func testOmissionProvenanceAndRepeatOrigins() throws {
        var query = try CDKRGroupQueryManipulator.toRGroupQuery(parse("C* |$;R1$,RG:_R1={O}|"))
        query.rGroupDefinitions[1]?.occurrence = "0,1"
        var configurations: [CDKRGroupConfiguration] = []
        try query.enumerateConfigurationsWithAssignments {
            configurations.append($0); return true
        }
        XCTAssertEqual(configurations.count, 2)
        XCTAssertNil(configurations[0].assignments.first?.alternativeIndex)
        XCTAssertTrue(configurations[0].assignments[0].memberAtomIDByFragmentAtomID.isEmpty)
        XCTAssertEqual(configurations[1].assignments.first?.alternativeIndex, 0)
        let repeated = try CDKMarkushEnumerator.configurations(of: parse("FCCl |LN:1:1.3|"))
        for configuration in repeated {
            XCTAssertEqual(Set(configuration.sourceByMemberAtomID.keys), Set(configuration.molecule.atoms.map(\.id)))
            XCTAssertEqual(
                configuration.sourceByMemberAtomID.values.filter { $0.atomID == 2 }.count, configuration.repeatCounts[0]
            )
        }
    }

    func testVariableStereoIsClassifiedAsFamilyAndAbsoluteStereoIsConcrete() throws {
        for suffix in ["o1:1", "&1:1", "r"] {
            XCTAssertTrue(CDKMarkushStructure.requiresSpecifiedStructure(try parse("F[C@H](Cl)Br |\(suffix)|")))
        }
        XCTAssertFalse(CDKMarkushStructure.requiresSpecifiedStructure(try parse("F[C@H](Cl)Br |a:1|")))
    }

    func testRelationCancellationIsUnresolved() async throws {
        let graph = try parse("C* |$;R1$,RG:_R1={O},{N}|")
        let task = Task {
            withUnsafeCurrentTask { $0?.cancel() }
            return CDKMarkushRelations.compare(graph, graph, relation: .equivalent)
        }
        let result = await task.value
        XCTAssertEqual(result.verdict, .unresolved)
        XCTAssertEqual(result.reason, .cancelled)
        XCTAssertFalse(result.isExhaustive)
    }
}
