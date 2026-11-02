import XCTest
@testable import CDKSwiftNativePort

final class FiniteChemicalClassTests: XCTestCase {
    private func parse(_ text: String) throws -> Molecule {
        try CDKSmilesParser(generateCoordinates: false).parseSmiles(text)
    }
    private func definition(_ kind: CDKFiniteChemicalClass.Kind = .alkyl, min: Int = 1, max: Int = 6)
        -> CDKFiniteMarkushDefinition
    {
        .init(
            sourceIdentifier: "synthetic:finite-class", backbone: "F* |$;X$|", reviewed: true,
            variables: [
                .init(
                    label: "X", sourceText: "Explicit unsubstituted acyclic saturated class",
                    alternatives: [], chemicalClasses: [.init(kind: kind, minimumCarbons: min, maximumCarbons: max)])
            ])
    }

    func testAllConstitutionalAlkylAttachmentsAndIndependentIdentities() throws {
        // OEIS A000598, n=1...6. Attachment-rooted isomers, ignoring stereoisomers.
        // https://oeis.org/A000598 (reference counts only; no external runtime).
        let expected = [1, 1, 2, 4, 8, 17]
        for n in 1...6 {
            let compiled = try definition(min: n, max: n).compileWithEvidence()
            let members = try CDKMarkushEnumerator.configurations(of: compiled.molecule)
            XCTAssertEqual(members.count, expected[n - 1])
            let identities = try members.map { try CDKInChIGenerator(molecule: $0.molecule).getInchi() }
            XCTAssertEqual(Set(identities).count, members.count, "Different attachment isomers must remain distinct.")
            for member in members {
                XCTAssertEqual(member.molecule.atoms.filter { $0.element == "C" }.count, n)
                XCTAssertEqual(member.molecule.bonds.count, n)
                XCTAssertTrue(member.molecule.bonds.allSatisfy { $0.order == .single })
                let choices = try compiled.evidence.choices(for: member.assignments)
                XCTAssertEqual(choices.count, 1)
                XCTAssertEqual(choices[0].chemicalClass?.minimumCarbons, n)
                XCTAssertNil(choices[0].sourceAlternativeIndex)
                XCTAssertFalse(choices[0].omitted)
            }
        }
        let family = try definition(min: 4, max: 4).compile()
        let expectedIsomers = ["FCCCC", "FC(C)CC", "FCC(C)C", "FC(C)(C)C"]
        let actual = try CDKMarkushEnumerator.configurations(of: family).map {
            try CDKInChIGenerator(molecule: $0.molecule).getInchi()
        }
        XCTAssertEqual(
            Set(actual), try Set(expectedIsomers.map { try CDKInChIGenerator(molecule: parse($0)).getInchi() }))
    }

    func testMembershipRangeBranchingStereoAndExcludedChemistry() throws {
        let family = try definition(min: 1, max: 4).compile()
        for text in ["FC", "FCC", "FC(C)C", "FC(C)(C)C", "F[C@H](C)CC", "F[C@@H](C)CC"] {
            XCTAssertEqual(CDKMarkushMatcher.match(try parse(text), in: family).verdict, .match, text)
        }
        for text in ["F", "FCCCCCC", "FC1CC1", "FC=C", "FCCO", "FC(F)C", "F[13CH3]"] {
            let result = CDKMarkushMatcher.match(try parse(text), in: family)
            XCTAssertEqual(result.verdict, .noMatch, text)
            XCTAssertTrue(result.isExhaustive)
        }
        let limited = CDKMarkushMatcher.match(
            try parse("FCCCCCC"), in: family,
            options: .init(maximumConfigurations: 2))
        XCTAssertEqual(limited.verdict, .unresolved)
        XCTAssertEqual(limited.reason, .configurationLimit)
        let alkoxy = try definition(.alkoxy, min: 3, max: 3).compile()
        for text in ["FOCCC", "FOC(C)C"] {
            XCTAssertEqual(CDKMarkushMatcher.match(try parse(text), in: alkoxy).verdict, .match)
        }
        XCTAssertEqual(CDKMarkushMatcher.match(try parse("FCCOC"), in: alkoxy).verdict, .noMatch)
    }

    func testExplicitAlternativesClassesAndLimits() throws {
        var source = definition(min: 1, max: 2)
        source.variables[0].alternatives = ["O", "N"]
        source.variables[0].excludedAlternativeIndices = [1]
        let compiled = try source.compileWithEvidence()
        XCTAssertEqual(try CDKMarkushEnumerator.configurations(of: compiled.molecule).count, 3)
        XCTAssertEqual(compiled.evidence.variables[0].alternatives[0][0].sourceAlternativeIndex, 0)
        XCTAssertNil(compiled.evidence.variables[0].alternatives[1][0].sourceAlternativeIndex)
        XCTAssertThrowsError(try source.compile(maximumAlternatives: 2))
        XCTAssertThrowsError(try source.compile(maximumAtoms: 2))
        for range in [(0, 1), (3, 2), (1, 7), (1, Int.max)] {
            XCTAssertThrowsError(try definition(min: range.0, max: range.1).compile())
        }
        for backbone in ["C=* |$;X$|", "C*C |$;X;$|"] {
            source.backbone = backbone
            XCTAssertThrowsError(try source.compile())
        }
        source = definition()
        source.unresolvedConditions = ["optionally substituted; substitution rules not supplied"]
        XCTAssertThrowsError(try source.compile())
    }

    func testVersionedRoundTripLegacyReviewGateAndFutureSchemaRejection() throws {
        let source = definition()
        let data = try JSONEncoder().encode(source)
        XCTAssertEqual(try JSONDecoder().decode(CDKFiniteMarkushDefinition.self, from: data), source)
        var json = try XCTUnwrap(JSONSerialization.jsonObject(with: data) as? [String: Any])
        XCTAssertEqual(json["schemaVersion"] as? Int, 2)
        XCTAssertEqual(json["reviewed"] as? Bool, false, "Older readers must not execute unknown class semantics.")
        XCTAssertEqual(json["reviewedDefinition"] as? Bool, true)
        json["schemaVersion"] = 99
        XCTAssertThrowsError(
            try JSONDecoder().decode(
                CDKFiniteMarkushDefinition.self,
                from: JSONSerialization.data(withJSONObject: json)))
        json.removeValue(forKey: "schemaVersion")
        XCTAssertThrowsError(
            try JSONDecoder().decode(
                CDKFiniteMarkushDefinition.self,
                from: JSONSerialization.data(withJSONObject: json)))
        var legacy = source
        legacy.variables[0].chemicalClasses = nil
        legacy.variables[0].alternatives = ["C"]
        let oldData = try JSONEncoder().encode(legacy)
        XCTAssertEqual(try JSONDecoder().decode(CDKFiniteMarkushDefinition.self, from: oldData), legacy)
        XCTAssertNil((try JSONSerialization.jsonObject(with: oldData) as? [String: Any])?["schemaVersion"])
    }
}
