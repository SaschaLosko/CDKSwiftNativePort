import XCTest
@testable import CDKSwiftNativePort

final class FiniteMarkushDefinitionTests: XCTestCase {
    private func parse(_ text: String) throws -> Molecule {
        try CDKSmilesParser(flavor: .cdkDefault, generateCoordinates: false).parseSmiles(text)
    }
    func testReviewedM2SDefinitionsAndUnresolvedSourceTerms() throws {
        struct Record: Decodable {
            let id: Int
            let definition: CDKFiniteMarkushDefinition
            let expectedConfigurations: Int?
            let positiveSMILES: [String]
            let negativeSMILES: [String]
        }
        let url = try XCTUnwrap(
            Bundle.module.url(
                forResource: "m2s-definitions", withExtension: "json", subdirectory: "RGroupEnumerationFixtures"))
        let records = try JSONDecoder().decode([Record].self, from: Data(contentsOf: url))
        XCTAssertEqual(records.count, 103)
        XCTAssertEqual(records.filter { $0.definition.reviewed }.count, 16)
        XCTAssertEqual(records.filter { $0.expectedConfigurations != nil }.count, 10)
        for record in records {
            guard let expected = record.expectedConfigurations else {
                XCTAssertThrowsError(try record.definition.compile(), "M2S \(record.id)")
                continue
            }
            let family = try record.definition.compile()
            let members = try CDKMarkushEnumerator.configurations(of: family)
            XCTAssertEqual(members.count, expected, "M2S \(record.id)")
            for text in record.positiveSMILES {
                let result = CDKMarkushMatcher.match(try parse(text), in: family)
                XCTAssertEqual(result.verdict, .match, "M2S \(record.id): \(result.detail ?? "")")
                // Independent C InChI canonicalization checks the full graph,
                // isotope and stereo layers against manually specified members.
                // Do not replace this with formula/mass signatures or the matcher.
                let expectedInchi = try CDKInChIGenerator(molecule: parse(text)).getInchi()
                let witness = try XCTUnwrap(result.witness, "M2S \(record.id)")
                XCTAssertEqual(
                    try CDKInChIGenerator(molecule: witness.configuration.molecule).getInchi(), expectedInchi)
            }
            for text in record.negativeSMILES {
                XCTAssertEqual(
                    CDKMarkushMatcher.match(try parse(text), in: family).verdict, .noMatch, "M2S \(record.id)")
            }
        }
    }

    func testNestedDefinitionsExclusionsAndReviewGate() throws {
        var definition = CDKFiniteMarkushDefinition(
            sourceIdentifier: "synthetic:nested-1",
            backbone: "C* |$;X$|", reviewed: true,
            variables: [
                .init(
                    label: "X", sourceText: "O-Y or N; N excluded", alternatives: ["O* |$;Y$|", "N"],
                    excludedAlternativeIndices: [1]),
                .init(label: "Y", sourceText: "methyl or ethyl", alternatives: ["C", "CC"]),
            ])
        let family = try definition.compile()
        XCTAssertEqual(try CDKMarkushEnumerator.configurations(of: family).count, 2)
        for text in ["COC", "COCC"] {
            XCTAssertEqual(CDKMarkushMatcher.match(try parse(text), in: family).verdict, .match)
        }
        XCTAssertEqual(CDKMarkushMatcher.match(try parse("CN"), in: family).verdict, .noMatch)
        definition.reviewed = false
        XCTAssertThrowsError(try definition.compile())
        definition.reviewed = true
        definition.unresolvedConditions = ["unless another unspecified substituent is present"]
        XCTAssertThrowsError(try definition.compile())
        definition.unresolvedConditions = []
        definition.variables[1].alternatives = ["C* |$;X$|"]
        XCTAssertThrowsError(try definition.compile())
    }

    func testNeutralImineStereo() throws {
        let family = try parse("C/N=C(/Cl)Br")
        XCTAssertEqual(CDKMarkushMatcher.match(try parse("Cl/C(Br)=N/C"), in: family).verdict, .match)
        XCTAssertEqual(CDKMarkushMatcher.match(try parse("C/N=C(\\Cl)Br"), in: family).verdict, .noMatch)
        XCTAssertEqual(CDKMarkushMatcher.match(try parse("CN=C(Cl)Br"), in: family).reason, .unspecifiedStereochemistry)
    }

    func testExplicitPositionVariationAndAmbiguousCarbonPlaceholder() throws {
        let family = try parse("O*.Fc1ccccc1 |m:1:4.5.6|")
        let members = try CDKMarkushEnumerator.configurations(of: family)
        XCTAssertEqual(members.count, 3)
        for member in members {
            XCTAssertEqual(member.positionVariationTargets.count, 1)
            XCTAssertEqual(CDKMarkushMatcher.match(member.molecule, in: family).verdict, .match)
        }
        XCTAssertEqual(CDKMarkushMatcher.match(try parse("CO"), in: family).verdict, .noMatch)
        XCTAssertThrowsError(try CDKMarkushEnumerator.configurations(of: parse("OC.c1ccccc1 |m:1:2.3|")))
    }

    func testAbsoluteAndRelativeStereoPreserveGroupCorrelations() throws {
        let family = try parse("F[C@H](Cl)[C@H](Br)I |o1:1,3|")
        let members = try CDKMarkushEnumerator.configurations(of: family)
        XCTAssertEqual(members.count, 2)
        XCTAssertEqual(members.map(\.stereoInversions), [[], [1]])
        for text in ["F[C@H](Cl)[C@H](Br)I", "F[C@@H](Cl)[C@@H](Br)I"] {
            XCTAssertEqual(CDKMarkushMatcher.match(try parse(text), in: family).verdict, .match)
        }
        XCTAssertEqual(CDKMarkushMatcher.match(try parse("F[C@H](Cl)[C@@H](Br)I"), in: family).verdict, .noMatch)
        XCTAssertEqual(CDKMarkushMatcher.match(try parse("FC(Cl)C(Br)I"), in: family).verdict, .unresolved)
        XCTAssertEqual(
            CDKMarkushMatcher.match(try parse("F[C@H](Cl)Br |a:1|"), in: try parse("F[C@H](Cl)Br")).verdict, .match)
        let mixture = try parse("F[C@H](Cl)Br |&1:1|")
        XCTAssertThrowsError(try CDKMarkushEnumerator.configurations(of: mixture))
        XCTAssertEqual(CDKMarkushRelations.compare(mixture, mixture, relation: .equivalent).verdict, .unresolved)
        XCTAssertThrowsError(try CDKMarkushEnumerator.configurations(of: family, maximumConfigurations: 1))
    }
}
