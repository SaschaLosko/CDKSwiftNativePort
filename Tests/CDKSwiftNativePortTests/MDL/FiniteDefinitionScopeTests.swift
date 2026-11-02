import XCTest
@testable import CDKSwiftNativePort

final class FiniteDefinitionScopeTests: XCTestCase {
    private func parse(_ text: String) throws -> Molecule {
        try CDKSmilesParser(generateCoordinates: false).parseSmiles(text)
    }
    private var source: CDKFiniteMarkushDefinition {
        .init(
            sourceIdentifier: "synthetic:local-definitions", backbone: "F* |$;X$|", reviewed: true,
            variables: [
                .init(
                    label: "X", sourceText: "O-Y or N-Y; each Y has a local definition",
                    alternatives: ["O* |$;Y$|", "N* |$;Y$|"],
                    localDefinitions: [
                        .init(
                            alternativeIndex: 0,
                            variables: [.init(label: "Y", sourceText: "methyl", alternatives: ["C"])]),
                        .init(
                            alternativeIndex: 1,
                            variables: [.init(label: "Y", sourceText: "ethyl", alternatives: ["CC"])]),
                    ]),
                .init(label: "Y", sourceText: "Outer Y is chlorine", alternatives: ["Cl"]),
            ])
    }

    func testLocalScopesDoNotLeakAndEvidenceExplainsNestedChoices() throws {
        let compilation = try source.compileWithEvidence()
        let family = compilation.molecule
        XCTAssertEqual(try CDKMarkushEnumerator.configurations(of: family).count, 2)
        for (index, text) in ["FOC", "FNCC"].enumerated() {
            let result = CDKMarkushMatcher.match(try parse(text), in: family)
            XCTAssertEqual(result.verdict, .match)
            let witness = try XCTUnwrap(result.witness)
            XCTAssertEqual(
                try CDKInChIGenerator(molecule: witness.configuration.molecule).getInchi(),
                try CDKInChIGenerator(molecule: parse(text)).getInchi())
            let choices = try compilation.evidence.choices(for: witness.configuration.assignments)
            XCTAssertEqual(choices.map(\.variableLabel), ["X", "Y"])
            XCTAssertEqual(choices[0].sourceAlternativeIndex, index)
            XCTAssertEqual(choices[1].definitionPath, "/variables/0/alternatives/\(index)/0")
            XCTAssertEqual(choices[1].sourceText, index == 0 ? "methyl" : "ethyl")
        }
        for text in ["FOCC", "FNC", "FOCl", "FNCl"] {
            XCTAssertEqual(CDKMarkushMatcher.match(try parse(text), in: family).verdict, .noMatch)
        }
        let encoded = try JSONEncoder().encode(source)
        XCTAssertEqual(try JSONDecoder().decode(CDKFiniteMarkushDefinition.self, from: encoded), source)
        XCTAssertEqual(
            try JSONDecoder().decode(
                CDKFiniteDefinitionEvidence.self,
                from: JSONEncoder().encode(compilation.evidence)), compilation.evidence)
    }

    func testShadowedSameLabelIsNotACycleButRealCyclesAreRejected() throws {
        var definition = source
        definition.variables = [
            .init(
                label: "X", sourceText: "outer X", alternatives: ["O* |$;X$|"],
                localDefinitions: [
                    .init(
                        alternativeIndex: 0,
                        variables: [.init(label: "X", sourceText: "inner X", alternatives: ["C"])])
                ])
        ]
        XCTAssertEqual(CDKMarkushMatcher.match(try parse("FOC"), in: try definition.compile()).verdict, .match)
        definition.variables[0].localDefinitions![0].variables[0].alternatives = ["C* |$;X$|"]
        XCTAssertThrowsError(try definition.compile())
    }

    func testInheritedDefinitionsKeepTheirDeclarationScope() throws {
        var definition = source
        definition.variables[0].alternatives = ["O* |$;Z$|"]
        definition.variables[0].localDefinitions = [
            .init(
                alternativeIndex: 0,
                variables: [.init(label: "Y", sourceText: "Local Y", alternatives: ["C"])])
        ]
        definition.variables.append(.init(label: "Z", sourceText: "Z uses outer Y", alternatives: ["N* |$;Y$|"]))
        let family = try definition.compile()
        XCTAssertEqual(CDKMarkushMatcher.match(try parse("FONCl"), in: family).verdict, .match)
        XCTAssertEqual(CDKMarkushMatcher.match(try parse("FONC"), in: family).verdict, .noMatch)
    }

    func testLocalConditionsExclusionsAndOmissionArePreserved() throws {
        var definition = source
        definition.variables = [
            .init(
                label: "X", sourceText: "one local branch", alternatives: ["C(*)(*)N* |$;A;A;;B$|"],
                localDefinitions: [
                    .init(
                        alternativeIndex: 0,
                        variables: [
                            .init(
                                label: "A", sourceText: "O or Cl; Cl excluded", alternatives: ["O", "Cl"],
                                excludedAlternativeIndices: [1], occurrence: "0,1", requiresLabel: "B", restH: true),
                            .init(label: "B", sourceText: "methyl", alternatives: ["C"], occurrence: "0,1"),
                        ])
                ])
        ]
        let compilation = try definition.compileWithEvidence()
        let members = try CDKMarkushEnumerator.configurations(of: compilation.molecule)
        XCTAssertEqual(members.count, 4)
        XCTAssertEqual(CDKMarkushMatcher.match(try parse("FC(O)NC"), in: compilation.molecule).verdict, .match)
        XCTAssertEqual(CDKMarkushMatcher.match(try parse("FC(O)N"), in: compilation.molecule).verdict, .noMatch)
        XCTAssertEqual(CDKMarkushMatcher.match(try parse("FC(Cl)NC"), in: compilation.molecule).verdict, .noMatch)
        let omitted = try members.flatMap { try compilation.evidence.choices(for: $0.assignments) }.filter(\.omitted)
        XCTAssertFalse(omitted.isEmpty)
        XCTAssertTrue(omitted.allSatisfy { $0.siteAtomID != nil })
    }

    func testMalformedScopesAndCrossScaffoldDependenciesFailClosed() throws {
        var definition = source
        definition.variables[0].localDefinitions![0].alternativeIndex = 9
        XCTAssertThrowsError(try definition.compile())
        definition = source
        definition.variables[0].localDefinitions!.append(definition.variables[0].localDefinitions![0])
        XCTAssertThrowsError(try definition.compile())
        definition = source
        definition.variables[0].localDefinitions![0].variables[0].requiresLabel = "X"
        XCTAssertThrowsError(try definition.compile())
        definition = source
        definition.variables[0].localDefinitions![0].variables.append(
            definition.variables[0].localDefinitions![0].variables[0])
        XCTAssertThrowsError(try definition.compile())
    }
}
