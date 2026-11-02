import XCTest

@testable import CDKSwiftNativePort

final class MarkushMatchingTests: XCTestCase {
    private func parse(_ text: String) throws -> Molecule {
        try CDKSmilesParser(flavor: .cdkDefault).parseSmiles(text)
    }
    private func match(
        _ compound: String, _ family: String,
        options: CDKMarkushMatchOptions = .init()
    ) throws -> CDKMarkushMatchResult {
        CDKMarkushMatcher.match(try parse(compound), in: try parse(family), options: options)
    }
    private let simpleFamily = "C* |$;R1$,RG:_R1={O},{N}|"

    func testFiniteMembershipHasConcreteWitnessesAndRequiresWholeGraph() throws {
        let source = try parse(simpleFamily)
        let snapshot = source
        for text in ["CO", "CN", "OC"] {
            let compound = try parse(text)
            let result = CDKMarkushMatcher.match(compound, in: source)
            XCTAssertEqual(result.verdict, .match, "\(text): \(result.detail ?? "")")
            let witness = try XCTUnwrap(result.witness)
            XCTAssertEqual(witness.compoundAtomIDByMemberAtomID.count, 2)
            for (memberID, compoundID) in witness.compoundAtomIDByMemberAtomID {
                XCTAssertEqual(
                    witness.configuration.molecule.atom(id: memberID)?.element,
                    compound.atom(id: compoundID)?.element)
            }
            XCTAssertFalse(result.isExhaustive)
            XCTAssertEqual(result.variableNumberBySourceLabel, ["R1": 1])
        }
        for text in ["C", "CC", "COO", "CO.[Na+]", "[13CH3]O", "C[O-]"] {
            let result = try match(text, simpleFamily)
            XCTAssertEqual(result.verdict, .noMatch, "\(text): \(result.detail ?? "")")
            XCTAssertTrue(result.isExhaustive)
            XCTAssertEqual(result.configurationsChecked, 2)
            XCTAssertNil(result.witness)
        }
        XCTAssertEqual(source, snapshot)
    }

    func testOrdinaryHydrogenVerticesAndAromaticRepresentations() throws {
        for (compound, family) in [
            ("[H]OC", "CO"), ("[H]C([H])([H])O[H]", "CO"),
            ("c1ccccc1O", "OC1=CC=CC=C1"), ("[2H]O", "DO"),
        ] {
            let result = try match(compound, family)
            XCTAssertEqual(result.verdict, .match, "\(compound): \(result.detail ?? "")")
        }
        XCTAssertEqual(try match("[2H]OC", "CO").verdict, .noMatch)
        XCTAssertEqual(try match("C=O", "CO").verdict, .noMatch)
    }

    func testOriginalRepeatFamilyMatchesEveryMemberAndRejectsWrongRepeat() throws {
        let family = try parse("C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|")
        let members = try CDKMarkushEnumerator.configurations(of: family)
        XCTAssertEqual(members.count, 9)
        for member in members {
            let result = CDKMarkushMatcher.match(member.molecule, in: family)
            XCTAssertEqual(result.verdict, .match, result.detail ?? "")
            XCTAssertEqual(result.witness?.configuration.repeatCounts, member.repeatCounts)
        }
        let outside = try parse("C1CCCCNCC(Cl)C1")
        XCTAssertEqual(CDKMarkushMatcher.match(outside, in: family).verdict, .noMatch)
    }

    func testLimitsDoNotMasqueradeAsNegativesAndPositiveMayStopEarly() throws {
        let limited = CDKMarkushMatchOptions(maximumConfigurations: 1)
        XCTAssertEqual(try match("CO", simpleFamily, options: limited).verdict, .match)
        let result = try match("CN", simpleFamily, options: limited)
        XCTAssertEqual(result.verdict, .unresolved)
        XCTAssertEqual(result.reason, .configurationLimit)
        XCTAssertEqual(result.configurationsChecked, 1)
        XCTAssertFalse(result.isExhaustive)
        XCTAssertEqual(try match("CC", simpleFamily, options: .init(maximumSearchSteps: 1)).reason, .searchLimit)
        XCTAssertEqual(try match("CO", simpleFamily, options: .init(maximumAtoms: 1)).reason, .atomLimit)
        XCTAssertEqual(
            try match("CO", simpleFamily, options: .init(timeoutSeconds: .leastNonzeroMagnitude)).reason, .timedOut)
        XCTAssertEqual(try match("CO", simpleFamily, options: .init(maximumConfigurations: 0)).reason, .invalidInput)
    }

    func testNumericAndNamedVariableLabelsHaveExplicitStableMappings() throws {
        for label in ["R", "R42", "R100", "X", "Y", "B", "Rd", "A13"] {
            let family = try parse("C* |$;\(label)$,RG:_\(label)={O},{N}|")
            let snapshot = family
            let result = CDKMarkushMatcher.match(try parse("CO"), in: family)
            XCTAssertEqual(result.verdict, .match, "\(label): \(result.detail ?? "")")
            let expected = label == "R" ? 0 : (label == "R42" ? 42 : (label == "R100" ? 100 : 1))
            XCTAssertEqual(result.variableNumberBySourceLabel, [label: expected])
            XCTAssertEqual(family, snapshot)
            XCTAssertEqual(try CDKMarkushEnumerator.configurations(of: family).count, 2)
        }
        let numeric = try parse("C* |$;R42$,RG:_R42={O}|")
        let query = try CDKRGroupQueryManipulator.toRGroupQuery(numeric)
        XCTAssertNotNil(query.rGroupDefinitions[42])
        XCTAssertFalse(CDKRGFileWriter.canWrite(query))
        XCTAssertThrowsError(try CDKRGFileWriter.write(query))
        let encoded = try CDKSmilesGenerator(flavor: .cxsmiles).createPreservingMarkush(numeric)
        XCTAssertEqual(try match("CO", encoded).verdict, .match)
        XCTAssertEqual(try match("B", "B").verdict, .match)
        XCTAssertEqual(try match("[Y]", "[Y]").verdict, .match)
    }

    func testUndefinedLabelsAndNamedVariableCollisionsRemainUnresolved() throws {
        for label in ["X", "B", "Y", "R42", "R100"] {
            let result = try match("CO", "C* |$;\(label)$|")
            XCTAssertEqual(result.verdict, .unresolved)
            XCTAssertNil(result.witness)
        }
        let family = try parse("*C* |$X;;R1$,RG:_X={O},_R1={N}|")
        let result = CDKMarkushMatcher.match(try parse("OCN"), in: family)
        XCTAssertEqual(result.verdict, .match, result.detail ?? "")
        XCTAssertEqual(result.variableNumberBySourceLabel, ["X": 2, "R1": 1])
    }

    func testM2SHighNumberBackbonesAdaptButStillNeedDefinitions() throws {
        struct Record: Decodable { let id: Int; let cxsmiles: String }
        let url = try XCTUnwrap(
            Bundle.module.url(
                forResource: "m2s-backbones", withExtension: "json",
                subdirectory: "RGroupEnumerationFixtures"))
        let records = try JSONDecoder().decode([Record].self, from: Data(contentsOf: url))
        for record in records where [35, 47].contains(record.id) {
            let family = try parse(record.cxsmiles)
            XCTAssertNoThrow(try CDKRGroupQueryManipulator.toRGroupQuery(family))
            XCTAssertEqual(CDKMarkushMatcher.match(try parse("CO"), in: family).verdict, .unresolved)
        }
    }

    func testTetrahedralPermutationEnantiomersAndMissingStereo() throws {
        // OpenSMILES tetrahedral-order examples, plus implicit/explicit hydrogen equivalence.
        let family = "N[C@](Br)(O)C"
        for text in [family, "N[C@@](Br)(C)O", "Br[C@](O)(N)C"] {
            let result = try match(text, family)
            XCTAssertEqual(result.verdict, .match, "\(text): \(result.detail ?? "")")
        }
        XCTAssertEqual(try match("N[C@@](Br)(O)C", family).verdict, .noMatch)
        XCTAssertEqual(try match("NC(Br)(O)C", family).reason, .unspecifiedStereochemistry)
        XCTAssertEqual(try match(family, "NC(Br)(O)C").verdict, .match)
        XCTAssertEqual(try match("[H][C@@](F)(Cl)Br", "F[C@H](Cl)Br").verdict, .match)
    }

    func testStereoSurvivesSubstitutionNextToRootCenterAndInAlternative() throws {
        let root = "F[C@H](*)Cl |$;;R1$,RG:_R1={Br}|"
        XCTAssertEqual(try match("F[C@H](Br)Cl", root).verdict, .match)
        XCTAssertEqual(try match("F[C@@H](Br)Cl", root).verdict, .noMatch)
        let alternative = "N* |$;R1$,RG:_R1={*[C@H](F)Cl |$_AP1$|}|"
        let good = try match("N[C@H](F)Cl", alternative)
        XCTAssertEqual(good.verdict, .match, good.detail ?? "")
        XCTAssertEqual(try match("N[C@@H](F)Cl", alternative).verdict, .noMatch)
    }

    func testAlkeneDirectionsReversalAndOppositeGeometry() throws {
        let family = "F/C=C/Cl"
        for text in [family, "Cl/C=C/F", "F\\C=C\\Cl"] {
            let result = try match(text, family)
            XCTAssertEqual(result.verdict, .match, "\(text): \(result.detail ?? "")")
        }
        XCTAssertEqual(try match("F/C=C\\Cl", family).verdict, .noMatch)
        XCTAssertEqual(try match("FC=CCl", family).reason, .unspecifiedStereochemistry)
        XCTAssertEqual(try match(family, "FC=CCl").verdict, .match)
        XCTAssertEqual(try match("F/C(Br)=C(Cl)/I", "Br/C(F)=C(I)/Cl").verdict, .match)
    }

    func testTetrahedralWedgesDoNotInventAlkeneStereo() throws {
        let unspecified = "F[C@H](C)C=C[C@H](O)Br"
        for concrete in ["F[C@H](C)/C=C/[C@H](O)Br", "F[C@H](C)/C=C\\[C@H](O)Br"] {
            let result = try match(concrete, unspecified)
            XCTAssertEqual(result.verdict, .match, result.detail ?? "")
        }
    }

    func testExplicitAlkeneReferencesAgreeWithDirectionalSMILES() throws {
        var family = try parse("FC=CCl")
        let bondIndex = try XCTUnwrap(family.bonds.firstIndex { $0.order == .double })
        family.bonds[bondIndex].doubleBondStereo = .trans
        family.bonds[bondIndex].stereoReferenceAtomIDs = family.atoms.map(\.id)
        XCTAssertEqual(CDKMarkushMatcher.match(try parse("F/C=C/Cl"), in: family).verdict, .match)
        XCTAssertEqual(CDKMarkushMatcher.match(try parse("F/C=C\\Cl"), in: family).verdict, .noMatch)
        family.bonds[bondIndex].stereoReferenceAtomIDs = Array(family.atoms.map(\.id).reversed())
        XCTAssertEqual(CDKMarkushMatcher.match(try parse("F/C=C/Cl"), in: family).verdict, .match)
    }

    func testRestHAndOccurrenceDependenciesUseClosedFamilySemantics() throws {
        let family = "C* |$;R1$,RG:_R1={O},LOG={_R1:;H;0,1}|"
        XCTAssertEqual(try match("C", family).verdict, .match)
        XCTAssertEqual(try match("CO", family).verdict, .match)
        XCTAssertEqual(try match("CN", family).verdict, .noMatch)
        let noRestH = family.replacingOccurrences(of: ":;H;", with: ":;;")
        XCTAssertEqual(try match("CN", noRestH).verdict, .noMatch)
        let dependent = "*C* |$R1;;R2$,RG:_R1={O},_R2={N},LOG={_R1:_R2;;0,1._R2:;;0,1}|"
        XCTAssertEqual(try match("CO", dependent).verdict, .noMatch)
        XCTAssertEqual(try match("OCN", dependent).verdict, .match)
    }

    func testMalformedOrUnsupportedGraphsNeverProduceConfidentResults() throws {
        var duplicate = try parse("CO")
        duplicate.atoms.append(duplicate.atoms[0])
        XCTAssertEqual(CDKMarkushMatcher.match(duplicate, in: try parse("CO")).verdict, .unresolved)
        var enhanced = try parse("F[C@H](Cl)Br")
        enhanced.atoms[1].cxStereoGroup = 1
        XCTAssertEqual(CDKMarkushMatcher.match(enhanced, in: enhanced).reason, .unsupported)
        var incomplete = try parse("F[C@H](Cl)Br")
        incomplete.atoms[1].ligandOrderingAtomIDs = nil
        XCTAssertEqual(CDKMarkushMatcher.match(incomplete, in: incomplete).reason, .unsupported)
        var parallel = try parse("CO")
        parallel.bonds.append(Bond(id: 100, a1: parallel.atoms[0].id, a2: parallel.atoms[1].id, order: .single))
        XCTAssertEqual(CDKMarkushMatcher.match(parallel, in: parallel).reason, .invalidInput)
    }

    func testReorderedIDsAndStereoSymmetrySearch() throws {
        let source = try parse("F[C@H](Cl)Br")
        var reordered = source
        let idMap = Dictionary(uniqueKeysWithValues: source.atoms.enumerated().map { ($1.id, 100 - $0 * 7) })
        reordered.atoms = source.atoms.reversed().map { atom in
            var copy = RGroupEnumerationPlan.copyAtom(atom, id: idMap[atom.id]!)
            copy.ligandOrderingAtomIDs = atom.ligandOrderingAtomIDs?.map { idMap[$0]! }
            return copy
        }
        reordered.bonds = source.bonds.reversed().map { bond in
            RGroupEnumerationPlan.copyBond(
                bond, id: bond.id + 100,
                a1: idMap[bond.a1]!, a2: idMap[bond.a2]!, mapReference: { idMap[$0] })
        }
        let result = CDKMarkushMatcher.match(reordered, in: source)
        XCTAssertEqual(result.verdict, .match, result.detail ?? "")
        XCTAssertEqual(result.witness?.compoundAtomIDByMemberAtomID, idMap)
        // Identical disconnected components must permit a later stereo-compatible mapping.
        XCTAssertEqual(
            try match(
                "F[C@H](Cl)Br.F[C@@H](Cl)Br",
                "F[C@@H](Cl)Br.F[C@H](Cl)Br"
            ).verdict, .match)
    }

    func testUnknownStereoDoesNotHideKnownMismatchAtAnotherCenter() throws {
        XCTAssertEqual(try match("F[C@@H](Cl)C(Br)I", "F[C@H](Cl)[C@H](Br)I").verdict, .noMatch)
        let family = "N* |$;R1$,RG:_R1={*[C@H](F)Cl |$_AP1$|},{*C(F)Cl |$_AP1$|}|"
        XCTAssertEqual(try match("NC(F)Cl", family).verdict, .match)
    }

    func testEmptyFiniteFamilyCanProduceExhaustiveNegative() throws {
        let result = try match("CO", "C* |$;R1$,RG:_R1={O},LOG={_R1:;;2}|")
        XCTAssertEqual(result.verdict, .noMatch)
        XCTAssertTrue(result.isExhaustive)
        XCTAssertEqual(result.configurationsChecked, 0)
    }

    func testUnsupportedStereoAndMalformedFamilyAreUnresolved() throws {
        for text in ["F/C=CCl", "F[C@H](Cl)Br |r|", "N* |$;R1$,RG:_R1={*[C@H](F)Cl |$_AP1$,r|}|"] {
            let result = try match("NC(F)Cl", text)
            XCTAssertEqual(result.verdict, .unresolved, text)
        }
        var duplicate = try parse(simpleFamily)
        duplicate.atoms.append(duplicate.atoms[0])
        XCTAssertEqual(CDKMarkushMatcher.match(try parse("CO"), in: duplicate).verdict, .unresolved)
        var invalidLabel = try parse(simpleFamily)
        invalidLabel.atoms[1].rGroupLabel = -1
        XCTAssertEqual(CDKMarkushMatcher.match(try parse("CO"), in: invalidLabel).verdict, .unresolved)
        var crossComponent = try parse(simpleFamily)
        let fragmentID = try XCTUnwrap(crossComponent.atoms.first { $0.rGroupMembership != nil }?.id)
        crossComponent.bonds.append(Bond(id: 999, a1: crossComponent.atoms[0].id, a2: fragmentID, order: .single))
        XCTAssertEqual(CDKMarkushMatcher.match(try parse("CO"), in: crossComponent).verdict, .unresolved)
        var dangling = try parse(simpleFamily)
        dangling.bonds.append(Bond(id: 999, a1: dangling.atoms[0].id, a2: -1000, order: .single))
        XCTAssertEqual(CDKMarkushMatcher.match(try parse("CO"), in: dangling).verdict, .unresolved)
    }

    func testCancellationIsUnresolved() async throws {
        let source = try parse(simpleFamily)
        let compound = try parse("CO")
        let task = Task {
            withUnsafeCurrentTask { $0?.cancel() }
            return CDKMarkushMatcher.match(compound, in: source)
        }
        let result = await task.value
        XCTAssertEqual(result.verdict, .unresolved)
        XCTAssertEqual(result.reason, .cancelled)
    }
}
