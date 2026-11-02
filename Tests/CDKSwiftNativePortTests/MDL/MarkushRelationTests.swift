import XCTest
@testable import CDKSwiftNativePort

final class MarkushRelationTests: XCTestCase {
    private func parse(_ text: String) throws -> Molecule {
        try CDKSmilesParser(flavor: .cdkDefault).parseSmiles(text)
    }
    private let family = "C* |$;R1$,RG:_R1={O},{N}|"

    func testRelationsRequireExhaustionAndProvideCounterexamples() throws {
        let a = try parse(family)
        let b = try parse("C* |$;X$,RG:_X={N},{F}|")
        let overlap = CDKMarkushRelations.compare(a, b, relation: .overlaps)
        XCTAssertEqual(overlap.verdict, .holds)
        let witness = try XCTUnwrap(overlap.witness)
        XCTAssertEqual(CDKMarkushMatcher.match(witness.molecule, in: try parse("CN")).verdict, .match)
        XCTAssertEqual(CDKMarkushRelations.compare(a, b, relation: .disjoint).verdict, .doesNotHold)
        let subset = CDKMarkushRelations.compare(a, b, relation: .containedIn)
        XCTAssertEqual(subset.verdict, .doesNotHold)
        XCTAssertEqual(CDKMarkushMatcher.match(try XCTUnwrap(subset.witness).molecule, in: b).verdict, .noMatch)
        XCTAssertEqual(CDKMarkushRelations.compare(try parse("CN"), a, relation: .containedIn).verdict, .holds)
        let same = CDKMarkushRelations.compare(a, try parse("C* |$;Y$,RG:_Y={N},{O},{O}|"), relation: .equivalent)
        XCTAssertEqual(same.verdict, .holds)
        XCTAssertTrue(same.isExhaustive)
        let different = CDKMarkushRelations.compare(try parse("CN"), a, relation: .equivalent)
        XCTAssertEqual(different.verdict, .doesNotHold)
        XCTAssertEqual(different.direction, .rightToLeft)
        XCTAssertEqual(CDKMarkushRelations.compare(a, try parse("CC"), relation: .disjoint).verdict, .holds)
    }

    func testLimitsAndUnsupportedFamiliesStayUnresolved() throws {
        let a = try parse(family)
        for relation in [CDKMarkushRelation.overlaps, .containedIn, .equivalent, .disjoint] {
            let limited = CDKMarkushRelations.compare(
                a, try parse("CC"), relation: relation,
                options: .init(maximumSearchSteps: 2))
            XCTAssertEqual(limited.verdict, .unresolved)
            XCTAssertEqual(limited.reason, .searchLimit)
            XCTAssertFalse(limited.isExhaustive)
            XCTAssertEqual(
                CDKMarkushRelations.compare(try parse("C* |$;X$|"), a, relation: relation).verdict, .unresolved)
        }
        XCTAssertEqual(
            CDKMarkushRelations.compare(
                a, try parse("CN"), relation: .overlaps,
                options: .init(maximumConfigurations: 1)
            ).reason, .configurationLimit)
    }

    func testStereoOverlapTriesBothDirectionsWithoutClaimingContainment() throws {
        let unspecified = try parse("FC(Cl)Br")
        let specified = try parse("F[C@H](Cl)Br")
        XCTAssertEqual(CDKMarkushRelations.compare(unspecified, specified, relation: .overlaps).verdict, .holds)
        XCTAssertEqual(CDKMarkushRelations.compare(unspecified, specified, relation: .containedIn).verdict, .unresolved)
        XCTAssertEqual(CDKMarkushRelations.compare(specified, unspecified, relation: .containedIn).verdict, .holds)
        XCTAssertEqual(
            CDKMarkushRelations.compare(specified, try parse("F[C@@H](Cl)Br"), relation: .overlaps).verdict,
            .doesNotHold)
    }

    func testWitnessMapsOriginalSitesAndSelectedFragmentAtoms() throws {
        let source = try parse("C* |$;X$,RG:_X={O},{N}|")
        let original = source
        let result = CDKMarkushMatcher.match(try parse("NC"), in: source)
        let member = try XCTUnwrap(result.witness?.configuration)
        let assignment = try XCTUnwrap(member.assignments.first)
        XCTAssertEqual(assignment.alternativeIndex, 1)
        XCTAssertEqual(assignment.siteAtomID, source.atoms.first { $0.element == "X" || $0.aliasLabel == "X" }?.id)
        XCTAssertEqual(assignment.scaffoldAttachmentAtomIDs.count, 1)
        XCTAssertEqual(assignment.fragmentAttachmentAtomIDs.count, 1)
        for (fragment, memberID) in assignment.memberAtomIDByFragmentAtomID {
            XCTAssertEqual(source.atom(id: fragment)?.element, "N")
            XCTAssertEqual(member.sourceByMemberAtomID[memberID]?.atomID, fragment)
            XCTAssertEqual(member.sourceByMemberAtomID[memberID]?.siteAtomID, assignment.siteAtomID)
            XCTAssertEqual(member.sourceByMemberAtomID[memberID]?.alternativeIndex, 1)
        }
        XCTAssertEqual(member.sourceByMemberAtomID.count, member.molecule.atoms.count)
        XCTAssertEqual(source, original)
        let repeated = try CDKMarkushEnumerator.configurations(of: parse("FCCl |LN:1:1.3|"))
        XCTAssertEqual(repeated.last?.sourceByMemberAtomID.values.filter { $0.atomID == 2 }.count, 3)
    }
}
