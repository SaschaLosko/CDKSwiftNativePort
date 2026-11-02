import XCTest
@testable import CDKSwiftNativePort

final class MarkushRepeatRangeTests: XCTestCase {
    private let parser = CDKSmilesParser(flavor: .cdkDefault, generateCoordinates: false)

    func testZeroRepeatJoinsBoundariesAndPreservesSource() throws {
        let source = try parser.parseSmiles("FCCCl |Sg:n:1,2:0-2:ht|")
        let snapshot = source
        let members = try CDKMarkushEnumerator.configurations(of: source)
        let expected = try ["FCl", "FCCCl", "FCCCCCl"].map { try parser.parseSmiles($0) }
        XCTAssertEqual(members.map(\.repeatCounts), [[0], [1], [2]])
        for (member, expectedGraph) in zip(members, expected) {
            XCTAssertEqual(
                try CDKInChIGenerator(molecule: member.molecule).getInchi(),
                try CDKInChIGenerator(molecule: expectedGraph).getInchi())
            XCTAssertEqual(Set(member.sourceByMemberAtomID.keys), Set(member.molecule.atoms.map(\.id)))
        }
        XCTAssertEqual(source, snapshot)
    }

    func testZeroRepeatRejectsSelfOrDuplicateEdges() throws {
        for text in ["C1CC1 |Sg:n:1,2:0:ht|", "C1CCC1 |Sg:n:1,2:0:ht|"] {
            XCTAssertThrowsError(try CDKMarkushEnumerator.configurations(of: parser.parseSmiles(text)))
        }
    }

    func testLargerExplicitRepeatKeepsConfigurationAndAtomBudgets() throws {
        let family = try parser.parseSmiles("FCOCl |Sg:n:1,2:65-70:ht|")
        let members = try CDKMarkushEnumerator.configurations(of: family)
        XCTAssertEqual(members.count, 6)
        XCTAssertEqual(members.last?.molecule.atoms.count, 142)
        XCTAssertThrowsError(try CDKMarkushEnumerator.configurations(of: family, maximumConfigurations: 5))
        XCTAssertThrowsError(try CDKMarkushEnumerator.configurations(of: family, maximumAtoms: 141))
        XCTAssertThrowsError(
            try CDKMarkushEnumerator.configurations(of: parser.parseSmiles("FCCCl |Sg:n:1,2:1001:ht|")))
    }

    func testGraphOnlyParsingPreservesChemicalIdentityAndExplicitCoordinates() throws {
        for text in ["[13CH3][C@H](O)C(=O)[O-]", "F/C=C/Cl", "C* |$;R1$,RG:_R1={O},{N}|"] {
            let graph = try parser.parseSmiles(text)
            let drawn = try CDKSmilesParser(flavor: .cdkDefault).parseSmiles(text)
            XCTAssertEqual(graph.atoms.map(\.element), drawn.atoms.map(\.element))
            if !CDKMarkushStructure.requiresSpecifiedStructure(graph) {
                XCTAssertEqual(
                    try CDKInChIGenerator(molecule: graph).getInchi(),
                    try CDKInChIGenerator(molecule: drawn).getInchi())
            } else {
                XCTAssertEqual(
                    CDKMarkushRelations.compare(graph, drawn, relation: .equivalent).verdict, .holds)
            }
        }
        let coordinates = try parser.parseSmiles("CC |(1,2,;3,4,)|")
        XCTAssertEqual(coordinates.atoms[0].position.x, 1)
        XCTAssertEqual(coordinates.atoms[1].position.y, 4)
    }
}
