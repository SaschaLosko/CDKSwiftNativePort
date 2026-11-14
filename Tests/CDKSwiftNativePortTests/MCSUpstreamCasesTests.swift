// Selected numerical fixtures from CDK 2.14.1 MCSTest (LGPL-2.1-or-later).
// https://github.com/cdk/cdk/blob/cdk-2.14.1/tool/smarts/src/test/java/org/openscience/cdk/isomorphism/MCSTest.java
import XCTest
@testable import CDKSwiftNativePort

final class MCSUpstreamCasesTests: XCTestCase {
    func testUpstreamStereoScoresAndMappingMultiplicities() throws {
        let alanine = "N[C@@H](C)C(=O)O"
        let threonine = "C[C@H]([C@@H](C(=O)O)N)O"
        let allo = "C[C@@H]([C@@H](C(=O)O)N)O"
        let mirror = "C[C@@H]([C@H](C(=O)O)N)O"
        let cases: [(String, String, Int, Int, Int)] = [
            (alanine, alanine, 6, 5, 1),
            (alanine, "N[C@H](C)C(=O)O", 5, 4, 2),
            (alanine, "N[C@@H](CS)C(=O)O", 6, 5, 1),
            (alanine, "N[C@H](CS)C(=O)O", 5, 4, 2),
            (alanine, "N[C@@]([H])(C)C(=O)O", 6, 5, 1),
            (alanine, "N[C@]([H])(C)C(=O)O", 5, 4, 2),
            ("N[C@]([H])(C)C(=O)O", alanine, 5, 4, 2),
            (alanine, "N[C@H](C)C(=O)OC", 5, 4, 2),
            ("C/C=C/C", "C/C=C\\C", 3, 2, 4),
            ("C/C=C/C", "C(\\C)=C/C", 4, 3, 2),
            ("F/C(Cl)=C/F", "Cl\\C(F)=C/F", 5, 4, 1),
            ("F/C(Cl)=C/F", "F/C(Cl)=C\\F", 4, 3, 1),
            ("C[C@@H](O)[C@@H](C)O", "C[C@@H](O)[C@@H](C)O", 6, 5, 2),
            ("C[C@@H](O)[C@H](C)O", "C[C@@H](O)[C@H](C)O", 6, 5, 1),
            ("C[C@@H](O)[C@@H](C)O", "C[C@@H](O)[C@H](C)O", 5, 4, 4),
            ("C[C@](N)(O)CC", "C[C@@](N)(O)Cl", 4, 3, 1),
            (threonine, allo, 7, 6, 2),
            (threonine, mirror, 6, 5, 3),
            (threonine, mirror + " |&1:1,2|", 8, 7, 1),
            (threonine, allo + " |&1:1,2|", 7, 6, 3),
            (threonine, allo + " |o1:1,2|", 7, 6, 3),
            (threonine, allo + " |&1:1,&2:2|", 8, 7, 1),
            ("C[C@@H]1CCCC[C@H]1C", "C[C@@H]1CCCC[C@@H]1C", 8, 7, 2),
        ]
        for (query, target, atoms, bonds, count) in cases {
            try check(
                query, target, options: .init(matchStereochemistry: true, maximumSearchSteps: 5_000_000),
                expected: [atoms, bonds, count])
        }
    }

    func testUpstreamDisconnectedScoringAndFragmentBounds() throws {
        let cases: [(String, String, Int, Int, Int, Int)] = [
            ("CCOCC", "CCSCC", 1, 4, 2, 8),
            ("CCCC.CC", "CCCCC", 1, 5, 3, 16),
            ("CCCC.CC", "CCCCC", 2, 4, 3, 4),
            ("C1CC1.CC", "CC1CC1", 1, 3, 3, 6),
            ("CCOC(C)=O", "CC(=O)O.OCC", 1, 6, 4, 3),
            ("CCOC(C)=O", "CC(=O)O.OCC", 2, 6, 4, 1),
            ("CC.O", "OCC", 1, 2, 1, 2),
            ("CC(C)=O", "CC(C)O", 1, 3, 2, 2),
        ]
        for (query, target, minimum, atoms, bonds, count) in cases {
            try check(query, target, options: .init(minimumBondsPerFragment: minimum), expected: [atoms, bonds, count])
        }
    }

    func testUpstreamCombinedRingAndStereoConstraints() throws {
        try check(
            "C[C@@H]1C[C@H]1C", "C[C@@H]1C[C@@H]1C", options: .init(completeRings: true, matchStereochemistry: true),
            expected: [4, 4, 4])
        try check("C1CCCCC1", "CC1CCCCC1", options: .init(completeRings: true), expected: [6, 6, 12])
        try check("Cc1ccccc1", "CC1CCCCC1", options: .init(completeRings: true), expected: [1, 0, 1])
        try check(
            "[C;R]C(=O)[N,O]", "NC(=O)C1CCCCC1", options: .init(completeRings: true), expected: [3, 2, 1], smarts: true)
    }

    private func check(_ query: String, _ target: String, options: CDKMCSOptions, expected: [Int], smarts: Bool = false)
        throws
    {
        let parser = CDKSmilesParser(generateCoordinates: false)
        let search =
            smarts
            ? try CDKMCS(query: CDKSmartsParser.parse(query), options: options)
            : try CDKMCS(query: parser.parseSmiles(query), options: options)
        let result = try search.search(in: parser.parseSmiles(target))
        XCTAssertTrue(result.isMaximumProven, "\(query) / \(target): \(result.termination)")
        XCTAssertEqual(
            [result.mappedAtomCount, result.mappedBondCount, result.mappings.count], expected, "\(query) / \(target)")
    }
}
