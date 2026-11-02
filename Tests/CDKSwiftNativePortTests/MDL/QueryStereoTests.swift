import XCTest
@testable import CDKSwiftNativePort

final class QueryStereoTests: XCTestCase {
    private func match(_ query: String, _ target: String) throws -> CDKQueryMatchResult {
        CDKQueryMatcher.match(try CDKSmartsParser.parse(query), in: try CDKQueryStructure.parseSmiles(target))
    }

    func testTetrahedralOrderingAndMissingStereo() throws {
        XCTAssertEqual(try match("F[C@H](Cl)Br", "F[C@H](Cl)Br").verdict, .match)
        XCTAssertEqual(try match("F[C@H](Cl)Br", "Cl[C@@H](F)Br").verdict, .match)
        XCTAssertEqual(try match("F[C@H](Cl)Br", "F[C@@H](Cl)Br").verdict, .noMatch)
        XCTAssertEqual(try match("F[C@H](Cl)Br", "FC(Cl)Br").verdict, .unresolved)
        XCTAssertEqual(try match("F[C@?H](Cl)Br", "FC(Cl)Br").verdict, .match)
        XCTAssertEqual(try match("FC(Cl)Br", "F[C@@H](Cl)Br").verdict, .match)
        XCTAssertEqual(try match("F[C@H](Cl)Br", "[H][C@@](F)(Cl)Br").verdict, .match)
    }

    // Examples from CDK 2.13 SquarePlanar and Octahedral class documentation.
    func testCoordinationStereoReadAndMatch() throws {
        XCTAssertEqual(try match("[NH3][Pt@SP1]([NH3])(Cl)Cl", "[NH3][Pt@SP2](Cl)([NH3])Cl").verdict, .match)
        XCTAssertEqual(try match("[NH3][Pt@SP1]([NH3])(Cl)Cl", "[NH3][Pt@SP2]([NH3])(Cl)Cl").verdict, .noMatch)
        XCTAssertEqual(try match("C[Co@OH8](F)(Br)(Cl)(I)S", "C[Co@OH1](F)(Cl)(Br)(I)S").verdict, .match)
        XCTAssertEqual(try match("C[Co@OH1](F)(Cl)(Br)(I)S", "C[Co@OH2](F)(Cl)(Br)(I)S").verdict, .noMatch)
        XCTAssertEqual(try match("F[P@TB1](Cl)(Br)(I)N", "F[P@TB2](Cl)(Br)(I)N").verdict, .noMatch)
        XCTAssertEqual(try match("F[P@TB1](Cl)(Br)(I)N", "F[P](Cl)(Br)(I)N").verdict, .unresolved)
        let structure = try CDKQueryStructure.parseSmiles("C[Co@OH8](F)(Br)(Cl)(I)S")
        XCTAssertEqual(
            try JSONDecoder().decode(CDKQueryStructure.self, from: JSONEncoder().encode(structure)), structure)
        XCTAssertThrowsError(try CDKQueryStructure.parseSmiles("[C,N]"))
        XCTAssertThrowsError(try CDKQueryStructure.parseSmiles("[C;H3]"))
        XCTAssertThrowsError(try CDKQueryStructure.parseSmiles("*"))
    }

    func testEveryUpstreamCoordinationPermutation() throws {
        let url = try XCTUnwrap(
            Bundle.module.url(
                forResource: "query-stereo-permutations", withExtension: "json",
                subdirectory: "RGroupEnumerationFixtures"))
        let data = try JSONDecoder().decode([String: [[Int]]].self, from: Data(contentsOf: url))
        for geometry in [CDKQueryAtomStereo.Geometry.squarePlanar, .trigonalBipyramidal, .octahedral] {
            let rows = try XCTUnwrap(data[geometry.rawValue]); let n = geometry.carrierCount
            XCTAssertEqual(rows.count, geometry.configurationCount)
            var covered: Set<[Int]> = []
            for (index, row) in rows.enumerated() {
                let source = CDKQueryAtomStereo(
                    focus: 100, carriers: Array(0..<n), geometry: geometry, configuration: index + 1)
                for offset in stride(from: 0, to: row.count, by: n) {
                    let carriers = Array(row[offset..<offset + n])
                    covered.insert(carriers)
                    let normalized = CDKQueryAtomStereo(
                        focus: 100, carriers: carriers, geometry: geometry, configuration: 1)
                    let map = Dictionary(uniqueKeysWithValues: (0..<n).map { ($0, $0) })
                    XCTAssertTrue(
                        QueryStereoSupport.equivalent(source, normalized, mapping: map),
                        "\(geometry) order \(index + 1)")
                    for other in 1...geometry.configurationCount where other != index + 1 {
                        var different = source; different.configuration = other
                        XCTAssertFalse(QueryStereoSupport.equivalent(different, normalized, mapping: map))
                    }
                }
            }
            XCTAssertEqual(covered.count, (1...n).reduce(1, *))
        }
    }

    func testRelativeGroupsMustInvertTogether() throws {
        var query = try CDKSmartsParser.parse("F[C@H](Cl)[C@H](Br)I")
        for index in query.stereo.indices { query.stereo[index].relativeGroup = 1 }
        XCTAssertEqual(
            CDKQueryMatcher.match(query, in: try CDKQueryStructure.parseSmiles("F[C@@H](Cl)[C@@H](Br)I")).verdict,
            .match)
        XCTAssertEqual(
            CDKQueryMatcher.match(query, in: try CDKQueryStructure.parseSmiles("F[C@@H](Cl)[C@H](Br)I")).verdict,
            .noMatch)
    }

    func testPartialCoordinationQueryMapsUnspecifiedLigandsAndCrossGeometry() throws {
        var query = try CDKSmartsParser.parse("Cl[Pt](Br)I")
        let focus = query.molecule.atoms[1].id
        query.stereo = [.init(focus: focus, carriers: [0, 2, focus, 3], geometry: .squarePlanar, configuration: 1)]
        let result = CDKQueryMatcher.match(query, in: try CDKQueryStructure.parseSmiles("Cl[Pt@SP1](Br)(F)I"))
        XCTAssertEqual(result.verdict, .match, result.detail ?? "")
        var tbQuery = try CDKSmartsParser.parse("F[Co](Cl)Br")
        tbQuery.stereo = [.init(focus: 1, carriers: [0, 2, 1, 1, 3], geometry: .trigonalBipyramidal, configuration: 1)]
        let octahedral = try CDKQueryStructure.parseSmiles("F[Co@OH1](N)(Cl)(I)(O)Br")
        XCTAssertEqual(CDKQueryMatcher.match(tbQuery, in: octahedral).verdict, .match)
    }

    func testGeometricStereoAndRecursiveUnknownNegation() throws {
        XCTAssertEqual(try match("F/C=C/Cl", "F/C=C/Cl").verdict, .match)
        XCTAssertEqual(try match("F/C=C/Cl", "F/C=C\\Cl").verdict, .noMatch)
        XCTAssertEqual(try match("F/C=C/Cl", "FC=CCl").verdict, .unresolved)
        // Other atoms definitely fail the anchored query.
        XCTAssertEqual(try match("[!$(F[C@H](Cl)Br)]", "FC(Cl)Br").verdict, .match)
        let result = try match("[F;!$(F[C@H](Cl)Br)]", "FC(Cl)Br")
        XCTAssertEqual(result.verdict, .unresolved)
    }
}
