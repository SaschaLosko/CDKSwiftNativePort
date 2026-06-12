import XCTest
@testable import CDKSwiftNativePort

final class SubstructureSearchTests: XCTestCase {
    private let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    func testFindsCarbonylSubstructure() throws {
        let target = try parser.parseSmiles("CC(=O)N")
        let carbonyl = try parser.parseSmiles("C=O")
        let nitroso = try parser.parseSmiles("N=O")

        XCTAssertTrue(target.containsSubstructure(carbonyl))
        XCTAssertFalse(target.containsSubstructure(nitroso))
    }

    func testFindsBenzeneInNaphthaleneButNotReverse() throws {
        let benzene = try parser.parseSmiles("c1ccccc1")
        let naphthalene = try parser.parseSmiles("c1ccc2ccccc2c1")

        XCTAssertTrue(naphthalene.containsSubstructure(benzene))
        XCTAssertFalse(benzene.containsSubstructure(naphthalene))
        XCTAssertFalse(naphthalene.substructureMatches(of: benzene).isEmpty)
    }

    func testChargeSensitiveSubstructureMatching() throws {
        let anionTarget = try parser.parseSmiles("CC[C-]")
        let neutralTarget = try parser.parseSmiles("CCC")
        let carbanion = try parser.parseSmiles("[C-]")

        XCTAssertTrue(anionTarget.containsSubstructure(carbanion))
        XCTAssertFalse(neutralTarget.containsSubstructure(carbanion))
    }

    func testWildcardQueryAtomMatchesAnyElement() throws {
        let methanol = try parser.parseSmiles("CO")
        let wildcardAlcohol = try parser.parseSmiles("*O")

        XCTAssertTrue(methanol.containsSubstructure(wildcardAlcohol))
    }

    func testReturnsAllMappedOccurrencesForRepeatedGroup() throws {
        let propane = try parser.parseSmiles("CCC")
        let ethyl = try parser.parseSmiles("CC")
        let matches = propane.substructureMatches(of: ethyl)

        let uniqueTargetAtomSets = Set(matches.map { Set($0.targetAtomIDs) })
        XCTAssertEqual(uniqueTargetAtomSets, Set([Set([1, 2]), Set([2, 3])]))
        XCTAssertGreaterThanOrEqual(matches.count, 2)
    }
}
