import XCTest
@testable import CDKSwiftNativePort

final class QueryExpressionTests: XCTestCase {
    private func verdict(_ query: String, _ target: String, mode: CDKQueryMatcher.Mode = .substructure) throws
        -> CDKQueryMatchResult
    {
        CDKQueryMatcher.match(
            try CDKSmartsParser.parse(query), in: try CDKSmilesParser(generateCoordinates: false).parseSmiles(target),
            mode: mode)
    }

    // Independently specified examples from Daylight's SMARTS theory, sections 4.1–4.4.
    func testAtomPrimitivesAndBooleanPrecedence() throws {
        let cases: [(String, String, Bool)] = [
            ("[C,N]", "N", true), ("[C,N]", "O", false), ("[!C;R]", "C1CO1", true),
            ("[!C;R]", "CCO", false), ("[c,n;H1]", "n1ccccc1", true),
            ("[n;H1]", "n1ccccc1", false), ("[n;H1]", "c1cc[nH]c1", true),
            ("[X3&H0]", "CC(=O)C", true), ("[D3]", "CC(C)C", true),
            ("[D4]", "C", false), ("[X4]", "C", true), ("[v4]", "C", true),
            ("[35Cl]", "[35Cl-]", true), ("[35Cl]", "Cl", false),
            ("[#7;+1]", "[NH4+]", true), ("[#7;+0]", "[NH4+]", false),
            ("[F,Cl,Br,I]", "Br", true), ("[F,Cl,Br,I]", "C", false),
            ("[He]", "[He]", true), ("[Xe]", "[Xe]", true),
            ("[H]", "[H][H]", true), ("[2H]", "[2H]O", true),
            ("[h]", "C", true), ("[h0]", "C", false),
            ("[r3]", "C1CC1", true), ("[r4]", "C1CC1", false),
            ("[R0]", "CC", true), ("[R1]", "C1CC1", true), ("[x2]", "C1CC1", true),
            ("[nH;v3]", "c1cc[nH]c1", true), ("[c;v4]", "c1ccc2ccccc2c1", true),
        ]
        for (query, target, expected) in cases {
            let result = try verdict(query, target)
            XCTAssertEqual(
                result.verdict, expected ? .match : .noMatch, "\(query) in \(target): \(result.detail ?? "")")
            if !expected { XCTAssertTrue(result.isExhaustive) }
        }
        XCTAssertEqual(try verdict("[C,N;H3]", "N").verdict, .match)
        XCTAssertEqual(try verdict("[C,N;H3]", "C").verdict, .noMatch)
        XCTAssertEqual(try verdict("[C,N&H3]", "C").verdict, .match)
    }

    func testRecursiveQueriesAreAnchoredAndNegationPreservesChemistry() throws {
        XCTAssertEqual(try verdict("[C;$(*=O)]", "CC=O").verdict, .match)
        XCTAssertEqual(try verdict("[N;$(*=O)]", "NC=O").verdict, .noMatch)
        XCTAssertEqual(try verdict("[C;!$(*=O)]", "C=O").verdict, .noMatch)
        XCTAssertEqual(try verdict("[$(*C);$(*CC)]", "CCC").verdict, .match)
        XCTAssertEqual(try verdict("[$(*C);$(*CC)]", "CC").verdict, .noMatch)
        XCTAssertEqual(try verdict("[C;$([C;$(*=O)])]", "C=O").verdict, .match)
    }

    func testBondQueriesRingsComponentsAndWholeGraph() throws {
        let biphenyl = try CDKSmilesParser(generateCoordinates: false).parseSmiles("c1ccccc1-c2ccccc2")
        XCTAssertEqual(
            biphenyl.bonds.filter { $0.order == .single }.count, 1, "The biaryl connecting bond must remain single.")
        XCTAssertEqual(try CDKSmartsParser.parse("c-c").bonds.values.first, .order(.single))
        let cases: [(String, String, Bool)] = [
            ("C~C", "C=C", true), ("C-,=C", "C=C", true), ("C!#C", "C#C", false),
            ("C@C", "C1CC1", true), ("C!@C", "C1CC1", false),
            ("c:c", "C1=CC=CC=C1", true), ("c-c", "c1ccccc1-c2ccccc2", true),
            ("c-c", "c1ccccc1", false), ("C.C", "CCC", true), ("C.C", "C", false),
        ]
        for (query, target, expected) in cases {
            XCTAssertEqual(try verdict(query, target).verdict, expected ? .match : .noMatch, query)
        }
        XCTAssertEqual(try verdict("CO", "CCO").verdict, .match)
        XCTAssertEqual(try verdict("CO", "CCO", mode: .wholeGraph).verdict, .noMatch)
        XCTAssertEqual(try verdict("C.C", "CC", mode: .wholeGraph).verdict, .noMatch)
        XCTAssertEqual(try verdict("CO", "CO", mode: .wholeGraph).verdict, .match)
        XCTAssertEqual(try verdict("C1CC1", "C1CC1").verdict, .match)
    }

    func testLimitsMalformedQueriesAndCodingRoundTrip() throws {
        let query = try CDKSmartsParser.parse("[C;$(*C)]~[O,N]")
        XCTAssertEqual(try JSONDecoder().decode(CDKQueryPattern.self, from: JSONEncoder().encode(query)), query)
        let target = try CDKSmilesParser(generateCoordinates: false).parseSmiles("CCO")
        let limited = CDKQueryMatcher.match(query, in: target, options: .init(maximumSearchSteps: 1))
        XCTAssertEqual(limited.verdict, .unresolved); XCTAssertEqual(limited.reason, .searchLimit)
        XCTAssertFalse(limited.isExhaustive)
        for source in [
            "", "[]", "C(", "C()", "C1CC", "C11", "C..C", "C-", "[C,]", "[$(C]", "[C^3]", "C>>O",
            "[Pt@SP4](F)(Cl)(Br)I",
        ] {
            XCTAssertThrowsError(try CDKSmartsParser.parse(source), source)
        }
        var invalid = query
        invalid.atoms.removeValue(forKey: invalid.molecule.atoms[0].id)
        XCTAssertEqual(CDKQueryMatcher.match(invalid, in: target).reason, .invalidInput)
        XCTAssertEqual(try verdict("*", "*").reason, .unsupported)
    }

    func testImportedAtomListsAndRestHConstraints() throws {
        var molecule = try CDKSmilesParser(generateCoordinates: false).parseSmiles("CO")
        molecule.atoms[1].atomList = ["O", "N"]
        var query = try CDKQueryPattern(molecule: molecule)
        let target = try CDKSmilesParser(generateCoordinates: false).parseSmiles("CNC")
        XCTAssertEqual(CDKQueryMatcher.match(query, in: target).verdict, .match)
        molecule.atoms[1].rGroupRestH = true
        query = try CDKQueryPattern(molecule: molecule)
        XCTAssertEqual(CDKQueryMatcher.match(query, in: target).verdict, .noMatch)
        XCTAssertEqual(
            CDKQueryMatcher.match(query, in: try CDKSmilesParser(generateCoordinates: false).parseSmiles("CN")).verdict,
            .match)
        let kekule = try CDKQueryPattern(
            molecule: CDKSmilesParser(generateCoordinates: false).parseSmiles("C1=CC=CC=C1"))
        XCTAssertEqual(
            CDKQueryMatcher.match(kekule, in: try CDKSmilesParser(generateCoordinates: false).parseSmiles("c1ccccc1"))
                .verdict, .match)
    }

    func testCancellationNeverBecomesNoMatch() async throws {
        let query = try CDKSmartsParser.parse("CC")
        let target = try CDKSmilesParser(generateCoordinates: false).parseSmiles("CCC")
        let task = Task.detached {
            while !Task.isCancelled { await Task.yield() }
            return CDKQueryMatcher.match(query, in: target)
        }
        task.cancel()
        let result = await task.value
        XCTAssertEqual(result.verdict, .unresolved); XCTAssertEqual(result.reason, .cancelled)
        XCTAssertFalse(result.isExhaustive)
    }
}
