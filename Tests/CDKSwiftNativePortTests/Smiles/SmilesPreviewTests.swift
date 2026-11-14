import XCTest
@testable import CDKSwiftNativePort

final class SmilesPreviewTests: XCTestCase {
    func testUnfinishedPrefixesNeverBecomeCompleteMolecules() throws {
        let cases: [(String, CDKSmilesPreview.Issue, Int, Int)] = [
            ("C0", .ringClosure, 1, 0), ("C1CCC", .ringClosure, 4, 3),
            ("C(C", .branch, 2, 1), ("C(", .branch, 1, 0),
            ("C=", .bond, 1, 0), ("CC.", .fragment, 2, 1),
            ("C[NH", .atom, 1, 0), ("[13", .atom, 0, 0),
            ("C[NH2:", .atom, 1, 0), ("C%1", .ringIndex, 1, 0),
        ]
        for (text, issue, atoms, bonds) in cases {
            let preview = try CDKSmilesPreview.evaluate(text)
            XCTAssertEqual(preview.status, .incomplete, text)
            XCTAssertTrue(preview.issues.contains(issue), text)
            XCTAssertEqual(preview.moleculeForDepiction?.atoms.count ?? 0, atoms, text)
            XCTAssertEqual(preview.moleculeForDepiction?.bonds.count ?? 0, bonds, text)
            XCTAssertThrowsError(try CDKSmilesParser().parseSmiles(text), text)
        }
    }

    func testCompletionReparsesOriginalTextIncludingStereo() throws {
        for text in [
            "CO", "C1CCCCC1", "C(C)O", "[13CH4]", "C%12CCCCC%12", "N[C@@H](C)C(=O)O", "F/C=C/F", "[Na+].[Cl-]",
        ] {
            let result = try CDKSmilesPreview.evaluate(text)
            XCTAssertEqual(result.status, .complete, text)
            XCTAssertTrue(result.issues.isEmpty)
            let strict = try CDKSmilesParser(generateCoordinates: false).parseSmiles(text)
            XCTAssertEqual(result.moleculeForDepiction?.atoms, strict.atoms, text)
            XCTAssertEqual(result.moleculeForDepiction?.bonds, strict.bonds, text)
        }
        XCTAssertEqual(try CDKSmilesPreview.evaluate("CO methanol").status, .complete)
    }

    func testMalformedInteriorIsNotRecoveredAsAUsablePreview() throws {
        for text in ["C)O", "C()", "C..O", "C==O", "C(=)", "C11", "C1C1", "C%X", "C[N!", "C[Q", "Q"] {
            let result = try CDKSmilesPreview.evaluate(text)
            XCTAssertEqual(result.status, .invalid, text)
            XCTAssertNil(result.moleculeForDepiction, text)
        }
    }

    func testStrictParserRejectsIncompleteFragmentAndMalformedBondSequences() {
        for text in ["C.", ".C", "C..O", "C==O", "C()", "C(=)", "C11", "C1C1"] {
            XCTAssertThrowsError(try CDKSmilesParser().parseSmiles(text), text)
        }
    }

    func testPartialStereoIsNotInventedFromMissingNeighbors() throws {
        let result = try CDKSmilesPreview.evaluate("N[C@@H](C)C1")
        XCTAssertEqual(result.status, .incomplete)
        XCTAssertTrue(try XCTUnwrap(result.moleculeForDepiction).atoms.allSatisfy { $0.chirality == .none })
        XCTAssertTrue(try XCTUnwrap(result.moleculeForDepiction).atoms.allSatisfy { $0.ligandOrderingAtomIDs == nil })
    }

    func testBoundsUnsupportedFormatsAndEmptyInputAreExplicit() throws {
        XCTAssertEqual(try CDKSmilesPreview.evaluate("  ").status, .empty)
        for text in ["CC>O>CO", "C* |$;R1$|", "CC\nCO", "InChI=1S/CH4/h1H4", "<cml>"] {
            XCTAssertEqual(try CDKSmilesPreview.evaluate(text).status, .unsupported)
        }
        for text in [String(repeating: "C", count: 97), String(repeating: "C", count: 1025)] {
            let result = try CDKSmilesPreview.evaluate(text)
            XCTAssertEqual(result.status, .limited)
            XCTAssertNil(result.moleculeForDepiction)
        }
    }

    func testCancellationIsNotReportedAsInvalidChemistry() async throws {
        let task = Task {
            withUnsafeCurrentTask { $0?.cancel() }
            return try CDKSmilesPreview.evaluate("C1CCC")
        }
        do { _ = try await task.value; XCTFail("Expected cancellation") } catch is CancellationError {} catch {
            XCTFail("Unexpected error: \(error)")
        }
    }
}
