import XCTest
@testable import CDKSwiftNativePort

/// Independent identities and ordered-graph assertions complement the CDK
/// 2.14 #1340/#1370 examples; cached source text is never the oracle.
final class CDK214StereoTests: XCTestCase {
    private func identity(_ molecule: Molecule) throws -> String {
        try CDKInChIOfficialLibraryGenerator.generate(for: molecule).inchi
    }

    func testInChIDefinedAndUnknownAlkenesStayDistinctWithoutSourceCache() throws {
        var identities: [String] = []
        for sign in ["+", "-", "?"] {
            let source = "InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3\(sign)"
            let reader = CDKInChIToStructure(inchi: source)
            var molecule = try reader.getAtomContainer()
            XCTAssertEqual(reader.getStatus(), .success)
            molecule.dataFields = [:]
            let double = try XCTUnwrap(molecule.bonds.first { $0.order == .double })
            if sign == "?" {
                XCTAssertNil(double.doubleBondStereo)
                XCTAssertFalse(try identity(molecule).contains("/b"))
            } else {
                XCTAssertEqual(double.doubleBondStereo, sign == "+" ? .trans : .cis)
                XCTAssertEqual(double.stereoReferenceAtomIDs?.count, 4)
                XCTAssertEqual(try identity(molecule), source)
            }
            identities.append(try identity(molecule))
        }
        XCTAssertEqual(Set(identities).count, 3)
    }

    func testUnknownAlkeneRetainsTetrahedralCenter() throws {
        let source = "InChI=1S/C5H10O/c1-3-4-5(2)6/h3-6H,1-2H3/b4-3?/t5-/m1/s1"
        let molecule = try CDKInChIToStructure(inchi: source).getAtomContainer()
        XCTAssertNil(molecule.bonds.first { $0.order == .double }?.doubleBondStereo)
        XCTAssertTrue(try identity(molecule).contains("/t5-/m1/s1"))
        XCTAssertFalse(try identity(molecule).contains("/b"))
    }

    func testConjugatedAlkeneUsesReferenceConnectivityAndHydrogens() throws {
        let source = "InChI=1S/C6H8/c1-3-5-6-4-2/h3-6H,1-2H2/b6-5+"
        let molecule = try CDKInChIToStructure(inchi: source).getAtomContainer()
        XCTAssertEqual(molecule.bonds.filter { $0.order == .double }.count, 3)
        XCTAssertEqual(try identity(molecule), source)
    }

    func testChangingOrderedStereoInvalidatesIdentifierCache() throws {
        let source = "InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+"
        var molecule = try CDKInChIToStructure(inchi: source).getAtomContainer()
        let index = try XCTUnwrap(molecule.bonds.firstIndex { $0.doubleBondStereo != nil })
        molecule.bonds[index].doubleBondStereo = .cis
        let output = try CDKInChIGenerator(molecule: molecule).getInchi()
        XCTAssertNotEqual(output, source)
        XCTAssertTrue(output.hasSuffix("/b4-3-"))
    }

    func testCTfileIdentityMatrixWithZeroCoordinatesAndReorderedArrays() throws {
        let inputs = [
            "N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O", "[H][C@](F)(Cl)Br",
            "F[C@](Cl)(Br)I", "F/C=C/F", "F/C=C\\F", "FC=CF", "C/C=C/C=C\\C", "F/C=C/[C@H](Br)Cl", "FC=C[C@H](Br)Cl",
        ]
        for input in inputs {
            var source = try CDKSmilesParser().parseSmiles(input)
            source.atoms.reverse(); source.bonds.reverse()
            for index in source.atoms.indices { source.atoms[index].position = .zero }
            for v3000 in [false, true] {
                let text = try v3000 ? CDKMDLV3000Writer.write(source) : CDKMDLV2000Writer.write(source)
                let copy = try CDKMDLV2000Reader.read(text: text)
                XCTAssertEqual(try identity(copy), try identity(source), "\(input) V3000=\(v3000)")
            }
        }
    }

    func testUndefinedAlkeneKeepsUnknownGeometryAfterMOLAndCML() throws {
        let source = try CDKSmilesParser().parseSmiles("FC=CF")
        let mol = try CDKMDLV2000Reader.read(text: CDKMDLV2000Writer.write(source))
        let cml = try XCTUnwrap(CDKCMLReader.read(text: CDKCMLWriter.write([mol])).first)
        XCTAssertEqual(try identity(source), try identity(mol))
        XCTAssertEqual(try identity(source), try identity(cml))
    }

    func testMalformedAlkeneCannotBeExportedSilently() throws {
        var molecule = try CDKSmilesParser().parseSmiles("F/C=C/F")
        let index = try XCTUnwrap(molecule.bonds.firstIndex { $0.doubleBondStereo != nil })
        molecule.bonds[index].stereoReferenceAtomIDs = [999, 2, 3, 4]
        XCTAssertThrowsError(try CDKMDLV2000Writer.write(molecule))
        XCTAssertThrowsError(try CDKMDLV3000Writer.write(molecule))
    }
}
