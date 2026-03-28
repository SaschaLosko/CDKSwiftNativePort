import XCTest
@testable import CDKSwiftNativePort

final class RGFileReaderWriterTests: XCTestCase {
    func testRGroupListOccurrenceValidationMatchesCDKRules() throws {
        XCTAssertTrue(CDKRGroupList.isValidOccurrenceSyntax("0,1-3"))
        XCTAssertTrue(CDKRGroupList.isValidOccurrenceSyntax("0,2"))
        XCTAssertTrue(CDKRGroupList.isValidOccurrenceSyntax(">0"))
        XCTAssertTrue(CDKRGroupList.isValidOccurrenceSyntax("<2"))
        XCTAssertFalse(CDKRGroupList.isValidOccurrenceSyntax("<0"))
        XCTAssertFalse(CDKRGroupList.isValidOccurrenceSyntax("3-1"))
        XCTAssertFalse(CDKRGroupList.isValidOccurrenceSyntax("abc"))

        let list = try CDKRGroupList(rGroupNumber: 1, occurrence: "0,2")
        XCTAssertEqual(list.matchOccurrence(maxAttachments: 3), [0, 2])
    }

    func testManipulatorRoundTripsFlattenedMarkushWithLogic() throws {
        let flattened = try CDKRGFileReader.readFlattenedMolecule(text: RGFileFixtures.simpleQuery)

        let query = try CDKRGroupQueryManipulator.toRGroupQuery(flattened)
        let rebuilt = CDKRGroupQueryManipulator.toFlatMolecule(query)

        XCTAssertEqual(query.rootStructure.atomCount, 7)
        XCTAssertEqual(query.rGroupDefinitions[1]?.occurrence, "0,1-3")
        XCTAssertEqual(query.rGroupDefinitions[1]?.restH, true)
        XCTAssertEqual(query.rGroupDefinitions[1]?.rGroups.count, 3)
        XCTAssertEqual(rebuilt.rGroupLogicDefinitions[1]?.occurrence, "0,1-3")
        XCTAssertEqual(rebuilt.rGroupLogicDefinitions[1]?.restH, true)
        XCTAssertEqual(rebuilt.atoms.filter { $0.rGroupMembership == "R1" }.count,
                       flattened.atoms.filter { $0.rGroupMembership == "R1" }.count)
    }

    func testRGFileReaderParsesUpstreamSimpleFixture() throws {
        let query = try CDKRGFileReader.read(text: RGFileFixtures.simpleQuery)

        XCTAssertEqual(query.rGroupDefinitions.count, 1)
        XCTAssertEqual(query.rootStructure.atomCount, 7)
        XCTAssertEqual(query.rGroupQueryAtoms(1).count, 1)

        let definition = try XCTUnwrap(query.rGroupDefinitions[1])
        XCTAssertEqual(definition.occurrence, "0,1-3")
        XCTAssertEqual(definition.requiredRGroupNumber, 0)
        XCTAssertTrue(definition.restH)
        XCTAssertEqual(definition.rGroups.count, 3)

        XCTAssertEqual(attachmentSymbol(for: definition.rGroups[0], using: \.firstAttachmentPointAtomID), "N")
        XCTAssertEqual(attachmentSymbol(for: definition.rGroups[1], using: \.firstAttachmentPointAtomID), "O")
        XCTAssertEqual(attachmentSymbol(for: definition.rGroups[2], using: \.firstAttachmentPointAtomID), "S")
        XCTAssertTrue(definition.rGroups.allSatisfy { $0.secondAttachmentPointAtomID == nil })
    }

    func testRGFileReaderParsesUpstreamElaborateFixture() throws {
        let query = try CDKRGFileReader.read(text: RGFileFixtures.elaborateQuery)

        XCTAssertEqual(query.rGroupDefinitions.count, 3)
        XCTAssertEqual(query.rootStructure.atomCount, 14)
        XCTAssertEqual(query.rGroupQueryAtoms().count, 4)
        XCTAssertEqual(query.rGroupQueryAtoms(1).count, 1)

        let rootR1 = try XCTUnwrap(query.rGroupQueryAtoms(1).first)
        let rootR1Order = try XCTUnwrap(rootR1.ligandOrderingAtomIDs)
        XCTAssertEqual(rootR1Order.count, 2)
        XCTAssertEqual(query.rootStructure.atom(id: rootR1Order[0])?.symbolToDraw, "N")
        XCTAssertEqual(query.rootStructure.atom(id: rootR1Order[1])?.symbolToDraw, "C")

        let r1 = try XCTUnwrap(query.rGroupDefinitions[1])
        XCTAssertEqual(r1.rGroups.count, 2)
        XCTAssertFalse(r1.restH)
        XCTAssertEqual(r1.rGroups.map { attachmentSymbol(for: $0, using: \.secondAttachmentPointAtomID) }, ["O", "O"])

        let r2 = try XCTUnwrap(query.rGroupDefinitions[2])
        XCTAssertEqual(r2.rGroups.count, 2)
        XCTAssertEqual(r2.occurrence, "0,2")
        XCTAssertEqual(r2.requiredRGroupNumber, 11)
        XCTAssertFalse(r2.restH)

        let r11 = try XCTUnwrap(query.rGroupDefinitions[11])
        XCTAssertEqual(r11.rGroups.count, 1)
        XCTAssertEqual(r11.requiredRGroupNumber, 0)
        XCTAssertTrue(r11.restH)
        XCTAssertEqual(attachmentSymbol(for: r11.rGroups[0], using: \.firstAttachmentPointAtomID), "Pt")
        XCTAssertNil(r11.rGroups[0].secondAttachmentPointAtomID)
    }

    func testRGFileReaderHandlesDetachedRGroupWithoutAttachmentPoints() throws {
        let query = try CDKRGFileReader.read(text: RGFileFixtures.detachedQuery)

        XCTAssertEqual(query.rGroupDefinitions.count, 1)
        XCTAssertEqual(query.rootStructure.atomCount, 6)
        XCTAssertEqual(query.rGroupQueryAtoms().count, 1)
        XCTAssertTrue(query.areSubstituentsDefined())

        let definition = try XCTUnwrap(query.rGroupDefinitions[1])
        XCTAssertEqual(definition.rGroups.count, 2)
        XCTAssertEqual(definition.requiredRGroupNumber, 0)
        XCTAssertFalse(definition.restH)
        XCTAssertTrue(definition.rGroups.allSatisfy { $0.firstAttachmentPointAtomID == nil })
        XCTAssertTrue(definition.rGroups.allSatisfy { $0.secondAttachmentPointAtomID == nil })
    }

    func testRGFileReaderAcceptsMissingRGroupBlocksLikeCDKLenientMode() throws {
        let query = try CDKRGFileReader.read(text: RGFileFixtures.missingDefinitionsQuery)

        XCTAssertEqual(query.rGroupDefinitions.count, 3)
        XCTAssertEqual(query.rootStructure.atomCount, 14)
        XCTAssertFalse(query.areSubstituentsDefined())
        XCTAssertEqual(query.rGroupDefinitions[11]?.rGroups.count, 1)
        XCTAssertEqual(query.rGroupDefinitions[1]?.rGroups.count, 0)
        XCTAssertEqual(query.rGroupDefinitions[2]?.rGroups.count, 0)
    }

    func testRGFileWriterRoundTripsSimpleFixture() throws {
        let query = try CDKRGFileReader.read(text: RGFileFixtures.simpleQuery)

        let text = try CDKRGFileWriter.write(query)
        let reparsed = try CDKRGFileReader.read(text: text)

        XCTAssertTrue(text.contains("$MDL"))
        XCTAssertTrue(text.contains("$RGP"))
        XCTAssertTrue(text.contains("M  LOG  1   1   0   1   0,1-3"))
        XCTAssertEqual(reparsed.rGroupDefinitions[1]?.occurrence, "0,1-3")
        XCTAssertEqual(reparsed.rGroupDefinitions[1]?.restH, true)
        XCTAssertEqual(reparsed.rGroupDefinitions[1]?.rGroups.count, 3)
    }

    func testRGFileWriterRejectsNonQueryMolecule() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("CCO")

        XCTAssertFalse(CDKRGFileWriter.canWrite(molecule))
        XCTAssertThrowsError(try CDKRGFileWriter.write(molecule))
    }

    private func attachmentSymbol(for group: CDKRGroup,
                                  using keyPath: KeyPath<CDKRGroup, Int?>) -> String? {
        guard let atomID = group[keyPath: keyPath] else { return nil }
        return group.group.atom(id: atomID)?.symbolToDraw
    }
}

enum RGFileFixtures {
    static let simpleQuery = """
    $MDL  REV  1 0118101730
    $MOL
    $HDR
      RGroup query unit test: simple query test.
      Marvin  01181017302D

    $END HDR
    $CTAB
      7  7  0  0  0  0            999 V2000
        4.5003    0.0208    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        4.5003   -0.8043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        5.2124   -1.2126    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        5.9244   -0.8043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        5.9244    0.0208    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
        5.2124    0.4375    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        6.9938    1.3152    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
      5  7  1  0  0  0  0
      1  2  1  0  0  0  0
      1  6  1  0  0  0  0
      2  3  1  0  0  0  0
      3  4  1  0  0  0  0
      4  5  1  0  0  0  0
      5  6  1  0  0  0  0
    M  LOG  1   1   0   1   0,1-3
    M  RGP  1   7   1
    M  END
    $END CTAB
    $RGP
      1
    $CTAB
      2  1  0  0  0  0            999 V2000
        3.5545   -5.0670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        4.2689   -4.6544    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
      1  2  1  0  0  0  0
    M  APO  1   2   1
    M  END
    $END CTAB
    $CTAB
      3  2  0  0  0  0            999 V2000
        1.4874   -4.0503    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        2.2019   -3.6378    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
        1.4874   -4.8753    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
      1  2  1  0  0  0  0
      1  3  1  0  0  0  0
    M  APO  1   3   1
    M  END
    $END CTAB
    $CTAB
      3  2  0  0  0  0            999 V2000
       -0.5544   -4.0378    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        0.1600   -3.6253    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       -0.5438   -4.8627    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
      1  2  1  0  0  0  0
      1  3  1  0  0  0  0
    M  APO  1   3   1
    M  END
    $END CTAB
    $END RGP
    $END MOL
    """

    static let elaborateQuery = """
    $MDL  REV  1   0112101424
    $MOL
    $HDR
      RGroup query unit test: more elaborate query, R1 double bound.
      SMMXDraw01121015182D

    $END HDR
    $CTAB
     14 15  0  0  0  0  0  0  0  0999 V2000
       10.6592  -11.8997    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
       14.2529   -9.5247    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
        8.2529  -10.0559    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
       11.7372   -8.5950    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
       12.3850  -10.2192    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
       12.1300  -11.0038    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       11.3050  -11.0038    0.0000 Pb  0  0  3  0  0  0  0  0  0  0  0  0
       11.0501  -10.2192    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
       11.7568   -9.4200    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
       12.4516   -6.9122    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       12.4516   -6.0872    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       11.0227   -6.0872    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       11.0227   -6.9122    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       11.7372   -7.3247    0.0000 N   0  0  3  0  0  0  0  0  0  0  0  0
     14 13  1  0  0  0  0
     14 10  1  0  0  0  0
     13 12  2  0  0  0  0
     12 11  1  0  0  0  0
     11 10  2  0  0  0  0
      9  8  1  0  0  0  0
      9  5  1  0  0  0  0
      8  7  1  0  0  0  0
      7  6  1  0  0  0  0
      6  5  1  0  0  0  0
      9  4  1  0  0  0  0
      4 14  1  0  0  0  0
      8  3  1  0  0  0  0
      5  2  1  0  0  0  0
      7  1  1  0  0  0  0
    M  AAL   4  2  14   1   9   2
    M  RGP  4   1  11   2   2   3   2   4   1
    M  LOG  1  11   0   1   0,>0
    M  LOG  1   2  11   0   0,2
    M  LOG  1   1   0   0   >0
    M  END
    $END CTAB
    $RGP
      11
    $CTAB
      2  1  0  0  0  0  0  0  0  0999 V2000
       11.2380  -21.5306    0.0000 Pt  0  0  0  0  0  0  0  0  0  0  0  0
       11.6522  -20.8151    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
      1  2  1  0  0  0  0
    M  APO  1   1   1
    M  END
    $END CTAB
    $END RGP
    $RGP
       2
    $CTAB
      4  4  0  0  0  0  0  0  0  0999 V2000
       11.7388  -18.7390    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
       12.5797  -18.7390    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
       12.5797  -19.5659    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       11.7388  -19.5659    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
      4  1  1  0  0  0  0
      3  2  1  0  0  0  0
      3  4  1  0  0  0  0
      2  1  1  0  0  0  0
    M  APO  1   1   1
    M  END
    $END CTAB
    $CTAB
      3  3  0  0  0  0  0  0  0  0999 V2000
       15.5109  -18.8699    0.0000 Si  0  0  3  0  0  0  0  0  0  0  0  0
       16.3392  -18.8699    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       15.9251  -18.1544    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
      2  3  1  0  0  0  0
      1  3  1  0  0  0  0
      1  2  1  0  0  0  0
    M  APO  1   1   1
    M  END
    $END CTAB
    $END RGP
    $RGP
       1
    $CTAB
      3  2  0  0  0  0  0  0  0  0999 V2000
       12.3749  -16.4596    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       11.8446  -14.9667    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
       11.8446  -13.7882    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
      1  2  1  0  0  0  0
      2  3  1  0  0  0  0
    M  APO  2   1   1   3   2
    M  END
    $END CTAB
    $CTAB
      3  2  0  0  0  0  0  0  0  0999 V2000
       15.8985  -15.0059    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
       16.7252  -15.0059    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
       15.4851  -14.2899    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
      1  3  1  0  0  0  0
      1  2  1  0  0  0  0
    M  APO  2   1   1   2   2
    M  END
    $END CTAB
    $END RGP
    $END MOL
    """

    static let detachedQuery = """
    $MDL  REV  1   0120101650
    $MOL
    $HDR
      RGroup query unit test: test detached R-group (not bound)
      SMMXDraw01201016502D

    $END HDR
    $CTAB
      6  5  0  0  0  0  0  0  0  0999 V2000
       10.7188   -7.4161    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       11.1343   -8.6951    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       11.3911   -7.9093    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       10.3032   -8.6951    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       10.0463   -7.9093    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        7.1563   -4.7813    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
      5  4  2  0  0  0  0
      4  2  1  0  0  0  0
      3  1  1  0  0  0  0
      2  3  2  0  0  0  0
      1  5  1  0  0  0  0
    M  RGP  1   6   1
    M  LOG  1   1   0   0  >0
    M  END
    $END CTAB
    $RGP
       1
    $CTAB
      2  1  0  0  0  0  0  0  0  0999 V2000
        3.4759  -11.2839    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
        4.3026  -11.2839    0.0000 C   0  5  0  0  0  0  0  0  0  0  0  0
      1  2  3  0  0  0  0
    M  CHG  2   1   1   2  -1
    M  END
    $END CTAB
    $CTAB
      3  2  0  0  0  0  0  0  0  0999 V2000
        6.1589  -11.4906    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
        6.8749  -11.0772    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
        5.4428  -11.0772    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
      1  2  2  0  0  0  0
      1  3  2  0  0  0  0
    M  APO  1   1  0
    M  END
    $END CTAB
    $END RGP
    $END MOL
    """

    static let missingDefinitionsQuery = """
    $MDL  REV  1   0112101424
    $MOL
    $HDR
      RGroup query unit test: missing $RGP blocks (for R1 and R2)
      SMMXDraw01121015182D

    $END HDR
    $CTAB
     14 15  0  0  0  0  0  0  0  0999 V2000
       10.6592  -11.8997    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
       14.2529   -9.5247    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
        8.2529  -10.0559    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
       11.7372   -8.5950    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
       12.3850  -10.2192    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
       12.1300  -11.0038    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       11.3050  -11.0038    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
       11.0501  -10.2192    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
       11.7568   -9.4200    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
       12.4516   -6.9122    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       12.4516   -6.0872    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       11.0227   -6.0872    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       11.0227   -6.9122    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       11.7372   -7.3247    0.0000 N   0  0  3  0  0  0  0  0  0  0  0  0
     14 13  1  0  0  0  0
     14 10  1  0  0  0  0
     13 12  2  0  0  0  0
     12 11  1  0  0  0  0
     11 10  2  0  0  0  0
      9  8  1  0  0  0  0
      9  5  1  0  0  0  0
      8  7  1  0  0  0  0
      7  6  1  0  0  0  0
      6  5  1  0  0  0  0
      9  4  1  0  0  0  0
      4 14  1  0  0  0  0
      8  3  1  0  0  0  0
      5  2  1  0  0  0  0
      7  1  1  0  0  0  0
    M  AAL   4  2  14   1   9   2
    M  RGP  4   1  11   2   2   3   2   4   1
    M  LOG  1  11   0   1   >0
    M  LOG  1   2  11   0   >1
    M  LOG  1   1  11   0   >0
    M  END
    $END CTAB
    $RGP
      11
    $CTAB
      2  1  0  0  0  0  0  0  0  0999 V2000
       11.2380  -21.5306    0.0000 Pt  0  0  0  0  0  0  0  0  0  0  0  0
       11.6522  -20.8151    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
      1  2  1  0  0  0  0
    M  APO  1   1   1
    M  END
    $END CTAB
    $END RGP
    $END MOL
    """
}
