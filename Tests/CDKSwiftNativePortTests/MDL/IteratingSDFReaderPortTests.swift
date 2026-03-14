import XCTest
@testable import CDKSwiftNativePort

final class IteratingSDFReaderPortTests: XCTestCase {

    func testReadsMultipleRecordsAndSkipsInvalidBlocks() throws {
        let sdf = validMol(name: "First") + "\n$$$$\n" + "not a molfile\n$$$$\n" + validMol(name: "Second") + "\n$$$$\n"

        let molecules = try CDKIteratingSDFReader.read(text: sdf)

        XCTAssertEqual(molecules.count, 2)
        XCTAssertEqual(molecules.map(\.name), ["First", "Second"])
    }

    func testReadsSingleMolfileWithoutRecordDelimiter() throws {
        let molecules = try CDKIteratingSDFReader.read(text: validMol(name: "Single"))

        XCTAssertEqual(molecules.count, 1)
        XCTAssertEqual(molecules[0].name, "Single")
        XCTAssertEqual(molecules[0].atomCount, 2)
    }

    func testPreservesSDFDataFieldsForValidRecords() throws {
        let sdf = (
            annotatedMol(name: "First",
                         fields: [
                            ("ID", ["001"]),
                            ("Tags", ["alpha", "beta"])
                         ])
            + "\n$$$$\n"
            + "not a molfile\n$$$$\n"
            + annotatedMol(name: "Second",
                           fields: [
                            ("Tags", ["gamma"]),
                            ("Quoted", ["he said \"hi\""])
                           ])
            + "\n$$$$\n"
        )

        let molecules = try CDKIteratingSDFReader.read(text: sdf)

        XCTAssertEqual(molecules.count, 2)
        XCTAssertEqual(molecules[0].orderedDataFieldNames, ["ID", "Tags"])
        XCTAssertEqual(molecules[0].dataFieldValues(named: "Tags"), ["alpha", "beta"])
        XCTAssertEqual(molecules[1].orderedDataFieldNames, ["Tags", "Quoted"])
        XCTAssertEqual(molecules[1].dataFieldValues(named: "Quoted"), ["he said \"hi\""])
    }

    func testReadsTitleAndDataItems() throws {
        let sdf = annotatedMol(name: "2-methylbenzo-1,4-quinone",
                               fields: [("E_NSC", ["1"]), ("E_CAS", ["553-97-9"])]) + "\n$$$$\n"

        let molecules = try CDKIteratingSDFReader.read(text: sdf)

        XCTAssertEqual(molecules.count, 1)
        XCTAssertEqual(molecules[0].name, "2-methylbenzo-1,4-quinone")
        XCTAssertEqual(molecules[0].dataFieldValues(named: "E_NSC"), ["1"])
        XCTAssertEqual(molecules[0].dataFieldValues(named: "E_CAS"), ["553-97-9"])
    }

    func testReadsFieldValuesStartingWithGreaterThan() throws {
        let sdf = annotatedMol(name: "Tagged", fields: [("IC50_uM", [">1"])]) + "\n$$$$\n"

        let molecules = try CDKIteratingSDFReader.read(text: sdf)

        XCTAssertEqual(molecules.count, 1)
        XCTAssertEqual(molecules[0].dataFieldValues(named: "IC50_uM"), [">1"])
    }

    func testReadsV3000RecordWithinSDF() throws {
        let sdf = """
V3000Record
CDKSwiftNativePort

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 O 1.2 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
> <chembl_id>
CHEMBL1

$$$$
"""

        let molecules = try CDKIteratingSDFReader.read(text: sdf)

        XCTAssertEqual(molecules.count, 1)
        XCTAssertEqual(molecules[0].atomCount, 2)
        XCTAssertEqual(molecules[0].dataFieldValues(named: "chembl_id"), ["CHEMBL1"])
    }

    func testSkipOptionCanStopOnBrokenRecord() throws {
        let sdf = validMol(name: "First") + "\n$$$$\nnot a molfile\n$$$$\n" + validMol(name: "Second") + "\n$$$$\n"

        XCTAssertThrowsError(try CDKIteratingSDFReader.read(text: sdf,
                                                            options: CDKIteratingSDFReaderOptions(skipInvalidRecords: false)))
    }

    func testReadsEmptyStructureEntryAndPreservesFields() throws {
        let sdf = """
First
CDKSwiftNativePort

  1  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <Species>
human

$$$$
Empty
CDKSwiftNativePort

  0  0  0  0  0  0  0  0  0  0  0 V2000
M  END
> <Species>
rat

$$$$
"""

        let molecules = try CDKIteratingSDFReader.read(text: sdf)

        XCTAssertEqual(molecules.count, 2)
        XCTAssertEqual(molecules[1].atomCount, 0)
        XCTAssertEqual(molecules[1].dataFieldValues(named: "Species"), ["rat"])
    }

    func testAcceptsExtraSpacesAroundDelimitersAndHeaders() throws {
        let sdf = """
Spaced
CDKSwiftNativePort

  2  1  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
   >   <chembl_id>
CHEMBL564829

  $$$$  
"""

        let molecules = try CDKIteratingSDFReader.read(text: sdf)

        XCTAssertEqual(molecules.count, 1)
        XCTAssertEqual(molecules[0].dataFieldValues(named: "chembl_id"), ["CHEMBL564829"])
    }

    func testReadsFromFile() throws {
        let url = FileManager.default.temporaryDirectory
            .appendingPathComponent("cdk_iterating_sdf_\(UUID().uuidString).sdf")
        defer { try? FileManager.default.removeItem(at: url) }

        try validMol(name: "FromFile").write(to: url, atomically: true, encoding: .utf8)

        let molecules = try CDKIteratingSDFReader.readFile(url: url)
        XCTAssertEqual(molecules.count, 1)
        XCTAssertEqual(molecules[0].name, "FromFile")
    }

    func testThrowsOnEmptyInput() {
        XCTAssertThrowsError(try CDKIteratingSDFReader.read(text: "   \n\n"))
    }

    func testThrowsWhenNoValidRecordsExist() {
        XCTAssertThrowsError(try CDKIteratingSDFReader.read(text: "junk\n$$$$\nmore junk\n"))
    }

    private func validMol(name: String) -> String {
        """
\(name)
CDKSwiftNativePort

  2  1  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
"""
    }

    private func annotatedMol(name: String, fields: [(String, [String])]) -> String {
        var lines = validMol(name: name).components(separatedBy: "\n")

        if lines.last?.isEmpty == true {
            lines.removeLast()
        }

        for (fieldName, values) in fields {
            lines.append(">  <\(fieldName)>")
            lines.append(contentsOf: values)
            lines.append("")
        }

        return lines.joined(separator: "\n")
    }
}
