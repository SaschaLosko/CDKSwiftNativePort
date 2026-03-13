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
