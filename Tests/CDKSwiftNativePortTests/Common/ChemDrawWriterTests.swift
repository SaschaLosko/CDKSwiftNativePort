import Foundation
#if canImport(FoundationXML)
    import FoundationXML
#endif
import XCTest
@testable import CDKSwiftNativePort

final class ChemDrawWriterTests: XCTestCase {
    private let parser = CDKSmilesParser()

    func testMoleculeInterchangeFixtures() throws {
        let smiles = [
            "CC(=O)Oc1ccccc1C(=O)O", "[13CH3][NH3+]", "[Na+].[Cl-]", "c1cc[nH]c1",
            "N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O", "F/C=C/F", "F/C=C\\F", "C#N", "[2H]O[2H]", "F[C@](Cl)(Br)I",
            "F[C@@](Cl)(Br)I", "C[C@H](O)[C@@H](C)O", "CC=CC", "C/C=C/C=C/C", "C/C=C\\C=C/C", "[C@@H]1(O)CCC[C@H]1O",
        ]
        for (index, smiles) in smiles.enumerated() {
            var molecule = try parser.parseSmiles(smiles)
            molecule.name = "Müller & <test> \"β\""
            let xml = try CDKFileExporter.writeData(molecule: molecule, as: .cdxml)
            XCTAssertTrue(XMLParser(data: xml).parse())
            let binary = try CDKFileExporter.writeData(molecule: molecule, as: .cdx)
            XCTAssertEqual(binary.prefix(8), Data("VjCD0100".utf8))
            let records = try CDXTestDecoder.decode(binary)
            XCTAssertEqual(records.filter { $0.tag == 0x8004 }.count, molecule.atomCount)
            XCTAssertEqual(records.filter { $0.tag == 0x8005 }.count, molecule.bondCount)
            XCTAssertEqual(Set(records.map(\.id)).count, records.count)
            try evidence("molecule-\(index)", xml: xml, binary: binary, smiles: smiles)
        }
    }

    func testReactionRolesAndAtomMapReferences() throws {
        let reaction = CDKReaction(
            reactants: [try parser.parseSmiles("[CH3:1][OH:2]")],
            agents: [try parser.parseSmiles("O")],
            products: [try parser.parseSmiles("[CH2:1]=[O:2]")])
        let xml = try CDKFileExporter.writeData(reaction: reaction, as: .cdxml)
        XCTAssertTrue(XMLParser(data: xml).parse())
        let binary = try CDKFileExporter.writeData(reaction: reaction, as: .cdx)
        let records = try CDXTestDecoder.decode(binary)
        let step = try XCTUnwrap(records.first { $0.tag == 0x800E })
        let fragments = Set(records.filter { $0.tag == 0x8003 }.map(\.id))
        for tag: UInt16 in [0x0C01, 0x0C02, 0x0C05] {
            let reference = try XCTUnwrap(step.properties[tag])
            XCTAssertEqual(reference.count, 4)
            XCTAssertTrue(fragments.contains(CDXTestDecoder.uint32(reference, 0)))
        }
        XCTAssertEqual(step.properties[0x0C00]?.count, 16)
        try evidence("reaction", xml: xml, binary: binary, smiles: "[CH3:1][OH:2]>O>[CH2:1]=[O:2]")
    }

    func testBinaryURLWriterAndTextAPIGuard() throws {
        let molecule = try parser.parseSmiles("CO")
        XCTAssertThrowsError(try CDKFileExporter.write(molecule: molecule, as: .cdx))
        let url = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString + ".cdx")
        defer { try? FileManager.default.removeItem(at: url) }
        try CDKFileExporter.write(molecule: molecule, to: url)
        XCTAssertEqual(try Data(contentsOf: url), try CDKFileExporter.writeData(molecule: molecule, as: .cdx))
    }

    func testInvalidOrUnsupportedChemistryFailsInsteadOfLosingInformation() throws {
        var molecule = try parser.parseSmiles("CO")
        molecule.atoms.append(molecule.atoms[0])
        XCTAssertThrowsError(try CDKChemDrawWriter.cdx(molecules: [molecule]))
        molecule = try parser.parseSmiles("CO")
        molecule.atoms[0].queryType = .anyAtom
        XCTAssertThrowsError(try CDKChemDrawWriter.cdxml(molecules: [molecule]))
        molecule.atoms[0].queryType = nil
        molecule.atoms[0].position.x = .infinity
        XCTAssertThrowsError(try CDKChemDrawWriter.cdx(molecules: [molecule]))
        XCTAssertThrowsError(try CDKChemDrawWriter.cdx(molecules: []))
    }

    func testExtendedPropertyLengthAndXMLNames() throws {
        var molecule = try parser.parseSmiles("CO")
        molecule.name = String(repeating: "β", count: 40000)
        let records = try CDXTestDecoder.decode(CDKChemDrawWriter.cdx(molecules: [molecule]))
        XCTAssertEqual(records.first { $0.tag == 0x8003 }?.properties[0x0008]?.count, 80012)
        XCTAssertTrue(XMLParser(data: Data(try CDKChemDrawWriter.cdxml(molecules: [molecule]).utf8)).parse())
    }

    private func evidence(_ name: String, xml: Data, binary: Data, smiles: String) throws {
        guard let directory = ProcessInfo.processInfo.environment["CDK_CHEMDRAW_EVIDENCE_DIR"] else { return }
        let root = URL(fileURLWithPath: directory)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: true)
        try xml.write(to: root.appendingPathComponent(name + ".cdxml"))
        try binary.write(to: root.appendingPathComponent(name + ".cdx"))
        try smiles.write(to: root.appendingPathComponent(name + ".smi"), atomically: true, encoding: .utf8)
    }
}

// Bounds-checked, test-only decoder independent of the production serializer.
private enum CDXTestDecoder {
    struct Record { var tag: UInt16; var id: UInt32; var properties: [UInt16: Data] = [:] }
    static func uint32(_ data: Data, _ offset: Int) -> UInt32 {
        (0..<4).reduce(0) { $0 | UInt32(data[data.startIndex + offset + $1]) << (8 * $1) }
    }
    static func decode(_ data: Data) throws -> [Record] {
        var offset = 28
        var records: [Record] = [Record(tag: 0x8000, id: 0)]
        var stack: [Int] = [0]
        func read(_ count: Int) throws -> Data {
            guard count >= 0, offset + count <= data.count else { throw ChemError.emptyInput }
            defer { offset += count }
            return data.subdata(in: offset..<(offset + count))
        }
        func short() throws -> UInt16 {
            let bytes = try read(2)
            return UInt16(bytes[0]) | UInt16(bytes[1]) << 8
        }
        while offset < data.count {
            let tag = try short()
            if tag == 0 { if !stack.isEmpty { stack.removeLast() }; continue }
            if tag & 0x8000 != 0 {
                records.append(Record(tag: tag, id: uint32(try read(4), 0)))
                stack.append(records.count - 1)
            } else {
                let length = try short()
                let count = length == .max ? Int(uint32(try read(4), 0)) : Int(length)
                let index = try XCTUnwrap(stack.last)
                records[index].properties[tag] = try read(count)
            }
        }
        XCTAssertTrue(stack.isEmpty)
        return records
    }
}
