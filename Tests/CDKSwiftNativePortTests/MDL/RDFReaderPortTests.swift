import XCTest
@testable import CDKSwiftNativePort

final class RDFReaderPortTests: XCTestCase {
    func testExtractsEmbeddedReactionRecords() throws {
        let text = """
        $RDFILE 1
        $DATM      01/01/2026 00:00
        $RFMT
        $RXN
        RDFReaction
          CDKSwiftNativePort

          1  0
        $MOL
        Reactant
          CDK

          1  0  0  0  0  0            999 V2000
            0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
        M  END
        $DTYPE NAME
        $DATUM Example
        """

        let molecules = try CDKRDFReader.read(text: text)
        XCTAssertEqual(molecules.count, 1)
        XCTAssertEqual(molecules[0].atoms.first?.element.uppercased(), "N")
    }

    func testRejectsRdfWithoutReactionBlocks() {
        XCTAssertThrowsError(try CDKRDFReader.read(text: "$RDFILE 1\n"))
    }
}
