import XCTest
@testable import CDKSwiftNativePort

final class RXNReaderPortTests: XCTestCase {
    func testParsesBasicV2000ReactionBlock() throws {
        let text = """
        $RXN
        DemoReaction
          CDKSwiftNativePort

          1  1
        $MOL
        Reactant
          CDK

          1  0  0  0  0  0            999 V2000
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        M  END
        $MOL
        Product
          CDK

          1  0  0  0  0  0            999 V2000
            0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
        M  END
        """

        let molecules = try CDKRXNReader.read(text: text)
        XCTAssertEqual(molecules.count, 2)
        XCTAssertEqual(molecules[0].atoms.first?.element.uppercased(), "C")
        XCTAssertEqual(molecules[1].atoms.first?.element.uppercased(), "O")
    }

    func testRejectsTextWithoutReactionHeader() {
        XCTAssertThrowsError(try CDKRXNReader.read(text: "not an rxn"))
    }
}
