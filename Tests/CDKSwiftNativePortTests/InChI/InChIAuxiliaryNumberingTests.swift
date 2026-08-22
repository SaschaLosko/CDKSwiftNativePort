import XCTest
@testable import CDKSwiftNativePort

final class InChIAuxiliaryNumberingTests: XCTestCase {
    func testCDK213FixedHLayerUsesIndependentBaseComponentCursor() throws {
        let molecule = Molecule(atoms: (1...14).map { id in
            Atom(id: id, element: "C", position: .zero)
        })
        let auxInfo = "AuxInfo=1/1/N:1,2;3,4;5,6;7,8;9,10;11,12;13,14/F:4m;10,9;2m"

        let numbers = try CDKInChIAuxiliaryNumbering.universalSmilesNumbers(from: auxInfo,
                                                                            molecule: molecule)

        XCTAssertEqual(numbers.sorted(), Array(1...14))
        XCTAssertLessThan(numbers[9], numbers[8], "The explicit /F: order must override the /N: order.")
        XCTAssertEqual(numbers[10...13], [11, 12, 13, 14])
    }
}
