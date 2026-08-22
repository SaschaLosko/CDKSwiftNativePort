import XCTest
@testable import CDKSwiftNativePort

final class MoleculeRingTests: XCTestCase {
    private let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    func testCDK213SmallestRingSizesForFusedThreeAndSixMemberRings() throws {
        let molecule = try parser.parseSmiles("C1CC2OC2CC1")

        XCTAssertEqual(molecule.atoms.map { molecule.smallestRingSize(containingAtomID: $0.id) },
                       [6, 6, 3, 3, 3, 6, 6])
        XCTAssertEqual(molecule.atoms.map {
            molecule.smallestRingSize(containingAtomID: $0.id, maxSize: 6)
        },
        [6, 6, 3, 3, 3, 6, 6])
        XCTAssertEqual(molecule.atoms.map {
            molecule.smallestRingSize(containingAtomID: $0.id, maxSize: 5)
        },
        [nil, nil, 3, 3, 3, nil, nil])
    }
}
