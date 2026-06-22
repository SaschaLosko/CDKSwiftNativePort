import Foundation
#if canImport(CoreGraphics)
    import CoreGraphics
#endif
import XCTest
@testable import CDKSwiftNativePort

final class CircularFingerprinterTests: XCTestCase {
    private let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    func testGetBitFingerprintMatchesCDKReferenceForTrivialAmide() throws {
        let molecule = try parser.parseSmiles("CCC(=O)N")
        let fingerprint = CDKCircularFingerprinter().bitFingerprint(for: molecule)

        XCTAssertEqual(fingerprint.size, 1024)
        XCTAssertEqual(
            fingerprint.sortedBitIndices,
            [19, 152, 293, 340, 439, 480, 507, 726, 762, 947, 993]
        )
    }

    func testGetCountFingerprintMatchesCDKReferenceForTrivialAmide() throws {
        let molecule = try parser.parseSmiles("CCC(=O)N")
        let fingerprint = CDKCircularFingerprinter().countFingerprint(for: molecule)

        let expected: [Int32: Int] = [
            -414_937_772: 1,
            -1_027_418_143: 1,
            1_627_608_083: 1,
            -868_007_456: 1,
            -1_006_701_866: 1,
            -1_059_145_289: 1,
            -801_752_141: 1,
            790_592_664: 1,
            -289_109_509: 1,
            -1_650_154_758: 1,
            1_286_833_445: 1,
        ]

        XCTAssertEqual(fingerprint.countsByHash, expected)
        XCTAssertEqual(fingerprint.populatedBinCount, expected.count)
    }

    func testAromaticBenzeneMatchesKekuleBenzeneFingerprint() throws {
        let aromaticBenzene = try parser.parseSmiles("c1ccccc1")
        let kekuleBenzene = try parser.parseSmiles("C1=CC=CC=C1")
        let fingerprinter = CDKCircularFingerprinter()

        XCTAssertEqual(
            fingerprinter.bitFingerprint(for: aromaticBenzene),
            fingerprinter.bitFingerprint(for: kekuleBenzene))
        XCTAssertEqual(
            fingerprinter.countFingerprint(for: aromaticBenzene),
            fingerprinter.countFingerprint(for: kekuleBenzene))
    }

    func testCircularFingerprintKeepsAtomEnvironmentGroups() throws {
        let molecule = try parser.parseSmiles("CC")
        let features = CDKCircularFingerprinter().calculate(molecule)
        let atomGroups = Set(features.map { $0.atomIDs })

        XCTAssertEqual(atomGroups, Set([[1], [2], [1, 2]]))
        XCTAssertEqual(features.filter { $0.atomIDs.count == 1 }.count, 2)
    }

    func testFunctionalFingerprintHandlesHydrogenOnlyMolecule() throws {
        let proton = try parser.parseSmiles("[H+]")
        let fingerprint = CDKCircularFingerprinter(fingerprintClass: .fcfp2).bitFingerprint(for: proton)

        XCTAssertEqual(fingerprint.cardinality, 0)
    }

    func testFunctionalFingerprintHandlesImineRingRegression() {
        let molecule = Molecule(
            name: "Pyrazole-like regression",
            atoms: [
                Atom(id: 1, element: "H", position: .zero),
                Atom(id: 2, element: "N", position: .zero),
                Atom(id: 3, element: "C", position: .zero),
                Atom(id: 4, element: "C", position: .zero),
                Atom(id: 5, element: "C", position: .zero),
                Atom(id: 6, element: "N", position: .zero),
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single),
                Bond(id: 2, a1: 2, a2: 3, order: .single),
                Bond(id: 3, a1: 3, a2: 4, order: .double),
                Bond(id: 4, a1: 4, a2: 5, order: .single),
                Bond(id: 5, a1: 5, a2: 6, order: .double),
                Bond(id: 6, a1: 2, a2: 6, order: .single),
            ]
        )

        let fingerprint = CDKCircularFingerprinter(fingerprintClass: .fcfp2).bitFingerprint(for: molecule)
        XCTAssertGreaterThanOrEqual(fingerprint.cardinality, 0)
    }
}
