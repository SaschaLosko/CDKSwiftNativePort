import Foundation
import XCTest
@testable import CDKSwiftNativePort

/// Small, reproducible chemistry laws complement the attributed reference corpus.
/// Identity comparisons use Standard InChI; serialized spelling is not the oracle.
final class ChemistryInvariantTests: XCTestCase {
    private let fixtures = [
        "CCO", "CC(=O)[O-]", "[13CH3]CO", "[Na+].[Cl-]", "c1ccccc1O",
        "F[C@H](Cl)Br", "F[C@@H](Cl)Br", "F/C=C/Cl", "F/C=C\\Cl",
    ]

    func testIdentityAndFingerprintSurviveSeededGraphRenumbering() throws {
        for smiles in fixtures {
            let source = try parse(smiles)
            let expected = try identity(source)
            let fingerprinter = CDKCircularFingerprinter(fingerprintClass: .ecfp6, foldedSize: 2048)
            let bits = fingerprinter.bitFingerprint(for: source).sortedBitIndices
            for seed in UInt64(1)...8 {
                let permuted = try renumber(source, seed: seed)
                XCTAssertEqual(try identity(permuted), expected, "\(smiles), seed \(seed)")
                XCTAssertEqual(
                    fingerprinter.bitFingerprint(for: permuted).sortedBitIndices, bits,
                    "\(smiles), seed \(seed)")
            }
        }
    }

    func testIsomericSmilesRoundTripPreservesSpecifiedChemistry() throws {
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor: [
            .useAromaticSymbols, .isomeric, .strict,
        ])
        for smiles in fixtures {
            let source = try parse(smiles)
            let restored = try parse(generator.create(source))
            XCTAssertEqual(try identity(restored), try identity(source), smiles)
            XCTAssertEqual(
                restored.atoms.reduce(0) { $0 + $1.charge }, source.atoms.reduce(0) { $0 + $1.charge }, smiles)
            XCTAssertEqual(
                restored.atoms.compactMap(\.isotopeMassNumber).sorted(),
                source.atoms.compactMap(\.isotopeMassNumber).sorted(), smiles)
        }
    }

    func testPersistenceRoundTripPreservesStereoAndCoordinateProvenance() throws {
        for smiles in fixtures {
            let source = try parse(smiles)
            let restored = try JSONDecoder().decode(Molecule.self, from: JSONEncoder().encode(source))
            XCTAssertEqual(restored, source, smiles)
            XCTAssertEqual(try identity(restored), try identity(source), smiles)
            XCTAssertEqual(restored.coordinatesAreGenerated, source.coordinatesAreGenerated, smiles)
        }
    }

    func testSpecifiedStereoisomersStayDistinct() throws {
        for (left, right) in [("F[C@H](Cl)Br", "F[C@@H](Cl)Br"), ("F/C=C/Cl", "F/C=C\\Cl")] {
            XCTAssertNotEqual(try identity(parse(left)), try identity(parse(right)))
        }
    }

    func testDirectionalStereoHasStableExplicitReferencesWithoutCoordinates() throws {
        let parser = CDKSmilesParser(generateCoordinates: false)
        for (smiles, stereo) in [("F/C=C/Cl", DoubleBondStereo.trans), ("F/C=C\\Cl", .cis)] {
            let molecule = try parser.parseSmiles(smiles)
            let bond = try XCTUnwrap(molecule.bonds.first { $0.order == .double })
            XCTAssertEqual(bond.doubleBondStereo, stereo)
            XCTAssertEqual(bond.stereoReferenceAtomIDs?.count, 4)
            XCTAssertTrue(try identity(molecule).contains("/b"), smiles)
        }
        XCTAssertThrowsError(try parser.parseSmiles("F/C(\\Cl)=C/F"))
    }

    private func parse(_ smiles: String) throws -> Molecule {
        try CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault).parseSmiles(smiles)
    }

    private func identity(_ molecule: Molecule) throws -> String {
        try CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule).getInchi()
    }

    private func renumber(_ source: Molecule, seed: UInt64) throws -> Molecule {
        var generator = InvariantRandom(seed: seed)
        let ids = source.atoms.map(\.id).shuffled(using: &generator)
        let mapping = Dictionary(uniqueKeysWithValues: ids.enumerated().map { ($0.element, 101 + $0.offset * 7) })
        // Change identity fields through Codable to retain every other atom/bond
        // property. Stereo references follow the same bijection, never array order.
        var object = try XCTUnwrap(JSONSerialization.jsonObject(with: JSONEncoder().encode(source)) as? [String: Any])
        var atoms = try XCTUnwrap(object["atoms"] as? [[String: Any]])
        for i in atoms.indices {
            let old = try XCTUnwrap(atoms[i]["id"] as? Int)
            atoms[i]["id"] = try XCTUnwrap(mapping[old])
            if let order = atoms[i]["ligandOrderingAtomIDs"] as? [Int] {
                atoms[i]["ligandOrderingAtomIDs"] = order.map { mapping[$0] ?? $0 }
            }
        }
        var bonds = try XCTUnwrap(object["bonds"] as? [[String: Any]])
        for i in bonds.indices {
            bonds[i]["id"] = 503 + i * 11
            for key in ["a1", "a2", "coordinateBondReferenceAtomID"] {
                if let old = bonds[i][key] as? Int { bonds[i][key] = try XCTUnwrap(mapping[old]) }
            }
            if let order = bonds[i]["stereoReferenceAtomIDs"] as? [Int] {
                bonds[i]["stereoReferenceAtomIDs"] = order.map { mapping[$0] ?? $0 }
            }
        }
        atoms.shuffle(using: &generator)
        bonds.shuffle(using: &generator)
        object["atoms"] = atoms
        object["bonds"] = bonds
        return try JSONDecoder().decode(Molecule.self, from: JSONSerialization.data(withJSONObject: object))
    }
}

private struct InvariantRandom: RandomNumberGenerator {
    var seed: UInt64
    mutating func next() -> UInt64 {
        seed = seed &* 6364136223846793005 &+ 1442695040888963407
        return seed
    }
}
