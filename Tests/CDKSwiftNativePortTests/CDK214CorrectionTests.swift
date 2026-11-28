import Foundation
import XCTest
@testable import CDKSwiftNativePort

/// Reference cases: cdk/cdk #1347, #1364, #1350, #1369, #1379, #1317, #1322, #1367,
/// #1335 and #1301 (2.14.1).
final class CDK214CorrectionTests: XCTestCase {
    func testTPSAUpstreamCorrectionsAndControls() throws {
        for (smiles, expected) in [
            ("C1CS1", 25.30), ("CP1CC1", 13.59), ("C1CN1", 21.94),
            ("C1CO1", 12.53), ("CP=C", 34.14), ("CP=NC", 46.50), ("CP(C)C", 13.59),
        ] {
            let molecule = try CDKSmilesParser().parseSmiles(smiles)
            XCTAssertEqual(CDKTPSADescriptor.calculate(for: molecule), expected, accuracy: 0.001, smiles)
        }
    }

    func testCTfileMassReferenceAndAbsoluteOverride() throws {
        for (symbol, reference) in [("Cu", 64), ("Ni", 59), ("Zn", 65), ("Sg", 266), ("R", 0)] {
            let text = """
                Mass


                  1  0  0  0  0  0            999 V2000
                    0.0000    0.0000    0.0000 \(symbol.padding(toLength: 3, withPad: " ", startingAt: 0)) 1  0  0  0  0  0  0  0  0  0  0  0
                M  END
                """
            let molecule = try CDKMDLV2000Reader.read(text: text)
            XCTAssertEqual(molecule.atoms[0].isotopeMassNumber, reference + 1, symbol)
            let explicit = text.replacingOccurrences(of: "M  END", with: "M  ISO  1   1  70\nM  END")
            XCTAssertEqual(try CDKMDLV2000Reader.read(text: explicit).atoms[0].isotopeMassNumber, 70)
        }
        for (mass, delta) in [(61, -3), (64, 0), (68, 4), (69, 0)] {
            let molecule = Molecule(
                name: "Cu", atoms: [Atom(id: 1, element: "Cu", position: .zero, isotopeMassNumber: mass)], bonds: [])
            let text = try CDKMDLV2000Writer.write(molecule)
            let line = text.components(separatedBy: "\n")[4]
            XCTAssertEqual(Int(line.dropFirst(34).prefix(2).trimmingCharacters(in: .whitespaces)), delta)
            XCTAssertEqual(try CDKMDLV2000Reader.read(text: text).atoms[0].isotopeMassNumber, mass)
        }
    }

    func testReactionFlagsBecomeParticipantStereoAndSurviveExport() throws {
        let parser = CDKSmilesParser()
        let writer = CDKSmilesGenerator(flavor: [.isomeric, .cxsmiles, .cxEnhancedStereo])
        for flag in ["r", "r:0", "&1:1,6"] {
            let reaction = try parser.parseReactionSmiles("C[C@H](O)CN>>C[C@H](O)CN(=O)=O |\(flag)|")
            let first = reaction.reactants[0].atoms.first { $0.chirality != .none }!
            let second = reaction.products[0].atoms.first { $0.chirality != .none }!
            XCTAssertEqual(first.cxStereoGroup, CDKCxSmilesParser.encodeStereoGroup(kind: "and", number: 1))
            XCTAssertEqual(second.cxStereoGroup != nil, flag != "r:0")
            let copy = try parser.parseReactionSmiles(writer.create(reaction))
            XCTAssertEqual(copy.reactants[0].atoms.first { $0.chirality != .none }?.cxStereoGroup, first.cxStereoGroup)
            if flag == "r:0" {
                XCTAssertEqual(
                    CDKCxSmilesParser.decodeStereoGroup(
                        copy.products[0].atoms.first { $0.chirality != .none }!.cxStereoGroup!)?.kind, "abs")
            } else {
                XCTAssertEqual(
                    copy.products[0].atoms.first { $0.chirality != .none }?.cxStereoGroup, second.cxStereoGroup)
            }
        }
    }

    func testNestedAbbreviationsExpansionHighlightAndNonmutation() throws {
        var molecule = try CDKSmilesParser().parseSmiles("c1ccccc1OP(=O)(Oc1ccccc1)Oc1ccccc1")
        func group(_ count: Int, _ label: String) -> MoleculeSgroup {
            let ids = Set(molecule.atoms.prefix(count).map(\.id))
            return MoleculeSgroup(
                kind: .generic, keyword: "SUP", atomIDs: ids.sorted(),
                crossingBondIDs: molecule.bonds.filter { ids.contains($0.a1) != ids.contains($0.a2) }.map(\.id),
                subscriptText: label)
        }
        let inner = group(6, "Ph")
        var outer = group(7, "OPh")
        for groups in [[inner, outer], [outer, inner]] {
            molecule.sgroups = groups
            let original = molecule
            for _ in 0..<2 {
                let rendered = CDKDepictionPreprocessor.prepareForRendering(molecule: molecule, style: RenderStyle())
                XCTAssertEqual(rendered.atoms.compactMap(\.aliasLabel), ["OPh"])
            }
            XCTAssertEqual(molecule, original)
        }
        outer.expanded = true
        molecule.sgroups = [outer, inner]
        XCTAssertEqual(
            CDKDepictionPreprocessor.prepareForRendering(molecule: molecule, style: RenderStyle()).atoms.compactMap(
                \.aliasLabel), ["Ph"])
        molecule.highlightedAtomIDs = [molecule.atoms[0].id]
        XCTAssertTrue(
            CDKDepictionPreprocessor.prepareForRendering(molecule: molecule, style: RenderStyle()).atoms.compactMap(
                \.aliasLabel
            ).isEmpty)
    }

    func testDisconnectedAbbreviationRemainsVisible() throws {
        var molecule = try CDKSmilesParser().parseSmiles("NN")
        molecule.sgroups = [
            MoleculeSgroup(kind: .generic, keyword: "SUP", atomIDs: molecule.atoms.map(\.id), subscriptText: "N2H4")
        ]
        let rendered = CDKDepictionPreprocessor.prepareForRendering(molecule: molecule, style: RenderStyle())
        XCTAssertEqual(rendered.atoms.compactMap(\.aliasLabel), ["N2H4"])
        XCTAssertTrue(rendered.bonds.isEmpty)
    }
    func testSMARTSRingValueZeroMatchesChainAtoms() throws {
        let parser = CDKSmilesParser(flavor: .cdkDefault, generateCoordinates: false)
        func verdict(_ smarts: String, _ smiles: String) throws -> CDKMarkushMatchResult.Verdict {
            CDKQueryMatcher.match(try CDKSmartsParser.parse(smarts), in: try parser.parseSmiles(smiles)).verdict
        }
        for (smarts, smiles, expected) in [
            ("[r0]", "CC", CDKMarkushMatchResult.Verdict.match),
            ("[r0;R]", "C1CC1", .noMatch),
            ("[R0]", "CC", .match),
            ("[x0]", "CC", .match),
            ("[r0]", "C1CCCCC1", .noMatch),
            ("[r0]C1CC1", "CC1CC1", .match),
            ("[!r0]", "C1CC1", .match),
            ("[!r0]", "CC", .noMatch),
            // Spiro[2.5]octane: no five-membered ring, the spiro atom is in a three-ring.
            ("[r5]", "C1CC11CCCCC1", .noMatch),
            ("[r3]", "C1CC11CCCCC1", .match),
            ("[R2;r3]", "C1C2CCC12", .match),
            ("[R2;r4]", "C1C2CCC12", .noMatch),
        ] {
            XCTAssertEqual(try verdict(smarts, smiles), expected, "\(smarts) in \(smiles)")
        }
    }

    func testPDBElementSymbolsAreRightJustifiedInColumns77And78() throws {
        let molecule = try CDKSmilesParser().parseSmiles("ClC")
        let atomLines = try CDKPDBWriter.write([molecule]).split(separator: "\n").filter { $0.hasPrefix("ATOM") }
        let elements = atomLines.map { line -> String in
            let characters = Array(line)
            XCTAssertEqual(characters.count, 78, String(line))
            return String(characters[76...77])
        }
        XCTAssertEqual(elements, ["CL", " C"])
        let parsed = try XCTUnwrap(CDKPDBReader.read(text: try CDKPDBWriter.write([molecule])).first)
        XCTAssertEqual(parsed.atoms.map(\.element).sorted(), ["C", "Cl"])
    }
    func testCongestedSubstituentsDoNotOverlap() throws {
        let molecule = try CDKSmilesParser().parseSmiles(
            "COC(=O)C=C(C(=O)O)C(CC1=C(Br)C(O)=C(O)C(Br)=C1Br)C(=O)O")
        let bondLengths = molecule.bonds.map { bond -> CGFloat in
            let p1 = molecule.atom(id: bond.a1)!.position
            let p2 = molecule.atom(id: bond.a2)!.position
            return hypot(p2.x - p1.x, p2.y - p1.y)
        }.sorted()
        let bondLength = bondLengths[bondLengths.count / 2]
        for (index, first) in molecule.atoms.enumerated() {
            for second in molecule.atoms[(index + 1)...] {
                let distance = hypot(first.position.x - second.position.x, first.position.y - second.position.y)
                XCTAssertGreaterThan(distance, bondLength / 4, "atoms \(first.id) and \(second.id)")
            }
        }
    }

    /// Hidden abbreviation atoms are removed before label placement, so a contracted
    /// nitro group is aligned against its one visible bond.
    func testContractedGroupKeepsOnlyVisibleNeighbors() throws {
        var molecule = try CDKSmilesParser().parseSmiles("c1ccccc1[N+](=O)[O-]")
        let ids = Set(molecule.atoms.filter { $0.element != "C" }.map(\.id))
        molecule.sgroups = [
            MoleculeSgroup(
                kind: .generic, keyword: "SUP", atomIDs: ids.sorted(),
                crossingBondIDs: molecule.bonds.filter { ids.contains($0.a1) != ids.contains($0.a2) }.map(\.id),
                subscriptText: "NO2")
        ]
        let rendered = CDKDepictionPreprocessor.prepareForRendering(molecule: molecule, style: RenderStyle())
        let label = try XCTUnwrap(rendered.atoms.first { $0.aliasLabel == "NO2" })
        XCTAssertEqual(rendered.atoms.count, 7)
        XCTAssertEqual(rendered.neighbors(of: label.id).count, 1)
    }
}
