import Foundation
import XCTest
@testable import CDKSwiftNativePort

/// Reference cases: cdk/cdk #1347, #1364, #1350, #1369 and #1379 (2.14.1).
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
}
