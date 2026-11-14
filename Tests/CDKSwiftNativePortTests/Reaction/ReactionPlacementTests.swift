import Foundation
import XCTest
@testable import CDKSwiftNativePort

final class ReactionPlacementTests: XCTestCase {
    private func reaction() throws -> CDKReaction {
        let parser = CDKSmilesParser(generateCoordinates: false)
        return CDKReaction(
            reactantParticipants: [.init(molecule: try parser.parseSmiles("CC"), role: .reactant)],
            agentParticipants: [
                .init(molecule: try parser.parseSmiles("O"), role: .agent, arrowPosition: .below),
                .init(molecule: try parser.parseSmiles("N"), role: .agent, arrowPosition: .above),
            ],
            productParticipants: [.init(molecule: try parser.parseSmiles("C=C"), role: .product)])
    }

    func testPlacementDoesNotChangeChemicalRolesOrSMILES() throws {
        let placed = try reaction()
        var ordinary = placed
        for index in ordinary.agentParticipants.indices { ordinary.agentParticipants[index].arrowPosition = nil }
        XCTAssertEqual(CDKSmilesGenerator().create(placed), CDKSmilesGenerator().create(ordinary))
        XCTAssertEqual(placed.agents, ordinary.agents)
        XCTAssertEqual(placed.agentParticipants.map(\.role), [.agent, .agent])
    }

    func testParticipantCodingRetainsPositionAndReadsOlderDocuments() throws {
        var participant = try reaction().agentParticipants[0]
        let encoder = JSONEncoder()
        XCTAssertEqual(
            try JSONDecoder().decode(CDKReactionParticipant.self, from: encoder.encode(participant)), participant)
        participant.arrowPosition = nil
        let legacy = try encoder.encode(participant)
        XCTAssertFalse(String(decoding: legacy, as: UTF8.self).contains("arrowPosition"))
        XCTAssertNil(try JSONDecoder().decode(CDKReactionParticipant.self, from: legacy).arrowPosition)
    }

    func testCMLPreservesParticipantPositionAndCoefficient() throws {
        var source = try reaction()
        source.agentParticipants[0].stoichiometry = 2
        let text = try CDKCMLReactionWriter.write(source)
        XCTAssertTrue(text.contains("cdk:arrowPosition=\"below\""))
        let restored = try CDKCMLReactionReader.readReaction(text: text)
        XCTAssertEqual(restored.agentParticipants.map(\.arrowPosition), [.below, .above])
        XCTAssertEqual(restored.agentParticipants[0].stoichiometry, 2)
        XCTAssertEqual(
            restored.participants.map { $0.molecule.atoms.map(\.element) },
            source.participants.map { $0.molecule.atoms.map(\.element) })
    }

    func testInvalidPresentationAndCoefficientDoNotSilentlyDisappear() throws {
        var source = try reaction()
        let text = try CDKCMLReactionWriter.write(source)
        XCTAssertThrowsError(
            try CDKCMLReactionReader.readReaction(text: text.replacingOccurrences(of: "below", with: "future")))
        source.agentParticipants[0].stoichiometry = .nan
        XCTAssertThrowsError(try CDKCMLReactionWriter.write(source))
    }

    func testCMLCoefficientPrecisionIsIndependentOfCoordinateFormatting() throws {
        for count in [1.2345678901234567, 1.0 / 3.0, 1.0e-12, 1.0e12] {
            var source = try reaction()
            source.reactantParticipants[0].stoichiometry = count
            source.agentParticipants[0].stoichiometry = count
            source.productParticipants[0].stoichiometry = count
            let restored = try CDKCMLReactionReader.readReaction(text: CDKCMLReactionWriter.write(source))
            XCTAssertEqual(restored.reactantParticipants[0].stoichiometry, count)
            XCTAssertEqual(restored.agentParticipants[0].stoichiometry, count)
            XCTAssertEqual(restored.productParticipants[0].stoichiometry, count)
        }
    }

    func testCMLSharedMoleculeReferenceDoesNotShareParticipantPlacement() throws {
        let molecule = try reaction().agents[0]
        let source = CDKReaction(
            reactantParticipants: [],
            agentParticipants: [
                .init(molecule: molecule, role: .agent, arrowPosition: .above),
                .init(molecule: molecule, role: .agent, arrowPosition: .below),
            ], productParticipants: [])
        let restored = try CDKCMLReactionReader.readReaction(text: CDKCMLReactionWriter.write(source))
        XCTAssertEqual(restored.agentParticipants.map(\.arrowPosition), [.above, .below])
    }

    func testRenderingAndHitTestingKeepOriginalParticipantIndices() throws {
        let source = try reaction()
        let style = RenderStyle(
            showCarbons: true, showImplicitHydrogens: false, showAtomIDs: false,
            bondWidth: 2, fontSize: 16, padding: 20)
        let canvas = CGRect(x: 0, y: 0, width: 1000, height: 600)
        for zoom: CGFloat in [1, 1.2] {
            let pan = CGSize(width: 12, height: -5)
            let scene = CDKMetalReactionDepictionSceneBuilder.build(
                reaction: source, style: style, canvasRect: canvas, zoom: zoom, pan: pan)
            let oxygen = try XCTUnwrap(scene.labels.first { $0.text == "O" })
            let nitrogen = try XCTUnwrap(scene.labels.first { $0.text == "N" })
            XCTAssertGreaterThan(oxygen.position.y, canvas.midY + pan.height)
            XCTAssertLessThan(nitrogen.position.y, canvas.midY + pan.height)
            for (label, index) in [(oxygen, 0), (nitrogen, 1)] {
                XCTAssertEqual(
                    CDKMetalReactionDepictionSceneBuilder.participant(
                        at: label.position, in: source, style: style, canvasRect: canvas, zoom: zoom, pan: pan),
                    CDKReactionParticipantSelection(role: .agent, index: index))
            }
        }
    }
}
