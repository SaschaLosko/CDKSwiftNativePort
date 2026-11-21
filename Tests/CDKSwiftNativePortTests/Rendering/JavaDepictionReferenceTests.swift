import Foundation
import XCTest
@testable import CDKSwiftNativePort

final class JavaDepictionReferenceTests: XCTestCase {
    private struct Reference: Decodable {
        struct ZoomSample: Decodable {
            let zoom: CGFloat
            let firstBondLength: CGFloat
            let stroke: CGFloat
        }

        struct ArrowSample: Decodable {
            let direction: CDKReactionDirection
            let shafts: Int
            let leftHead: Bool
            let rightHead: Bool
            let crossStrokes: Int
        }

        struct Isotope: Decodable {
            let atomNumber: Int
            let element: String
            let mass: Int
            let visible: Bool
        }

        let version: String
        let revision: String
        let zoomSamples: [ZoomSample]
        let arrowSamples: [ArrowSample]
        let isotope: Isotope
    }

    func testReactionZoomMatchesJavaScalingAtDocumentAndCompactSizes() throws {
        let reference = try loadReference()
        let reaction = try CDKRXNReader.readReaction(text: fixture("participants", extension: "rxn"))
        let javaBase = try XCTUnwrap(reference.zoomSamples.first)
        for canvas in [
            CGSize(width: 960, height: 720), CGSize(width: 560, height: 418), CGSize(width: 200, height: 100),
        ] {
            let base = scene(reaction, canvas: canvas)
            let baseBond = try XCTUnwrap(base.bondSegments.first)
            let baseLabel = try XCTUnwrap(base.labels.first { $0.text == "OH" })
            for sample in reference.zoomSamples {
                let scaled = scene(reaction, canvas: canvas, zoom: sample.zoom)
                let bond = try XCTUnwrap(scaled.bondSegments.first)
                let label = try XCTUnwrap(scaled.labels.first { $0.text == baseLabel.text })
                let javaGeometryScale = sample.firstBondLength / javaBase.firstBondLength
                XCTAssertEqual(length(of: bond) / length(of: baseBond), javaGeometryScale, accuracy: 0.002)
                XCTAssertEqual(bond.width / baseBond.width, sample.stroke / javaBase.stroke, accuracy: 0.004)
                XCTAssertEqual(label.fontSize / baseLabel.fontSize, javaGeometryScale, accuracy: 0.002)
            }
        }

        let document = scene(reaction)
        XCTAssertGreaterThanOrEqual(try XCTUnwrap(document.labels.first { $0.text == "OH" }).fontSize, 12)
        XCTAssertGreaterThanOrEqual(try XCTUnwrap(document.bondSegments.first).width, 1.5)
    }

    func testAutomaticReactionArrowMeaningMatchesJava() throws {
        let reference = try loadReference()
        XCTAssertEqual(Set(reference.arrowSamples.map(\.direction)), Set(CDKReactionDirection.allCases))
        var reaction = try CDKRXNReader.readReaction(text: fixture("participants", extension: "rxn"))
        reaction.direction = .undirected
        let undirected = try XCTUnwrap(reference.arrowSamples.first { $0.direction == .undirected })
        let participants = scene(reaction).bondSegments.dropLast(undirected.shafts)
        for sample in reference.arrowSamples {
            reaction.direction = sample.direction
            let segments = scene(reaction).bondSegments
            XCTAssertEqual(Array(segments.prefix(participants.count)), Array(participants), "\(sample.direction)")
            let arrow = segments.dropFirst(participants.count)
            let positions = arrow.flatMap { [$0.from.x, $0.to.x] }
            let left = try XCTUnwrap(positions.min())
            let right = try XCTUnwrap(positions.max())
            let horizontal = arrow.filter { abs($0.from.y - $0.to.y) < 0.001 }
            let diagonals = arrow.filter { abs($0.from.y - $0.to.y) >= 0.001 }
            let leftHead = diagonals.contains { abs($0.to.x - left) < 0.001 }
            let rightHead = diagonals.contains { abs($0.to.x - right) < 0.001 }
            let crossStrokes = diagonals.count { $0.to.x > left + 0.001 && $0.to.x < right - 0.001 }
            XCTAssertEqual(horizontal.count, sample.shafts, "\(sample.direction)")
            XCTAssertEqual(leftHead, sample.leftHead, "\(sample.direction)")
            XCTAssertEqual(rightHead, sample.rightHead, "\(sample.direction)")
            XCTAssertEqual(crossStrokes, sample.crossStrokes, "\(sample.direction)")
            if sample.direction == .bidirectional {
                let leftTip = try XCTUnwrap(diagonals.first { abs($0.to.x - left) < 0.001 })
                let rightTip = try XCTUnwrap(diagonals.first { abs($0.to.x - right) < 0.001 })
                XCTAssertLessThan(rightTip.to.y, leftTip.to.y)
            }
        }
    }

    func testIsotopeLabelRemainsVisibleWithImplicitCarbonsAsInJava() throws {
        let reference = try loadReference().isotope
        let molecule = try CDKMDLV2000Reader.read(text: fixture("carbon-13", extension: "mol"))
        let atom = try XCTUnwrap(molecule.atoms.first { $0.id == reference.atomNumber })
        XCTAssertEqual(atom.element, reference.element)
        XCTAssertEqual(atom.isotopeMassNumber, reference.mass)

        let style = RenderStyle()
        XCTAssertFalse(style.showCarbons)
        let moleculeScene = CDKMetalDepictionSceneBuilder.build(
            molecule: molecule, style: style, canvasRect: CGRect(x: 0, y: 0, width: 560, height: 418),
            zoom: 1, pan: .zero)
        let label = "¹³C"
        XCTAssertEqual(moleculeScene.labels.contains { $0.text.hasPrefix(label) }, reference.visible)
        let reaction = CDKReaction(reactants: [molecule], agents: [], products: [molecule])
        XCTAssertEqual(scene(reaction).labels.count { $0.text.hasPrefix(label) }, reference.visible ? 2 : 0)
    }

    private func scene(
        _ reaction: CDKReaction, canvas: CGSize = CGSize(width: 560, height: 418), zoom: CGFloat = 1
    ) -> CDKMetalDepictionScene {
        CDKMetalReactionDepictionSceneBuilder.build(
            reaction: reaction, style: RenderStyle(), canvasRect: CGRect(origin: .zero, size: canvas),
            zoom: zoom, pan: .zero)
    }

    private func length(of segment: CDKMetalDepictionScene.LineSegment) -> CGFloat {
        hypot(segment.to.x - segment.from.x, segment.to.y - segment.from.y)
    }

    private func loadReference() throws -> Reference {
        let reference = try JSONDecoder().decode(
            Reference.self, from: Data(fixture("expectations", extension: "json").utf8))
        XCTAssertEqual(reference.version, "2.13")
        XCTAssertEqual(reference.revision, "3bc6efe84119e5da0d00804bd4cf5f3bdd9c3d6d")
        return reference
    }

    private func fixture(_ name: String, extension fileExtension: String) throws -> String {
        let url = try XCTUnwrap(
            Bundle.module.url(
                forResource: name, withExtension: fileExtension, subdirectory: "JavaReference"))
        return try String(contentsOf: url, encoding: .utf8)
    }
}
