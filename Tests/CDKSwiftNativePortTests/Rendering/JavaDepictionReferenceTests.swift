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

        struct LegibilitySample: Decodable {
            let width: CGFloat
            let height: CGFloat
            let singleStroke: CGFloat
            let singleBondLength: CGFloat
            let oxygenInkHeight: CGFloat
            let oxygenClearance: CGFloat
            let carbon13InkHeight: CGFloat
            let chargedNitrogenInkHeight: CGFloat
            let doubleSeparation: CGFloat
            let doubleStroke: CGFloat
            let aromaticSeparation: CGFloat
            let aromaticStroke: CGFloat
            let stereoStroke: CGFloat
            let stereoMaxHash: CGFloat

            var canvas: CGSize { CGSize(width: width, height: height) }
        }

        let version: String
        let revision: String
        let zoomSamples: [ZoomSample]
        let arrowSamples: [ArrowSample]
        let isotope: Isotope
        let legibility: [LegibilitySample]
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

    func testFittedMultipleAndStereoBondsRetainJavaReferenceSeparation() throws {
        let charged = try molecule("carbon-13")
        let stereo = try molecule("stereo")
        let reaction = try CDKRXNReader.readReaction(text: fixture("participants", extension: "rxn"))
        let samples = try loadReference().legibility
        XCTAssertEqual(samples.map(\.width), [200, 560, 960])
        for sample in samples {
            let chargedScene = scene(charged, canvas: sample.canvas)
            // The first three bonds are singles. The fourth is the carbonyl pair.
            let carbonyl = Array(chargedScene.bondSegments.dropFirst(3).prefix(2))
            try assertParallelStrokes(
                carbonyl, javaSeparation: sample.doubleSeparation, javaStroke: sample.doubleStroke)

            let reactionScene = scene(reaction, canvas: sample.canvas)
            // Use the first aromatic pair in participant two. A compact layout may
            // clip the preceding labelled bonds, so locate the first shorter inner line.
            let bonds = reactionScene.bondSegments
            let inner = try XCTUnwrap(
                bonds.indices.dropFirst().first { index in
                    let primary = bonds[index - 1]
                    let secondary = bonds[index]
                    return abs(directionalCosine(primary, secondary)) > 0.99
                        && (0.55...0.90).contains(length(of: secondary) / length(of: primary))
                        && perpendicularDistance(secondary.from, from: primary) > primary.width
                })
            try assertParallelStrokes(
                [bonds[inner - 1], bonds[inner]],
                javaSeparation: sample.aromaticSeparation, javaStroke: sample.aromaticStroke)

            let stereoScene = scene(stereo, canvas: sample.canvas)
            let single = try XCTUnwrap(stereoScene.bondSegments.first)
            let hashes = stereoScene.bondSegments.dropFirst().dropLast()
            XCTAssertGreaterThanOrEqual(hashes.count, 3)
            XCTAssertTrue(hashes.allSatisfy { $0.width >= 0.75 })
            let widest = try XCTUnwrap(hashes.map { length(of: $0) }.max())
            let javaSpread = sample.stereoMaxHash / sample.stereoStroke
            XCTAssertGreaterThan(widest / single.width, javaSpread * 0.65)
            XCTAssertLessThan(widest / single.width, javaSpread * 1.55)
        }
    }

    #if canImport(CoreText)
        func testFittedGlyphInkAndBondClearanceAgainstJavaReference() throws {
            let ethanol = try molecule("ethanol")
            let charged = try molecule("carbon-13")
            for sample in try loadReference().legibility {
                let fitted = scene(ethanol, canvas: sample.canvas)
                let label = try XCTUnwrap(fitted.labels.first { $0.text == "OH" })
                let oxygen = try ink(for: label)
                let plain = try XCTUnwrap(fitted.bondSegments.first)
                let adjoining = try XCTUnwrap(fitted.bondSegments.last)
                let clearance = min(distance(adjoining.from, to: oxygen), distance(adjoining.to, to: oxygen))
                let javaClearance = sample.oxygenClearance / sample.singleStroke
                XCTAssertGreaterThan(clearance / adjoining.width, javaClearance * 0.30)
                XCTAssertLessThan(clearance / adjoining.width, javaClearance * 1.40)
                let javaStrokeToInk = sample.singleStroke / sample.oxygenInkHeight
                XCTAssertGreaterThan(plain.width / oxygen.height, javaStrokeToInk * 0.85)
                XCTAssertLessThan(plain.width / oxygen.height, javaStrokeToInk * 2.8)
                XCTAssertGreaterThanOrEqual(plain.width, 1.4)
                XCTAssertLessThan(plain.width / length(of: plain), sample.singleStroke / sample.singleBondLength)

                let labelled = scene(charged, canvas: sample.canvas)
                let isotope = try XCTUnwrap(labelled.labels.first { $0.text.hasPrefix("¹³C") })
                let nitrogen = try XCTUnwrap(
                    labelled.labels.first { $0.text.hasPrefix("NH") && $0.text.hasSuffix("+") })
                let isotopeHeight = try ink(for: isotope).height
                let nitrogenHeight = try ink(for: nitrogen).height
                XCTAssertGreaterThanOrEqual(isotopeHeight, 5.5)
                XCTAssertGreaterThanOrEqual(nitrogenHeight, 5.5)
                XCTAssertEqual(
                    isotopeHeight / nitrogenHeight, sample.carbon13InkHeight / sample.chargedNitrogenInkHeight,
                    accuracy: 0.32)
            }
        }

        private func ink(for label: CDKMetalDepictionScene.AtomLabel) throws -> CGRect {
            let layout = CDKDepictionLabelLayout(
                text: label.text, fontSize: label.fontSize, italicized: label.italicized)
            let bounds = try XCTUnwrap(layout.outline(centeredAt: label.position, yDown: true)).boundingBoxOfPath
            XCTAssertEqual(bounds.midX, label.position.x, accuracy: 0.05)
            XCTAssertEqual(bounds.midY, label.position.y, accuracy: 0.05)
            XCTAssertEqual(bounds.height, layout.inkBounds.height, accuracy: 0.05)
            return bounds
        }

        private func distance(_ point: CGPoint, to rect: CGRect) -> CGFloat {
            hypot(max(0, rect.minX - point.x, point.x - rect.maxX), max(0, rect.minY - point.y, point.y - rect.maxY))
        }
    #endif

    private func assertParallelStrokes(
        _ strokes: [CDKMetalDepictionScene.LineSegment], javaSeparation: CGFloat, javaStroke: CGFloat
    ) throws {
        XCTAssertEqual(strokes.count, 2)
        let first = try XCTUnwrap(strokes.first)
        let second = try XCTUnwrap(strokes.dropFirst().first)
        XCTAssertGreaterThan(abs(directionalCosine(first, second)), 0.99)
        let separation = perpendicularDistance(second.from, from: first)
        let stroke = max(first.width, second.width)
        XCTAssertGreaterThan(separation - stroke, 0.55)
        XCTAssertGreaterThan(separation / stroke, javaSeparation / javaStroke * 0.40)
        XCTAssertLessThan(separation / stroke, javaSeparation / javaStroke * 1.05)
    }

    private func directionalCosine(
        _ first: CDKMetalDepictionScene.LineSegment, _ second: CDKMetalDepictionScene.LineSegment
    ) -> CGFloat {
        let (dx, dy) = (first.to.x - first.from.x, first.to.y - first.from.y)
        let (ux, uy) = (second.to.x - second.from.x, second.to.y - second.from.y)
        return (dx * ux + dy * uy) / (length(of: first) * length(of: second))
    }

    private func perpendicularDistance(_ point: CGPoint, from line: CDKMetalDepictionScene.LineSegment) -> CGFloat {
        let (dx, dy) = (line.to.x - line.from.x, line.to.y - line.from.y)
        return abs(dx * (point.y - line.from.y) - dy * (point.x - line.from.x)) / length(of: line)
    }

    private func scene(_ molecule: Molecule, canvas: CGSize) -> CDKMetalDepictionScene {
        CDKMetalDepictionSceneBuilder.build(
            molecule: molecule, style: RenderStyle(), canvasRect: CGRect(origin: .zero, size: canvas), zoom: 1,
            pan: .zero)
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

    private func molecule(_ name: String) throws -> Molecule {
        try CDKMDLV2000Reader.read(text: fixture(name, extension: "mol"))
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
