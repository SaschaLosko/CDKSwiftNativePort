import Foundation
import XCTest
@testable import CDKSwiftNativePort

final class JavaStereoDepictionReferenceTests: XCTestCase {
    private struct Reference: Decodable {
        struct Canvas: Decodable {
            let width: CGFloat
            let height: CGFloat

            var size: CGSize { CGSize(width: width, height: height) }
        }

        struct Sample: Decodable {
            struct Participant: Decodable {
                enum AlkeneSides: String, Decodable { case same, opposite }
                enum StereoDisplay: String, Decodable { case solid, hashed }

                let labels: Int
                let bondLines: Int
                let solidWedges: Int
                let alkeneSides: AlkeneSides?
                let stereoDisplay: StereoDisplay?
            }

            let name: String
            let source: String
            let participants: [Participant]
        }

        let version: String
        let revision: String
        let canvases: [Canvas]
        let samples: [Sample]
    }

    func testEnantiomersMatchJavaWedgeAndHashConventions() throws {
        let reference = try loadReference()
        let samples = reference.samples.filter { $0.name == "chiral-at" || $0.name == "chiral-at-at" }
        XCTAssertEqual(samples.count, 2)
        for sample in samples {
            let molecule = try CDKSmilesParser().parseSmiles(sample.source)
            let expected = try XCTUnwrap(sample.participants.first)
            let stereoBonds = molecule.bonds.filter { $0.stereo != .none }
            XCTAssertEqual(stereoBonds.count, 1, sample.name)
            let stereo = try XCTUnwrap(stereoBonds.first).stereo
            let display: Reference.Sample.Participant.StereoDisplay =
                (stereo == .up || stereo == .upReversed) ? .solid : .hashed
            XCTAssertEqual(display, expected.stereoDisplay, sample.name)

            for canvas in reference.canvases {
                let scene = scene(molecule, canvas: canvas.size)
                XCTAssertEqual(scene.labels.count, expected.labels, sample.name)
                XCTAssertGreaterThan(scene.bondSegments.count, expected.bondLines, sample.name)
                let svg = CDKDepictionGenerator.toSVG(molecule: molecule, canvasSize: canvas.size)
                XCTAssertEqual(svg.components(separatedBy: "<polygon ").count - 1, expected.solidWedges, sample.name)
                if display == .hashed {
                    XCTAssertGreaterThanOrEqual(
                        svg.components(separatedBy: "<line ").count - 1, expected.bondLines - 3, sample.name)
                }
            }
        }
    }

    func testAlkenesMatchJavaSubstituentSidesAndRenderOrdinaryLines() throws {
        let reference = try loadReference()
        let alkenes = reference.samples.filter { $0.name == "alkene-e" || $0.name == "alkene-z" }
        XCTAssertEqual(alkenes.count, 2)
        for sample in alkenes {
            let molecule = try CDKSmilesParser().parseSmiles(sample.source)
            let participant = try XCTUnwrap(sample.participants.first)
            try assertAlkene(molecule, equals: participant, name: sample.name)
            for canvas in reference.canvases {
                let scene = scene(molecule, canvas: canvas.size)
                XCTAssertEqual(scene.bondSegments.count, participant.bondLines, sample.name)
                let doubleLine = try XCTUnwrap(scene.bondSegments.dropFirst().first)
                let labels = scene.labels.filter { $0.text == "F" }
                XCTAssertEqual(labels.count, participant.labels, sample.name)
                let first = try XCTUnwrap(labels.first).position
                let last = try XCTUnwrap(labels.last).position
                XCTAssertEqual(
                    substituentSides(first, last, axisFrom: doubleLine.from, axisTo: doubleLine.to),
                    participant.alkeneSides, sample.name)
                let svg = CDKDepictionGenerator.toSVG(molecule: molecule, canvasSize: canvas.size)
                XCTAssertFalse(svg.contains("<polygon "), sample.name)
                XCTAssertEqual(svg.components(separatedBy: "<line ").count - 1, participant.bondLines, sample.name)
            }
        }
    }

    func testStereoReactionsRetainParticipantGeometryAndArrow() throws {
        let reference = try loadReference()
        let reactions = reference.samples.filter { $0.source.contains(">>") }
        XCTAssertEqual(reactions.count, 2)
        for sample in reactions {
            let reaction = try CDKSmilesParser().parseReactionSmiles(sample.source)
            let molecules = reaction.reactants + reaction.products
            XCTAssertEqual(molecules.count, sample.participants.count, sample.name)
            for (molecule, participant) in zip(molecules, sample.participants) {
                if participant.alkeneSides != nil {
                    try assertAlkene(molecule, equals: participant, name: sample.name)
                }
            }
            for canvas in reference.canvases {
                let drawn = CDKMetalReactionDepictionSceneBuilder.build(
                    reaction: reaction, style: RenderStyle(), canvasRect: CGRect(origin: .zero, size: canvas.size),
                    zoom: 1, pan: .zero)
                XCTAssertEqual(drawn.filledPolygons.count, 1, sample.name)
                XCTAssertEqual(drawn.labels.count, sample.participants.map(\.labels).reduce(0, +), sample.name)
                // The scene stores the forward shaft alongside bond strokes.
                let bondLines = sample.participants.map(\.bondLines).reduce(0, +) + 1
                if sample.name == "alkene-reaction" {
                    XCTAssertEqual(drawn.bondSegments.count, bondLines, sample.name)
                } else {
                    XCTAssertGreaterThan(drawn.bondSegments.count, bondLines, sample.name)
                }
            }
        }
    }

    func testCrowdedFixtureKeepsEveryJavaLabelAtBothSizes() throws {
        let reference = try loadReference()
        let sample = try XCTUnwrap(reference.samples.first { $0.name == "crowded-labels" })
        let expected = try XCTUnwrap(sample.participants.first)
        let molecule = try CDKSmilesParser().parseSmiles(sample.source)
        var style = RenderStyle()
        style.padding = 32
        for canvas in reference.canvases {
            let drawn = CDKMetalDepictionSceneBuilder.build(
                molecule: molecule, style: style, canvasRect: CGRect(origin: .zero, size: canvas.size),
                zoom: 1, pan: .zero)
            XCTAssertEqual(drawn.labels.count, expected.labels)
            XCTAssertEqual(drawn.bondSegments.count, expected.bondLines)
            for prefix in ["¹³C", "NH", "Cl", "F", "O-"] {
                XCTAssertTrue(drawn.labels.contains { $0.text.hasPrefix(prefix) }, prefix)
            }
            let bounds = CGRect(origin: .zero, size: canvas.size)
            XCTAssertTrue(drawn.labels.allSatisfy { bounds.contains($0.position) })
            #if canImport(CoreText)
                let ink = try drawn.labels.map { label in
                    let glyphs = CDKDepictionLabelLayout(text: label.text, fontSize: label.fontSize)
                    return try XCTUnwrap(glyphs.outline(centeredAt: label.position, yDown: true)).boundingBoxOfPath
                }
                for (index, rectangle) in ink.enumerated() {
                    XCTAssertTrue(bounds.contains(rectangle), "\(drawn.labels[index].text) \(rectangle) in \(bounds)")
                    XCTAssertFalse(ink.dropFirst(index + 1).contains { rectangle.intersects($0) })
                }
            #endif
        }
    }

    private func assertAlkene(
        _ molecule: Molecule, equals expected: Reference.Sample.Participant, name: String
    ) throws {
        let bond = try XCTUnwrap(molecule.bonds.first { $0.doubleBondStereo != nil }, name)
        let refs = try XCTUnwrap(bond.stereoReferenceAtomIDs, name)
        let positions = try refs.map { id in try XCTUnwrap(molecule.atoms.first { $0.id == id }).position }
        guard positions.count == 4 else { return XCTFail("Expected four alkene references in \(name)") }
        XCTAssertEqual(bond.doubleBondStereo, expected.alkeneSides == .same ? .cis : .trans, name)
        XCTAssertEqual(
            substituentSides(positions[0], positions[3], axisFrom: positions[1], axisTo: positions[2]),
            expected.alkeneSides, name)
    }

    private func substituentSides(
        _ first: CGPoint, _ last: CGPoint, axisFrom start: CGPoint, axisTo end: CGPoint
    ) -> Reference.Sample.Participant.AlkeneSides? {
        func side(_ point: CGPoint) -> CGFloat {
            (end.x - start.x) * (point.y - start.y) - (end.y - start.y) * (point.x - start.x)
        }
        let sides = side(first) * side(last)
        guard abs(sides) > 1e-8 else { return nil }
        return sides > 0 ? .same : .opposite
    }

    private func scene(_ molecule: Molecule, canvas: CGSize) -> CDKMetalDepictionScene {
        CDKMetalDepictionSceneBuilder.build(
            molecule: molecule, style: RenderStyle(), canvasRect: CGRect(origin: .zero, size: canvas),
            zoom: 1, pan: .zero)
    }

    private func loadReference() throws -> Reference {
        let url = try XCTUnwrap(
            Bundle.module.url(forResource: "stereochemistry", withExtension: "json", subdirectory: "JavaReference"))
        let reference = try JSONDecoder().decode(Reference.self, from: Data(contentsOf: url))
        XCTAssertEqual(reference.version, "2.14.1")
        XCTAssertEqual(reference.revision, "4253aa225f474a2a5b26ad9ed3ecf73262355693")
        XCTAssertEqual(reference.canvases.map(\.width), [224, 560])
        return reference
    }
}
