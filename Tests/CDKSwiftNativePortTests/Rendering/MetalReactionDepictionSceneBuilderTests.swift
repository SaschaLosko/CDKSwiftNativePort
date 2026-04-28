import Foundation
#if canImport(CoreGraphics)
import CoreGraphics
#endif
import XCTest
@testable import CDKSwiftNativePort

final class MetalReactionDepictionSceneBuilderTests: XCTestCase {
    func testBuildsReactionSceneWithArrowAndPlusMarkers() {
        let reaction = CDKReaction(reactants: [ethane(name: "R1"), ethane(name: "R2")],
                                   agents: [],
                                   products: [ethane(name: "P1")])
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.2,
                                fontSize: 16,
                                padding: 20)

        let scene = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                style: style,
                                                                canvasRect: CGRect(x: 0, y: 0, width: 880, height: 420),
                                                                zoom: 1.0,
                                                                pan: .zero,
                                                                rotationDegrees: 0)

        XCTAssertFalse(scene.bondSegments.isEmpty)
        XCTAssertTrue(scene.labels.contains(where: { $0.text == "+" }),
                      "Expected plus markers between multiple reaction participants.")

        let horizontalArrowCandidates = scene.bondSegments.filter {
            abs($0.from.y - $0.to.y) < 1.5 && ($0.to.x - $0.from.x) > 36
        }
        XCTAssertFalse(horizontalArrowCandidates.isEmpty,
                       "Expected a reaction arrow baseline in the rendered reaction scene.")
    }

    func testThreeReactantReactionWrapsParticipantsWithoutFlatPlusRow() throws {
        let reaction = try brominationStyleReaction()
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.1,
                                fontSize: 15,
                                padding: 22)
        let canvas = CGRect(x: 0, y: 0, width: 1360, height: 620)

        let scene = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                style: style,
                                                                canvasRect: canvas,
                                                                zoom: 1.0,
                                                                pan: .zero,
                                                                rotationDegrees: 0)

        XCTAssertFalse(scene.labels.contains(where: { $0.text == "+" }),
                       "Wrapped multi-row reaction layouts should suppress plus glyphs and pack participants more like publication-style depictions.")

        let leftLabels = scene.labels.filter {
            $0.text != "+" && $0.position.x < (canvas.midX - 50)
        }
        let rowBands = Set(leftLabels.map { Int(($0.position.y / 72).rounded()) })
        XCTAssertGreaterThanOrEqual(rowBands.count, 2,
                                    "Expected multi-row packing on the reactant side for dense three-reactant reactions.")

        let rightLabels = scene.labels.filter {
            $0.text != "+" && $0.position.x > (canvas.midX + 70)
        }
        let rightWidth = (rightLabels.map(\.position.x).max() ?? 0) - (rightLabels.map(\.position.x).min() ?? 0)
        XCTAssertGreaterThan(rightWidth, 120,
                             "Single products should retain a reasonably wide lane instead of being squeezed into a narrow vertical strip.")

        let rightBondPoints = scene.bondSegments
            .filter { min($0.from.x, $0.to.x) > (canvas.midX + 90) }
            .flatMap { [$0.from, $0.to] }
        XCTAssertFalse(rightBondPoints.isEmpty)
        let rightBounds = CGRect(x: rightBondPoints.map(\.x).min() ?? 0,
                                 y: rightBondPoints.map(\.y).min() ?? 0,
                                 width: (rightBondPoints.map(\.x).max() ?? 0) - (rightBondPoints.map(\.x).min() ?? 0),
                                 height: (rightBondPoints.map(\.y).max() ?? 0) - (rightBondPoints.map(\.y).min() ?? 0))
        XCTAssertGreaterThan(rightBounds.width, rightBounds.height * 1.05,
                             "The bromination product should read horizontally on the product side instead of as a tall vertical participant.")
    }

    func testReactionSceneScalesLabelFontWithZoom() {
        let reaction = CDKReaction(reactants: [ethane(name: "R")], agents: [], products: [ethane(name: "P")])
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14,
                                padding: 20)

        let zoom1 = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                style: style,
                                                                canvasRect: CGRect(x: 0, y: 0, width: 760, height: 360),
                                                                zoom: 1.0,
                                                                pan: .zero,
                                                                rotationDegrees: 0)
        let zoom2 = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                style: style,
                                                                canvasRect: CGRect(x: 0, y: 0, width: 760, height: 360),
                                                                zoom: 2.0,
                                                                pan: .zero,
                                                                rotationDegrees: 0)

        let maxFont1 = zoom1.labels.map(\.fontSize).max() ?? 0
        let maxFont2 = zoom2.labels.map(\.fontSize).max() ?? 0
        XCTAssertGreaterThan(maxFont2, maxFont1 * 1.5)
    }

    func testReactionSceneScalesStrokeAndLabelSizingWithCanvasResize() {
        let reaction = CDKReaction(reactants: [ethane(name: "R"), ethane(name: "R2")],
                                   agents: [],
                                   products: [ethane(name: "P")])
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14,
                                padding: 20)

        let compact = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                  style: style,
                                                                  canvasRect: CGRect(x: 0, y: 0, width: 760, height: 360),
                                                                  zoom: 1.0,
                                                                  pan: .zero,
                                                                  rotationDegrees: 0)
        let expanded = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                   style: style,
                                                                   canvasRect: CGRect(x: 0, y: 0, width: 1320, height: 700),
                                                                   zoom: 1.0,
                                                                   pan: .zero,
                                                                   rotationDegrees: 0)

        let compactWidths = compact.bondSegments
            .map(\.width)
            .filter { $0 > 0.01 }
        let expandedWidths = expanded.bondSegments
            .map(\.width)
            .filter { $0 > 0.01 }
        let compactAvgWidth = compactWidths.reduce(0, +) / CGFloat(max(1, compactWidths.count))
        let expandedAvgWidth = expandedWidths.reduce(0, +) / CGFloat(max(1, expandedWidths.count))
        XCTAssertGreaterThan(expandedAvgWidth, compactAvgWidth * 1.15,
                             "Reaction line widths should grow when the canvas auto-fit grows with window size.")

        let compactMaxFont = compact.labels.map(\.fontSize).max() ?? 0
        let expandedMaxFont = expanded.labels.map(\.fontSize).max() ?? 0
        XCTAssertGreaterThan(expandedMaxFont, compactMaxFont * 1.15,
                             "Reaction label sizes should scale up with window-size auto-fit, not stay static.")
    }

    func testReactionSceneKeepsShrinkingOnVeryCompactCanvas() {
        let reaction = CDKReaction(reactants: [moleculeWithHeteroLabel(name: "R")],
                                   agents: [],
                                   products: [])
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14,
                                padding: 20)

        let compact = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                  style: style,
                                                                  canvasRect: CGRect(x: 0, y: 0, width: 760, height: 420),
                                                                  zoom: 1.0,
                                                                  pan: .zero,
                                                                  rotationDegrees: 0)
        let extraCompact = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                       style: style,
                                                                       canvasRect: CGRect(x: 0, y: 0, width: 540, height: 340),
                                                                       zoom: 1.0,
                                                                       pan: .zero,
                                                                       rotationDegrees: 0)

        let compactWidths = compact.bondSegments
            .map(\.width)
            .filter { $0 > 0.01 }
        let extraCompactWidths = extraCompact.bondSegments
            .map(\.width)
            .filter { $0 > 0.01 }
        let compactAvgWidth = compactWidths.reduce(0, +) / CGFloat(max(1, compactWidths.count))
        let extraCompactAvgWidth = extraCompactWidths.reduce(0, +) / CGFloat(max(1, extraCompactWidths.count))
        XCTAssertLessThan(extraCompactAvgWidth, compactAvgWidth * 0.72,
                          "Expected bond widths to keep shrinking on very compact reaction canvases.")

        let compactMaxFont = compact.labels.map(\.fontSize).max() ?? 0
        let extraCompactMaxFont = extraCompact.labels.map(\.fontSize).max() ?? 0
        XCTAssertLessThanOrEqual(extraCompactMaxFont, compactMaxFont,
                                 "Expected reaction labels to never grow on very compact canvases.")
    }

    func testHighlightedParticipantAddsOuterGlowForBondsAndLabels() {
        let reaction = CDKReaction(reactants: [ethane(name: "R1"), ethane(name: "R2")],
                                   agents: [],
                                   products: [ethane(name: "P1")])
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14,
                                padding: 20)
        let canvas = CGRect(x: 0, y: 0, width: 920, height: 420)

        let baseline = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                   style: style,
                                                                   canvasRect: canvas,
                                                                   zoom: 1.0,
                                                                   pan: .zero,
                                                                   rotationDegrees: 0)
        let highlighted = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                      style: style,
                                                                      canvasRect: canvas,
                                                                      zoom: 1.0,
                                                                      pan: .zero,
                                                                      rotationDegrees: 0,
                                                                      highlightedParticipant: CDKReactionParticipantSelection(role: .product,
                                                                                                                              index: 0))
        let baselineGlowSegments = baseline.bondSegments.filter(isOuterGlowSegment)
        let highlightedGlowSegments = highlighted.bondSegments.filter(isOuterGlowSegment)
        XCTAssertTrue(baselineGlowSegments.isEmpty,
                      "Without selection there should be no outer-glow bond overlays.")
        XCTAssertFalse(highlightedGlowSegments.isEmpty,
                       "Selected participants should render outer-glow bond overlays.")

        XCTAssertEqual(highlighted.labels.filter(\.drawsBackground).count,
                       baseline.labels.count,
                       "Outer glow should add overlays, not replace normal atom labels.")
        XCTAssertTrue(highlighted.labels.contains(where: { !$0.drawsBackground }),
                      "Selected participants should render backgroundless outer-glow label overlays.")
        let baseFonts = highlighted.labels.filter(\.drawsBackground).map(\.fontSize)
        let glowFonts = highlighted.labels.filter { !$0.drawsBackground }.map(\.fontSize)
        if let maxBaseFont = baseFonts.max(), let maxGlowFont = glowFonts.max() {
            XCTAssertLessThanOrEqual(maxGlowFont, maxBaseFont * 1.10,
                                     "Outer-glow labels should hug glyph edges and not be heavily oversized.")
        }
    }

    func testOuterGlowHighlightFollowsRoleAndIndexAcrossReactionSides() {
        let reaction = CDKReaction(reactants: [moleculeWithHeteroLabel(name: "R1")],
                                   agents: [],
                                   products: [moleculeWithHeteroLabel(name: "P1")])
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14,
                                padding: 20)
        let canvas = CGRect(x: 0, y: 0, width: 960, height: 420)

        let reactantHighlighted = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                              style: style,
                                                                              canvasRect: canvas,
                                                                              zoom: 1.0,
                                                                              pan: .zero,
                                                                              rotationDegrees: 0,
                                                                              highlightedParticipant: CDKReactionParticipantSelection(role: .reactant,
                                                                                                                                      index: 0))
        let productHighlighted = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                             style: style,
                                                                             canvasRect: canvas,
                                                                             zoom: 1.0,
                                                                             pan: .zero,
                                                                             rotationDegrees: 0,
                                                                             highlightedParticipant: CDKReactionParticipantSelection(role: .product,
                                                                                                                                     index: 0))

        let leftMaxX = canvas.midX - 24
        let rightMinX = canvas.midX + 24

        let reactantLeftGlow = glowSegmentCount(in: reactantHighlighted, minX: 0, maxX: leftMaxX)
        let reactantRightGlow = glowSegmentCount(in: reactantHighlighted, minX: rightMinX, maxX: canvas.maxX)
        let productLeftGlow = glowSegmentCount(in: productHighlighted, minX: 0, maxX: leftMaxX)
        let productRightGlow = glowSegmentCount(in: productHighlighted, minX: rightMinX, maxX: canvas.maxX)

        XCTAssertGreaterThan(reactantLeftGlow, reactantRightGlow,
                             "Reactant selection should emphasize left-side reactant geometry.")
        XCTAssertGreaterThan(productRightGlow, productLeftGlow,
                             "Product selection should emphasize right-side product geometry.")
    }

    func testOuterGlowCanBeDisabledViaBuildParameter() {
        let reaction = CDKReaction(reactants: [ethane(name: "R1")],
                                   agents: [],
                                   products: [ethane(name: "P1")])
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14,
                                padding: 20)
        let canvas = CGRect(x: 0, y: 0, width: 900, height: 380)

        let highlightedWithoutGlow = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                                  style: style,
                                                                                  canvasRect: canvas,
                                                                                  zoom: 1.0,
                                                                                  pan: .zero,
                                                                                  rotationDegrees: 0,
                                                                                  highlightedParticipant: CDKReactionParticipantSelection(role: .reactant,
                                                                                                                                          index: 0),
                                                                                  withOuterGlowHighlight: false)

        XCTAssertFalse(highlightedWithoutGlow.labels.contains(where: { !$0.drawsBackground }))
        XCTAssertFalse(highlightedWithoutGlow.bondSegments.contains(where: isOuterGlowSegment))
    }

    func testOutOfRangeHighlightedParticipantFallsBackToUnhighlightedRendering() {
        let reaction = CDKReaction(reactants: [ethane(name: "R1")],
                                   agents: [],
                                   products: [ethane(name: "P1")])
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14,
                                padding: 20)
        let canvas = CGRect(x: 0, y: 0, width: 900, height: 380)

        let baseline = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                   style: style,
                                                                   canvasRect: canvas,
                                                                   zoom: 1.0,
                                                                   pan: .zero,
                                                                   rotationDegrees: 0)
        let invalid = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                  style: style,
                                                                  canvasRect: canvas,
                                                                  zoom: 1.0,
                                                                  pan: .zero,
                                                                  rotationDegrees: 0,
                                                                  highlightedParticipant: CDKReactionParticipantSelection(role: .agent,
                                                                                                                          index: 0))

        XCTAssertEqual(invalid, baseline,
                       "Invalid participant selections should be ignored instead of dimming all participants.")
    }

    func testParticipantHitTestingFindsReactantAndProductFromRenderedScene() {
        let reaction = CDKReaction(reactants: [moleculeWithHeteroLabel(name: "R1")],
                                   agents: [],
                                   products: [moleculeWithHeteroLabel(name: "P1")])
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14,
                                padding: 20)
        let canvas = CGRect(x: 0, y: 0, width: 960, height: 420)

        let scene = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                style: style,
                                                                canvasRect: canvas,
                                                                zoom: 1.0,
                                                                pan: .zero,
                                                                rotationDegrees: 0)
        let atomLabels = scene.labels.filter { $0.text != "+" }
        guard let leftPoint = atomLabels.min(by: { $0.position.x < $1.position.x })?.position,
              let rightPoint = atomLabels.max(by: { $0.position.x < $1.position.x })?.position else {
            XCTFail("Expected atom labels to derive hit-test sample points.")
            return
        }

        let leftHit = CDKMetalReactionDepictionSceneBuilder.participant(at: leftPoint,
                                                                        in: reaction,
                                                                        style: style,
                                                                        canvasRect: canvas,
                                                                        zoom: 1.0,
                                                                        pan: .zero,
                                                                        rotationDegrees: 0)
        let rightHit = CDKMetalReactionDepictionSceneBuilder.participant(at: rightPoint,
                                                                         in: reaction,
                                                                         style: style,
                                                                         canvasRect: canvas,
                                                                         zoom: 1.0,
                                                                         pan: .zero,
                                                                         rotationDegrees: 0)

        XCTAssertEqual(leftHit, CDKReactionParticipantSelection(role: .reactant, index: 0))
        XCTAssertEqual(rightHit, CDKReactionParticipantSelection(role: .product, index: 0))
    }

    func testParticipantHitTestingSupportsZoomPanRotationAndMissesArrowGap() {
        let reaction = CDKReaction(reactants: [moleculeWithHeteroLabel(name: "R1")],
                                   agents: [],
                                   products: [moleculeWithHeteroLabel(name: "P1")])
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14,
                                padding: 20)
        let canvas = CGRect(x: 0, y: 0, width: 980, height: 420)
        let zoom: CGFloat = 1.48
        let pan = CGSize(width: 46, height: -18)
        let rotation: CGFloat = 17

        let scene = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                style: style,
                                                                canvasRect: canvas,
                                                                zoom: zoom,
                                                                pan: pan,
                                                                rotationDegrees: rotation)
        let atomLabels = scene.labels.filter { $0.text != "+" }
        guard let leftPoint = atomLabels.min(by: { $0.position.x < $1.position.x })?.position,
              let rightPoint = atomLabels.max(by: { $0.position.x < $1.position.x })?.position else {
            XCTFail("Expected atom labels to derive transformed hit-test sample points.")
            return
        }

        let leftHit = CDKMetalReactionDepictionSceneBuilder.participant(at: leftPoint,
                                                                        in: reaction,
                                                                        style: style,
                                                                        canvasRect: canvas,
                                                                        zoom: zoom,
                                                                        pan: pan,
                                                                        rotationDegrees: rotation)
        let rightHit = CDKMetalReactionDepictionSceneBuilder.participant(at: rightPoint,
                                                                         in: reaction,
                                                                         style: style,
                                                                         canvasRect: canvas,
                                                                         zoom: zoom,
                                                                         pan: pan,
                                                                         rotationDegrees: rotation)

        XCTAssertEqual(leftHit, CDKReactionParticipantSelection(role: .reactant, index: 0))
        XCTAssertEqual(rightHit, CDKReactionParticipantSelection(role: .product, index: 0))

        let gapCenter = viewportTransformed(CGPoint(x: canvas.midX, y: canvas.midY),
                                            canvasRect: canvas,
                                            zoom: zoom,
                                            pan: pan,
                                            rotationDegrees: rotation)
        let miss = CDKMetalReactionDepictionSceneBuilder.participant(at: gapCenter,
                                                                     in: reaction,
                                                                     style: style,
                                                                     canvasRect: canvas,
                                                                     zoom: zoom,
                                                                     pan: pan,
                                                                     rotationDegrees: rotation)
        XCTAssertNil(miss, "Point in the reaction arrow gap should not select any participant.")
    }

    func testDenseReactionClampsLabelFontForCDKLikeReadability() {
        let hetero = moleculeWithHeteroLabel(name: "Het")
        let reaction = CDKReaction(reactants: [hetero, hetero, hetero],
                                   agents: [hetero],
                                   products: [hetero, hetero])
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.35,
                                fontSize: 28,
                                padding: 24)

        let scene = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                style: style,
                                                                canvasRect: CGRect(x: 0, y: 0, width: 1440, height: 720),
                                                                zoom: 1.0,
                                                                pan: .zero,
                                                                rotationDegrees: 0)

        let maxAtomLabel = scene.labels.filter { $0.text != "+" }.map(\.fontSize).max() ?? 0
        XCTAssertLessThanOrEqual(maxAtomLabel, 16.0,
                                 "Reaction atom labels should be compact and not explode with a large base font style.")
    }

    func testLowQualityParticipantCoordinatesTriggerPerParticipantRelayout() throws {
        var distorted = try parseSmiles("NC1=NC=NC2=N1N=CN2")
        for idx in distorted.atoms.indices {
            let x = CGFloat(idx % 4) * 0.08
            let y = CGFloat(idx / 4) * 0.06
            distorted.atoms[idx].position = CGPoint(x: x, y: y)
        }

        let originalPenalty = CDKReactionParticipantLayoutRefiner.qualityPenalty(for: distorted)
        XCTAssertTrue(CDKReactionParticipantLayoutRefiner.shouldRecompute2DLayout(for: distorted))

        let refined = CDKReactionParticipantLayoutRefiner.recomputeIfLowQuality(distorted)
        let refinedPenalty = CDKReactionParticipantLayoutRefiner.qualityPenalty(for: refined)

        XCTAssertLessThan(refinedPenalty, originalPenalty,
                          "Low-quality participant coordinates should be improved by SDG relayout.")
        XCTAssertNotEqual(refined.boundingBox(), distorted.boundingBox(),
                          "Refined coordinates should not keep the same highly compressed participant box.")
    }

    func testAlreadyGoodParticipantCoordinatesAreNotRecomputed() throws {
        let base = try parseSmiles("CC(=O)OC1=CC=CC=C1C(=O)O")
        let good = Depiction2DGenerator.generate(for: base)

        XCTAssertFalse(CDKReactionParticipantLayoutRefiner.shouldRecompute2DLayout(for: good))

        let refined = CDKReactionParticipantLayoutRefiner.recomputeIfLowQuality(good)
        let goodPenalty = CDKReactionParticipantLayoutRefiner.qualityPenalty(for: good)
        let refinedPenalty = CDKReactionParticipantLayoutRefiner.qualityPenalty(for: refined)
        XCTAssertLessThanOrEqual(refinedPenalty, goodPenalty + 0.001,
                                 "Participants with good coordinates should stay layout-stable and not degrade during reaction refinement.")

        if let goodBox = good.boundingBox(), let refinedBox = refined.boundingBox() {
            let goodRatio = goodBox.width / max(0.0001, goodBox.height)
            let refinedRatio = refinedBox.width / max(0.0001, refinedBox.height)
            XCTAssertGreaterThanOrEqual(refinedRatio + 0.01, goodRatio,
                                        "Reaction refinement may reorient a good participant, but it should not make it less readable horizontally.")
        }
    }

    func testUnmappedElongatedReactionParticipantsPreferHorizontalOrientation() throws {
        let reaction = try brominationStyleReaction()
        let elongatedParticipants = [reaction.reactants[0], reaction.products[0]]

        for molecule in elongatedParticipants {
            guard let originalBox = molecule.boundingBox() else {
                XCTFail("Expected bounding box for elongated participant.")
                continue
            }
            let refined = CDKReactionParticipantLayoutRefiner.recomputeIfLowQuality(molecule)
            guard let refinedBox = refined.boundingBox() else {
                XCTFail("Expected bounding box after reaction-layout refinement.")
                continue
            }

            let originalRatio = originalBox.width / max(0.0001, originalBox.height)
            let refinedRatio = refinedBox.width / max(0.0001, refinedBox.height)
            XCTAssertGreaterThanOrEqual(refinedRatio + 0.01, originalRatio,
                                        "Reaction-layout refinement should not make an elongated unmapped participant less horizontal.")
            XCTAssertGreaterThan(refinedRatio, 1.05,
                                 "Elongated unmapped reaction participants should prefer a horizontal presentation for reaction readability.")
        }
    }

    func testReactionMapHighlightAppliesConsistentColorAcrossSides() {
        let reaction = mappedSiliconReaction()
        var style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                atomColoringMode: .atomMapHighlight,
                                colorBondsByAtom: true,
                                bondWidth: 2.0,
                                fontSize: 13,
                                padding: 20)
        style.showAtomMapNumbers = false

        let scene = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                style: style,
                                                                canvasRect: CGRect(x: 0, y: 0, width: 980, height: 380),
                                                                zoom: 1.0,
                                                                pan: .zero,
                                                                rotationDegrees: 0)

        let mappedSiliconLabels = scene.labels.filter { $0.text == "Si" }
        XCTAssertEqual(mappedSiliconLabels.count, 2, "Expected mapped Si labels on both sides of the reaction.")
        XCTAssertNotEqual(mappedSiliconLabels[0].color, CDKRenderColor.ink)
        XCTAssertEqual(mappedSiliconLabels[0].color, mappedSiliconLabels[1].color,
                       "Mapped atoms should share the same reaction-map highlight color.")
    }

    func testReactionMapNumbersAreRenderedWhenEnabled() {
        let reaction = mappedSiliconReaction()
        var style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                atomColoringMode: .monochrome,
                                colorBondsByAtom: false,
                                bondWidth: 2.0,
                                fontSize: 13,
                                padding: 20)
        style.showAtomMapNumbers = true

        let scene = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                style: style,
                                                                canvasRect: CGRect(x: 0, y: 0, width: 980, height: 380),
                                                                zoom: 1.0,
                                                                pan: .zero,
                                                                rotationDegrees: 0)

        XCTAssertTrue(scene.labels.contains(where: { $0.text == "Si:7" }),
                      "Expected map-number annotations in reaction labels (CDK-like atom map numbers).")
    }

    func testReactionMapHighlightShowsMappedCarbonsWhenCarbonsHidden() {
        let reaction = mappedCarbonReaction()
        var style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                atomColoringMode: .atomMapHighlight,
                                colorBondsByAtom: true,
                                bondWidth: 2.0,
                                fontSize: 13,
                                padding: 20)
        style.showAtomMapNumbers = false

        let scene = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                style: style,
                                                                canvasRect: CGRect(x: 0, y: 0, width: 980, height: 380),
                                                                zoom: 1.0,
                                                                pan: .zero,
                                                                rotationDegrees: 0)

        let mappedCarbons = scene.labels.filter { $0.text == "C" }
        XCTAssertEqual(mappedCarbons.count, 2, "Expected mapped carbon labels on both reaction sides when map highlight is enabled.")
        XCTAssertNotEqual(mappedCarbons[0].color, CDKRenderColor.ink)
        XCTAssertEqual(mappedCarbons[0].color, mappedCarbons[1].color,
                       "Mapped carbons should share the same atom-map highlight color.")
    }

    func testMappedReactionSmilesParticipantsStayReasonablyCompact() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let rsmi = "[CH3:1][C:2](=[O:3])[CH2:4][S:5][CH2:6][CH:7]=[CH2:8].[CH3:9][CH:10]([NH2:11])[CH:12]([OH:13])[CH:14]=[CH2:15]>>[CH3:1][C:2](=[O:3])[CH2:4][S:5][CH2:6][CH:7]=[CH:14][CH:12]([CH:10]([CH3:9])[NH2:11])[OH:13]"
        let reaction = try parser.parseReactionSmiles(rsmi)
        let participants = reaction.reactants + reaction.agents + reaction.products
        XCTAssertEqual(participants.count, 3)

        for molecule in participants {
            let refined = CDKReactionParticipantLayoutRefiner.recomputeIfLowQuality(molecule)
            let penalty = CDKReactionParticipantLayoutRefiner.qualityPenalty(for: refined)
            XCTAssertLessThan(penalty, 120, "Mapped reaction participant should not retain very poor 2D geometry after refinement.")

            if let box = refined.boundingBox() {
                let minSide = max(0.0001, min(box.width, box.height))
                let aspect = max(box.width, box.height) / minSide
                XCTAssertLessThan(aspect, 7.5, "Mapped reaction participant should avoid extreme stretched aspect ratios.")
                XCTAssertGreaterThan(box.width / max(0.0001, box.height), 0.62,
                                     "Reaction participant should prefer a horizontal presentation for readability.")
            }

            if let correlation = mapNumberXCorrelation(refined) {
                XCTAssertGreaterThan(correlation, 0.15,
                                     "Mapped reaction participant should generally progress left-to-right by atom-map number.")
            }
        }
    }

    func testMappedReactionSmilesReactantsAlignWithProductTemplate() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let rsmi = "[CH3:1][C:2](=[O:3])[CH2:4][S:5][CH2:6][CH:7]=[CH2:8].[CH3:9][CH:10]([NH2:11])[CH:12]([OH:13])[CH:14]=[CH2:15]>>[CH3:1][C:2](=[O:3])[CH2:4][S:5][CH2:6][CH:7]=[CH:14][CH:12]([CH:10]([CH3:9])[NH2:11])[OH:13]"
        let reaction = try parser.parseReactionSmiles(rsmi)

        let refinedProducts = reaction.products.map { CDKReactionParticipantLayoutRefiner.recomputeIfLowQuality($0) }
        let template = mapXTemplate(from: refinedProducts)
        XCTAssertFalse(template.isEmpty)

        for reactant in reaction.reactants {
            let refined = CDKReactionParticipantLayoutRefiner.recomputeIfLowQuality(reactant)
            let before = mapTemplateCorrelation(refined, template: template)
            let aligned = CDKReactionParticipantLayoutRefiner.reorientForMapTemplate(refined, mapXTemplate: template)
            let after = mapTemplateCorrelation(aligned, template: template)

            if let before, let after {
                XCTAssertGreaterThanOrEqual(after + 0.001, before,
                                            "Mapped reactant orientation should not diverge from product-side map template.")
            }
        }
    }

    func testMappedReactionSmilesParticipantsKeepReasonableBondLengthSpread() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let rsmi = "[CH3:1][C:2](=[O:3])[CH2:4][S:5][CH2:6][CH:7]=[CH2:8].[CH3:9][CH:10]([NH2:11])[CH:12]([OH:13])[CH:14]=[CH2:15]>>[CH3:1][C:2](=[O:3])[CH2:4][S:5][CH2:6][CH:7]=[CH:14][CH:12]([CH:10]([CH3:9])[NH2:11])[OH:13]"
        let reaction = try parser.parseReactionSmiles(rsmi)
        let participants = reaction.reactants + reaction.agents + reaction.products
        XCTAssertEqual(participants.count, 3)

        for molecule in participants {
            let refined = CDKReactionParticipantLayoutRefiner.recomputeIfLowQuality(molecule)
            let lengths = refined.bonds.compactMap { bond -> CGFloat? in
                guard let a1 = refined.atom(id: bond.a1), let a2 = refined.atom(id: bond.a2) else { return nil }
                let d = a1.position.distance(to: a2.position)
                return d.isFinite && d > 0.0001 ? d : nil
            }
            guard let minLen = lengths.min(), let maxLen = lengths.max(), !lengths.isEmpty else {
                XCTFail("No bond lengths available for participant.")
                continue
            }
            let spread = maxLen / max(0.0001, minLen)
            XCTAssertLessThan(spread, 2.15, "Mapped reaction participant has excessive bond-length spread (layout likely distorted).")
        }
    }

    func testMappedReactionSmilesSceneStaysWithinCanvasBoundsAtDefaultView() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let rsmi = "[CH3:1][C:2](=[O:3])[CH2:4][S:5][CH2:6][CH:7]=[CH2:8].[CH3:9][CH:10]([NH2:11])[CH:12]([OH:13])[CH:14]=[CH2:15]>>[CH3:1][C:2](=[O:3])[CH2:4][S:5][CH2:6][CH:7]=[CH:14][CH:12]([CH:10]([CH3:9])[NH2:11])[OH:13]"
        let reaction = try parser.parseReactionSmiles(rsmi)
        let canvas = CGRect(x: 0, y: 0, width: 1500, height: 700)
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 1.6,
                                fontSize: 14,
                                padding: 24)
        let scene = CDKMetalReactionDepictionSceneBuilder.build(reaction: reaction,
                                                                style: style,
                                                                canvasRect: canvas,
                                                                zoom: 1.0,
                                                                pan: .zero,
                                                                rotationDegrees: 0)

        let points = scene.bondSegments.flatMap { [$0.from, $0.to] } + scene.labels.map(\.position)
        guard !points.isEmpty else {
            XCTFail("Expected non-empty reaction scene.")
            return
        }
        let minX = points.map(\.x).min() ?? 0
        let maxX = points.map(\.x).max() ?? 0
        let minY = points.map(\.y).min() ?? 0
        let maxY = points.map(\.y).max() ?? 0

        XCTAssertGreaterThanOrEqual(minX, canvas.minX - 12, "Reaction scene should not render heavily off-canvas on the left.")
        XCTAssertLessThanOrEqual(maxX, canvas.maxX + 12, "Reaction scene should not render heavily off-canvas on the right.")
        XCTAssertGreaterThanOrEqual(minY, canvas.minY - 12, "Reaction scene should not render heavily off-canvas on the bottom.")
        XCTAssertLessThanOrEqual(maxY, canvas.maxY + 12, "Reaction scene should not render heavily off-canvas on the top.")
    }

    private func ethane(name: String) -> Molecule {
        Molecule(name: name,
                 atoms: [
                    Atom(id: 1, element: "C", position: CGPoint(x: 0, y: 0)),
                    Atom(id: 2, element: "C", position: CGPoint(x: 1.2, y: 0))
                 ],
                 bonds: [
                    Bond(id: 1, a1: 1, a2: 2, order: .single)
                 ])
    }

    private func moleculeWithHeteroLabel(name: String) -> Molecule {
        Molecule(name: name,
                 atoms: [
                    Atom(id: 1, element: "C", position: CGPoint(x: 0, y: 0)),
                    Atom(id: 2, element: "O", position: CGPoint(x: 1.2, y: 0))
                 ],
                 bonds: [
                    Bond(id: 1, a1: 1, a2: 2, order: .single)
                 ])
    }

    private func mappedSiliconReaction() -> CDKReaction {
        let reactant = Molecule(name: "R",
                                atoms: [
                                    Atom(id: 1, element: "Si", position: CGPoint(x: 0.0, y: 0.0), atomMapNumber: 7),
                                    Atom(id: 2, element: "O", position: CGPoint(x: 1.3, y: 0.0))
                                ],
                                bonds: [
                                    Bond(id: 1, a1: 1, a2: 2, order: .single)
                                ])
        let product = Molecule(name: "P",
                               atoms: [
                                   Atom(id: 10, element: "Si", position: CGPoint(x: 0.0, y: 0.0), atomMapNumber: 7),
                                   Atom(id: 11, element: "N", position: CGPoint(x: 1.3, y: 0.0))
                               ],
                               bonds: [
                                   Bond(id: 10, a1: 10, a2: 11, order: .single)
                               ])
        return CDKReaction(reactants: [reactant], agents: [], products: [product])
    }

    private func mappedCarbonReaction() -> CDKReaction {
        let reactant = Molecule(name: "R",
                                atoms: [
                                    Atom(id: 1, element: "C", position: CGPoint(x: 0.0, y: 0.0), atomMapNumber: 5),
                                    Atom(id: 2, element: "O", position: CGPoint(x: 1.2, y: 0.0))
                                ],
                                bonds: [
                                    Bond(id: 1, a1: 1, a2: 2, order: .single)
                                ])
        let product = Molecule(name: "P",
                               atoms: [
                                   Atom(id: 10, element: "C", position: CGPoint(x: 0.0, y: 0.0), atomMapNumber: 5),
                                   Atom(id: 11, element: "N", position: CGPoint(x: 1.2, y: 0.0))
                               ],
                               bonds: [
                                   Bond(id: 10, a1: 10, a2: 11, order: .single)
                               ])
        return CDKReaction(reactants: [reactant], agents: [], products: [product])
    }

    private func brominationStyleReaction() throws -> CDKReaction {
        let sorbate = try parseSmiles("COC(=O)C=CC=CC")
        let nbs = try parseSmiles("BrN1C(=O)CCC1=O")
        let peroxide = try parseSmiles("O=C(OOC(=O)c1ccccc1)c1ccccc1")
        let benzene = try parseSmiles("c1ccccc1")
        let product = try parseSmiles("COC(=O)C=CC(Br)CC")
        return CDKReaction(reactants: [sorbate, nbs, peroxide],
                           agents: [benzene],
                           products: [product])
    }

    private func parseSmiles(_ smiles: String) throws -> Molecule {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        return try parser.parseSmiles(smiles)
    }

    private func mapNumberXCorrelation(_ molecule: Molecule) -> CGFloat? {
        let mapped = molecule.atoms.compactMap { atom -> (CGFloat, CGFloat)? in
            guard let map = atom.atomMapNumber, map > 0 else { return nil }
            return (CGFloat(map), atom.position.x)
        }
        guard mapped.count >= 3 else { return nil }

        let n = CGFloat(mapped.count)
        let meanMap = mapped.reduce(0) { $0 + $1.0 } / n
        let meanX = mapped.reduce(0) { $0 + $1.1 } / n
        var covariance: CGFloat = 0
        var varMap: CGFloat = 0
        var varX: CGFloat = 0
        for pair in mapped {
            let dm = pair.0 - meanMap
            let dx = pair.1 - meanX
            covariance += dm * dx
            varMap += dm * dm
            varX += dx * dx
        }
        let denom = sqrt(max(0.000_000_1, varMap * varX))
        guard denom > 0 else { return nil }
        return covariance / denom
    }

    private func mapXTemplate(from molecules: [Molecule]) -> [Int: CGFloat] {
        var values: [Int: [CGFloat]] = [:]
        for molecule in molecules {
            guard let box = molecule.boundingBox() else { continue }
            let denom = max(0.0001, box.width)
            for atom in molecule.atoms {
                guard let map = atom.atomMapNumber, map > 0 else { continue }
                values[map, default: []].append((atom.position.x - box.minX) / denom)
            }
        }

        var out: [Int: CGFloat] = [:]
        for (map, xs) in values where !xs.isEmpty {
            out[map] = xs.reduce(0, +) / CGFloat(xs.count)
        }
        return out
    }

    private func mapTemplateCorrelation(_ molecule: Molecule, template: [Int: CGFloat]) -> CGFloat? {
        guard let box = molecule.boundingBox() else { return nil }
        let denom = max(0.0001, box.width)
        let mapped = molecule.atoms.compactMap { atom -> (CGFloat, CGFloat)? in
            guard let map = atom.atomMapNumber, let target = template[map] else { return nil }
            return ((atom.position.x - box.minX) / denom, target)
        }
        guard mapped.count >= 3 else { return nil }

        let n = CGFloat(mapped.count)
        let meanX = mapped.reduce(0) { $0 + $1.0 } / n
        let meanY = mapped.reduce(0) { $0 + $1.1 } / n
        var covariance: CGFloat = 0
        var varX: CGFloat = 0
        var varY: CGFloat = 0
        for pair in mapped {
            let dx = pair.0 - meanX
            let dy = pair.1 - meanY
            covariance += dx * dy
            varX += dx * dx
            varY += dy * dy
        }
        let denomCorr = sqrt(max(0.000_000_1, varX * varY))
        guard denomCorr > 0 else { return nil }
        return covariance / denomCorr
    }

    private func glowSegmentCount(in scene: CDKMetalDepictionScene, minX: CGFloat, maxX: CGFloat) -> Int {
        scene.bondSegments.filter { segment in
            let midpointX = (segment.from.x + segment.to.x) * 0.5
            return midpointX >= minX &&
                   midpointX <= maxX &&
                   isOuterGlowSegment(segment)
        }.count
    }

    private func isOuterGlowSegment(_ segment: CDKMetalDepictionScene.LineSegment) -> Bool {
        let glow = CDKRenderColor.outerGlowHighlight
        return abs(segment.color.red - glow.red) < 0.000_1 &&
               abs(segment.color.green - glow.green) < 0.000_1 &&
               abs(segment.color.blue - glow.blue) < 0.000_1
    }

    private func viewportTransformed(_ point: CGPoint,
                                     canvasRect: CGRect,
                                     zoom: CGFloat,
                                     pan: CGSize,
                                     rotationDegrees: CGFloat) -> CGPoint {
        let viewportCenter = CGPoint(x: canvasRect.midX, y: canvasRect.midY)
        let radians = rotationDegrees * (.pi / 180)
        let c = cos(radians)
        let s = sin(radians)
        let rx = point.x - viewportCenter.x
        let ry = point.y - viewportCenter.y
        let rotated = CGPoint(x: (rx * c) - (ry * s) + viewportCenter.x,
                              y: (rx * s) + (ry * c) + viewportCenter.y)
        return CGPoint(x: ((rotated.x - viewportCenter.x) * zoom) + viewportCenter.x + pan.width,
                       y: ((rotated.y - viewportCenter.y) * zoom) + viewportCenter.y + pan.height)
    }
}
