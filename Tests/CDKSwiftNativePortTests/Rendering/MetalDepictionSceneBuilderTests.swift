import Foundation
#if canImport(CoreGraphics)
import CoreGraphics
#endif
import XCTest
@testable import CDKSwiftNativePort

final class MetalDepictionSceneBuilderTests: XCTestCase {
    private let smilesParser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    func testAlternatingRingRetainsExplicitDoubleBondDepiction() throws {
        let molecule = alternatingSixRing
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 12.0,
                                padding: 24.0)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: CGRect(x: 0, y: 0, width: 480, height: 360),
                                                        zoom: 1.0,
                                                        pan: .zero)

        XCTAssertTrue(scene.gridSegments.isEmpty,
                      "Live metal depictions should render on a plain background without the old decorative grid.")

        let labelPositionByAtomID = Dictionary(uniqueKeysWithValues: scene.labels.map { ($0.id, $0.position) })
        XCTAssertEqual(labelPositionByAtomID.count, molecule.atomCount)

        for bond in molecule.bonds where bond.order == .double {
            let from = try XCTUnwrap(labelPositionByAtomID[bond.a1], "Missing label position for atom \(bond.a1)")
            let to = try XCTUnwrap(labelPositionByAtomID[bond.a2], "Missing label position for atom \(bond.a2)")

            let supportingSegments = countSupportingSegments(forBondFrom: from, to: to, segments: scene.bondSegments)
            XCTAssertGreaterThanOrEqual(
                supportingSegments,
                2,
                "Expected ring double bond \(bond.id) to render with at least two strokes, got \(supportingSegments)."
            )
        }
    }

    func testStereoWedgesUsePronouncedGeometryInMetalScene() {
        let molecule = Molecule(
            name: "StereoWedges",
            atoms: [
                Atom(id: 1, element: "C", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "C", position: CGPoint(x: 1.1, y: 0.0)),
                Atom(id: 3, element: "C", position: CGPoint(x: 0.0, y: 1.0)),
                Atom(id: 4, element: "C", position: CGPoint(x: -1.2, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single, stereo: .up),
                Bond(id: 2, a1: 1, a2: 3, order: .single, stereo: .down),
                Bond(id: 3, a1: 1, a2: 4, order: .single, stereo: .none)
            ]
        )
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 12.0,
                                padding: 24.0)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: CGRect(x: 0, y: 0, width: 480, height: 360),
                                                        zoom: 1.0,
                                                        pan: .zero)

        let emphasizedWedgeSegments = scene.bondSegments.filter { $0.opacity >= 0.969 }
        XCTAssertGreaterThanOrEqual(
            emphasizedWedgeSegments.count,
            16,
            "Expected stereo up wedge to render as a dense filled wedge, but only got \(emphasizedWedgeSegments.count) emphasized slices."
        )
    }

    func testStandardAtomLabelsDoNotRequestOpaqueBackgrounds() {
        let molecule = Molecule(
            name: "LabelBackgrounds",
            atoms: [
                Atom(id: 1, element: "O", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "C", position: CGPoint(x: 1.2, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single)
            ]
        )
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 12.0,
                                padding: 24.0)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: CGRect(x: 0, y: 0, width: 360, height: 240),
                                                        zoom: 1.0,
                                                        pan: .zero)

        XCTAssertFalse(scene.labels.contains(where: \.drawsBackground),
                       "Upstream-style atom labels should render without opaque backdrop rectangles.")
    }

    func testStereoTerminalLabelsKeepExtraClearanceFromWedgeTips() throws {
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 24.0,
                                padding: 24.0)

        let solidWedge = Molecule(
            name: "TerminalSolidWedge",
            atoms: [
                Atom(id: 1, element: "C", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "Cl", position: CGPoint(x: 0.0, y: 1.35))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single, stereo: .up)
            ]
        )
        let solidScene = CDKMetalDepictionSceneBuilder.build(molecule: solidWedge,
                                                             style: style,
                                                             canvasRect: CGRect(x: 0, y: 0, width: 420, height: 320),
                                                             zoom: 1.0,
                                                             pan: .zero)
        let solidLabel = try XCTUnwrap(solidScene.labels.first(where: { $0.id == 2 })?.position)
        let solidNearest = nearestEndpointDistance(in: solidScene.bondSegments, to: solidLabel)
        XCTAssertGreaterThan(solidNearest, 7.0,
                             "Solid stereo wedges should leave extra white space before terminal labels.")

        let hashedWedge = Molecule(
            name: "TerminalHashedWedge",
            atoms: [
                Atom(id: 1, element: "C", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "Cl", position: CGPoint(x: 0.0, y: 1.35))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single, stereo: .down)
            ]
        )
        let hashedScene = CDKMetalDepictionSceneBuilder.build(molecule: hashedWedge,
                                                              style: style,
                                                              canvasRect: CGRect(x: 0, y: 0, width: 420, height: 320),
                                                              zoom: 1.0,
                                                              pan: .zero)
        let hashedLabel = try XCTUnwrap(hashedScene.labels.first(where: { $0.id == 2 })?.position)
        let hashedNearest = nearestEndpointDistance(in: hashedScene.bondSegments, to: hashedLabel)
        XCTAssertGreaterThan(hashedNearest, 7.0,
                             "Hashed stereo wedges should leave extra white space before terminal labels.")
    }

    func testReversedStereoWedgesBroadenAwayFromStereocenter() throws {
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 12.0,
                                padding: 24.0)

        let solidWedge = Molecule(
            name: "ReversedSolidWedge",
            atoms: [
                Atom(id: 1, element: "N", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "Cl", position: CGPoint(x: 1.35, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 2, a2: 1, order: .single, stereo: .upReversed)
            ]
        )
        let solidScene = CDKMetalDepictionSceneBuilder.build(molecule: solidWedge,
                                                             style: style,
                                                             canvasRect: CGRect(x: 0, y: 0, width: 420, height: 240),
                                                             zoom: 1.0,
                                                             pan: .zero)
        let solidStart = try XCTUnwrap(solidScene.labels.first(where: { $0.id == 1 })?.position)
        let solidEnd = try XCTUnwrap(solidScene.labels.first(where: { $0.id == 2 })?.position)
        let solidNearStart = averageStereoStripeLength(in: solidScene.bondSegments,
                                                       from: solidStart,
                                                       to: solidEnd,
                                                       tRange: 0.12...0.38)
        let solidNearEnd = averageStereoStripeLength(in: solidScene.bondSegments,
                                                     from: solidStart,
                                                     to: solidEnd,
                                                     tRange: 0.62...0.88)
        XCTAssertGreaterThan(solidNearEnd,
                             solidNearStart * 1.7,
                             "Expected reversed solid wedge to broaden away from the stereocenter.")

        let hashedWedge = Molecule(
            name: "ReversedHashedWedge",
            atoms: [
                Atom(id: 1, element: "N", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "Cl", position: CGPoint(x: 1.35, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 2, a2: 1, order: .single, stereo: .downReversed)
            ]
        )
        let hashedScene = CDKMetalDepictionSceneBuilder.build(molecule: hashedWedge,
                                                              style: style,
                                                              canvasRect: CGRect(x: 0, y: 0, width: 420, height: 240),
                                                              zoom: 1.0,
                                                              pan: .zero)
        let hashedStart = try XCTUnwrap(hashedScene.labels.first(where: { $0.id == 1 })?.position)
        let hashedEnd = try XCTUnwrap(hashedScene.labels.first(where: { $0.id == 2 })?.position)
        let hashedNearStart = averageStereoStripeLength(in: hashedScene.bondSegments,
                                                        from: hashedStart,
                                                        to: hashedEnd,
                                                        tRange: 0.12...0.38)
        let hashedNearEnd = averageStereoStripeLength(in: hashedScene.bondSegments,
                                                      from: hashedStart,
                                                      to: hashedEnd,
                                                      tRange: 0.62...0.88)
        XCTAssertGreaterThan(hashedNearEnd,
                             hashedNearStart * 1.7,
                             "Expected reversed hashed wedge to broaden away from the stereocenter.")
    }

    func testBondSegmentsAreClippedAwayFromLabelCenters() {
        let molecule = Molecule(
            name: "OO",
            atoms: [
                Atom(id: 1, element: "O", position: CGPoint(x: 0, y: 0)),
                Atom(id: 2, element: "O", position: CGPoint(x: 1.2, y: 0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single)
            ]
        )
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 24.0,
                                padding: 24.0)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: CGRect(x: 0, y: 0, width: 520, height: 360),
                                                        zoom: 1.0,
                                                        pan: .zero)

        let center1 = try? XCTUnwrap(scene.labels.first(where: { $0.id == 1 })?.position)
        let center2 = try? XCTUnwrap(scene.labels.first(where: { $0.id == 2 })?.position)
        let segment = scene.bondSegments.max(by: { hypot($0.to.x - $0.from.x, $0.to.y - $0.from.y) < hypot($1.to.x - $1.from.x, $1.to.y - $1.from.y) })

        guard let c1 = center1, let c2 = center2, let segment else {
            XCTFail("Expected labels and at least one bond segment.")
            return
        }

        let d11 = hypot(segment.from.x - c1.x, segment.from.y - c1.y)
        let d12 = hypot(segment.to.x - c1.x, segment.to.y - c1.y)
        let d21 = hypot(segment.from.x - c2.x, segment.from.y - c2.y)
        let d22 = hypot(segment.to.x - c2.x, segment.to.y - c2.y)

        XCTAssertGreaterThan(min(d11, d12), 4.0, "Bond endpoint should be clipped away from label 1 center.")
        XCTAssertGreaterThan(min(d21, d22), 4.0, "Bond endpoint should be clipped away from label 2 center.")
    }

    func testSuppressesSimpleExplicitHydrogensByDefault() {
        let molecule = Molecule(
            name: "AlcoholFragment",
            atoms: [
                Atom(id: 1, element: "O", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "H", position: CGPoint(x: 0.0, y: 1.0)),
                Atom(id: 3, element: "C", position: CGPoint(x: -1.1, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single),
                Bond(id: 2, a1: 1, a2: 3, order: .single)
            ]
        )
        let style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: true,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14.0,
                                padding: 24.0)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: CGRect(x: 0, y: 0, width: 480, height: 320),
                                                        zoom: 1.0,
                                                        pan: .zero)

        XCTAssertNil(scene.labels.first(where: { $0.id == 2 || $0.text == "H" }))
        XCTAssertNotNil(scene.labels.first(where: { $0.id == 1 && $0.text.hasPrefix("OH") }))
    }

    func testMetalSceneColorsCxHighlightedAtomsAndBondsWithCDKSelectionColor() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("CN |ha:1,hb:0|")
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14.0,
                                padding: 24.0)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: CGRect(x: 0, y: 0, width: 480, height: 320),
                                                        zoom: 1.0,
                                                        pan: .zero)

        let carbonLabel = try XCTUnwrap(scene.labels.first(where: { $0.id == molecule.atoms[0].id }))
        let nitrogenLabel = try XCTUnwrap(scene.labels.first(where: { $0.id == molecule.atoms[1].id }))

        XCTAssertEqual(carbonLabel.color, .ink)
        XCTAssertEqual(nitrogenLabel.color, .selectionHighlight)
        XCTAssertFalse(scene.bondSegments.isEmpty)
        XCTAssertTrue(scene.bondSegments.allSatisfy { $0.color == .selectionHighlight },
                      "Expected selected CX bond to use CDK's selection highlight color.")
    }

    func testMetalSceneShowsDisconnectedSelectedCarbonLabelLikeCDKSelectionVisibility() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("CC |ha:0|")
        let style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14.0,
                                padding: 24.0)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: CGRect(x: 0, y: 0, width: 480, height: 320),
                                                        zoom: 1.0,
                                                        pan: .zero,
                                                        includeTerminalCarbonLabelsWhenCarbonsHidden: false)

        let selectedLabel = try XCTUnwrap(scene.labels.first(where: { $0.id == molecule.atoms[0].id }))
        XCTAssertEqual(selectedLabel.text, "C")
        XCTAssertEqual(selectedLabel.color, .selectionHighlight)
    }

    func testMetalSceneKeepsSelectedCarbonLabelHiddenWhenAdjacentBondIsSelected() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("CC |ha:0,hb:0|")
        let style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14.0,
                                padding: 24.0)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: CGRect(x: 0, y: 0, width: 480, height: 320),
                                                        zoom: 1.0,
                                                        pan: .zero,
                                                        includeTerminalCarbonLabelsWhenCarbonsHidden: false)

        XCTAssertFalse(scene.labels.contains(where: { $0.id == molecule.atoms[0].id }),
                       "Expected selected carbon label to remain hidden when a connected selected bond is present.")
    }

    func testMetalSceneUsesGlowHighlightOverlayInsteadOfSelectionColorInOuterGlowMode() throws {
        let molecule = try smilesParser.parseSmiles("CN |ha:1,hb:0|")
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                atomColoringMode: .monochrome,
                                highlightStyle: .outerGlow,
                                bondWidth: 2.0,
                                fontSize: 14.0,
                                padding: 24.0)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: CGRect(x: 0, y: 0, width: 480, height: 320),
                                                        zoom: 1.0,
                                                        pan: .zero)

        let baseNitrogen = try XCTUnwrap(scene.labels.first(where: { $0.id == 2 }))
        XCTAssertEqual(baseNitrogen.color, CDKRenderColor.ink)
        let glowNitrogen = try XCTUnwrap(scene.labels.first(where: { $0.usesGlowOverlay && $0.text == "N" }))
        XCTAssertEqual(glowNitrogen.color, CDKRenderColor.outerGlowHighlight.withAlpha(0.92))
        XCTAssertTrue(glowNitrogen.suppressesMatchedBackground)
        XCTAssertTrue(scene.bondSegments.contains(where: { $0.color == .outerGlowHighlight }))
        XCTAssertTrue(scene.bondSegments.contains(where: { $0.color == .ink }))
    }

    func testMetalScenePreservesBaseLabelBackgroundInOuterGlowWhiteEdgeMode() throws {
        let molecule = try smilesParser.parseSmiles("CN |ha:1,hb:0|")
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                atomColoringMode: .monochrome,
                                highlightStyle: .outerGlowWhiteEdge,
                                bondWidth: 2.0,
                                fontSize: 14.0,
                                padding: 24.0)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: CGRect(x: 0, y: 0, width: 480, height: 320),
                                                        zoom: 1.0,
                                                        pan: .zero)

        let glowNitrogen = try XCTUnwrap(scene.labels.first(where: { $0.usesGlowOverlay && $0.text == "N" }))
        XCTAssertFalse(glowNitrogen.suppressesMatchedBackground)
    }

    func testMetalSceneRendersQueryAtomListLabels() {
        let molecule = Molecule(name: "QueryAtom",
                                atoms: [
                                    Atom(id: 1,
                                         element: "L",
                                         position: CGPoint(x: 0, y: 0),
                                         queryType: .anyAtom,
                                         atomList: ["C", "N"]),
                                    Atom(id: 2, element: "O", position: CGPoint(x: 1.2, y: 0))
                                ],
                                bonds: [
                                    Bond(id: 1, a1: 1, a2: 2, order: .single)
                                ])
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14.0,
                                padding: 24.0)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: CGRect(x: 0, y: 0, width: 520, height: 320),
                                                        zoom: 1.0,
                                                        pan: .zero)

        XCTAssertTrue(scene.labels.contains(where: { $0.id == 1 && $0.text == "[C,N]" }))
    }

    func testMetalSceneRendersQueryBondSemantics() {
        let molecule = Molecule(name: "QueryBond",
                                atoms: [
                                    Atom(id: 1, element: "C", position: CGPoint(x: 0, y: 0)),
                                    Atom(id: 2, element: "C", position: CGPoint(x: 1.4, y: 0))
                                ],
                                bonds: [
                                    Bond(id: 1, a1: 1, a2: 2, order: .single, queryType: .singleOrDouble)
                                ])
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14.0,
                                padding: 24.0)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: CGRect(x: 0, y: 0, width: 520, height: 320),
                                                        zoom: 1.0,
                                                        pan: .zero)

        XCTAssertGreaterThanOrEqual(scene.bondSegments.count, 2)
        let widths = scene.bondSegments.map(\.width)
        XCTAssertTrue(widths.contains(where: { $0 < style.bondWidth }),
                      "Expected query secondary line to render narrower than the primary bond stroke.")
    }

    func testMetalSceneRendersPolymerSgroupBracketAnnotations() {
        var molecule = Molecule(name: "PolymerBracket",
                                atoms: [
                                    Atom(id: 1, element: "C", position: CGPoint(x: 0, y: 0)),
                                    Atom(id: 2, element: "C", position: CGPoint(x: 1.4, y: 0))
                                ],
                                bonds: [
                                    Bond(id: 1, a1: 1, a2: 2, order: .single)
                                ])
        molecule.sgroups = [
            MoleculeSgroup(kind: .polymer,
                           keyword: "COP",
                           atomIDs: [1, 2],
                           subtype: "RAN")
        ]
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14.0,
                                padding: 24.0)

        let plainScene = CDKMetalDepictionSceneBuilder.build(molecule: Molecule(name: molecule.name,
                                                                                atoms: molecule.atoms,
                                                                                bonds: molecule.bonds),
                                                             style: style,
                                                             canvasRect: CGRect(x: 0, y: 0, width: 520, height: 320),
                                                             zoom: 1.0,
                                                             pan: .zero)
        let annotatedScene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                                 style: style,
                                                                 canvasRect: CGRect(x: 0, y: 0, width: 520, height: 320),
                                                                 zoom: 1.0,
                                                                 pan: .zero)

        XCTAssertTrue(annotatedScene.labels.contains(where: { $0.text == "ran" }))
        XCTAssertGreaterThan(annotatedScene.bondSegments.count, plainScene.bondSegments.count)
    }

    func testGenericSgroupRendererSkipsRoundBracketLinkNodeRepeatUnits() throws {
        let molecule = try smilesParser.parseSmiles("C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|")
        let positionsByAtomID = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0.position) })

        let annotations = CDKSgroupRendering.bracketAnnotations(molecule: molecule,
                                                                positionsByAtomID: positionsByAtomID,
                                                                fontSize: 16,
                                                                bondWidth: 2.0)

        XCTAssertTrue(annotations.isEmpty,
                      "Round-bracket LN repeat units should stay on the dedicated Markush annotation path and not also emit generic square-bracket Sgroup graphics.")
    }

    func testGenericSgroupRendererStillDrawsSquareBracketRepeatUnits() {
        var molecule = Molecule(
            name: "SquareBracketSRU",
            atoms: [
                Atom(id: 1, element: "C", position: CGPoint(x: 0, y: 0)),
                Atom(id: 2, element: "C", position: CGPoint(x: 1.4, y: 0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single)
            ]
        )
        molecule.sgroups = [
            MoleculeSgroup(kind: .structureRepeatUnit,
                           keyword: "SRU",
                           atomIDs: [1, 2],
                           crossingBondIDs: [1],
                           subscriptText: "n",
                           roundBrackets: false)
        ]
        let positionsByAtomID = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0.position) })

        let annotations = CDKSgroupRendering.bracketAnnotations(molecule: molecule,
                                                                positionsByAtomID: positionsByAtomID,
                                                                fontSize: 16,
                                                                bondWidth: 2.0)

        XCTAssertEqual(annotations.count, 1)
        XCTAssertFalse(annotations[0].segments.isEmpty)
        XCTAssertEqual(annotations[0].subscriptText, "n")
    }

    func testMarkushLinkRendererSkipsSquareBracketRepeatUnits() {
        var molecule = Molecule(
            name: "SquareBracketSRU",
            atoms: [
                Atom(id: 1, element: "C", position: CGPoint(x: 0, y: 0)),
                Atom(id: 2, element: "C", position: CGPoint(x: 1.4, y: 0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single)
            ]
        )
        molecule.sgroups = [
            MoleculeSgroup(kind: .structureRepeatUnit,
                           keyword: "SRU",
                           atomIDs: [1, 2],
                           crossingBondIDs: [1],
                           subscriptText: "n",
                           roundBrackets: false)
        ]
        let positionsByAtomID = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0.position) })

        let annotations = CDKMarkushRendering.linkAnnotations(molecule: molecule,
                                                              positionsByAtomID: positionsByAtomID,
                                                              fontSize: 16)

        XCTAssertTrue(annotations.isEmpty)
    }

    func testMetalSceneContractsDisplayShortcutsForAbbreviations() {
        var molecule = Molecule(name: "Abbreviation",
                                atoms: [
                                    Atom(id: 1, element: "C", position: CGPoint(x: 0, y: 0)),
                                    Atom(id: 2, element: "C", position: CGPoint(x: 1.0, y: 0)),
                                    Atom(id: 3, element: "C", position: CGPoint(x: 2.2, y: 0))
                                ],
                                bonds: [
                                    Bond(id: 1, a1: 1, a2: 2, order: .single),
                                    Bond(id: 2, a1: 2, a2: 3, order: .single)
                                ])
        molecule.sgroups = [
            MoleculeSgroup(kind: .generic,
                           keyword: "SUP",
                           atomIDs: [1, 2],
                           crossingBondIDs: [2],
                           subscriptText: "Et")
        ]
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14.0,
                                padding: 24.0)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: CGRect(x: 0, y: 0, width: 520, height: 320),
                                                        zoom: 1.0,
                                                        pan: .zero)

        XCTAssertTrue(scene.labels.contains(where: { $0.text == "Et" }))
        XCTAssertFalse(scene.labels.contains(where: { $0.id == 1 && $0.text == "C" }))
    }

    func testMetalSceneIncludesMarkushBackgroundBoxAndLinkNodeLabels() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|")
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 16.0,
                                padding: 24.0)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: CGRect(x: 0, y: 0, width: 720, height: 520),
                                                        zoom: 1.0,
                                                        pan: .zero)

        XCTAssertEqual(scene.backgroundBoxes.count, 1)
        XCTAssertTrue(scene.labels.contains(where: { $0.text == "R¹ is" }))
        XCTAssertTrue(scene.labels.contains(where: { $0.text == "1-3" }))
        XCTAssertTrue(scene.labels.contains(where: { $0.text == "R¹" && $0.italicized }))
        XCTAssertTrue(scene.labels.contains(where: { $0.text == "R¹ is" && $0.italicized }))

        let rootR = try XCTUnwrap(scene.labels.first(where: { $0.text == "R¹" }))
        let legendPrefix = try XCTUnwrap(scene.labels.first(where: { $0.text == "R¹ is" }))
        XCTAssertEqual(legendPrefix.fontSize,
                       rootR.fontSize,
                       accuracy: 0.0001,
                       "Expected Markush legend prefix to use the same dynamic font size as the scaffold R-group label.")

        let markushLabelIDs = Set(molecule.atoms.filter { $0.rGroupMembership != nil }.map(\.id))
        XCTAssertTrue(scene.labels.filter { markushLabelIDs.contains($0.id) }.allSatisfy { !$0.drawsBackground })
        XCTAssertTrue(scene.labels.filter { ["R¹ is", "(", ")", "1-3"].contains($0.text) }.allSatisfy { !$0.drawsBackground })
    }

    func testMarkushBackgroundBoxContainsVisibleLabelsAndAttachmentWaves() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|")
        let style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: true,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 16.0,
                                padding: 24.0)

        let positionsByAtomID = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0.position) })
        let box = try XCTUnwrap(
            CDKMarkushRendering.rGroupBoxes(molecule: molecule,
                                            positionsByAtomID: positionsByAtomID,
                                            style: style,
                                            padding: max(1.4, style.bondWidth * 0.55),
                                            fontSize: style.fontSize,
                                            bondWidth: style.bondWidth)
                .first
        )

        for atom in molecule.atoms where atom.rGroupMembership == "R1" && atom.attachmentPoint == nil {
            let degree = molecule.neighbors(of: atom.id).count
            guard CDKLabelText.shouldDrawLabel(atom: atom, degree: degree, style: style),
                  let position = positionsByAtomID[atom.id] else {
                continue
            }

            let implicitH = style.showImplicitHydrogens ? molecule.implicitHydrogenCount(for: atom.id) : 0
            let text = CDKLabelText.build(atom: atom, style: style, implicitHydrogenCount: implicitH)
            let offset = CDKLabelText.centerOffset(atom: atom,
                                                   style: style,
                                                   implicitHydrogenCount: implicitH,
                                                   fontSize: style.fontSize)
            let center = position.offsetBy(dx: offset.dx, dy: offset.dy)
            let rect = CDKLabelClipping.makeLabelRect(center: center,
                                                      estimatedTextSize: CDKLabelClipping.estimateLabelSize(text: text,
                                                                                                            fontSize: style.fontSize))
            XCTAssertTrue(box.boxRect.contains(rect),
                          "Expected Markush box to contain visible member label \(text).")
        }

        for bond in molecule.bonds {
            guard let atom1 = molecule.atom(id: bond.a1),
                  let atom2 = molecule.atom(id: bond.a2),
                  atom1.rGroupMembership == "R1",
                  atom2.rGroupMembership == "R1",
                  let p1 = positionsByAtomID[bond.a1],
                  let p2 = positionsByAtomID[bond.a2] else {
                continue
            }

            if atom1.attachmentPoint != nil {
                for segment in CDKMarkushRendering.attachmentPointSegments(center: p1, other: p2, style: style) {
                    XCTAssertTrue(box.boxRect.insetBy(dx: -1, dy: -1).contains(segment.from))
                    XCTAssertTrue(box.boxRect.insetBy(dx: -1, dy: -1).contains(segment.to))
                }
            }

            if atom2.attachmentPoint != nil {
                for segment in CDKMarkushRendering.attachmentPointSegments(center: p2, other: p1, style: style) {
                    XCTAssertTrue(box.boxRect.insetBy(dx: -1, dy: -1).contains(segment.from))
                    XCTAssertTrue(box.boxRect.insetBy(dx: -1, dy: -1).contains(segment.to))
                }
            }
        }
    }

    func testMetalSceneDrawsAttachmentPointWaveSegmentsForMarkushDefinitions() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C* |$;R1$,RG:_R1={Cl* |$;_AP1$|}|")
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14.0,
                                padding: 24.0)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: CGRect(x: 0, y: 0, width: 640, height: 420),
                                                        zoom: 1.0,
                                                        pan: .zero)

        XCTAssertGreaterThan(scene.bondSegments.count, molecule.bondCount + 4)
        XCTAssertFalse(scene.labels.contains(where: { $0.text.contains("_AP") }))
    }

    func testMetalSceneHidesMarkushTerminalCarbonLabelsWhenCarbonsAreHidden() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|")
        let style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: true,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 16.0,
                                padding: 24.0)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: CGRect(x: 0, y: 0, width: 720, height: 520),
                                                        zoom: 1.0,
                                                        pan: .zero)

        let r1AtomIDs = Set(molecule.atoms.filter { $0.rGroupMembership == "R1" }.map(\.id))
        XCTAssertFalse(scene.labels.contains(where: { r1AtomIDs.contains($0.id) && $0.text == "C" }))
        XCTAssertTrue(scene.labels.contains(where: { $0.text == "O" }))
        XCTAssertTrue(scene.labels.contains(where: { $0.text == "Cl" }))
        XCTAssertTrue(scene.labels.contains(where: { $0.text == "N" }))
    }

    func testCanKeepExplicitHydrogensWhenEnabled() {
        let molecule = Molecule(
            name: "AlcoholFragment",
            atoms: [
                Atom(id: 1, element: "O", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "H", position: CGPoint(x: 0.0, y: 1.0)),
                Atom(id: 3, element: "C", position: CGPoint(x: -1.1, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single),
                Bond(id: 2, a1: 1, a2: 3, order: .single)
            ]
        )
        var style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: true,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14.0,
                                padding: 24.0)
        style.showExplicitHydrogens = true

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: CGRect(x: 0, y: 0, width: 480, height: 320),
                                                        zoom: 1.0,
                                                        pan: .zero)

        XCTAssertNotNil(scene.labels.first(where: { $0.id == 2 && $0.text == "H" }))
    }

    func testSuffixHydrogenKeepsAtomSymbolNearAnchor() {
        let molecule = Molecule(
            name: "AlcoholFragment",
            atoms: [
                Atom(id: 1, element: "O", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "H", position: CGPoint(x: 0.0, y: 1.0)),
                Atom(id: 3, element: "C", position: CGPoint(x: -1.1, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single),
                Bond(id: 2, a1: 1, a2: 3, order: .single)
            ]
        )

        let styleWithHydrogen = RenderStyle(showCarbons: false,
                                            showImplicitHydrogens: true,
                                            showAtomIDs: false,
                                            bondWidth: 2.0,
                                            fontSize: 14.0,
                                            padding: 24.0)
        let styleWithoutHydrogen = RenderStyle(showCarbons: false,
                                               showImplicitHydrogens: false,
                                               showAtomIDs: false,
                                               bondWidth: 2.0,
                                               fontSize: 14.0,
                                               padding: 24.0)

        let sceneWithHydrogen = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                                     style: styleWithHydrogen,
                                                                     canvasRect: CGRect(x: 0, y: 0, width: 480, height: 320),
                                                                     zoom: 1.0,
                                                                     pan: .zero)
        let sceneWithoutHydrogen = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                                        style: styleWithoutHydrogen,
                                                                        canvasRect: CGRect(x: 0, y: 0, width: 480, height: 320),
                                                                        zoom: 1.0,
                                                                        pan: .zero)

        let oxygenWithHydrogen = try? XCTUnwrap(sceneWithHydrogen.labels.first(where: { $0.id == 1 && $0.text.hasPrefix("OH") }))
        let oxygenWithoutHydrogen = try? XCTUnwrap(sceneWithoutHydrogen.labels.first(where: { $0.id == 1 && $0.text == "O" }))

        XCTAssertNotNil(oxygenWithHydrogen)
        XCTAssertNotNil(oxygenWithoutHydrogen)
        if let oxygenWithHydrogen, let oxygenWithoutHydrogen {
            XCTAssertGreaterThan(oxygenWithHydrogen.position.x, oxygenWithoutHydrogen.position.x + 1.0)
        }
    }

    func testRingOxygenLabelRemainsAtAtomAnchor() {
        let molecule = Molecule(
            name: "Oxacycle",
            atoms: [
                Atom(id: 1, element: "O", position: CGPoint(x: 1.0, y: 0.0)),
                Atom(id: 2, element: "C", position: CGPoint(x: 0.5, y: 0.8660254)),
                Atom(id: 3, element: "C", position: CGPoint(x: -0.5, y: 0.8660254)),
                Atom(id: 4, element: "C", position: CGPoint(x: -1.0, y: 0.0)),
                Atom(id: 5, element: "C", position: CGPoint(x: -0.5, y: -0.8660254)),
                Atom(id: 6, element: "C", position: CGPoint(x: 0.5, y: -0.8660254))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single),
                Bond(id: 2, a1: 2, a2: 3, order: .single),
                Bond(id: 3, a1: 3, a2: 4, order: .single),
                Bond(id: 4, a1: 4, a2: 5, order: .single),
                Bond(id: 5, a1: 5, a2: 6, order: .single),
                Bond(id: 6, a1: 6, a2: 1, order: .single)
            ]
        )
        let style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14.0,
                                padding: 24.0)
        let canvasRect = CGRect(x: 0, y: 0, width: 480, height: 320)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: canvasRect,
                                                        zoom: 1.0,
                                                        pan: .zero)

        guard let oxygenLabel = scene.labels.first(where: { $0.id == 1 }) else {
            XCTFail("Expected oxygen label for atom 1.")
            return
        }

        guard let box = molecule.boundingBox(),
              let oxygen = molecule.atoms.first(where: { $0.id == 1 }) else {
            XCTFail("Expected oxygen atom and molecule bounds.")
            return
        }

        let pad = style.padding
        let available = CGRect(x: canvasRect.minX + pad,
                               y: canvasRect.minY + pad,
                               width: max(1, canvasRect.width - 2 * pad),
                               height: max(1, canvasRect.height - 2 * pad))
        let scaleX = available.width / max(0.0001, box.width)
        let scaleY = available.height / max(0.0001, box.height)
        let scale = min(scaleX, scaleY)
        let center = CGPoint(x: available.midX, y: available.midY)
        let transform = CGAffineTransform.identity
            .translatedBy(x: center.x, y: center.y)
            .scaledBy(x: scale, y: -scale)
            .translatedBy(x: -box.midX, y: -box.midY)
        let expected = oxygen.position.applying(transform)
        let distance = hypot(expected.x - oxygenLabel.position.x, expected.y - oxygenLabel.position.y)

        XCTAssertLessThan(distance, 1.5, "Expected ring oxygen label to remain anchored at atom coordinate.")
    }

    func testRingOxygenLabelAnchorRespectsRotationPanAndZoom() {
        let molecule = Molecule(
            name: "Oxacycle",
            atoms: [
                Atom(id: 1, element: "O", position: CGPoint(x: 1.0, y: 0.0)),
                Atom(id: 2, element: "C", position: CGPoint(x: 0.5, y: 0.8660254)),
                Atom(id: 3, element: "C", position: CGPoint(x: -0.5, y: 0.8660254)),
                Atom(id: 4, element: "C", position: CGPoint(x: -1.0, y: 0.0)),
                Atom(id: 5, element: "C", position: CGPoint(x: -0.5, y: -0.8660254)),
                Atom(id: 6, element: "C", position: CGPoint(x: 0.5, y: -0.8660254))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single),
                Bond(id: 2, a1: 2, a2: 3, order: .single),
                Bond(id: 3, a1: 3, a2: 4, order: .single),
                Bond(id: 4, a1: 4, a2: 5, order: .single),
                Bond(id: 5, a1: 5, a2: 6, order: .single),
                Bond(id: 6, a1: 6, a2: 1, order: .single)
            ]
        )
        let style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14.0,
                                padding: 24.0)
        let canvasRect = CGRect(x: 0, y: 0, width: 720, height: 520)
        let zoom: CGFloat = 1.65
        let pan = CGSize(width: 32, height: -26)
        let rotationDegrees: CGFloat = 73

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: canvasRect,
                                                        zoom: zoom,
                                                        pan: pan,
                                                        rotationDegrees: rotationDegrees)

        guard let oxygenLabel = scene.labels.first(where: { $0.id == 1 }) else {
            XCTFail("Expected oxygen label for atom 1.")
            return
        }

        guard let box = molecule.boundingBox(),
              let oxygen = molecule.atoms.first(where: { $0.id == 1 }) else {
            XCTFail("Expected oxygen atom and molecule bounds.")
            return
        }

        let pad = style.padding
        let available = CGRect(x: canvasRect.minX + pad,
                               y: canvasRect.minY + pad,
                               width: max(1, canvasRect.width - 2 * pad),
                               height: max(1, canvasRect.height - 2 * pad))
        let scaleX = available.width / max(0.0001, box.width)
        let scaleY = available.height / max(0.0001, box.height)
        let scale = min(scaleX, scaleY)
        let center = CGPoint(x: available.midX, y: available.midY)
        let transform = CGAffineTransform.identity
            .translatedBy(x: center.x, y: center.y)
            .scaledBy(x: scale, y: -scale)
            .translatedBy(x: -box.midX, y: -box.midY)
        let unrotated = oxygen.position.applying(transform)
        let expected = applyViewportTransform(unrotated,
                                              in: canvasRect,
                                              zoom: zoom,
                                              pan: pan,
                                              rotationDegrees: rotationDegrees)
        let distance = hypot(expected.x - oxygenLabel.position.x, expected.y - oxygenLabel.position.y)

        XCTAssertLessThan(distance, 1.5, "Expected oxygen label to stay anchored during rotation/pan/zoom.")
    }

    func testLabelPlacementIsStableAcrossZoomLevels() {
        let molecule = Molecule(
            name: "PhosphateCluster",
            atoms: [
                Atom(id: 1, element: "P", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "O", position: CGPoint(x: 1.0, y: 0.1)),
                Atom(id: 3, element: "O", position: CGPoint(x: -0.95, y: 0.25)),
                Atom(id: 4, element: "O", position: CGPoint(x: 0.15, y: 1.1)),
                Atom(id: 5, element: "N", position: CGPoint(x: -0.25, y: -1.0)),
                Atom(id: 6, element: "O", position: CGPoint(x: 1.9, y: 0.55)),
                Atom(id: 7, element: "C", position: CGPoint(x: -1.35, y: -0.85))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .double),
                Bond(id: 2, a1: 1, a2: 3, order: .single),
                Bond(id: 3, a1: 1, a2: 4, order: .single),
                Bond(id: 4, a1: 1, a2: 5, order: .single),
                Bond(id: 5, a1: 2, a2: 6, order: .single),
                Bond(id: 6, a1: 5, a2: 7, order: .single)
            ]
        )

        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: true,
                                showAtomIDs: false,
                                bondWidth: 2.35,
                                fontSize: 28.0,
                                padding: 24.0)
        let canvas = CGRect(x: 0, y: 0, width: 900, height: 620)
        let pan = CGSize(width: 36, height: -24)

        let lowZoomScene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                                style: style,
                                                                canvasRect: canvas,
                                                                zoom: 0.45,
                                                                pan: pan)
        let highZoomScene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                                 style: style,
                                                                 canvasRect: canvas,
                                                                 zoom: 1.85,
                                                                 pan: pan)

        let lowByID = Dictionary(uniqueKeysWithValues: lowZoomScene.labels.map { ($0.id, $0.position) })
        let highByID = Dictionary(uniqueKeysWithValues: highZoomScene.labels.map { ($0.id, $0.position) })
        let commonIDs = Set(lowByID.keys).intersection(highByID.keys)
        XCTAssertFalse(commonIDs.isEmpty, "Expected at least one common label between zoom levels.")

        for atomID in commonIDs {
            guard let low = lowByID[atomID], let high = highByID[atomID] else {
                XCTFail("Missing label coordinates for atom \(atomID).")
                continue
            }
            let lowCanonical = canonicalizeViewportPoint(low, in: canvas, zoom: 0.45, pan: pan)
            let highCanonical = canonicalizeViewportPoint(high, in: canvas, zoom: 1.85, pan: pan)
            let drift = hypot(lowCanonical.x - highCanonical.x, lowCanonical.y - highCanonical.y)
            XCTAssertLessThan(drift, 0.75, "Label for atom \(atomID) drifted across zoom levels (\(drift)).")
        }
    }

    func testLabelFontScalesWithCanvasResizeAutoFit() {
        let molecule = Molecule(
            name: "FormylFragment",
            atoms: [
                Atom(id: 1, element: "O", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "C", position: CGPoint(x: 1.22, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single)
            ]
        )
        let style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.35,
                                fontSize: 24.0,
                                padding: 24.0)

        let smallCanvas = CGRect(x: 0, y: 0, width: 700, height: 500)
        let largeCanvas = CGRect(x: 0, y: 0, width: 1000, height: 700)
        let smallScene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                              style: style,
                                                              canvasRect: smallCanvas,
                                                              zoom: 1.0,
                                                              pan: .zero)
        let largeScene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                              style: style,
                                                              canvasRect: largeCanvas,
                                                              zoom: 1.0,
                                                              pan: .zero)

        guard let smallLabel = smallScene.labels.first(where: { $0.id == 1 }),
              let largeLabel = largeScene.labels.first(where: { $0.id == 1 }) else {
            XCTFail("Expected oxygen label in both scenes.")
            return
        }

        XCTAssertGreaterThan(largeLabel.fontSize,
                             smallLabel.fontSize * 1.15,
                             "Expected label font to grow when viewport auto-fit enlarges structure.")

        let smallBondLength = longestBondLength(in: smallScene.bondSegments)
        let largeBondLength = longestBondLength(in: largeScene.bondSegments)
        let smallRatio = smallLabel.fontSize / max(1.0, smallBondLength)
        let largeRatio = largeLabel.fontSize / max(1.0, largeBondLength)
        XCTAssertLessThan(abs(smallRatio - largeRatio),
                          0.08,
                          "Expected label-to-bond scale to stay consistent across viewport resize auto-fit.")

        let smallWidth = averageBondWidth(in: smallScene.bondSegments)
        let largeWidth = averageBondWidth(in: largeScene.bondSegments)
        XCTAssertGreaterThan(largeWidth,
                             smallWidth * 1.15,
                             "Expected bond stroke width to grow when viewport auto-fit enlarges structure.")

        let smallWidthRatio = smallWidth / max(1.0, smallBondLength)
        let largeWidthRatio = largeWidth / max(1.0, largeBondLength)
        XCTAssertLessThan(abs(smallWidthRatio - largeWidthRatio),
                          0.02,
                          "Expected bond-width-to-length ratio to stay consistent across viewport resize auto-fit.")
    }

    func testLabelAndBondSizingKeepShrinkingOnVeryCompactCanvas() {
        let molecule = Molecule(
            name: "FormylFragment",
            atoms: [
                Atom(id: 1, element: "O", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "C", position: CGPoint(x: 1.22, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single)
            ]
        )
        let style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.35,
                                fontSize: 24.0,
                                padding: 24.0)

        let compactCanvas = CGRect(x: 0, y: 0, width: 620, height: 420)
        let extraCompactCanvas = CGRect(x: 0, y: 0, width: 500, height: 340)
        let compactScene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                                style: style,
                                                                canvasRect: compactCanvas,
                                                                zoom: 1.0,
                                                                pan: .zero)
        let extraCompactScene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                                     style: style,
                                                                     canvasRect: extraCompactCanvas,
                                                                     zoom: 1.0,
                                                                     pan: .zero)

        guard let compactLabel = compactScene.labels.first(where: { $0.id == 1 }),
              let extraCompactLabel = extraCompactScene.labels.first(where: { $0.id == 1 }) else {
            XCTFail("Expected oxygen label in both scenes.")
            return
        }

        XCTAssertLessThan(extraCompactLabel.fontSize,
                          compactLabel.fontSize * 0.92,
                          "Expected labels to continue shrinking when the canvas becomes very compact.")

        let compactBondLength = longestBondLength(in: compactScene.bondSegments)
        let extraCompactBondLength = longestBondLength(in: extraCompactScene.bondSegments)
        let compactRatio = compactLabel.fontSize / max(1.0, compactBondLength)
        let extraCompactRatio = extraCompactLabel.fontSize / max(1.0, extraCompactBondLength)
        XCTAssertLessThan(abs(compactRatio - extraCompactRatio),
                          0.08,
                          "Expected label-to-bond scale to remain stable on compact canvases.")

        let compactBondWidth = averageBondWidth(in: compactScene.bondSegments)
        let extraCompactBondWidth = averageBondWidth(in: extraCompactScene.bondSegments)
        XCTAssertLessThan(extraCompactBondWidth,
                          compactBondWidth * 0.92,
                          "Expected bond widths to keep shrinking on very compact canvases.")
    }

    func testSparseSceneGetsBoundedLabelBoostAboveViewportBaseline() {
        let molecule = Molecule(
            name: "FormylFragment",
            atoms: [
                Atom(id: 1, element: "O", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "C", position: CGPoint(x: 1.22, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single)
            ]
        )
        let style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.35,
                                fontSize: 24.0,
                                padding: 24.0)
        let canvas = CGRect(x: 0, y: 0, width: 760, height: 540)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: canvas,
                                                        zoom: 1.0,
                                                        pan: .zero)

        guard let oxygenLabel = scene.labels.first(where: { $0.id == 1 }) else {
            XCTFail("Expected oxygen label in sparse scene.")
            return
        }

        let viewportBaseline = viewportFitBaselineFontSize(for: molecule, style: style, canvasRect: canvas)
        XCTAssertGreaterThan(oxygenLabel.fontSize,
                             viewportBaseline * 1.18,
                             "Expected sparse scenes to receive a visible label boost over viewport-only sizing.")
        XCTAssertLessThanOrEqual(oxygenLabel.fontSize,
                                 viewportBaseline * 1.50 + 0.0001,
                                 "Expected sparse-scene label boost to remain bounded.")
    }

    func testDenseSceneDoesNotReceiveSparseProminenceBoost() {
        let atoms = (0..<10).map { index in
            Atom(id: index + 1, element: "N", position: CGPoint(x: CGFloat(index) * 1.18, y: 0.0))
        }
        let bonds = (0..<9).map { index in
            Bond(id: index + 1, a1: index + 1, a2: index + 2, order: .single)
        }
        let molecule = Molecule(name: "DenseNitrogenChain", atoms: atoms, bonds: bonds)
        let style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.35,
                                fontSize: 24.0,
                                padding: 24.0)
        let canvas = CGRect(x: 0, y: 0, width: 760, height: 360)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: canvas,
                                                        zoom: 1.0,
                                                        pan: .zero)

        guard let nitrogenLabel = scene.labels.first(where: { $0.id == 1 }) else {
            XCTFail("Expected visible nitrogen labels in dense scene.")
            return
        }

        let viewportBaseline = viewportFitBaselineFontSize(for: molecule, style: style, canvasRect: canvas)
        XCTAssertLessThanOrEqual(nitrogenLabel.fontSize,
                                 viewportBaseline * 1.06 + 0.0001,
                                 "Expected dense scenes to stay close to viewport-only label sizing.")
    }

    func testSparseSceneGetsVisibleBondBoostAboveViewportBaseline() {
        let molecule = Molecule(
            name: "FormylFragment",
            atoms: [
                Atom(id: 1, element: "O", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "C", position: CGPoint(x: 1.22, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single)
            ]
        )
        let style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.35,
                                fontSize: 24.0,
                                padding: 24.0)
        let canvas = CGRect(x: 0, y: 0, width: 760, height: 540)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: canvas,
                                                        zoom: 1.0,
                                                        pan: .zero)

        let averageWidth = averageBondWidth(in: scene.bondSegments)
        let viewportBaseline = viewportFitBaselineBondWidth(for: molecule, style: style, canvasRect: canvas)
        XCTAssertGreaterThan(averageWidth,
                             viewportBaseline * 1.08,
                             "Expected sparse scenes to receive a visible bond-width boost over viewport-only sizing.")
        XCTAssertLessThanOrEqual(averageWidth,
                                 viewportBaseline * 1.34 + 0.0001,
                                 "Expected sparse-scene bond boost to remain bounded.")
    }

    func testZoomStrengthensSparseSceneLabelsAndBondsBeyondLinearScaling() {
        let molecule = Molecule(
            name: "FormylFragment",
            atoms: [
                Atom(id: 1, element: "O", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "C", position: CGPoint(x: 1.22, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single)
            ]
        )
        let style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.35,
                                fontSize: 24.0,
                                padding: 24.0)
        let canvas = CGRect(x: 0, y: 0, width: 760, height: 540)

        let baseScene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                            style: style,
                                                            canvasRect: canvas,
                                                            zoom: 1.0,
                                                            pan: .zero)
        let zoomedScene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                              style: style,
                                                              canvasRect: canvas,
                                                              zoom: 1.5,
                                                              pan: .zero)

        guard let baseLabel = baseScene.labels.first(where: { $0.id == 1 }),
              let zoomedLabel = zoomedScene.labels.first(where: { $0.id == 1 }) else {
            XCTFail("Expected oxygen label in both zoom scenes.")
            return
        }

        let baseBondWidth = averageBondWidth(in: baseScene.bondSegments)
        let zoomedBondWidth = averageBondWidth(in: zoomedScene.bondSegments)
        XCTAssertGreaterThan(zoomedLabel.fontSize,
                             baseLabel.fontSize * 1.55,
                             "Expected sparse-scene label size to scale more strongly than the previous linear zoom baseline.")
        XCTAssertGreaterThan(zoomedBondWidth,
                             baseBondWidth * 1.52,
                             "Expected sparse-scene bond width to scale more strongly than the previous linear zoom baseline.")
    }

    func testZoomedOutSparseSceneRetainsLabelAndBondReadabilityAboveLinearShrink() {
        let molecule = Molecule(
            name: "FormylFragment",
            atoms: [
                Atom(id: 1, element: "O", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "C", position: CGPoint(x: 1.22, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single)
            ]
        )
        let style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.35,
                                fontSize: 24.0,
                                padding: 24.0)
        let canvas = CGRect(x: 0, y: 0, width: 760, height: 540)

        let baseScene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                            style: style,
                                                            canvasRect: canvas,
                                                            zoom: 1.0,
                                                            pan: .zero)
        let zoomedOutScene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                                 style: style,
                                                                 canvasRect: canvas,
                                                                 zoom: 0.75,
                                                                 pan: .zero)

        guard let baseLabel = baseScene.labels.first(where: { $0.id == 1 }),
              let zoomedOutLabel = zoomedOutScene.labels.first(where: { $0.id == 1 }) else {
            XCTFail("Expected oxygen label in both zoom scenes.")
            return
        }

        let baseBondWidth = averageBondWidth(in: baseScene.bondSegments)
        let zoomedOutBondWidth = averageBondWidth(in: zoomedOutScene.bondSegments)
        XCTAssertGreaterThan(zoomedOutLabel.fontSize,
                             baseLabel.fontSize * 0.80,
                             "Expected zoomed-out sparse-scene labels to shrink more gently than linear zoom.")
        XCTAssertGreaterThan(zoomedOutBondWidth,
                             baseBondWidth * 0.79,
                             "Expected zoomed-out sparse-scene bond widths to shrink more gently than linear zoom.")
    }

    func testCanLowerMinimumLabelFontSizeForTinyPreviewCanvases() {
        let molecule = Molecule(
            name: "FormylFragment",
            atoms: [
                Atom(id: 1, element: "O", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "C", position: CGPoint(x: 1.22, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single)
            ]
        )
        let style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 1.95,
                                fontSize: 8.0,
                                padding: 8.0)
        let tinyCanvas = CGRect(x: 0, y: 0, width: 94, height: 54)

        let defaultScene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                                style: style,
                                                                canvasRect: tinyCanvas,
                                                                zoom: 1.0,
                                                                pan: .zero,
                                                                rotationDegrees: 0)
        let loweredMinimumScene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                                       style: style,
                                                                       canvasRect: tinyCanvas,
                                                                       zoom: 1.0,
                                                                       pan: .zero,
                                                                       rotationDegrees: 0,
                                                                       minimumLabelFontSize: 4.8)

        guard let defaultLabel = defaultScene.labels.first(where: { $0.id == 1 }),
              let loweredLabel = loweredMinimumScene.labels.first(where: { $0.id == 1 }) else {
            XCTFail("Expected oxygen label in both tiny-preview scenes.")
            return
        }

        XCTAssertEqual(defaultLabel.fontSize, 8.0, accuracy: 0.0001)
        XCTAssertEqual(loweredLabel.fontSize, 4.8, accuracy: 0.0001)
    }

    func testMappedCarbonsRemainVisibleWhenMapNumbersEnabledAndCarbonsHidden() {
        let molecule = Molecule(
            name: "MappedCarbon",
            atoms: [
                Atom(id: 1, element: "C", position: CGPoint(x: 0.0, y: 0.0), atomMapNumber: 7),
                Atom(id: 2, element: "O", position: CGPoint(x: 1.2, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single)
            ]
        )
        var style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14.0,
                                padding: 24.0)
        style.showAtomMapNumbers = true

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: CGRect(x: 0, y: 0, width: 480, height: 320),
                                                        zoom: 1.0,
                                                        pan: .zero)
        XCTAssertTrue(scene.labels.contains(where: { $0.id == 1 && $0.text == "C:7" }),
                      "Mapped carbons should remain labeled when map-number rendering is enabled.")
    }

    func testCanOptIntoAromaticCarbonLabelsWhenCarbonsHidden() {
        let style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14.0,
                                padding: 24.0)
        let canvas = CGRect(x: 0, y: 0, width: 420, height: 300)

        let defaultScene = CDKMetalDepictionSceneBuilder.build(molecule: aromaticSixRing,
                                                                style: style,
                                                                canvasRect: canvas,
                                                                zoom: 1.0,
                                                                pan: .zero)
        let optInScene = CDKMetalDepictionSceneBuilder.build(molecule: aromaticSixRing,
                                                              style: style,
                                                              canvasRect: canvas,
                                                              zoom: 1.0,
                                                              pan: .zero,
                                                              includeAromaticCarbonLabelsWhenCarbonsHidden: true)

        XCTAssertFalse(defaultScene.labels.contains(where: { $0.text == "C" }),
                       "Default rendering should suppress aromatic carbon labels when carbons are hidden.")
        XCTAssertTrue(optInScene.labels.contains(where: { $0.text == "C" }),
                      "Callers can opt into aromatic carbon labels when needed.")
        XCTAssertGreaterThan(defaultScene.bondSegments.count, 0,
                             "Suppressing aromatic carbon labels must not remove bond rendering.")
    }

    func testCanOptIntoTerminalCarbonLabelsWhenCarbonsHidden() {
        let molecule = Molecule(
            name: "TerminalCarbon",
            atoms: [
                Atom(id: 1, element: "C", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "O", position: CGPoint(x: 1.22, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single)
            ]
        )
        let style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 14.0,
                                padding: 24.0)
        let canvas = CGRect(x: 0, y: 0, width: 420, height: 300)

        let defaultScene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                                style: style,
                                                                canvasRect: canvas,
                                                                zoom: 1.0,
                                                                pan: .zero)
        let optInScene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                              style: style,
                                                              canvasRect: canvas,
                                                              zoom: 1.0,
                                                              pan: .zero,
                                                              includeTerminalCarbonLabelsWhenCarbonsHidden: true)

        XCTAssertFalse(defaultScene.labels.contains(where: { $0.id == 1 && $0.text == "C" }),
                       "Default rendering should suppress terminal carbon labels when carbons are hidden.")
        XCTAssertTrue(optInScene.labels.contains(where: { $0.id == 1 && $0.text == "C" }),
                      "Callers can opt into terminal carbon labels when needed.")
        XCTAssertTrue(defaultScene.labels.contains(where: { $0.id == 2 && $0.text.hasPrefix("O") }))
    }

    func testMarkushLegendStaysBelowRotatedRootDuringInteractiveRotation() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|")
        let style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: true,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 16.0,
                                padding: 24.0)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: CGRect(x: 0, y: 0, width: 900, height: 620),
                                                        zoom: 1.35,
                                                        pan: CGSize(width: 26, height: -18),
                                                        rotationDegrees: 112)

        let box = try XCTUnwrap(scene.backgroundBoxes.first?.rect)
        let nh = try XCTUnwrap(scene.labels.first(where: { $0.text == "NH" }))
        let repeatCount = try XCTUnwrap(scene.labels.first(where: { $0.text == "1-3" }))
        let rootR = try XCTUnwrap(scene.labels.first(where: { $0.text == "R¹" }))

        XCTAssertGreaterThan(box.minY, nh.position.y + (nh.fontSize * 0.35))
        XCTAssertGreaterThan(box.minY, repeatCount.position.y + (repeatCount.fontSize * 0.20))
        XCTAssertGreaterThan(box.minY, rootR.position.y + (rootR.fontSize * 0.15))

        let expandedBox = box.insetBy(dx: -2, dy: -2)
        let collidingRootSegments = scene.bondSegments.filter { segment in
            let bounds = segmentBounds(for: segment)
            guard bounds.intersects(box) else { return false }
            return !(expandedBox.contains(segment.from) && expandedBox.contains(segment.to))
        }
        XCTAssertTrue(collidingRootSegments.isEmpty,
                      "Expected the fixed legend box to stay clear of rotated scaffold bonds.")
    }

    func testMarkushLegendFragmentsRotateInsideFixedLegendBox() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|")
        let style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: true,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 16.0,
                                padding: 24.0)
        let canvas = CGRect(x: 0, y: 0, width: 900, height: 620)

        let scene0 = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                         style: style,
                                                         canvasRect: canvas,
                                                         zoom: 1.0,
                                                         pan: .zero,
                                                         rotationDegrees: 0)
        let scene90 = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                          style: style,
                                                          canvasRect: canvas,
                                                          zoom: 1.0,
                                                          pan: .zero,
                                                          rotationDegrees: 90)

        let box0 = try XCTUnwrap(scene0.backgroundBoxes.first?.rect)
        let box90 = try XCTUnwrap(scene90.backgroundBoxes.first?.rect)
        let oxygen0 = try XCTUnwrap(scene0.labels.first(where: { $0.text == "O" }))
        let oxygen90 = try XCTUnwrap(scene90.labels.first(where: { $0.text == "O" }))
        let chloride90 = try XCTUnwrap(scene90.labels.first(where: { $0.text == "Cl" }))
        let nitrile90 = try XCTUnwrap(scene90.labels.first(where: { $0.text == "N" }))

        let offset0 = CGPoint(x: oxygen0.position.x - box0.midX, y: oxygen0.position.y - box0.midY)
        let offset90 = CGPoint(x: oxygen90.position.x - box90.midX, y: oxygen90.position.y - box90.midY)
        XCTAssertGreaterThan(hypot(offset0.x - offset90.x, offset0.y - offset90.y),
                             8.0,
                             "Expected legend fragment contents to rotate relative to the fixed box.")

        for label in [oxygen90, chloride90, nitrile90] {
            XCTAssertTrue(box90.insetBy(dx: 2, dy: 2).contains(label.position),
                          "Expected rotated legend label \(label.text) to stay inside the fixed legend box.")
        }

        let expandedBox = box90.insetBy(dx: -2, dy: -2)
        XCTAssertTrue(scene90.bondSegments.filter { segmentBounds(for: $0).intersects(box90) }
            .allSatisfy { expandedBox.contains($0.from) && expandedBox.contains($0.to) },
                      "Expected any segment inside the fixed legend box to remain fully contained after rotation.")
    }

    private func countSupportingSegments(forBondFrom from: CGPoint,
                                         to: CGPoint,
                                         segments: [CDKMetalDepictionScene.LineSegment]) -> Int {
        let dx = to.x - from.x
        let dy = to.y - from.y
        let bondLength = hypot(dx, dy)
        guard bondLength > 0.0001 else { return 0 }

        let ux = dx / bondLength
        let uy = dy / bondLength
        let nx = -uy
        let ny = ux
        let bondMid = CGPoint(x: (from.x + to.x) * 0.5, y: (from.y + to.y) * 0.5)

        return segments.filter { segment in
            let sx = segment.to.x - segment.from.x
            let sy = segment.to.y - segment.from.y
            let segLength = hypot(sx, sy)
            guard segLength > bondLength * 0.55 else { return false }

            let sux = sx / segLength
            let suy = sy / segLength
            let alignment = abs((sux * ux) + (suy * uy))
            guard alignment >= 0.985 else { return false }

            let segMid = CGPoint(x: (segment.from.x + segment.to.x) * 0.5,
                                 y: (segment.from.y + segment.to.y) * 0.5)
            let rx = segMid.x - bondMid.x
            let ry = segMid.y - bondMid.y

            let alongDistance = abs((rx * ux) + (ry * uy))
            let perpendicularDistance = abs((rx * nx) + (ry * ny))

            return alongDistance <= bondLength * 0.30 && perpendicularDistance <= bondLength * 0.40
        }.count
    }

    private func averageStereoStripeLength(in segments: [CDKMetalDepictionScene.LineSegment],
                                           from start: CGPoint,
                                           to end: CGPoint,
                                           tRange: ClosedRange<CGFloat>) -> CGFloat {
        let dx = end.x - start.x
        let dy = end.y - start.y
        let bondLength = hypot(dx, dy)
        guard bondLength > 0.0001 else { return 0 }

        let ux = dx / bondLength
        let uy = dy / bondLength
        let nx = -uy
        let ny = ux

        let lengths = segments.compactMap { segment -> CGFloat? in
            let sx = segment.to.x - segment.from.x
            let sy = segment.to.y - segment.from.y
            let segLength = hypot(sx, sy)
            guard segLength > 0.0001 else { return nil }

            let sux = sx / segLength
            let suy = sy / segLength
            let perpendicularAlignment = abs((sux * nx) + (suy * ny))
            guard perpendicularAlignment >= 0.94 else { return nil }

            let mid = CGPoint(x: (segment.from.x + segment.to.x) * 0.5,
                              y: (segment.from.y + segment.to.y) * 0.5)
            let rx = mid.x - start.x
            let ry = mid.y - start.y
            let t = ((rx * ux) + (ry * uy)) / bondLength
            guard tRange.contains(t) else { return nil }

            return segLength
        }

        guard !lengths.isEmpty else { return 0 }
        return lengths.reduce(0, +) / CGFloat(lengths.count)
    }

    private func nearestEndpointDistance(in segments: [CDKMetalDepictionScene.LineSegment],
                                         to point: CGPoint) -> CGFloat {
        segments
            .flatMap { [hypot($0.from.x - point.x, $0.from.y - point.y),
                        hypot($0.to.x - point.x, $0.to.y - point.y)] }
            .min() ?? .greatestFiniteMagnitude
    }

    private func longestBondLength(in segments: [CDKMetalDepictionScene.LineSegment]) -> CGFloat {
        segments.map { hypot($0.to.x - $0.from.x, $0.to.y - $0.from.y) }.max() ?? 0
    }

    private func averageBondWidth(in segments: [CDKMetalDepictionScene.LineSegment]) -> CGFloat {
        guard !segments.isEmpty else { return 0 }
        let total = segments.reduce(CGFloat.zero) { $0 + $1.width }
        return total / CGFloat(segments.count)
    }

    private func segmentBounds(for segment: CDKMetalDepictionScene.LineSegment) -> CGRect {
        CGRect(x: min(segment.from.x, segment.to.x),
               y: min(segment.from.y, segment.to.y),
               width: abs(segment.to.x - segment.from.x),
               height: abs(segment.to.y - segment.from.y))
    }

    private func viewportFitBaselineFontSize(for molecule: Molecule,
                                             style: RenderStyle,
                                             canvasRect: CGRect) -> CGFloat {
        let depictionMolecule = CDKDepictionPreprocessor.prepareForRendering(molecule: molecule, style: style)
        guard let box = depictionMolecule.boundingBox() else { return style.fontSize }

        let available = CGRect(x: canvasRect.minX + style.padding,
                               y: canvasRect.minY + style.padding,
                               width: max(1, canvasRect.width - (2 * style.padding)),
                               height: max(1, canvasRect.height - (2 * style.padding)))
        let scaleX = available.width / max(0.0001, box.width)
        let scaleY = available.height / max(0.0001, box.height)
        let scale = min(scaleX, scaleY)

        let referenceViewport = CGSize(width: 960, height: 720)
        let referenceAvailable = CGRect(x: style.padding,
                                        y: style.padding,
                                        width: max(1, referenceViewport.width - (2 * style.padding)),
                                        height: max(1, referenceViewport.height - (2 * style.padding)))
        let referenceScaleX = referenceAvailable.width / max(0.0001, box.width)
        let referenceScaleY = referenceAvailable.height / max(0.0001, box.height)
        let referenceFitScale = min(referenceScaleX, referenceScaleY)
        let viewportAutoFitScale = min(1.90, scale / max(0.0001, referenceFitScale))
        return max(8.0, style.fontSize * viewportAutoFitScale)
    }

    private func viewportFitBaselineBondWidth(for molecule: Molecule,
                                              style: RenderStyle,
                                              canvasRect: CGRect) -> CGFloat {
        let depictionMolecule = CDKDepictionPreprocessor.prepareForRendering(molecule: molecule, style: style)
        guard let box = depictionMolecule.boundingBox() else { return style.bondWidth }

        let available = CGRect(x: canvasRect.minX + style.padding,
                               y: canvasRect.minY + style.padding,
                               width: max(1, canvasRect.width - (2 * style.padding)),
                               height: max(1, canvasRect.height - (2 * style.padding)))
        let scaleX = available.width / max(0.0001, box.width)
        let scaleY = available.height / max(0.0001, box.height)
        let scale = min(scaleX, scaleY)

        let referenceViewport = CGSize(width: 960, height: 720)
        let referenceAvailable = CGRect(x: style.padding,
                                        y: style.padding,
                                        width: max(1, referenceViewport.width - (2 * style.padding)),
                                        height: max(1, referenceViewport.height - (2 * style.padding)))
        let referenceScaleX = referenceAvailable.width / max(0.0001, box.width)
        let referenceScaleY = referenceAvailable.height / max(0.0001, box.height)
        let referenceFitScale = min(referenceScaleX, referenceScaleY)
        let viewportAutoFitScale = min(1.90, scale / max(0.0001, referenceFitScale))
        return max(1.0, style.bondWidth * viewportAutoFitScale)
    }

    private var alternatingSixRing: Molecule {
        let r = CGFloat(1.0)
        let atoms = [
            Atom(id: 1, element: "C", position: CGPoint(x: r, y: 0)),
            Atom(id: 2, element: "C", position: CGPoint(x: r * 0.5, y: r * 0.8660254)),
            Atom(id: 3, element: "C", position: CGPoint(x: -r * 0.5, y: r * 0.8660254)),
            Atom(id: 4, element: "C", position: CGPoint(x: -r, y: 0)),
            Atom(id: 5, element: "C", position: CGPoint(x: -r * 0.5, y: -r * 0.8660254)),
            Atom(id: 6, element: "C", position: CGPoint(x: r * 0.5, y: -r * 0.8660254))
        ]
        let bonds = [
            Bond(id: 1, a1: 1, a2: 2, order: .single),
            Bond(id: 2, a1: 2, a2: 3, order: .double),
            Bond(id: 3, a1: 3, a2: 4, order: .single),
            Bond(id: 4, a1: 4, a2: 5, order: .double),
            Bond(id: 5, a1: 5, a2: 6, order: .single),
            Bond(id: 6, a1: 6, a2: 1, order: .double)
        ]
        return Molecule(name: "AlternatingSixRing", atoms: atoms, bonds: bonds)
    }

    private var aromaticSixRing: Molecule {
        let r = CGFloat(1.0)
        let atoms = [
            Atom(id: 1, element: "C", position: CGPoint(x: r, y: 0), aromatic: true),
            Atom(id: 2, element: "C", position: CGPoint(x: r * 0.5, y: r * 0.8660254), aromatic: true),
            Atom(id: 3, element: "C", position: CGPoint(x: -r * 0.5, y: r * 0.8660254), aromatic: true),
            Atom(id: 4, element: "C", position: CGPoint(x: -r, y: 0), aromatic: true),
            Atom(id: 5, element: "C", position: CGPoint(x: -r * 0.5, y: -r * 0.8660254), aromatic: true),
            Atom(id: 6, element: "C", position: CGPoint(x: r * 0.5, y: -r * 0.8660254), aromatic: true)
        ]
        let bonds = [
            Bond(id: 1, a1: 1, a2: 2, order: .aromatic),
            Bond(id: 2, a1: 2, a2: 3, order: .aromatic),
            Bond(id: 3, a1: 3, a2: 4, order: .aromatic),
            Bond(id: 4, a1: 4, a2: 5, order: .aromatic),
            Bond(id: 5, a1: 5, a2: 6, order: .aromatic),
            Bond(id: 6, a1: 6, a2: 1, order: .aromatic)
        ]
        return Molecule(name: "AromaticSixRing", atoms: atoms, bonds: bonds)
    }

    private func canonicalizeViewportPoint(_ point: CGPoint,
                                           in canvas: CGRect,
                                           zoom: CGFloat,
                                           pan: CGSize) -> CGPoint {
        CGPoint(x: ((point.x - pan.width - canvas.midX) / zoom) + canvas.midX,
                y: ((point.y - pan.height - canvas.midY) / zoom) + canvas.midY)
    }

    private func applyViewportTransform(_ point: CGPoint,
                                        in canvas: CGRect,
                                        zoom: CGFloat,
                                        pan: CGSize,
                                        rotationDegrees: CGFloat) -> CGPoint {
        let center = CGPoint(x: canvas.midX, y: canvas.midY)
        let radians = rotationDegrees * (.pi / 180)
        let c = cos(radians)
        let s = sin(radians)
        let rx = point.x - center.x
        let ry = point.y - center.y
        let rotated = CGPoint(x: (rx * c) - (ry * s) + center.x,
                              y: (rx * s) + (ry * c) + center.y)
        return CGPoint(x: ((rotated.x - center.x) * zoom) + center.x + pan.width,
                       y: ((rotated.y - center.y) * zoom) + center.y + pan.height)
    }
}
