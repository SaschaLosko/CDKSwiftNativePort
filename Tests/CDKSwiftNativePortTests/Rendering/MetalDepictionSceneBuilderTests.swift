import CoreGraphics
import XCTest
@testable import CDKSwiftNativePort

final class MetalDepictionSceneBuilderTests: XCTestCase {

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
            8,
            "Expected stereo up wedge to emit multiple emphasized segments, but only got \(emphasizedWedgeSegments.count)."
        )
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
