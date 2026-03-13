import CoreGraphics
import XCTest
@testable import CDKSwiftNativePort

final class DepictionGeneratorTests: XCTestCase {
    private let smilesParser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    func testGeneratesSVGDocumentForSimpleMolecule() throws {
        let molecule = try smilesParser.parseSmiles("CCO")
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.35,
                                fontSize: 16,
                                padding: 24)

        let svg = CDKDepictionGenerator.toSVG(molecule: molecule,
                                              style: style,
                                              canvasSize: CGSize(width: 800, height: 500))

        XCTAssertTrue(svg.contains("<svg "))
        XCTAssertTrue(svg.contains("</svg>"))
        XCTAssertTrue(svg.contains("<line "))
        XCTAssertTrue(svg.contains("<text "))
        XCTAssertTrue(svg.contains("viewBox=\"0 0 800 500\""))
    }

    func testRendersExplicitDoubleBondWithTwoLines() throws {
        let molecule = try smilesParser.parseSmiles("C=C")
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.35,
                                fontSize: 16,
                                padding: 24)

        let svg = CDKDepictionGenerator.toSVG(molecule: molecule,
                                              style: style,
                                              canvasSize: CGSize(width: 600, height: 400))

        let lineCount = svg.components(separatedBy: "<line ").count - 1
        XCTAssertGreaterThanOrEqual(lineCount, 2, "Expected at least two line segments for a double bond.")
    }

    func testRendersStereoUpBondAsWedgePolygon() {
        let molecule = Molecule(
            name: "StereoUp",
            atoms: [
                Atom(id: 1, element: "C", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "C", position: CGPoint(x: 1.0, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single, stereo: .up)
            ]
        )
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.35,
                                fontSize: 16,
                                padding: 24)

        let svg = CDKDepictionGenerator.toSVG(molecule: molecule,
                                              style: style,
                                              canvasSize: CGSize(width: 400, height: 280))
        XCTAssertTrue(svg.contains("<polygon "), "Expected a filled wedge polygon for stereo up bond.")
    }

    func testSVGClipsBondEndpointsAwayFromLabelCenters() {
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
                                fontSize: 24,
                                padding: 24)

        let svg = CDKDepictionGenerator.toSVG(molecule: molecule,
                                              style: style,
                                              canvasSize: CGSize(width: 520, height: 360))

        let line = try? XCTUnwrap(firstLineCoordinates(in: svg))
        let labels = textCoordinates(in: svg)
        XCTAssertNotNil(line)
        XCTAssertGreaterThanOrEqual(labels.count, 2)
        guard let line else { return }

        let firstTwoLabels = Array(labels.prefix(2))
        for labelCenter in firstTwoLabels {
            let d1 = hypot(line.0.x - labelCenter.x, line.0.y - labelCenter.y)
            let d2 = hypot(line.1.x - labelCenter.x, line.1.y - labelCenter.y)
            XCTAssertGreaterThan(min(d1, d2), 4.0, "SVG bond endpoint should be clipped away from label center.")
        }
    }

    func testSVGSuppressesSimpleExplicitHydrogensByDefault() {
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
                                fontSize: 16,
                                padding: 24)

        let svg = CDKDepictionGenerator.toSVG(molecule: molecule,
                                              style: style,
                                              canvasSize: CGSize(width: 520, height: 360))

        XCTAssertFalse(svg.contains(">H</text>"))
        XCTAssertTrue(svg.contains(">OH</text>"))
    }

    func testSVGHidesAttachmentPointLabelsInMarkushDefinitions() throws {
        let molecule = try smilesParser.parseSmiles("C* |$;R1$,RG:_R1={Cl* |$;_AP1$|}|")
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 16,
                                padding: 24)

        let svg = CDKDepictionGenerator.toSVG(molecule: molecule,
                                              style: style,
                                              canvasSize: CGSize(width: 640, height: 420))

        XCTAssertTrue(svg.contains(">R¹</text>"))
        XCTAssertTrue(svg.contains(">Cl</text>"))
        XCTAssertFalse(svg.contains("_AP1"))
    }

    func testSVGDrawsMarkushBoxAndLinkNodeAnnotations() throws {
        let molecule = try smilesParser.parseSmiles("C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|")
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 16,
                                padding: 24)

        let svg = CDKDepictionGenerator.toSVG(molecule: molecule,
                                              style: style,
                                              canvasSize: CGSize(width: 720, height: 520))

        XCTAssertTrue(svg.contains(">R¹ is</text>"))
        XCTAssertTrue(svg.contains(">1-3</text>"))
        XCTAssertTrue(svg.contains("rx=\""))
    }

    func testSVGItalicizesRGroupLabels() throws {
        let molecule = try smilesParser.parseSmiles("C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|")
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 16,
                                padding: 24)

        let svg = CDKDepictionGenerator.toSVG(molecule: molecule,
                                              style: style,
                                              canvasSize: CGSize(width: 720, height: 520))

        XCTAssertTrue(svg.contains("font-family=\"Times New Roman, Georgia, serif\""))
        XCTAssertTrue(svg.contains("font-style=\"italic\""))
        XCTAssertTrue(svg.contains(">R¹</text>"))
    }

    func testSVGHidesMarkushTerminalCarbonLabelsWhenCarbonsAreHidden() throws {
        let molecule = try smilesParser.parseSmiles("C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|")
        let style = RenderStyle(showCarbons: false,
                                showImplicitHydrogens: true,
                                showAtomIDs: false,
                                bondWidth: 2.0,
                                fontSize: 16,
                                padding: 24)

        let svg = CDKDepictionGenerator.toSVG(molecule: molecule,
                                              style: style,
                                              canvasSize: CGSize(width: 720, height: 520))

        XCTAssertFalse(svg.contains(">C</text>"))
        XCTAssertTrue(svg.contains(">O</text>"))
        XCTAssertTrue(svg.contains(">Cl</text>"))
        XCTAssertTrue(svg.contains(">N</text>"))
    }

    private func firstLineCoordinates(in svg: String) -> (CGPoint, CGPoint)? {
        let pattern = #"<line x1="([\-0-9.]+)" y1="([\-0-9.]+)" x2="([\-0-9.]+)" y2="([\-0-9.]+)""#
        guard let regex = try? NSRegularExpression(pattern: pattern),
              let match = regex.firstMatch(in: svg, range: NSRange(svg.startIndex..., in: svg)),
              match.numberOfRanges == 5 else {
            return nil
        }

        func scalar(at idx: Int) -> CGFloat? {
            guard let range = Range(match.range(at: idx), in: svg) else { return nil }
            return CGFloat(Double(svg[range]) ?? .nan)
        }

        guard let x1 = scalar(at: 1), let y1 = scalar(at: 2), let x2 = scalar(at: 3), let y2 = scalar(at: 4),
              x1.isFinite, y1.isFinite, x2.isFinite, y2.isFinite else {
            return nil
        }
        return (CGPoint(x: x1, y: y1), CGPoint(x: x2, y: y2))
    }

    private func textCoordinates(in svg: String) -> [CGPoint] {
        let pattern = #"<text x="([\-0-9.]+)" y="([\-0-9.]+)""#
        guard let regex = try? NSRegularExpression(pattern: pattern) else { return [] }
        let range = NSRange(svg.startIndex..., in: svg)
        return regex.matches(in: svg, range: range).compactMap { match in
            guard match.numberOfRanges == 3,
                  let xRange = Range(match.range(at: 1), in: svg),
                  let yRange = Range(match.range(at: 2), in: svg),
                  let x = Double(svg[xRange]),
                  let y = Double(svg[yRange]) else {
                return nil
            }
            return CGPoint(x: x, y: y)
        }
    }
}
