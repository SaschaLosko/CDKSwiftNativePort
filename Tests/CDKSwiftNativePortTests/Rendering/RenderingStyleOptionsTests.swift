import CoreGraphics
import XCTest
@testable import CDKSwiftNativePort

final class RenderingStyleOptionsTests: XCTestCase {
    private let smilesParser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    func testMetalSceneUsesElementColorsForLabels() throws {
        let molecule = try smilesParser.parseSmiles("CO")
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                atomColoringMode: .cdk2D,
                                colorBondsByAtom: false,
                                aromaticDisplayMode: .innerLine,
                                bondWidth: 2.0,
                                fontSize: 14,
                                padding: 24)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: CGRect(x: 0, y: 0, width: 480, height: 320),
                                                        zoom: 1.0,
                                                        pan: .zero)

        let oxygenLabel = try XCTUnwrap(scene.labels.first(where: { $0.text.hasPrefix("O") }))
        let carbonLabel = try XCTUnwrap(scene.labels.first(where: { $0.text.hasPrefix("C") }))
        XCTAssertGreaterThan(oxygenLabel.color.red, oxygenLabel.color.green)
        XCTAssertGreaterThan(oxygenLabel.color.red, oxygenLabel.color.blue)
        XCTAssertNotEqual(oxygenLabel.color, carbonLabel.color)
    }

    func testMetalSceneColorsBondsWhenEnabled() throws {
        let molecule = try smilesParser.parseSmiles("CO")
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                atomColoringMode: .cdk2D,
                                colorBondsByAtom: true,
                                aromaticDisplayMode: .innerLine,
                                bondWidth: 2.0,
                                fontSize: 14,
                                padding: 24)

        let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                        style: style,
                                                        canvasRect: CGRect(x: 0, y: 0, width: 480, height: 320),
                                                        zoom: 1.0,
                                                        pan: .zero)

        let anyBond = try XCTUnwrap(scene.bondSegments.first)
        XCTAssertNotEqual(anyBond.color, CDKRenderColor.ink)
    }

    func testAromaticCircleModeAddsAdditionalGeometry() throws {
        let molecule = try smilesParser.parseSmiles("c1ccccc1")
        let baseStyle = RenderStyle(showCarbons: true,
                                    showImplicitHydrogens: false,
                                    showAtomIDs: false,
                                    atomColoringMode: .monochrome,
                                    colorBondsByAtom: false,
                                    aromaticDisplayMode: .innerLine,
                                    bondWidth: 2.0,
                                    fontSize: 14,
                                    padding: 24)
        let circleStyle = RenderStyle(showCarbons: true,
                                      showImplicitHydrogens: false,
                                      showAtomIDs: false,
                                      atomColoringMode: .monochrome,
                                      colorBondsByAtom: false,
                                      aromaticDisplayMode: .circle,
                                      bondWidth: 2.0,
                                      fontSize: 14,
                                      padding: 24)

        let innerScene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                             style: baseStyle,
                                                             canvasRect: CGRect(x: 0, y: 0, width: 480, height: 320),
                                                             zoom: 1.0,
                                                             pan: .zero)
        let circleScene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                              style: circleStyle,
                                                              canvasRect: CGRect(x: 0, y: 0, width: 480, height: 320),
                                                              zoom: 1.0,
                                                              pan: .zero)

        XCTAssertGreaterThan(circleScene.bondSegments.count, innerScene.bondSegments.count + 20)
    }

    func testSVGContainsCircleAndElementColorWhenEnabled() throws {
        let molecule = try smilesParser.parseSmiles("c1ccccc1O")
        let style = RenderStyle(showCarbons: true,
                                showImplicitHydrogens: false,
                                showAtomIDs: false,
                                atomColoringMode: .cdk2D,
                                colorBondsByAtom: true,
                                aromaticDisplayMode: .circle,
                                bondWidth: 2.0,
                                fontSize: 14,
                                padding: 24)

        let svg = CDKDepictionGenerator.toSVG(molecule: molecule,
                                              style: style,
                                              canvasSize: CGSize(width: 600, height: 400))
        let oxygenColorHex = CDKRenderingStyleResolver
            .atomColor(for: Atom(id: -1, element: "O", position: .zero), style: style)
            .svgHexRGB()

        XCTAssertTrue(svg.contains("<circle "))
        XCTAssertTrue(svg.contains("fill=\"\(oxygenColorHex)\""))
    }
}
