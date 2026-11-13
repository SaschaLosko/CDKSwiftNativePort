import Foundation
import XCTest
@testable import CDKSwiftNativePort

final class DepictionPreparationIsolationTests: XCTestCase {
    func testPreparedCacheDoesNotRetainViewportOrOtherMoleculeStyling() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C1CNCC(*)C1 |$;;;;;R1$,RG:_R1={OC},{Cl}|")
        let style = RenderStyle(showCarbons: true, showImplicitHydrogens: false)
        let canvas = CGRect(x: 0, y: 0, width: 480, height: 360)
        let original = CDKMetalDepictionSceneBuilder.build(
            molecule: molecule, style: style, canvasRect: canvas, zoom: 1, pan: .zero)
        XCTAssertFalse(original.labels.isEmpty)
        XCTAssertFalse(original.bondSegments.isEmpty)
        let transformed = CDKMetalDepictionSceneBuilder.build(
            molecule: molecule, style: style, canvasRect: canvas,
            zoom: 1.5, pan: CGSize(width: 20, height: -10), rotationDegrees: 70)
        XCTAssertNotEqual(original, transformed)

        // Exceed the preparation cache capacity with independent structures/styles.
        for length in 2...20 {
            let other = try parser.parseSmiles(String(repeating: "C", count: length) + "O")
            _ = CDKMetalDepictionSceneBuilder.build(
                molecule: other, style: RenderStyle(showCarbons: false),
                canvasRect: canvas, zoom: 0.7, pan: .zero)
        }
        let rebuilt = CDKMetalDepictionSceneBuilder.build(
            molecule: molecule, style: style, canvasRect: canvas, zoom: 1, pan: .zero)
        XCTAssertEqual(original, rebuilt)
    }
}
