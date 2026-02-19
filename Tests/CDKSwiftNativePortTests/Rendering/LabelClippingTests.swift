import CoreGraphics
import XCTest
@testable import CDKSwiftNativePort

final class LabelClippingTests: XCTestCase {

    func testClipSegmentTrimsEndpointsInsideLabelRects() {
        let start = CGPoint(x: 0, y: 0)
        let end = CGPoint(x: 100, y: 0)
        let rects = [
            CGRect(x: -10, y: -8, width: 20, height: 16),
            CGRect(x: 90, y: -8, width: 20, height: 16)
        ]

        let clipped = CDKLabelClipping.clipSegmentEndpoints(start,
                                                            end,
                                                            labelRects: rects,
                                                            padding: 2.0)

        let trimmed = try? XCTUnwrap(clipped)
        XCTAssertNotNil(trimmed)
        if let trimmed {
            XCTAssertGreaterThan(trimmed.0.x, 10)
            XCTAssertLessThan(trimmed.1.x, 90)
        }
    }

    func testClipSegmentReturnsNilWhenFullyOccluded() {
        let start = CGPoint(x: 0, y: 0)
        let end = CGPoint(x: 10, y: 0)
        let rect = CGRect(x: -5, y: -5, width: 20, height: 10)

        let clipped = CDKLabelClipping.clipSegmentEndpoints(start,
                                                            end,
                                                            labelRects: [rect],
                                                            padding: 2.0)

        XCTAssertNil(clipped)
    }

    func testGlyphObstacleTrimsSegmentFromLabelInterior() {
        let obstacle = CDKLabelClipping.makeGlyphObstacle(text: "O",
                                                          center: CGPoint(x: 40, y: 20),
                                                          fontSize: 24,
                                                          padding: 2.0)

        let start = CGPoint(x: 40, y: 20)
        let end = CGPoint(x: 80, y: 20)
        let clipped = CDKLabelClipping.clipSegmentEndpoints(start,
                                                            end,
                                                            labelObstacles: [obstacle],
                                                            padding: 2.0)

        let trimmed = try? XCTUnwrap(clipped)
        XCTAssertNotNil(trimmed)
        if let trimmed {
            XCTAssertGreaterThan(trimmed.0.x, 42.0)
            XCTAssertLessThan(trimmed.1.x, 80.01)
        }
    }
}
