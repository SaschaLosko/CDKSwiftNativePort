#if canImport(CoreGraphics)
    import CoreGraphics
#endif
import Foundation

enum CDKReactionArrowGeometry {
    static func build(
        direction: CDKReactionDirection, start: CGPoint, end: CGPoint, height: CGFloat
    ) -> (strokes: [CDKMetalDepictionScene.LineSegment], fills: [CDKMetalDepictionScene.FilledPolygon]) {
        let thickness = height / 3
        let stroke = height / 14
        let centerY = start.y
        let left = start.x
        let right = end.x
        var strokes: [CDKMetalDepictionScene.LineSegment] = []
        var fills: [CDKMetalDepictionScene.FilledPolygon] = []

        func line(_ x1: CGFloat, _ y1: CGFloat, _ x2: CGFloat, _ y2: CGFloat) {
            strokes.append(
                .init(
                    from: CGPoint(x: x1, y: centerY + y1), to: CGPoint(x: x2, y: centerY + y2),
                    width: stroke, opacity: 1, color: .ink))
        }

        func fill(_ vertices: [(CGFloat, CGFloat)]) {
            fills.append(.init(points: vertices.map { CGPoint(x: $0.0, y: centerY + $0.1) }, color: .ink))
        }

        func fullHead(tip: CGFloat, inward: CGFloat) {
            fill([
                (tip, 0),
                (tip + inward * height, -thickness),
                (tip + inward * height * 0.8, 0),
                (tip + inward * height, thickness),
            ])
        }

        // CDK 2.13 uses closed, notched heads except for the outlined
        // retrosynthetic head, the half-heads, and the headless parallel lines.
        switch direction {
        case .forward, .backward, .resonance, .noGo:
            let onLeft = direction == .backward || direction == .resonance
            let onRight = direction != .backward
            line(left + (onLeft ? height / 2 : 0), 0, right - (onRight ? height / 2 : 0), 0)
            if onRight { fullHead(tip: right, inward: -1) }
            if onLeft { fullHead(tip: left, inward: 1) }
            if direction == .noGo {
                let middle = (left + right) / 2
                line(middle - thickness, thickness, middle + thickness, -thickness)
                line(middle - thickness, -thickness, middle + thickness, thickness)
            }
        case .bidirectional:
            let lane = thickness / 2
            let tipOffset = lane - stroke / 2
            line(left, -lane, right - height / 2, -lane)
            line(left + height / 2, lane, right, lane)
            fill([(right, -tipOffset), (right - height, -1.5 * thickness), (right - 0.8 * height, -tipOffset)])
            fill([(left, tipOffset), (left + height, 1.5 * thickness), (left + 0.8 * height, tipOffset)])
        case .retroSynthetic:
            line(left, thickness, right - height / 2, thickness)
            line(left, -thickness, right - height / 2, -thickness)
            line(right - height, 2 * thickness, right, 0)
            line(right, 0, right - height, -2 * thickness)
        case .undirected:
            let lane = thickness / 2
            line(left + thickness, lane, right - 2 * thickness, lane)
            line(left + thickness, -lane, right - 2 * thickness, -lane)
        }
        return (strokes, fills)
    }
}
