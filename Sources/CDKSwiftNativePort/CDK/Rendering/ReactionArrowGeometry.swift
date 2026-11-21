#if canImport(CoreGraphics)
    import CoreGraphics
#endif
import Foundation

enum CDKReactionArrowGeometry {
    static func segments(
        direction: CDKReactionDirection, start: CGPoint, end: CGPoint, stroke: CGFloat, fontSize: CGFloat
    ) -> [CDKMetalDepictionScene.LineSegment] {
        let length = end.x - start.x
        let head = min(length / 3, max(stroke * 3, fontSize * 0.65))
        let halfHeight = head * 0.36
        let lane = max(stroke, halfHeight * 0.58)
        var segments: [CDKMetalDepictionScene.LineSegment] = []

        func line(_ a: CGPoint, _ b: CGPoint) {
            segments.append(.init(from: a, to: b, width: stroke, opacity: 0.95, color: .ink))
        }

        func arrowhead(at tip: CGPoint, towardRight: Bool, side: CGFloat? = nil, height: CGFloat) {
            let baseX = tip.x + (towardRight ? -head : head)
            if let side {
                line(CGPoint(x: baseX, y: tip.y + side * height), tip)
            } else {
                line(CGPoint(x: baseX, y: tip.y - height), tip)
                line(CGPoint(x: baseX, y: tip.y + height), tip)
            }
        }

        func shaft(offsetY: CGFloat = 0, endInset: CGFloat = 0) {
            line(
                CGPoint(x: start.x, y: start.y + offsetY),
                CGPoint(x: end.x - endInset, y: end.y + offsetY))
        }

        switch direction {
        case .forward, .backward, .resonance, .noGo:
            shaft()
            if direction == .backward || direction == .resonance {
                arrowhead(at: start, towardRight: false, height: halfHeight)
            }
            if direction != .backward {
                arrowhead(at: end, towardRight: true, height: halfHeight)
            }
            if direction == .noGo {
                let x = (start.x + end.x) / 2
                let d = halfHeight * 1.30
                line(CGPoint(x: x - d, y: start.y - d), CGPoint(x: x + d, y: start.y + d))
                line(CGPoint(x: x - d, y: start.y + d), CGPoint(x: x + d, y: start.y - d))
            }
        case .bidirectional:
            shaft(offsetY: -lane)
            shaft(offsetY: lane)
            arrowhead(at: CGPoint(x: end.x, y: end.y - lane), towardRight: true, side: -1, height: halfHeight)
            arrowhead(at: CGPoint(x: start.x, y: start.y + lane), towardRight: false, side: 1, height: halfHeight)
        case .retroSynthetic:
            shaft(offsetY: -lane, endInset: head * 0.5)
            shaft(offsetY: lane, endInset: head * 0.5)
            arrowhead(at: end, towardRight: true, height: lane * 2)
        case .undirected:
            shaft(offsetY: -lane)
            shaft(offsetY: lane)
        }
        return segments
    }
}
