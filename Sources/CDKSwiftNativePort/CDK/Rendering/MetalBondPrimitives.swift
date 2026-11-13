#if canImport(CoreGraphics)
    import CoreGraphics
#endif
import Foundation

/// Immutable drawing metrics shared by bond, stereo and query glyph emission.
struct MetalBondPrimitives {
    let style: RenderStyle
    let viewportStrokeScale: CGFloat
    let zoomedStrokeScale: CGFloat
    let baseBondWidth: CGFloat
    let doubleBondSeparation: CGFloat
    let stereoSolidWedgeHalfWidth: CGFloat
    let stereoHashedWedgeHalfWidth: CGFloat
    let stereoScale: CGFloat
    let hashedStereoScale: CGFloat

    func appendSegment(
        _ from: CGPoint,
        _ to: CGPoint,
        width: CGFloat,
        opacity: CGFloat,
        color: CDKRenderColor,
        into segments: inout [CDKMetalDepictionScene.LineSegment]
    ) {
        segments.append(
            CDKMetalDepictionScene.LineSegment(
                from: from,
                to: to,
                width: max(0.75, width),
                opacity: opacity.clamped(to: 0.05...1.0),
                color: color))
    }

    private func normalize(dx: CGFloat, dy: CGFloat) -> (CGFloat, CGFloat)? {
        let l = hypot(dx, dy)
        guard l > 0.0001 else { return nil }
        return (dx / l, dy / l)
    }

    func appendAromaticInnerLine(
        from p1: CGPoint,
        to p2: CGPoint,
        center: CGPoint,
        color: CDKRenderColor,
        into segments: inout [CDKMetalDepictionScene.LineSegment]
    ) {
        let mid = CGPoint(x: (p1.x + p2.x) * 0.5, y: (p1.y + p2.y) * 0.5)
        guard let inward = normalize(dx: center.x - mid.x, dy: center.y - mid.y) else { return }
        let trim = CGFloat(0.16)
        let inset = max(1.6, style.bondWidth * 0.95 * viewportStrokeScale) * zoomedStrokeScale
        let a = CGPoint(
            x: p1.x + (p2.x - p1.x) * trim + inward.0 * inset,
            y: p1.y + (p2.y - p1.y) * trim + inward.1 * inset)
        let b = CGPoint(
            x: p2.x + (p1.x - p2.x) * trim + inward.0 * inset,
            y: p2.y + (p1.y - p2.y) * trim + inward.1 * inset)
        appendSegment(
            a,
            b,
            width: max(1.0, baseBondWidth * 0.8),
            opacity: 0.76,
            color: color,
            into: &segments)
    }

    func appendConjugatedDoubleInnerLine(
        from p1: CGPoint,
        to p2: CGPoint,
        px: CGFloat,
        py: CGFloat,
        centers: [CGPoint],
        paritySeed: Int,
        color: CDKRenderColor,
        into segments: inout [CDKMetalDepictionScene.LineSegment]
    ) {
        let mid = CGPoint(x: (p1.x + p2.x) * 0.5, y: (p1.y + p2.y) * 0.5)
        let summed = centers.reduce((dx: CGFloat(0), dy: CGFloat(0))) { partial, c in
            (dx: partial.dx + (c.x - mid.x), dy: partial.dy + (c.y - mid.y))
        }
        let preferred = normalize(dx: summed.dx, dy: summed.dy)
        let sign: CGFloat
        if let preferred {
            sign = (preferred.0 * px + preferred.1 * py) >= 0 ? 1 : -1
        } else {
            sign = (paritySeed % 2 == 0) ? 1 : -1
        }

        let trim = CGFloat(0.15)
        let inset = doubleBondSeparation
        let a = CGPoint(
            x: p1.x + (p2.x - p1.x) * trim + px * inset * sign,
            y: p1.y + (p2.y - p1.y) * trim + py * inset * sign)
        let b = CGPoint(
            x: p2.x + (p1.x - p2.x) * trim + px * inset * sign,
            y: p2.y + (p1.y - p2.y) * trim + py * inset * sign)
        appendSegment(
            a,
            b,
            width: max(1.05, baseBondWidth * 0.85),
            opacity: 0.88,
            color: color,
            into: &segments)
    }

    func appendSolidWedge(
        from start: CGPoint,
        to end: CGPoint,
        endInset: CGFloat,
        color: CDKRenderColor,
        into segments: inout [CDKMetalDepictionScene.LineSegment]
    ) {
        let rawDX = end.x - start.x
        let rawDY = end.y - start.y
        let rawLen = max(0.0001, hypot(rawDX, rawDY))
        let inset = min(endInset, rawLen * 0.32)
        let trimmedEnd = CGPoint(
            x: end.x - (rawDX / rawLen) * inset,
            y: end.y - (rawDY / rawLen) * inset)
        let dx = trimmedEnd.x - start.x
        let dy = trimmedEnd.y - start.y
        let len = max(0.0001, hypot(dx, dy))
        let px = -dy / len
        let py = dx / len

        let maxHalfWidth = min(stereoSolidWedgeHalfWidth, len * 0.42)
        let stripeCount = max(16, min(40, Int(len / max(1.6, baseBondWidth * 0.28))))
        let step = len / CGFloat(stripeCount + 1)
        for i in 1...stripeCount {
            let t = CGFloat(i) / CGFloat(stripeCount + 1)
            let cx = start.x + dx * t
            let cy = start.y + dy * t
            let halfTick = max(baseBondWidth * 0.10 * stereoScale, t * maxHalfWidth)
            let left = CGPoint(x: cx + px * halfTick, y: cy + py * halfTick)
            let right = CGPoint(x: cx - px * halfTick, y: cy - py * halfTick)
            appendSegment(
                left,
                right,
                width: max(step * 1.35, baseBondWidth * 0.95 * stereoScale),
                opacity: 0.985,
                color: color,
                into: &segments)
        }
    }

    func appendHashedWedge(
        from start: CGPoint,
        to end: CGPoint,
        endInset: CGFloat,
        color: CDKRenderColor,
        into segments: inout [CDKMetalDepictionScene.LineSegment]
    ) {
        let rawDX = end.x - start.x
        let rawDY = end.y - start.y
        let rawLen = max(0.0001, hypot(rawDX, rawDY))
        let inset = min(endInset, rawLen * 0.32)
        let trimmedEnd = CGPoint(
            x: end.x - (rawDX / rawLen) * inset,
            y: end.y - (rawDY / rawLen) * inset)
        let dx = trimmedEnd.x - start.x
        let dy = trimmedEnd.y - start.y
        let len = max(0.0001, hypot(dx, dy))
        let px = -dy / len
        let py = dx / len

        let maxHalfWidth = min(stereoHashedWedgeHalfWidth, len * 0.34)
        let tickCount = max(9, min(22, Int(len / max(4.8, baseBondWidth * 1.55))))
        for i in 1...tickCount {
            let t = CGFloat(i) / CGFloat(tickCount + 1)
            let cx = start.x + dx * t
            let cy = start.y + dy * t
            let halfTick = max(baseBondWidth * 0.12 * stereoScale * hashedStereoScale, t * maxHalfWidth)
            let left = CGPoint(x: cx + px * halfTick, y: cy + py * halfTick)
            let right = CGPoint(x: cx - px * halfTick, y: cy - py * halfTick)
            appendSegment(
                left,
                right,
                width: max(0.75, baseBondWidth * 0.52 * stereoScale * hashedStereoScale),
                opacity: 0.96,
                color: color,
                into: &segments)
        }
    }

    func appendAromaticCircle(
        points: [CGPoint],
        color: CDKRenderColor,
        width: CGFloat,
        opacity: CGFloat,
        into segments: inout [CDKMetalDepictionScene.LineSegment]
    ) {
        guard points.count >= 5 else { return }
        let center = points.reduce(CGPoint.zero) { acc, p in
            CGPoint(
                x: acc.x + p.x / CGFloat(points.count),
                y: acc.y + p.y / CGFloat(points.count))
        }
        let minRadius = points.map { hypot($0.x - center.x, $0.y - center.y) }.min() ?? 0
        guard minRadius > 1 else { return }
        let radius = minRadius * 0.53
        let segmentsCount = 40
        var previous: CGPoint?
        for i in 0...segmentsCount {
            let t = CGFloat(i) / CGFloat(segmentsCount)
            let angle = t * 2 * .pi
            let point = CGPoint(
                x: center.x + cos(angle) * radius,
                y: center.y + sin(angle) * radius)
            if let previous {
                appendSegment(previous, point, width: width, opacity: opacity, color: color, into: &segments)
            }
            previous = point
        }
    }

    func appendAttachmentPointGlyph(
        center: CGPoint,
        other: CGPoint,
        color: CDKRenderColor,
        into segments: inout [CDKMetalDepictionScene.LineSegment]
    ) {
        for segment in CDKMarkushRendering.attachmentPointSegments(
            center: center,
            other: other,
            style: style,
            scale: viewportStrokeScale * zoomedStrokeScale)
        {
            appendSegment(
                segment.from,
                segment.to,
                width: max(1.0, baseBondWidth * 0.92),
                opacity: 0.95,
                color: color,
                into: &segments)
        }
    }

    func appendDashedSegment(
        _ from: CGPoint,
        _ to: CGPoint,
        width: CGFloat,
        opacity: CGFloat,
        color: CDKRenderColor,
        into segments: inout [CDKMetalDepictionScene.LineSegment]
    ) {
        let dx = to.x - from.x
        let dy = to.y - from.y
        let length = hypot(dx, dy)
        guard length > 0.0001 else { return }
        let direction = CGVector(dx: dx / length, dy: dy / length)
        let dashLength = max(width * 3.6, 8.0)
        let gapLength = max(width * 1.8, 4.8)
        var traveled: CGFloat = 0
        while traveled < length {
            let start = traveled
            let end = min(length, traveled + dashLength)
            let p1 = CGPoint(x: from.x + direction.dx * start, y: from.y + direction.dy * start)
            let p2 = CGPoint(x: from.x + direction.dx * end, y: from.y + direction.dy * end)
            appendSegment(p1, p2, width: width, opacity: opacity, color: color, into: &segments)
            traveled += dashLength + gapLength
        }
    }

    func appendGlowSegment(
        _ from: CGPoint,
        _ to: CGPoint,
        width: CGFloat,
        into segments: inout [CDKMetalDepictionScene.LineSegment]
    ) {
        let glowWidth = max(width * 2.35, width + max(1.6, style.bondWidth * 1.8))
        appendSegment(
            from,
            to,
            width: glowWidth,
            opacity: 0.34,
            color: .outerGlowHighlight,
            into: &segments)
    }

    func appendQuerySecondaryLine(
        from p1: CGPoint,
        to p2: CGPoint,
        px: CGFloat,
        py: CGFloat,
        centers: [CGPoint],
        paritySeed: Int,
        color: CDKRenderColor,
        dashed: Bool,
        into segments: inout [CDKMetalDepictionScene.LineSegment]
    ) {
        let mid = CGPoint(x: (p1.x + p2.x) * 0.5, y: (p1.y + p2.y) * 0.5)
        let summed = centers.reduce((dx: CGFloat(0), dy: CGFloat(0))) { partial, c in
            (dx: partial.dx + (c.x - mid.x), dy: partial.dy + (c.y - mid.y))
        }
        let preferred = normalize(dx: summed.dx, dy: summed.dy)
        let sign: CGFloat
        if let preferred {
            sign = (preferred.0 * px + preferred.1 * py) >= 0 ? 1 : -1
        } else {
            sign = (paritySeed % 2 == 0) ? 1 : -1
        }

        let trim = CGFloat(0.15)
        let inset = doubleBondSeparation
        let a = CGPoint(
            x: p1.x + (p2.x - p1.x) * trim + px * inset * sign,
            y: p1.y + (p2.y - p1.y) * trim + py * inset * sign)
        let b = CGPoint(
            x: p2.x + (p1.x - p2.x) * trim + px * inset * sign,
            y: p2.y + (p1.y - p2.y) * trim + py * inset * sign)
        if dashed {
            appendDashedSegment(
                a,
                b,
                width: max(1.0, baseBondWidth * 0.82),
                opacity: 0.92,
                color: color,
                into: &segments)
        } else {
            appendSegment(
                a,
                b,
                width: max(1.0, baseBondWidth * 0.82),
                opacity: 0.92,
                color: color,
                into: &segments)
        }
    }

}

private extension Comparable {
    func clamped(to range: ClosedRange<Self>) -> Self {
        min(max(self, range.lowerBound), range.upperBound)
    }
}
