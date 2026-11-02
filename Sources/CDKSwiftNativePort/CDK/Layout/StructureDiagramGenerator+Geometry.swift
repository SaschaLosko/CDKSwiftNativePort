import Foundation

#if canImport(CoreGraphics)
    import CoreGraphics
#endif

extension CDKStructureDiagramGenerator {
    static func boundingBox(component: Set<Int>, positions: [Int: CGPoint]) -> CGRect? {
        let points = component.compactMap { positions[$0] }
        guard let first = points.first else { return nil }
        var minX = first.x
        var minY = first.y
        var maxX = first.x
        var maxY = first.y
        for p in points {
            minX = min(minX, p.x)
            minY = min(minY, p.y)
            maxX = max(maxX, p.x)
            maxY = max(maxY, p.y)
        }
        return CGRect(
            x: minX, y: minY, width: max(0.0001, maxX - minX), height: max(0.0001, maxY - minY))
    }

    static func fanDirections(count: Int, baseAngle: CGFloat, totalSpread: CGFloat)
        -> [CGVector]
    {
        guard count > 0 else { return [] }
        if count == 1 { return [unitVector(angle: baseAngle)] }

        let start = baseAngle - totalSpread * 0.5
        let step = totalSpread / CGFloat(count - 1)
        return (0..<count).map { i in
            unitVector(angle: start + CGFloat(i) * step)
        }
    }

    static func radialDirections(count: Int, baseAngle: CGFloat) -> [CGVector] {
        guard count > 0 else { return [] }
        let step = (2 * .pi) / CGFloat(count)
        return (0..<count).map { i in
            unitVector(angle: baseAngle + CGFloat(i) * step)
        }
    }

    static func unitVector(angle: CGFloat) -> CGVector {
        CGVector(dx: cos(angle), dy: sin(angle))
    }

    static func normalize(_ v: CGVector) -> CGVector? {
        let len = hypot(v.dx, v.dy)
        guard len > 0.0001 else { return nil }
        return CGVector(dx: v.dx / len, dy: v.dy / len)
    }

    static func rotate(_ v: CGVector, by angle: CGFloat) -> CGVector {
        let c = cos(angle)
        let s = sin(angle)
        return CGVector(
            dx: v.dx * c - v.dy * s,
            dy: v.dx * s + v.dy * c
        )
    }

    static func segmentsIntersect(_ p1: CGPoint, _ p2: CGPoint, _ p3: CGPoint, _ p4: CGPoint)
        -> Bool
    {
        let o1 = orientation(p1, p2, p3)
        let o2 = orientation(p1, p2, p4)
        let o3 = orientation(p3, p4, p1)
        let o4 = orientation(p3, p4, p2)
        return o1 * o2 < 0 && o3 * o4 < 0
    }

    static func orientation(_ a: CGPoint, _ b: CGPoint, _ c: CGPoint) -> CGFloat {
        (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
    }

    static func reflect(point p: CGPoint, acrossLineFrom a: CGPoint, to b: CGPoint) -> CGPoint {
        let vx = b.x - a.x
        let vy = b.y - a.y
        let len2 = vx * vx + vy * vy
        guard len2 > 0.0001 else { return p }
        let t = ((p.x - a.x) * vx + (p.y - a.y) * vy) / len2
        let proj = CGPoint(x: a.x + t * vx, y: a.y + t * vy)
        return CGPoint(x: 2 * proj.x - p.x, y: 2 * proj.y - p.y)
    }

    static func segmentDistance(_ p1: CGPoint, _ p2: CGPoint, _ p3: CGPoint, _ p4: CGPoint)
        -> CGFloat
    {
        if segmentsIntersect(p1, p2, p3, p4) { return 0 }
        return min(
            pointSegmentDistance(p1, a: p3, b: p4),
            pointSegmentDistance(p2, a: p3, b: p4),
            pointSegmentDistance(p3, a: p1, b: p2),
            pointSegmentDistance(p4, a: p1, b: p2)
        )
    }

    static func pointSegmentDistance(_ p: CGPoint, a: CGPoint, b: CGPoint) -> CGFloat {
        let vx = b.x - a.x
        let vy = b.y - a.y
        let len2 = vx * vx + vy * vy
        guard len2 > 0.0001 else { return p.distance(to: a) }
        let t = max(0, min(1, ((p.x - a.x) * vx + (p.y - a.y) * vy) / len2))
        let proj = CGPoint(x: a.x + t * vx, y: a.y + t * vy)
        return p.distance(to: proj)
    }
}
