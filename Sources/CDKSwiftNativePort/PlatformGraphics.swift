#if canImport(CoreGraphics)
@_exported import CoreGraphics
#else
import Foundation

public struct CGVector: Sendable, Hashable, Codable {
    public var dx: CGFloat
    public var dy: CGFloat

    public init(dx: CGFloat, dy: CGFloat) {
        self.dx = dx
        self.dy = dy
    }

    public static let zero = CGVector(dx: 0, dy: 0)
}

public struct CGAffineTransform: Sendable, Hashable, Codable {
    public var a: CGFloat
    public var b: CGFloat
    public var c: CGFloat
    public var d: CGFloat
    public var tx: CGFloat
    public var ty: CGFloat

    public init(a: CGFloat, b: CGFloat, c: CGFloat, d: CGFloat, tx: CGFloat, ty: CGFloat) {
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.tx = tx
        self.ty = ty
    }

    public init(translationX tx: CGFloat, y ty: CGFloat) {
        self.init(a: 1, b: 0, c: 0, d: 1, tx: tx, ty: ty)
    }

    public static let identity = CGAffineTransform(a: 1, b: 0, c: 0, d: 1, tx: 0, ty: 0)

    public func translatedBy(x: CGFloat, y: CGFloat) -> CGAffineTransform {
        concatenating(CGAffineTransform(translationX: x, y: y))
    }

    public func scaledBy(x: CGFloat, y: CGFloat) -> CGAffineTransform {
        concatenating(CGAffineTransform(a: x, b: 0, c: 0, d: y, tx: 0, ty: 0))
    }

    public func concatenating(_ other: CGAffineTransform) -> CGAffineTransform {
        CGAffineTransform(
            a: (a * other.a) + (c * other.b),
            b: (b * other.a) + (d * other.b),
            c: (a * other.c) + (c * other.d),
            d: (b * other.c) + (d * other.d),
            tx: (a * other.tx) + (c * other.ty) + tx,
            ty: (b * other.tx) + (d * other.ty) + ty
        )
    }

    fileprivate func inverted() -> CGAffineTransform? {
        let determinant = (a * d) - (b * c)
        guard abs(determinant) > 0.0000001 else { return nil }

        let inverseDeterminant = 1 / determinant
        return CGAffineTransform(
            a: d * inverseDeterminant,
            b: -b * inverseDeterminant,
            c: -c * inverseDeterminant,
            d: a * inverseDeterminant,
            tx: ((c * ty) - (d * tx)) * inverseDeterminant,
            ty: ((b * tx) - (a * ty)) * inverseDeterminant
        )
    }
}

public extension CGPoint {
    func applying(_ transform: CGAffineTransform) -> CGPoint {
        CGPoint(x: (transform.a * x) + (transform.c * y) + transform.tx,
                y: (transform.b * x) + (transform.d * y) + transform.ty)
    }
}

public extension CGRect {
    func intersection(_ other: CGRect) -> CGRect {
        guard intersects(other) else { return .null }
        let minX = Swift.max(self.minX, other.minX)
        let minY = Swift.max(self.minY, other.minY)
        let maxX = Swift.min(self.maxX, other.maxX)
        let maxY = Swift.min(self.maxY, other.maxY)
        return CGRect(x: minX, y: minY, width: max(0, maxX - minX), height: max(0, maxY - minY))
    }
}

enum CGPathFillRule {
    case winding
    case evenOdd
}

enum CGLineCap {
    case butt
    case round
    case square
}

enum CGLineJoin {
    case miter
    case round
    case bevel
}

struct CGPath {
    fileprivate let shapes: [CGRect]

    fileprivate init(shapes: [CGRect]) {
        self.shapes = shapes
    }

    init(rect: CGRect, transform: CGAffineTransform?) {
        let transformedRect = transform.map { rect.cdkApplying($0) } ?? rect
        self.init(shapes: [transformedRect])
    }

    init(roundedRect rect: CGRect,
         cornerWidth: CGFloat,
         cornerHeight: CGFloat,
         transform: CGAffineTransform?) {
        let transformedRect = transform.map { rect.cdkApplying($0) } ?? rect
        self.init(shapes: [transformedRect])
    }

    var boundingBoxOfPath: CGRect {
        shapes.reduce(.null) { partial, rect in
            partial.union(rect)
        }
    }

    func contains(_ point: CGPoint,
                  using rule: CGPathFillRule,
                  transform: CGAffineTransform) -> Bool {
        let probePoint: CGPoint
        if let inverse = transform.inverted() {
            probePoint = point.applying(inverse)
        } else {
            probePoint = point
        }

        return shapes.contains { $0.contains(probePoint) }
    }

    func copy(strokingWithWidth width: CGFloat,
              lineCap: CGLineCap,
              lineJoin: CGLineJoin,
              miterLimit: CGFloat) -> CGPath {
        let inset = max(0, width) * 0.5
        return CGPath(shapes: shapes.map { $0.insetBy(dx: -inset, dy: -inset) })
    }

    func copy() -> CGPath {
        self
    }
}

final class CGMutablePath {
    private var shapes: [CGRect] = []

    var isEmpty: Bool { shapes.isEmpty }

    func addPath(_ path: CGPath, transform: CGAffineTransform = .identity) {
        shapes.append(contentsOf: path.shapes.map { $0.cdkApplying(transform) })
    }

    func copy() -> CGPath {
        CGPath(shapes: shapes)
    }
}

private extension CGRect {
    func cdkApplying(_ transform: CGAffineTransform) -> CGRect {
        let points = [
            origin,
            CGPoint(x: maxX, y: minY),
            CGPoint(x: minX, y: maxY),
            CGPoint(x: maxX, y: maxY)
        ].map { $0.applying(transform) }

        let minX = points.map(\.x).min() ?? origin.x
        let minY = points.map(\.y).min() ?? origin.y
        let maxX = points.map(\.x).max() ?? self.maxX
        let maxY = points.map(\.y).max() ?? self.maxY
        return CGRect(x: minX, y: minY, width: maxX - minX, height: maxY - minY)
    }
}
#endif
