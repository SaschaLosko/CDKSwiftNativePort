import CoreGraphics
import Foundation

public struct CDKMetalDepictionScene: Hashable {
    public struct LineSegment: Hashable {
        public let from: CGPoint
        public let to: CGPoint
        public let width: CGFloat
        public let opacity: CGFloat
        public let color: CDKRenderColor

        public init(from: CGPoint, to: CGPoint, width: CGFloat, opacity: CGFloat, color: CDKRenderColor) {
            self.from = from
            self.to = to
            self.width = width
            self.opacity = opacity
            self.color = color
        }
    }

    public struct AtomLabel: Identifiable, Hashable {
        public let id: Int
        public let text: String
        public let position: CGPoint
        public let fontSize: CGFloat
        public let aromatic: Bool
        public let color: CDKRenderColor

        public init(id: Int, text: String, position: CGPoint, fontSize: CGFloat, aromatic: Bool, color: CDKRenderColor) {
            self.id = id
            self.text = text
            self.position = position
            self.fontSize = fontSize
            self.aromatic = aromatic
            self.color = color
        }
    }

    public let gridSegments: [LineSegment]
    public let bondSegments: [LineSegment]
    public let labels: [AtomLabel]

    public init(gridSegments: [LineSegment], bondSegments: [LineSegment], labels: [AtomLabel]) {
        self.gridSegments = gridSegments
        self.bondSegments = bondSegments
        self.labels = labels
    }
}

private struct MetalRenderEdgeKey: Hashable {
    let a: Int
    let b: Int

    init(_ u: Int, _ v: Int) {
        a = min(u, v)
        b = max(u, v)
    }
}

private struct MetalPreparedDepictionKey: Hashable {
    let molecule: Molecule
    let style: RenderStyle
}

private struct MetalPreparedDepictionData {
    let depictionMolecule: Molecule
    let boundingBox: CGRect?
    let degreeByAtomID: [Int: Int]
    let aromaticRings: [[Int]]
    let aromaticBondIDs: Set<Int>
    let conjugatedRings: [[Int]]
}

private enum MetalPreparedDepictionCache {
    private static let maxEntries = 16
    private static let lock = NSLock()
    private static var entries: [MetalPreparedDepictionKey: MetalPreparedDepictionData] = [:]
    private static var lruOrder: [MetalPreparedDepictionKey] = []

    static func preparedData(molecule: Molecule, style: RenderStyle) -> MetalPreparedDepictionData {
        let key = MetalPreparedDepictionKey(molecule: molecule, style: style)

        lock.lock()
        if let cached = entries[key] {
            touchLRU(key)
            lock.unlock()
            return cached
        }
        lock.unlock()

        let depictionMolecule = CDKDepictionPreprocessor.prepareForRendering(molecule: molecule, style: style)
        let aromaticRings = depictionMolecule.aromaticDisplayRings()
        let aromaticBondIDs = depictionMolecule.aromaticDisplayBondIDs()
        let conjugatedRings = depictionMolecule.simpleCycles(maxSize: 10).filter { ring in
            guard ring.count >= 5 else { return false }
            let ringBonds: [Bond] = (0..<ring.count).compactMap { idx in
                depictionMolecule.bond(between: ring[idx], and: ring[(idx + 1) % ring.count])
            }
            guard ringBonds.count == ring.count else { return false }
            let piCount = ringBonds.reduce(0) { partial, bond in
                partial + ((bond.order == .double || bond.order == .aromatic) ? 1 : 0)
            }
            return piCount >= 2 && ringBonds.allSatisfy { $0.order != .triple }
        }

        var degreeByAtomID: [Int: Int] = [:]
        for bond in depictionMolecule.bonds {
            degreeByAtomID[bond.a1, default: 0] += 1
            degreeByAtomID[bond.a2, default: 0] += 1
        }

        let prepared = MetalPreparedDepictionData(depictionMolecule: depictionMolecule,
                                                  boundingBox: depictionMolecule.boundingBox(),
                                                  degreeByAtomID: degreeByAtomID,
                                                  aromaticRings: aromaticRings,
                                                  aromaticBondIDs: aromaticBondIDs,
                                                  conjugatedRings: conjugatedRings)

        lock.lock()
        if let cached = entries[key] {
            touchLRU(key)
            lock.unlock()
            return cached
        }
        entries[key] = prepared
        lruOrder.append(key)
        while lruOrder.count > maxEntries {
            let stale = lruOrder.removeFirst()
            entries.removeValue(forKey: stale)
        }
        lock.unlock()
        return prepared
    }

    private static func touchLRU(_ key: MetalPreparedDepictionKey) {
        guard let index = lruOrder.firstIndex(of: key) else { return }
        if index == lruOrder.count - 1 { return }
        let moved = lruOrder.remove(at: index)
        lruOrder.append(moved)
    }
}

public enum CDKMetalDepictionSceneBuilder {
    public static func build(molecule: Molecule,
                             style: RenderStyle,
                             canvasRect: CGRect,
                             zoom: CGFloat,
                             pan: CGSize,
                             rotationDegrees: CGFloat = 0) -> CDKMetalDepictionScene {
        let prepared = MetalPreparedDepictionCache.preparedData(molecule: molecule, style: style)
        let depictionMolecule = prepared.depictionMolecule
        guard let box = prepared.boundingBox, canvasRect.width > 1, canvasRect.height > 1 else {
            return CDKMetalDepictionScene(gridSegments: [], bondSegments: [], labels: [])
        }

        let gridStep: CGFloat = 40
        var gridSegments: [CDKMetalDepictionScene.LineSegment] = []
        var x: CGFloat = 0
        while x <= canvasRect.maxX {
            gridSegments.append(CDKMetalDepictionScene.LineSegment(from: CGPoint(x: x, y: canvasRect.minY),
                                                                   to: CGPoint(x: x, y: canvasRect.maxY),
                                                                   width: 1,
                                                                   opacity: 0.05,
                                                                   color: .grid))
            x += gridStep
        }

        var y: CGFloat = 0
        while y <= canvasRect.maxY {
            gridSegments.append(CDKMetalDepictionScene.LineSegment(from: CGPoint(x: canvasRect.minX, y: y),
                                                                   to: CGPoint(x: canvasRect.maxX, y: y),
                                                                   width: 1,
                                                                   opacity: 0.05,
                                                                   color: .grid))
            y += gridStep
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
        let baseTransform = CGAffineTransform.identity
            .translatedBy(x: center.x, y: center.y)
            .scaledBy(x: scale, y: -scale)
            .translatedBy(x: -box.midX, y: -box.midY)

        let viewportCenter = CGPoint(x: canvasRect.midX, y: canvasRect.midY)
        let normalizedRotation = rotationDegrees.truncatingRemainder(dividingBy: 360)
        let rotationRadians = normalizedRotation * (.pi / 180)
        let rotationCos = cos(rotationRadians)
        let rotationSin = sin(rotationRadians)

        func applyViewportTransform(_ point: CGPoint) -> CGPoint {
            let rx = point.x - viewportCenter.x
            let ry = point.y - viewportCenter.y
            let rotated = CGPoint(x: (rx * rotationCos) - (ry * rotationSin) + viewportCenter.x,
                                  y: (rx * rotationSin) + (ry * rotationCos) + viewportCenter.y)
            let zx = (rotated.x - viewportCenter.x) * zoom + viewportCenter.x + pan.width
            let zy = (rotated.y - viewportCenter.y) * zoom + viewportCenter.y + pan.height
            return CGPoint(x: zx, y: zy)
        }

        var basePositionByAtomID: [Int: CGPoint] = [:]
        var positionByAtomID: [Int: CGPoint] = [:]
        for atom in depictionMolecule.atoms {
            let basePosition = atom.position.applying(baseTransform)
            basePositionByAtomID[atom.id] = basePosition
            positionByAtomID[atom.id] = applyViewportTransform(basePosition)
        }

        let degree = prepared.degreeByAtomID
        let aromaticRings = prepared.aromaticRings
        let aromaticBondIDs = prepared.aromaticBondIDs
        var aromaticEdgeCenters: [MetalRenderEdgeKey: [CGPoint]] = [:]
        var conjugatedDoubleEdgeCenters: [MetalRenderEdgeKey: [CGPoint]] = [:]
        let baseBondWidth = max(1.2, style.bondWidth * zoom)
        let stereoSolidWedgeHalfWidth = max(7.2, style.bondWidth * 3.9) * zoom
        let stereoHashedWedgeHalfWidth = max(7.6, style.bondWidth * 4.15) * zoom
        let doubleBondSeparation = max(3.4, style.bondWidth * 2.25) * zoom
        let doubleBondHalfSeparation = doubleBondSeparation * 0.5
        let tripleBondOffset = max(3.1, style.bondWidth * 2.2) * zoom

        for ring in aromaticRings where ring.count >= 3 {
            let ringPoints = ring.compactMap { positionByAtomID[$0] }
            guard ringPoints.count == ring.count else { continue }

            let center = ringPoints.reduce(CGPoint.zero) { acc, p in
                CGPoint(x: acc.x + p.x / CGFloat(ringPoints.count),
                        y: acc.y + p.y / CGFloat(ringPoints.count))
            }

            for i in 0..<ring.count {
                let a = ring[i]
                let b = ring[(i + 1) % ring.count]
                guard let bond = depictionMolecule.bond(between: a, and: b), aromaticBondIDs.contains(bond.id) else { continue }
                aromaticEdgeCenters[MetalRenderEdgeKey(a, b), default: []].append(center)
            }
        }

        let conjugatedRings = prepared.conjugatedRings
        for ring in conjugatedRings where ring.count >= 5 {
            let ringPoints = ring.compactMap { positionByAtomID[$0] }
            guard ringPoints.count == ring.count else { continue }

            let center = ringPoints.reduce(CGPoint.zero) { acc, p in
                CGPoint(x: acc.x + p.x / CGFloat(ringPoints.count),
                        y: acc.y + p.y / CGFloat(ringPoints.count))
            }

            for i in 0..<ring.count {
                let a = ring[i]
                let b = ring[(i + 1) % ring.count]
                guard let bond = depictionMolecule.bond(between: a, and: b) else { continue }
                if bond.order == .double || aromaticBondIDs.contains(bond.id) {
                    conjugatedDoubleEdgeCenters[MetalRenderEdgeKey(a, b), default: []].append(center)
                }
            }
        }

        func appendSegment(_ from: CGPoint,
                           _ to: CGPoint,
                           width: CGFloat,
                           opacity: CGFloat,
                           color: CDKRenderColor,
                           into segments: inout [CDKMetalDepictionScene.LineSegment]) {
            segments.append(CDKMetalDepictionScene.LineSegment(from: from,
                                                               to: to,
                                                               width: max(0.75, width),
                                                               opacity: opacity.clamped(to: 0.05...1.0),
                                                               color: color))
        }

        func normalize(dx: CGFloat, dy: CGFloat) -> (CGFloat, CGFloat)? {
            let l = hypot(dx, dy)
            guard l > 0.0001 else { return nil }
            return (dx / l, dy / l)
        }

        func appendAromaticInnerLine(from p1: CGPoint,
                                     to p2: CGPoint,
                                     center: CGPoint,
                                     color: CDKRenderColor,
                                     into segments: inout [CDKMetalDepictionScene.LineSegment]) {
            let mid = CGPoint(x: (p1.x + p2.x) * 0.5, y: (p1.y + p2.y) * 0.5)
            guard let inward = normalize(dx: center.x - mid.x, dy: center.y - mid.y) else { return }
            let trim = CGFloat(0.16)
            let inset = max(1.6, style.bondWidth * 0.95) * zoom
            let a = CGPoint(x: p1.x + (p2.x - p1.x) * trim + inward.0 * inset,
                            y: p1.y + (p2.y - p1.y) * trim + inward.1 * inset)
            let b = CGPoint(x: p2.x + (p1.x - p2.x) * trim + inward.0 * inset,
                            y: p2.y + (p1.y - p2.y) * trim + inward.1 * inset)
            appendSegment(a,
                          b,
                          width: max(1.0, baseBondWidth * 0.8),
                          opacity: 0.76,
                          color: color,
                          into: &segments)
        }

        func appendConjugatedDoubleInnerLine(from p1: CGPoint,
                                             to p2: CGPoint,
                                             px: CGFloat,
                                             py: CGFloat,
                                             centers: [CGPoint],
                                             paritySeed: Int,
                                             color: CDKRenderColor,
                                             into segments: inout [CDKMetalDepictionScene.LineSegment]) {
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
            let a = CGPoint(x: p1.x + (p2.x - p1.x) * trim + px * inset * sign,
                            y: p1.y + (p2.y - p1.y) * trim + py * inset * sign)
            let b = CGPoint(x: p2.x + (p1.x - p2.x) * trim + px * inset * sign,
                            y: p2.y + (p1.y - p2.y) * trim + py * inset * sign)
            appendSegment(a,
                          b,
                          width: max(1.05, baseBondWidth * 0.85),
                          opacity: 0.88,
                          color: color,
                          into: &segments)
        }

        func appendSolidWedge(from start: CGPoint,
                              to end: CGPoint,
                              color: CDKRenderColor,
                              into segments: inout [CDKMetalDepictionScene.LineSegment]) {
            let dx = end.x - start.x
            let dy = end.y - start.y
            let len = max(0.0001, hypot(dx, dy))
            let px = -dy / len
            let py = dx / len

            let stripeCount = max(12, Int(len / 8.0))
            for i in 1...stripeCount {
                let t = CGFloat(i) / CGFloat(stripeCount + 1)
                let cx = start.x + dx * t
                let cy = start.y + dy * t
                let halfTick = max(baseBondWidth * 0.18, t * stereoSolidWedgeHalfWidth)
                let left = CGPoint(x: cx + px * halfTick, y: cy + py * halfTick)
                let right = CGPoint(x: cx - px * halfTick, y: cy - py * halfTick)
                appendSegment(left,
                              right,
                              width: max(1.1, baseBondWidth * 0.94),
                              opacity: 0.97,
                              color: color,
                              into: &segments)
            }
        }

        func appendHashedWedge(from start: CGPoint,
                               to end: CGPoint,
                               color: CDKRenderColor,
                               into segments: inout [CDKMetalDepictionScene.LineSegment]) {
            let dx = end.x - start.x
            let dy = end.y - start.y
            let len = max(0.0001, hypot(dx, dy))
            let px = -dy / len
            let py = dx / len

            let tickCount = 9
            for i in 1...tickCount {
                let t = CGFloat(i) / CGFloat(tickCount + 1)
                let cx = start.x + dx * t
                let cy = start.y + dy * t
                let halfTick = max(baseBondWidth * 0.30, t * stereoHashedWedgeHalfWidth)
                let left = CGPoint(x: cx + px * halfTick, y: cy + py * halfTick)
                let right = CGPoint(x: cx - px * halfTick, y: cy - py * halfTick)
                appendSegment(left,
                              right,
                              width: max(1.15, baseBondWidth * 0.86),
                              opacity: 0.96,
                              color: color,
                              into: &segments)
            }
        }

        func appendAromaticCircle(points: [CGPoint],
                                  color: CDKRenderColor,
                                  width: CGFloat,
                                  opacity: CGFloat,
                                  into segments: inout [CDKMetalDepictionScene.LineSegment]) {
            guard points.count >= 5 else { return }
            let center = points.reduce(CGPoint.zero) { acc, p in
                CGPoint(x: acc.x + p.x / CGFloat(points.count),
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
                let point = CGPoint(x: center.x + cos(angle) * radius,
                                    y: center.y + sin(angle) * radius)
                if let previous {
                    appendSegment(previous, point, width: width, opacity: opacity, color: color, into: &segments)
                }
                previous = point
            }
        }

        var bondSegments: [CDKMetalDepictionScene.LineSegment] = []
        bondSegments.reserveCapacity(max(16, depictionMolecule.bonds.count * 3))

        for bond in depictionMolecule.bonds {
            guard let p1 = positionByAtomID[bond.a1], let p2 = positionByAtomID[bond.a2] else { continue }
            let bondColor = CDKRenderingStyleResolver.bondColor(for: bond, molecule: depictionMolecule, style: style)

            let dx = p2.x - p1.x
            let dy = p2.y - p1.y
            let len = max(0.0001, hypot(dx, dy))
            let px = -dy / len
            let py = dx / len
            let edgeKey = MetalRenderEdgeKey(bond.a1, bond.a2)
            let aromaticCenters = aromaticEdgeCenters[edgeKey] ?? []
            let conjugatedCenters = conjugatedDoubleEdgeCenters[edgeKey] ?? []

            if bond.order == .aromatic {
                appendSegment(p1,
                              p2,
                              width: max(1.1, baseBondWidth * 0.88),
                              opacity: 0.80,
                              color: bondColor,
                              into: &bondSegments)
                if style.aromaticDisplayMode == .innerLine,
                   aromaticCenters.count == 1,
                   let center = aromaticCenters.first {
                    appendAromaticInnerLine(from: p1,
                                            to: p2,
                                            center: center,
                                            color: bondColor,
                                            into: &bondSegments)
                }
                continue
            }

            switch bond.order {
            case .single:
                switch bond.stereo {
                case .up:
                    appendSolidWedge(from: p1, to: p2, color: bondColor, into: &bondSegments)
                case .upReversed:
                    appendSolidWedge(from: p2, to: p1, color: bondColor, into: &bondSegments)
                case .down:
                    appendHashedWedge(from: p1, to: p2, color: bondColor, into: &bondSegments)
                case .downReversed:
                    appendHashedWedge(from: p2, to: p1, color: bondColor, into: &bondSegments)
                case .either:
                    appendSegment(p1, p2, width: max(1.0, baseBondWidth * 0.72), opacity: 0.78, color: bondColor, into: &bondSegments)
                case .none:
                    appendSegment(p1, p2, width: baseBondWidth, opacity: 0.95, color: bondColor, into: &bondSegments)
                }
            case .double:
                if !conjugatedCenters.isEmpty {
                    appendSegment(p1, p2, width: baseBondWidth, opacity: 0.95, color: bondColor, into: &bondSegments)
                    appendConjugatedDoubleInnerLine(from: p1,
                                                    to: p2,
                                                    px: px,
                                                    py: py,
                                                    centers: conjugatedCenters,
                                                    paritySeed: bond.id,
                                                    color: bondColor,
                                                    into: &bondSegments)
                } else {
                    appendSegment(CGPoint(x: p1.x + px * doubleBondHalfSeparation, y: p1.y + py * doubleBondHalfSeparation),
                                  CGPoint(x: p2.x + px * doubleBondHalfSeparation, y: p2.y + py * doubleBondHalfSeparation),
                                  width: baseBondWidth,
                                  opacity: 0.95,
                                  color: bondColor,
                                  into: &bondSegments)
                    appendSegment(CGPoint(x: p1.x - px * doubleBondHalfSeparation, y: p1.y - py * doubleBondHalfSeparation),
                                  CGPoint(x: p2.x - px * doubleBondHalfSeparation, y: p2.y - py * doubleBondHalfSeparation),
                                  width: baseBondWidth,
                                  opacity: 0.95,
                                  color: bondColor,
                                  into: &bondSegments)
                }
            case .triple:
                appendSegment(p1, p2, width: baseBondWidth, opacity: 0.95, color: bondColor, into: &bondSegments)
                appendSegment(CGPoint(x: p1.x + px * tripleBondOffset, y: p1.y + py * tripleBondOffset),
                              CGPoint(x: p2.x + px * tripleBondOffset, y: p2.y + py * tripleBondOffset),
                              width: baseBondWidth,
                              opacity: 0.95,
                              color: bondColor,
                              into: &bondSegments)
                appendSegment(CGPoint(x: p1.x - px * tripleBondOffset, y: p1.y - py * tripleBondOffset),
                              CGPoint(x: p2.x - px * tripleBondOffset, y: p2.y - py * tripleBondOffset),
                              width: baseBondWidth,
                              opacity: 0.95,
                              color: bondColor,
                              into: &bondSegments)
            case .aromatic:
                appendSegment(p1,
                              p2,
                              width: max(1.1, baseBondWidth * 0.88),
                              opacity: 0.80,
                              color: bondColor,
                              into: &bondSegments)
            }
        }

        if style.aromaticDisplayMode == .circle {
            for ring in aromaticRings where ring.count >= 5 {
                let ringPoints = ring.compactMap { positionByAtomID[$0] }
                guard ringPoints.count == ring.count else { continue }
                let ringColor = CDKRenderingStyleResolver.aromaticRingColor(atomIDs: ring, molecule: depictionMolecule, style: style)
                appendAromaticCircle(points: ringPoints,
                                     color: ringColor,
                                     width: max(1.0, baseBondWidth * 0.84),
                                     opacity: 0.82,
                                     into: &bondSegments)
            }
        }

        struct PendingLabelPlacement {
            let atomID: Int
            let text: String
            let anchor: CGPoint
            let centerOffset: CGVector
            let fontSize: CGFloat
            let aromatic: Bool
            let color: CDKRenderColor
            let neighborPositions: [CGPoint]
            let estimatedSize: CGSize
        }

        typealias BondSegment = (a1: Int, a2: Int, p1: CGPoint, p2: CGPoint)

        func normalizeVector(dx: CGFloat, dy: CGFloat) -> CGVector? {
            let len = hypot(dx, dy)
            guard len > 0.0001 else { return nil }
            return CGVector(dx: dx / len, dy: dy / len)
        }

        func rotate(_ vector: CGVector, by angle: CGFloat) -> CGVector {
            let c = cos(angle)
            let s = sin(angle)
            return CGVector(dx: (vector.dx * c) - (vector.dy * s),
                            dy: (vector.dx * s) + (vector.dy * c))
        }

        func preferredDirection(for item: PendingLabelPlacement) -> CGVector {
            let neighborVectors: [CGVector] = item.neighborPositions.compactMap {
                normalizeVector(dx: $0.x - item.anchor.x, dy: $0.y - item.anchor.y)
            }
            guard !neighborVectors.isEmpty else {
                return CGVector(dx: 0, dy: -1)
            }

            if neighborVectors.count == 1 {
                let v = neighborVectors[0]
                return CGVector(dx: -v.dx, dy: -v.dy)
            }

            let sum = neighborVectors.reduce(CGVector.zero) { partial, v in
                CGVector(dx: partial.dx + v.dx, dy: partial.dy + v.dy)
            }
            if let outward = normalizeVector(dx: -sum.dx, dy: -sum.dy), hypot(sum.dx, sum.dy) > 0.15 {
                return outward
            }

            let angles = neighborVectors.map { atan2($0.dy, $0.dx) }.sorted()
            guard !angles.isEmpty else { return CGVector(dx: 0, dy: -1) }

            var bestGap: CGFloat = -.infinity
            var bestAngle = angles[0]
            for i in angles.indices {
                let a = angles[i]
                let b = (i == angles.count - 1) ? angles[0] + (.pi * 2) : angles[i + 1]
                let gap = b - a
                if gap > bestGap {
                    bestGap = gap
                    bestAngle = a + (gap * 0.5)
                }
            }
            return CGVector(dx: cos(bestAngle), dy: sin(bestAngle))
        }

        func estimateLabelSize(text: String, fontSize: CGFloat) -> CGSize {
            let glyphCount = max(1, text.count)
            let width = max(fontSize * 0.9, (CGFloat(glyphCount) * fontSize * 0.62) + 8)
            let height = max(fontSize * 1.05, fontSize + 4)
            return CGSize(width: width, height: height)
        }

        func makeLabelRect(center: CGPoint, size: CGSize) -> CGRect {
            let paddedWidth = size.width + 8
            let paddedHeight = size.height + 6
            return CGRect(x: center.x - paddedWidth * 0.5,
                          y: center.y - paddedHeight * 0.5,
                          width: paddedWidth,
                          height: paddedHeight)
        }

        func ccw(_ a: CGPoint, _ b: CGPoint, _ c: CGPoint) -> CGFloat {
            ((b.x - a.x) * (c.y - a.y)) - ((b.y - a.y) * (c.x - a.x))
        }

        func onSegment(_ a: CGPoint, _ b: CGPoint, _ p: CGPoint) -> Bool {
            let epsilon: CGFloat = 0.0001
            return p.x >= min(a.x, b.x) - epsilon &&
                   p.x <= max(a.x, b.x) + epsilon &&
                   p.y >= min(a.y, b.y) - epsilon &&
                   p.y <= max(a.y, b.y) + epsilon
        }

        func segmentsIntersect(_ p1: CGPoint, _ p2: CGPoint, _ q1: CGPoint, _ q2: CGPoint) -> Bool {
            let d1 = ccw(p1, p2, q1)
            let d2 = ccw(p1, p2, q2)
            let d3 = ccw(q1, q2, p1)
            let d4 = ccw(q1, q2, p2)
            let epsilon: CGFloat = 0.0001

            if ((d1 > epsilon && d2 < -epsilon) || (d1 < -epsilon && d2 > epsilon)) &&
                ((d3 > epsilon && d4 < -epsilon) || (d3 < -epsilon && d4 > epsilon)) {
                return true
            }

            if abs(d1) <= epsilon && onSegment(p1, p2, q1) { return true }
            if abs(d2) <= epsilon && onSegment(p1, p2, q2) { return true }
            if abs(d3) <= epsilon && onSegment(q1, q2, p1) { return true }
            if abs(d4) <= epsilon && onSegment(q1, q2, p2) { return true }
            return false
        }

        func segmentIntersectsRect(_ a: CGPoint, _ b: CGPoint, _ rect: CGRect) -> Bool {
            let expanded = rect.insetBy(dx: -0.75, dy: -0.75)
            if expanded.contains(a) || expanded.contains(b) { return true }
            let segmentBounds = CGRect(x: min(a.x, b.x),
                                       y: min(a.y, b.y),
                                       width: abs(a.x - b.x),
                                       height: abs(a.y - b.y))
            if !segmentBounds.intersects(expanded) { return false }

            let topLeft = CGPoint(x: expanded.minX, y: expanded.minY)
            let topRight = CGPoint(x: expanded.maxX, y: expanded.minY)
            let bottomLeft = CGPoint(x: expanded.minX, y: expanded.maxY)
            let bottomRight = CGPoint(x: expanded.maxX, y: expanded.maxY)

            return segmentsIntersect(a, b, topLeft, topRight) ||
                   segmentsIntersect(a, b, topRight, bottomRight) ||
                   segmentsIntersect(a, b, bottomRight, bottomLeft) ||
                   segmentsIntersect(a, b, bottomLeft, topLeft)
        }

        let bondSegmentsForAvoidance: [BondSegment] = depictionMolecule.bonds.compactMap { bond in
            guard let p1 = basePositionByAtomID[bond.a1], let p2 = basePositionByAtomID[bond.a2] else { return nil }
            return (a1: bond.a1, a2: bond.a2, p1: p1, p2: p2)
        }

        var pendingLabels: [PendingLabelPlacement] = []
        pendingLabels.reserveCapacity(depictionMolecule.atoms.count)
        let labelLayoutFontSize = max(9, style.fontSize)
        let renderedLabelFontSize = max(9, style.fontSize * zoom)

        for atom in depictionMolecule.atoms {
            let atomDegree = degree[atom.id] ?? 0
            guard CDKLabelText.shouldDrawLabel(atom: atom, degree: atomDegree, style: style),
                  let anchor = basePositionByAtomID[atom.id] else { continue }

            let implicitH = style.showImplicitHydrogens ? depictionMolecule.implicitHydrogenCount(for: atom.id) : 0
            let text = CDKLabelText.build(atom: atom, style: style, implicitHydrogenCount: implicitH)

            let neighborPositions = depictionMolecule.neighbors(of: atom.id).compactMap { basePositionByAtomID[$0] }
            let estimated = estimateLabelSize(text: text, fontSize: labelLayoutFontSize)
            let atomColor = CDKRenderingStyleResolver.atomColor(for: atom, style: style)
            let centerOffset = CDKLabelText.centerOffset(atom: atom,
                                                         style: style,
                                                         implicitHydrogenCount: implicitH,
                                                         fontSize: labelLayoutFontSize)

            pendingLabels.append(PendingLabelPlacement(atomID: atom.id,
                                                       text: text,
                                                       anchor: anchor,
                                                       centerOffset: centerOffset,
                                                       fontSize: labelLayoutFontSize,
                                                       aromatic: atom.aromatic,
                                                       color: atomColor,
                                                       neighborPositions: neighborPositions,
                                                       estimatedSize: estimated))
        }

        var resolvedCentersByAtomID: [Int: CGPoint] = [:]
        var placedRects: [CGRect] = []
        placedRects.reserveCapacity(pendingLabels.count)

        let placementOrder = pendingLabels.sorted { lhs, rhs in
            let lhsHydrogen = lhs.text.hasPrefix("H") ? 1 : 0
            let rhsHydrogen = rhs.text.hasPrefix("H") ? 1 : 0
            if lhsHydrogen != rhsHydrogen { return lhsHydrogen < rhsHydrogen }
            if lhs.text.count != rhs.text.count { return lhs.text.count > rhs.text.count }
            return lhs.atomID < rhs.atomID
        }

        let angleOffsets: [CGFloat] = [0, .pi / 8, -.pi / 8, .pi / 4, -.pi / 4, .pi / 2, -.pi / 2, .pi * 0.75, -.pi * 0.75, .pi]

        for label in placementOrder {
            let nominalCenter = CGPoint(x: label.anchor.x + label.centerOffset.dx,
                                        y: label.anchor.y + label.centerOffset.dy)
            var bestCenter = nominalCenter
            var bestRect = makeLabelRect(center: nominalCenter, size: label.estimatedSize)
            var bestScore = CGFloat.greatestFiniteMagnitude
            let preferred = preferredDirection(for: label)

            var candidates: [CGPoint] = [nominalCenter]
            let baseRadius = max(7, label.fontSize * 0.62) + min(label.estimatedSize.width, label.estimatedSize.height) * 0.18
            for angle in angleOffsets {
                let direction = rotate(preferred, by: angle)
                candidates.append(CGPoint(x: nominalCenter.x + direction.dx * baseRadius,
                                          y: nominalCenter.y + direction.dy * baseRadius))
            }

            for (idx, center) in candidates.enumerated() {
                let candidateRect = makeLabelRect(center: center, size: label.estimatedSize)
                var score = CGFloat(idx) * 3.0

                if idx > 0 {
                    score += hypot(center.x - nominalCenter.x, center.y - nominalCenter.y) * 0.12
                }

                for neighbor in label.neighborPositions where candidateRect.contains(neighbor) {
                    score += 90
                }

                for placed in placedRects where candidateRect.intersects(placed) {
                    let overlap = candidateRect.intersection(placed)
                    score += 1200 + (overlap.width * overlap.height * 20)
                }

                for segment in bondSegmentsForAvoidance {
                    if segment.a1 == label.atomID || segment.a2 == label.atomID {
                        continue
                    }
                    if segmentIntersectsRect(segment.p1, segment.p2, candidateRect) {
                        score += 220
                    }
                }

                if score < bestScore {
                    bestScore = score
                    bestCenter = center
                    bestRect = candidateRect
                }
            }

            resolvedCentersByAtomID[label.atomID] = bestCenter
            placedRects.append(bestRect)
        }

        let labelClipPadding = max(1.4, style.bondWidth * 0.55) * zoom
        let labelClipObstacles: [CDKLabelObstacle] = pendingLabels.compactMap { label in
            guard let center = resolvedCentersByAtomID[label.atomID] else { return nil }
            let renderedCenter = applyViewportTransform(center)
            return CDKLabelClipping.makeGlyphObstacle(text: label.text,
                                                      center: renderedCenter,
                                                      fontSize: renderedLabelFontSize,
                                                      padding: labelClipPadding)
        }
        bondSegments = bondSegments.compactMap { segment in
            guard let (start, end) = CDKLabelClipping.clipSegmentEndpoints(segment.from,
                                                                            segment.to,
                                                                            labelObstacles: labelClipObstacles,
                                                                            padding: labelClipPadding) else {
                return nil
            }
            return CDKMetalDepictionScene.LineSegment(from: start,
                                                      to: end,
                                                      width: segment.width,
                                                      opacity: segment.opacity,
                                                      color: segment.color)
        }

        var labels: [CDKMetalDepictionScene.AtomLabel] = []
        labels.reserveCapacity(pendingLabels.count)
        for label in pendingLabels {
            let resolvedBaseCenter = resolvedCentersByAtomID[label.atomID] ?? label.anchor
            let renderedCenter = applyViewportTransform(resolvedBaseCenter)
            labels.append(CDKMetalDepictionScene.AtomLabel(id: label.atomID,
                                                           text: label.text,
                                                           position: renderedCenter,
                                                           fontSize: renderedLabelFontSize,
                                                           aromatic: label.aromatic,
                                                           color: label.color))
        }

        return CDKMetalDepictionScene(gridSegments: gridSegments, bondSegments: bondSegments, labels: labels)
    }
}

private extension Comparable {
    func clamped(to range: ClosedRange<Self>) -> Self {
        min(max(self, range.lowerBound), range.upperBound)
    }
}
