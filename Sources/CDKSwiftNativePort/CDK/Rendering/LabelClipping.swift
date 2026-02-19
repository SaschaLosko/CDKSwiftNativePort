import CoreGraphics
import CoreText
import Foundation

enum CDKLabelText {
    static func shouldDrawLabel(atom: Atom, degree: Int, style: RenderStyle) -> Bool {
        if atom.element == "C" && !style.showCarbons && atom.charge == 0 && !atom.aromatic {
            return degree <= 1
        }
        return true
    }

    static func chargeText(for atom: Atom) -> String {
        guard atom.charge != 0 else { return "" }
        if atom.charge == 1 { return "+" }
        if atom.charge == -1 { return "−" }
        return atom.charge > 0 ? "+\(atom.charge)" : "−\(-atom.charge)"
    }

    static func hydrogenText(atom: Atom, style: RenderStyle, implicitHydrogenCount: Int) -> String {
        let showH = style.showImplicitHydrogens && (atom.element != "C" || style.showCarbons)
        guard showH, implicitHydrogenCount > 0 else { return "" }
        return implicitHydrogenCount == 1 ? "H" : "H\(implicitHydrogenCount)"
    }

    static func build(atom: Atom, style: RenderStyle, implicitHydrogenCount: Int) -> String {
        let chargeText = chargeText(for: atom)
        let hText = hydrogenText(atom: atom, style: style, implicitHydrogenCount: implicitHydrogenCount)

        var label = atom.symbolToDraw + hText + chargeText
        if style.showAtomIDs {
            label += " \(atom.id)"
        }
        return label
    }

    static func centerOffset(atom: Atom,
                             style: RenderStyle,
                             implicitHydrogenCount: Int,
                             fontSize: CGFloat) -> CGVector {
        // Keep the atom symbol centered on the atomic anchor even when
        // suffix text (e.g. "H" in "OH") is displayed.
        guard !style.showAtomIDs else { return .zero }
        let hText = hydrogenText(atom: atom, style: style, implicitHydrogenCount: implicitHydrogenCount)
        guard !hText.isEmpty else { return .zero }

        let fullLabel = build(atom: atom, style: style, implicitHydrogenCount: implicitHydrogenCount)
        let coreLabel = atom.symbolToDraw

        let fullWidth = CDKLabelClipping.estimateLabelSize(text: fullLabel, fontSize: fontSize).width
        let coreWidth = CDKLabelClipping.estimateLabelSize(text: coreLabel, fontSize: fontSize).width
        let dx = max(0, (fullWidth - coreWidth) * 0.5)
        return CGVector(dx: dx, dy: 0)
    }
}

struct CDKLabelObstacle {
    let fillPath: CGPath
    let clipPath: CGPath
    let bounds: CGRect

    init(fillPath: CGPath, clipPath: CGPath) {
        self.fillPath = fillPath
        self.clipPath = clipPath
        self.bounds = fillPath.boundingBoxOfPath.union(clipPath.boundingBoxOfPath)
    }

    init(rect: CGRect) {
        let path = CGPath(rect: rect, transform: nil)
        self.init(fillPath: path, clipPath: path)
    }

    func contains(_ point: CGPoint) -> Bool {
        clipPath.contains(point, using: .winding, transform: .identity)
            || fillPath.contains(point, using: .winding, transform: .identity)
    }
}

enum CDKLabelClipping {
    static func estimateLabelSize(text: String, fontSize: CGFloat) -> CGSize {
        let glyphCount = max(1, text.count)
        let width = max(fontSize * 0.9, (CGFloat(glyphCount) * fontSize * 0.62) + 8)
        let height = max(fontSize * 1.05, fontSize + 4)
        return CGSize(width: width, height: height)
    }

    static func makeLabelRect(center: CGPoint, estimatedTextSize: CGSize) -> CGRect {
        let width = estimatedTextSize.width + 8
        let height = estimatedTextSize.height + 6
        return CGRect(x: center.x - width * 0.5,
                      y: center.y - height * 0.5,
                      width: width,
                      height: height)
    }

    static func makeGlyphObstacle(text: String,
                                  center: CGPoint,
                                  fontSize: CGFloat,
                                  padding: CGFloat) -> CDKLabelObstacle {
        guard let glyphPath = makeGlyphPath(text: text, center: center, fontSize: fontSize) else {
            let estimate = estimateLabelSize(text: text, fontSize: fontSize)
            let rect = makeLabelRect(center: center, estimatedTextSize: estimate)
                .insetBy(dx: -padding, dy: -padding)
            return CDKLabelObstacle(rect: rect)
        }

        let carrierInset = max(1.0, padding * 0.35)
        let carrierRect = glyphPath.boundingBoxOfPath.insetBy(dx: -carrierInset, dy: -carrierInset)
        let corner = max(2.0, fontSize * 0.16)
        let carrierPath = CGPath(roundedRect: carrierRect,
                                 cornerWidth: corner,
                                 cornerHeight: corner,
                                 transform: nil)

        if padding <= 0.0001 {
            return CDKLabelObstacle(fillPath: carrierPath, clipPath: glyphPath)
        }

        let expanded = glyphPath.copy(strokingWithWidth: max(0.2, padding * 2),
                                      lineCap: .round,
                                      lineJoin: .round,
                                      miterLimit: 2.0)
        return CDKLabelObstacle(fillPath: carrierPath, clipPath: expanded)
    }

    static func makeGlyphObstacles(_ labels: [(text: String, center: CGPoint, fontSize: CGFloat)],
                                   padding: CGFloat) -> [CDKLabelObstacle] {
        labels.map { makeGlyphObstacle(text: $0.text, center: $0.center, fontSize: $0.fontSize, padding: padding) }
    }

    static func clipSegmentEndpoints(_ start: CGPoint,
                                     _ end: CGPoint,
                                     labelObstacles: [CDKLabelObstacle],
                                     padding: CGFloat) -> (CGPoint, CGPoint)? {
        let dx = end.x - start.x
        let dy = end.y - start.y
        let len = hypot(dx, dy)
        guard len > 0.0001 else { return nil }

        let tPadding = min(0.12, max(0, padding / len))
        var t0: CGFloat = 0
        var t1: CGFloat = 1

        @inline(__always)
        func point(at t: CGFloat) -> CGPoint {
            CGPoint(x: start.x + dx * t, y: start.y + dy * t)
        }

        for obstacle in labelObstacles {
            let currentStart = point(at: t0)
            let currentEnd = point(at: t1)
            let currentBounds = CGRect(x: min(currentStart.x, currentEnd.x),
                                       y: min(currentStart.y, currentEnd.y),
                                       width: abs(currentEnd.x - currentStart.x),
                                       height: abs(currentEnd.y - currentStart.y))
            if !currentBounds.intersects(obstacle.bounds.insetBy(dx: -0.5, dy: -0.5)) {
                continue
            }

            let startInside = obstacle.contains(currentStart)
            let endInside = obstacle.contains(currentEnd)
            if !startInside && !endInside {
                continue
            }

            if startInside {
                guard let exitT = findBoundaryT(from: t0,
                                                to: t1,
                                                startInside: true,
                                                pointAt: point,
                                                contains: obstacle.contains) else {
                    return nil
                }
                t0 = min(t1, exitT + tPadding)
            }

            if endInside {
                guard let enterT = findBoundaryT(from: t1,
                                                 to: t0,
                                                 startInside: true,
                                                 pointAt: point,
                                                 contains: obstacle.contains) else {
                    return nil
                }
                t1 = max(t0, enterT - tPadding)
            }

            if t1 - t0 <= 0.0005 {
                return nil
            }
        }

        return (point(at: t0), point(at: t1))
    }

    // Backward-compatible rect clipping entry used by tests and fallback callers.
    static func clipSegmentEndpoints(_ start: CGPoint,
                                     _ end: CGPoint,
                                     labelRects: [CGRect],
                                     padding: CGFloat) -> (CGPoint, CGPoint)? {
        let obstacles = labelRects.map { rect in
            CDKLabelObstacle(rect: rect.insetBy(dx: -padding * 0.35, dy: -padding * 0.35))
        }
        return clipSegmentEndpoints(start, end, labelObstacles: obstacles, padding: padding)
    }

    private static func findBoundaryT(from startT: CGFloat,
                                      to endT: CGFloat,
                                      startInside: Bool,
                                      pointAt: (CGFloat) -> CGPoint,
                                      contains: (CGPoint) -> Bool) -> CGFloat? {
        let samples = 84
        let step = (endT - startT) / CGFloat(samples)
        var prevT = startT
        var prevInside = startInside

        for i in 1...samples {
            let t = startT + CGFloat(i) * step
            let inside = contains(pointAt(t))
            if inside != prevInside {
                var low = prevT
                var high = t
                for _ in 0..<16 {
                    let mid = (low + high) * 0.5
                    let midInside = contains(pointAt(mid))
                    if midInside == prevInside {
                        low = mid
                    } else {
                        high = mid
                    }
                }
                return (low + high) * 0.5
            }
            prevT = t
            prevInside = inside
        }

        return nil
    }

    private static func makeGlyphPath(text: String, center: CGPoint, fontSize: CGFloat) -> CGPath? {
        guard !text.isEmpty else { return nil }

        let fontName = "Helvetica-Bold" as CFString
        let ctFont = CTFontCreateWithName(fontName, max(6, fontSize), nil)
        let attributes: [NSAttributedString.Key: Any] = [
            NSAttributedString.Key(rawValue: kCTFontAttributeName as String): ctFont
        ]
        let attributed = NSAttributedString(string: text, attributes: attributes)
        let line = CTLineCreateWithAttributedString(attributed)
        let lineBounds = CTLineGetBoundsWithOptions(line, [.useGlyphPathBounds, .useOpticalBounds])
        if lineBounds.isNull || lineBounds.isEmpty || !lineBounds.width.isFinite || !lineBounds.height.isFinite {
            return nil
        }

        let lineOrigin = CGPoint(x: center.x - lineBounds.midX, y: center.y - lineBounds.midY)
        let runs = CTLineGetGlyphRuns(line)
        let runCount = CFArrayGetCount(runs)
        let path = CGMutablePath()

        for runIndex in 0..<runCount {
            let runValue = CFArrayGetValueAtIndex(runs, runIndex)
            let run = unsafeBitCast(runValue, to: CTRun.self)
            let runAttributes = CTRunGetAttributes(run) as NSDictionary
            let runFont = (runAttributes[kCTFontAttributeName] as! CTFont?) ?? ctFont
            let glyphCount = CTRunGetGlyphCount(run)
            if glyphCount <= 0 { continue }

            var glyphs = Array(repeating: CGGlyph(), count: glyphCount)
            var positions = Array(repeating: CGPoint.zero, count: glyphCount)
            CTRunGetGlyphs(run, CFRangeMake(0, 0), &glyphs)
            CTRunGetPositions(run, CFRangeMake(0, 0), &positions)

            for i in 0..<glyphCount {
                guard let glyphPath = CTFontCreatePathForGlyph(runFont, glyphs[i], nil) else { continue }
                let p = positions[i]
                let transform = CGAffineTransform(translationX: lineOrigin.x + p.x, y: lineOrigin.y + p.y)
                path.addPath(glyphPath, transform: transform)
            }
        }

        guard !path.isEmpty else { return nil }
        return path.copy()
    }
}
