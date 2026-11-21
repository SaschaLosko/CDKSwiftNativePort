#if canImport(CoreText)
    import CoreGraphics
    import CoreText
    import Foundation

    /// The shaped label and its ink bounds in a baseline-relative, y-up coordinate system.
    public struct CDKDepictionLabelLayout {
        private let line: CTLine
        public let inkBounds: CGRect

        public init(text: String, fontSize: CGFloat, italicized: Bool = false) {
            let size = fontSize.isFinite ? max(0.5, fontSize) : 12
            let regular =
                CTFontCreateUIFontForLanguage(.system, size, nil)
                ?? CTFontCreateWithName("Helvetica" as CFString, size, nil)
            let font: CTFont
            if italicized {
                let serif = CTFontCreateWithName("TimesNewRomanPS-ItalicMT" as CFString, size, nil)
                if CTFontCopyPostScriptName(serif) as String == "TimesNewRomanPS-ItalicMT" {
                    font = serif
                } else {
                    font =
                        CTFontCreateCopyWithSymbolicTraits(regular, size, nil, .traitItalic, .traitItalic)
                        ?? regular
                }
            } else {
                font = regular
            }
            let attributes: [NSAttributedString.Key: Any] = [
                NSAttributedString.Key(kCTFontAttributeName as String): font,
                NSAttributedString.Key(kCTForegroundColorFromContextAttributeName as String): true,
            ]
            let shaped = CTLineCreateWithAttributedString(NSAttributedString(string: text, attributes: attributes))
            line = shaped
            let bounds = CTLineGetBoundsWithOptions(shaped, .useGlyphPathBounds)
            if !bounds.isNull && !bounds.isEmpty && bounds.width.isFinite && bounds.height.isFinite {
                inkBounds = bounds
            } else {
                var ascent: CGFloat = 0
                var descent: CGFloat = 0
                let width = CTLineGetTypographicBounds(shaped, &ascent, &descent, nil)
                inkBounds = CGRect(x: 0, y: -descent, width: width, height: ascent + descent)
            }
        }

        /// Draws the line around its ink center using the context fill color.
        public func draw(in context: CGContext, centeredAt center: CGPoint) {
            context.textMatrix = .identity
            context.textPosition = CGPoint(x: center.x - inkBounds.midX, y: center.y - inkBounds.midY)
            CTLineDraw(line, context)
        }

        /// Returns the font outlines around the ink center for clipping or vector export.
        public func outline(centeredAt center: CGPoint, yDown: Bool) -> CGPath? {
            let runs = CTLineGetGlyphRuns(line) as! [CTRun]
            let path = CGMutablePath()
            let yScale: CGFloat = yDown ? -1 : 1
            for run in runs {
                let attributes = CTRunGetAttributes(run) as NSDictionary
                guard let value = attributes[kCTFontAttributeName] else { continue }
                let font = value as! CTFont
                let count = CTRunGetGlyphCount(run)
                var glyphs = Array(repeating: CGGlyph(), count: count)
                var positions = Array(repeating: CGPoint.zero, count: count)
                CTRunGetGlyphs(run, CFRangeMake(0, 0), &glyphs)
                CTRunGetPositions(run, CFRangeMake(0, 0), &positions)
                for (glyph, position) in zip(glyphs, positions) {
                    guard let shape = CTFontCreatePathForGlyph(font, glyph, nil) else { continue }
                    path.addPath(
                        shape,
                        transform: CGAffineTransform(
                            a: 1, b: 0, c: 0, d: yScale,
                            tx: center.x + position.x - inkBounds.midX,
                            ty: center.y + yScale * (position.y - inkBounds.midY)))
                }
            }
            return path.isEmpty ? nil : path.copy()
        }

        func horizontalOffset(atUTF16Index index: Int) -> CGFloat {
            CTLineGetOffsetForStringIndex(line, index, nil)
        }
    }
#endif
