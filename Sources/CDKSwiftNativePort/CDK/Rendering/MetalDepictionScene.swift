#if canImport(CoreGraphics)
    import CoreGraphics
#endif
import Foundation

public struct CDKMetalDepictionScene: Hashable {
    public struct BackgroundBox: Hashable {
        public let rect: CGRect
        public let cornerRadius: CGFloat
        public let color: CDKRenderColor
        public let opacity: CGFloat

        public init(
            rect: CGRect,
            cornerRadius: CGFloat,
            color: CDKRenderColor,
            opacity: CGFloat = 1.0
        ) {
            self.rect = rect
            self.cornerRadius = cornerRadius
            self.color = color
            self.opacity = opacity
        }
    }

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

    /// A filled outline whose first vertex can anchor a triangle fan.
    public struct FilledPolygon: Hashable {
        public let points: [CGPoint]
        public let color: CDKRenderColor

        public init(points: [CGPoint], color: CDKRenderColor) {
            self.points = points
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
        public let italicized: Bool
        public let drawsBackground: Bool
        public let usesGlowOverlay: Bool
        public let suppressesMatchedBackground: Bool

        public init(
            id: Int,
            text: String,
            position: CGPoint,
            fontSize: CGFloat,
            aromatic: Bool,
            color: CDKRenderColor,
            italicized: Bool = false,
            drawsBackground: Bool = false,
            usesGlowOverlay: Bool = false,
            suppressesMatchedBackground: Bool = false
        ) {
            self.id = id
            self.text = text
            self.position = position
            self.fontSize = fontSize
            self.aromatic = aromatic
            self.color = color
            self.italicized = italicized
            self.drawsBackground = drawsBackground
            self.usesGlowOverlay = usesGlowOverlay
            self.suppressesMatchedBackground = suppressesMatchedBackground
        }
    }

    public let gridSegments: [LineSegment]
    public let backgroundBoxes: [BackgroundBox]
    public let bondSegments: [LineSegment]
    public let labels: [AtomLabel]
    public let filledPolygons: [FilledPolygon]

    public init(
        gridSegments: [LineSegment],
        backgroundBoxes: [BackgroundBox] = [],
        bondSegments: [LineSegment],
        labels: [AtomLabel],
        filledPolygons: [FilledPolygon] = []
    ) {
        self.gridSegments = gridSegments
        self.backgroundBoxes = backgroundBoxes
        self.bondSegments = bondSegments
        self.labels = labels
        self.filledPolygons = filledPolygons
    }
}
