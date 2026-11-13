#if canImport(CoreGraphics)
    import CoreGraphics
#endif
import Foundation

enum MetalDepictionMetrics {
    static func representativeBondLength(
        molecule: Molecule,
        positionsByAtomID: [Int: CGPoint]
    ) -> CGFloat? {
        let bondLengths = molecule.bonds.compactMap { bond -> CGFloat? in
            guard let p1 = positionsByAtomID[bond.a1],
                let p2 = positionsByAtomID[bond.a2]
            else {
                return nil
            }
            let length = hypot(p2.x - p1.x, p2.y - p1.y)
            return (length.isFinite && length > 0.0001) ? length : nil
        }.sorted()

        guard !bondLengths.isEmpty else { return nil }
        let mid = bondLengths.count / 2
        if bondLengths.count.isMultiple(of: 2) {
            return (bondLengths[mid - 1] + bondLengths[mid]) * 0.5
        }
        return bondLengths[mid]
    }

    private static func labelDensityScale(visibleAtomLabelCount: Int) -> CGFloat {
        let extraLabels = max(0, visibleAtomLabelCount - 2)
        return (1.0 - (CGFloat(extraLabels) * 0.026)).clamped(to: 0.80...1.0)
    }

    private static func sceneProminenceBoost(
        showCarbons: Bool,
        atomCount: Int,
        bondCount: Int,
        visibleAtomLabelCount: Int
    ) -> CGFloat {
        let extraVisibleLabels = max(0, visibleAtomLabelCount - 1)
        let hiddenCarbonBias =
            showCarbons
            ? CGFloat.zero
            : max(0, 0.12 - (CGFloat(max(0, visibleAtomLabelCount - 4)) * 0.03))
        let macroPenalty =
            (CGFloat(max(0, atomCount - 40)) * 0.004)
            + (CGFloat(max(0, bondCount - 42)) * 0.003)
        return
            (1.60
            + hiddenCarbonBias
            - (CGFloat(extraVisibleLabels) * 0.07)
            - macroPenalty)
            .clamped(to: 1.0...1.72)
    }

    static func zoomVisualScale(
        _ zoom: CGFloat,
        zoomOutExponent: CGFloat,
        zoomInExponent: CGFloat
    ) -> CGFloat {
        let safeZoom = max(0.2, zoom)
        if safeZoom >= 1.0 {
            return pow(safeZoom, zoomInExponent)
        }
        return pow(safeZoom, zoomOutExponent)
    }

    static func autoFitLabelFontSize(
        style: RenderStyle,
        baseViewportFontSize: CGFloat,
        representativeBondLength: CGFloat?,
        atomCount: Int,
        bondCount: Int,
        visibleAtomLabelCount: Int,
        minimumLabelFontSize: CGFloat
    ) -> CGFloat {
        let clampedMinimum = max(2.5, minimumLabelFontSize)
        guard let representativeBondLength,
            representativeBondLength.isFinite,
            representativeBondLength > 0.0001
        else {
            return max(clampedMinimum, baseViewportFontSize)
        }
        if baseViewportFontSize <= clampedMinimum + 0.0001 {
            return clampedMinimum
        }

        let densityScale = labelDensityScale(visibleAtomLabelCount: visibleAtomLabelCount)
        let prominenceBoost = sceneProminenceBoost(
            showCarbons: style.showCarbons,
            atomCount: atomCount,
            bondCount: bondCount,
            visibleAtomLabelCount: visibleAtomLabelCount)
        let preferredMinimumRatio = ((style.fontSize / 152.0) * densityScale).clamped(to: 0.126...0.166)
        let ratioFloorFontSize = representativeBondLength * preferredMinimumRatio
        let prominenceFontSize = baseViewportFontSize * prominenceBoost
        let boostedFontSize = max(baseViewportFontSize, max(ratioFloorFontSize, prominenceFontSize))
        let absoluteMaxFontSize = max(
            clampedMinimum,
            baseViewportFontSize * (prominenceBoost + 0.10).clamped(to: 1.16...1.50))
        return boostedFontSize.clamped(to: clampedMinimum...absoluteMaxFontSize)
    }

    static func autoFitBondWidth(
        style: RenderStyle,
        baseViewportBondWidth: CGFloat,
        representativeBondLength: CGFloat?,
        atomCount: Int,
        bondCount: Int,
        visibleAtomLabelCount: Int
    ) -> CGFloat {
        guard let representativeBondLength,
            representativeBondLength.isFinite,
            representativeBondLength > 0.0001
        else {
            return max(1.0, baseViewportBondWidth)
        }

        let densityScale = labelDensityScale(visibleAtomLabelCount: visibleAtomLabelCount)
        let prominenceBoost = sceneProminenceBoost(
            showCarbons: style.showCarbons,
            atomCount: atomCount,
            bondCount: bondCount,
            visibleAtomLabelCount: visibleAtomLabelCount)
        let strokeProminenceBoost = 1.0 + ((prominenceBoost - 1.0) * 0.72)
        let preferredMinimumRatio = ((style.bondWidth / 72.0) * densityScale).clamped(to: 0.016...0.026)
        let ratioFloorWidth = representativeBondLength * preferredMinimumRatio
        let prominenceWidth = baseViewportBondWidth * strokeProminenceBoost
        let boostedBondWidth = max(baseViewportBondWidth, max(ratioFloorWidth, prominenceWidth))
        let absoluteMaxBondWidth = max(
            1.0,
            baseViewportBondWidth * (strokeProminenceBoost + 0.06).clamped(to: 1.10...1.34))
        return boostedBondWidth.clamped(to: 1.0...absoluteMaxBondWidth)
    }

}

private extension Comparable {
    func clamped(to range: ClosedRange<Self>) -> Self {
        min(max(self, range.lowerBound), range.upperBound)
    }
}
