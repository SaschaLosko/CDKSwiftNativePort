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

        public init(rect: CGRect,
                    cornerRadius: CGFloat,
                    color: CDKRenderColor,
                    opacity: CGFloat = 1.0) {
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

        public init(id: Int,
                    text: String,
                    position: CGPoint,
                    fontSize: CGFloat,
                    aromatic: Bool,
                    color: CDKRenderColor,
                    italicized: Bool = false,
                    drawsBackground: Bool = true,
                    usesGlowOverlay: Bool = false,
                    suppressesMatchedBackground: Bool = false) {
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

    public init(gridSegments: [LineSegment],
                backgroundBoxes: [BackgroundBox] = [],
                bondSegments: [LineSegment],
                labels: [AtomLabel]) {
        self.gridSegments = gridSegments
        self.backgroundBoxes = backgroundBoxes
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
    private static func representativeBondLength(molecule: Molecule,
                                                 positionsByAtomID: [Int: CGPoint]) -> CGFloat? {
        let bondLengths = molecule.bonds.compactMap { bond -> CGFloat? in
            guard let p1 = positionsByAtomID[bond.a1],
                  let p2 = positionsByAtomID[bond.a2] else {
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

    private static func sceneProminenceBoost(showCarbons: Bool,
                                             atomCount: Int,
                                             bondCount: Int,
                                             visibleAtomLabelCount: Int) -> CGFloat {
        let extraVisibleLabels = max(0, visibleAtomLabelCount - 1)
        let hiddenCarbonBias = showCarbons
            ? CGFloat.zero
            : max(0, 0.12 - (CGFloat(max(0, visibleAtomLabelCount - 4)) * 0.03))
        let macroPenalty = (CGFloat(max(0, atomCount - 40)) * 0.004)
            + (CGFloat(max(0, bondCount - 42)) * 0.003)
        return (1.60
            + hiddenCarbonBias
            - (CGFloat(extraVisibleLabels) * 0.07)
            - macroPenalty)
            .clamped(to: 1.0...1.72)
    }

    private static func zoomVisualScale(_ zoom: CGFloat,
                                        zoomOutExponent: CGFloat,
                                        zoomInExponent: CGFloat) -> CGFloat {
        let safeZoom = max(0.2, zoom)
        if safeZoom >= 1.0 {
            return pow(safeZoom, zoomInExponent)
        }
        return pow(safeZoom, zoomOutExponent)
    }

    private static func autoFitLabelFontSize(style: RenderStyle,
                                             baseViewportFontSize: CGFloat,
                                             representativeBondLength: CGFloat?,
                                             atomCount: Int,
                                             bondCount: Int,
                                             visibleAtomLabelCount: Int,
                                             minimumLabelFontSize: CGFloat) -> CGFloat {
        let clampedMinimum = max(2.5, minimumLabelFontSize)
        guard let representativeBondLength,
              representativeBondLength.isFinite,
              representativeBondLength > 0.0001 else {
            return max(clampedMinimum, baseViewportFontSize)
        }
        if baseViewportFontSize <= clampedMinimum + 0.0001 {
            return clampedMinimum
        }

        let densityScale = labelDensityScale(visibleAtomLabelCount: visibleAtomLabelCount)
        let prominenceBoost = sceneProminenceBoost(showCarbons: style.showCarbons,
                                                   atomCount: atomCount,
                                                   bondCount: bondCount,
                                                   visibleAtomLabelCount: visibleAtomLabelCount)
        let preferredMinimumRatio = ((style.fontSize / 152.0) * densityScale).clamped(to: 0.126...0.166)
        let ratioFloorFontSize = representativeBondLength * preferredMinimumRatio
        let prominenceFontSize = baseViewportFontSize * prominenceBoost
        let boostedFontSize = max(baseViewportFontSize, max(ratioFloorFontSize, prominenceFontSize))
        let absoluteMaxFontSize = max(clampedMinimum,
                                      baseViewportFontSize * (prominenceBoost + 0.10).clamped(to: 1.16...1.50))
        return boostedFontSize.clamped(to: clampedMinimum...absoluteMaxFontSize)
    }

    private static func autoFitBondWidth(style: RenderStyle,
                                         baseViewportBondWidth: CGFloat,
                                         representativeBondLength: CGFloat?,
                                         atomCount: Int,
                                         bondCount: Int,
                                         visibleAtomLabelCount: Int) -> CGFloat {
        guard let representativeBondLength,
              representativeBondLength.isFinite,
              representativeBondLength > 0.0001 else {
            return max(1.0, baseViewportBondWidth)
        }

        let densityScale = labelDensityScale(visibleAtomLabelCount: visibleAtomLabelCount)
        let prominenceBoost = sceneProminenceBoost(showCarbons: style.showCarbons,
                                                   atomCount: atomCount,
                                                   bondCount: bondCount,
                                                   visibleAtomLabelCount: visibleAtomLabelCount)
        let strokeProminenceBoost = 1.0 + ((prominenceBoost - 1.0) * 0.72)
        let preferredMinimumRatio = ((style.bondWidth / 72.0) * densityScale).clamped(to: 0.016...0.026)
        let ratioFloorWidth = representativeBondLength * preferredMinimumRatio
        let prominenceWidth = baseViewportBondWidth * strokeProminenceBoost
        let boostedBondWidth = max(baseViewportBondWidth, max(ratioFloorWidth, prominenceWidth))
        let absoluteMaxBondWidth = max(1.0,
                                       baseViewportBondWidth * (strokeProminenceBoost + 0.06).clamped(to: 1.10...1.34))
        return boostedBondWidth.clamped(to: 1.0...absoluteMaxBondWidth)
    }

    private static func usesGlowHighlights(style: RenderStyle) -> Bool {
        switch style.highlightStyle {
        case .outerGlow, .outerGlowWhiteEdge:
            return true
        case .none, .colored:
            return false
        }
    }

    private static func preservesBackgroundDuringGlow(style: RenderStyle) -> Bool {
        style.highlightStyle == .outerGlowWhiteEdge
    }

    public static func build(molecule: Molecule,
                             style: RenderStyle,
                             canvasRect: CGRect,
                             zoom: CGFloat,
                             pan: CGSize,
                             rotationDegrees: CGFloat = 0,
                             minimumLabelFontSize: CGFloat = 8.0,
                             includeAromaticCarbonLabelsWhenCarbonsHidden: Bool = true,
                             includeTerminalCarbonLabelsWhenCarbonsHidden: Bool = true) -> CDKMetalDepictionScene {
        let prepared = MetalPreparedDepictionCache.preparedData(molecule: molecule, style: style)
        let depictionMolecule = prepared.depictionMolecule
        let highlightedAtomIDs = style.highlightStyle == .none ? Set<Int>() : Set(depictionMolecule.highlightedAtomIDs)
        let highlightedBondIDs = style.highlightStyle == .none ? Set<Int>() : Set(depictionMolecule.highlightedBondIDs)
        guard let box = prepared.boundingBox, canvasRect.width > 1, canvasRect.height > 1 else {
            return CDKMetalDepictionScene(gridSegments: [], backgroundBoxes: [], bondSegments: [], labels: [])
        }

        let gridSegments: [CDKMetalDepictionScene.LineSegment] = []

        let pad = style.padding
        let available = CGRect(x: canvasRect.minX + pad,
                               y: canvasRect.minY + pad,
                               width: max(1, canvasRect.width - 2 * pad),
                               height: max(1, canvasRect.height - 2 * pad))

        let scaleX = available.width / max(0.0001, box.width)
        let scaleY = available.height / max(0.0001, box.height)
        let scale = min(scaleX, scaleY)
        // Keep auto-fit viewport changes visually consistent with direct zoom:
        // line widths, stereo wedges, and labels should all scale together.
        let referenceViewport = CGSize(width: 960, height: 720)
        let referenceAvailable = CGRect(x: style.padding,
                                        y: style.padding,
                                        width: max(1, referenceViewport.width - (2 * style.padding)),
                                        height: max(1, referenceViewport.height - (2 * style.padding)))
        let referenceScaleX = referenceAvailable.width / max(0.0001, box.width)
        let referenceScaleY = referenceAvailable.height / max(0.0001, box.height)
        let referenceFitScale = min(referenceScaleX, referenceScaleY)
        // Do not clamp the lower bound: when the viewport gets narrow, a floor here
        // makes labels/bonds stop shrinking while geometry keeps shrinking, which
        // visually reads as labels getting bigger again.
        let viewportAutoFitScale = min(1.90, scale / max(0.0001, referenceFitScale))
        let viewportStrokeScale = viewportAutoFitScale

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

        func rotateAround(_ point: CGPoint, center: CGPoint) -> CGPoint {
            let rx = point.x - center.x
            let ry = point.y - center.y
            return CGPoint(x: (rx * rotationCos) - (ry * rotationSin) + center.x,
                           y: (rx * rotationSin) + (ry * rotationCos) + center.y)
        }

        func applyZoomAndPan(_ point: CGPoint) -> CGPoint {
            let zx = (point.x - viewportCenter.x) * zoom + viewportCenter.x + pan.width
            let zy = (point.y - viewportCenter.y) * zoom + viewportCenter.y + pan.height
            return CGPoint(x: zx, y: zy)
        }

        func applyViewportTransform(_ point: CGPoint) -> CGPoint {
            applyZoomAndPan(rotateAround(point, center: viewportCenter))
        }

        var basePositionByAtomID: [Int: CGPoint] = [:]
        for atom in depictionMolecule.atoms {
            let basePosition = atom.position.applying(baseTransform)
            basePositionByAtomID[atom.id] = basePosition
        }

        let visibleAtomLabelCount = depictionMolecule.atoms.reduce(into: 0) { partial, atom in
            let atomDegree = prepared.degreeByAtomID[atom.id] ?? 0
            if CDKLabelText.shouldDrawLabel(atom: atom,
                                            degree: atomDegree,
                                            style: style,
                                            includeAromaticCarbonLabelsWhenCarbonsHidden: includeAromaticCarbonLabelsWhenCarbonsHidden,
                                            includeTerminalCarbonLabelsWhenCarbonsHidden: includeTerminalCarbonLabelsWhenCarbonsHidden,
                                            molecule: depictionMolecule,
                                            highlightedAtomIDs: highlightedAtomIDs,
                                            highlightedBondIDs: highlightedBondIDs) {
                partial += 1
            }
        }
        let representativeBaseBondLength = representativeBondLength(molecule: depictionMolecule,
                                                                    positionsByAtomID: basePositionByAtomID)
        let clampedMinimumLabelFontSize = max(2.5, minimumLabelFontSize)
        let baseViewportLabelFontSize = max(clampedMinimumLabelFontSize, style.fontSize * viewportAutoFitScale)
        let labelLayoutFontSize = autoFitLabelFontSize(style: style,
                                                       baseViewportFontSize: baseViewportLabelFontSize,
                                                       representativeBondLength: representativeBaseBondLength,
                                                       atomCount: depictionMolecule.atomCount,
                                                       bondCount: depictionMolecule.bondCount,
                                                       visibleAtomLabelCount: visibleAtomLabelCount,
                                                       minimumLabelFontSize: clampedMinimumLabelFontSize)
        // When the auto-fit view zooms out, keep labels and strokes from shrinking
        // as aggressively as the scaffold so sparse scenes stay readable.
        let zoomedLabelScale = zoomVisualScale(zoom,
                                               zoomOutExponent: 0.62,
                                               zoomInExponent: 1.18)
        let zoomedStrokeScale = zoomVisualScale(zoom,
                                                zoomOutExponent: 0.68,
                                                zoomInExponent: 1.14)
        let renderedLabelFontSize = max(clampedMinimumLabelFontSize, labelLayoutFontSize * zoomedLabelScale)

        func pointBounds(_ points: some Sequence<CGPoint>) -> CGRect? {
            var iterator = points.makeIterator()
            guard let first = iterator.next() else { return nil }
            var bounds = CGRect(x: first.x, y: first.y, width: 0, height: 0)
            while let point = iterator.next() {
                bounds = bounds.union(CGRect(x: point.x, y: point.y, width: 0, height: 0))
            }
            return bounds
        }

        let legendAtomIDs = Set(depictionMolecule.atoms.compactMap { atom in
            atom.rGroupMembership == nil ? nil : atom.id
        })
        let rootAtomIDs = Set(depictionMolecule.atoms.map(\.id)).subtracting(legendAtomIDs)

        var legendComponentCentersByAtomID: [Int: CGPoint] = [:]
        if !legendAtomIDs.isEmpty {
            var visited: Set<Int> = []
            let legendNeighbors = Dictionary(uniqueKeysWithValues: legendAtomIDs.map { atomID in
                (atomID, depictionMolecule.neighbors(of: atomID).filter { legendAtomIDs.contains($0) })
            })

            for seed in legendAtomIDs.sorted() where !visited.contains(seed) {
                var component: Set<Int> = [seed]
                var stack: [Int] = [seed]
                visited.insert(seed)

                while let atomID = stack.popLast() {
                    for neighbor in legendNeighbors[atomID] ?? [] where !visited.contains(neighbor) {
                        visited.insert(neighbor)
                        component.insert(neighbor)
                        stack.append(neighbor)
                    }
                }

                guard let componentBounds = pointBounds(component.compactMap { basePositionByAtomID[$0] }) else {
                    continue
                }

                let center = CGPoint(x: componentBounds.midX, y: componentBounds.midY)
                for atomID in component {
                    legendComponentCentersByAtomID[atomID] = center
                }
            }
        }

        let baseRootBounds = pointBounds(rootAtomIDs.compactMap { basePositionByAtomID[$0] })
        let baseLegendBounds = pointBounds(legendAtomIDs.compactMap { basePositionByAtomID[$0] })
        let legendGap = max(labelLayoutFontSize * 1.15,
                            (baseLegendBounds?.minY ?? 0) - (baseRootBounds?.maxY ?? 0))

        let rotatedRootBounds = pointBounds(rootAtomIDs.compactMap { atomID in
            basePositionByAtomID[atomID].map { rotateAround($0, center: viewportCenter) }
        })

        var legendTranslation = CGSize.zero
        if !legendAtomIDs.isEmpty {
            let rotatedLegendBounds = pointBounds(legendAtomIDs.compactMap { atomID in
                guard let basePosition = basePositionByAtomID[atomID] else { return nil }
                let center = legendComponentCentersByAtomID[atomID] ?? basePosition
                return rotateAround(basePosition, center: center)
            })

            if let rotatedLegendBounds,
               let rotatedRootBounds {
                legendTranslation = CGSize(width: rotatedRootBounds.midX - rotatedLegendBounds.midX,
                                           height: (rotatedRootBounds.maxY + legendGap) - rotatedLegendBounds.minY)
            }
        }

        func renderedPosition(for atomID: Int, basePosition: CGPoint) -> CGPoint {
            if legendAtomIDs.contains(atomID) {
                let center = legendComponentCentersByAtomID[atomID] ?? basePosition
                let rotated = rotateAround(basePosition, center: center)
                return applyZoomAndPan(rotated.offsetBy(dx: legendTranslation.width, dy: legendTranslation.height))
            }
            return applyViewportTransform(basePosition)
        }

        func renderedPoint(for rawPoint: CGPoint) -> CGPoint {
            applyViewportTransform(rawPoint.applying(baseTransform))
        }

        var positionByAtomID: [Int: CGPoint] = [:]
        for atom in depictionMolecule.atoms {
            guard let basePosition = basePositionByAtomID[atom.id] else { continue }
            positionByAtomID[atom.id] = renderedPosition(for: atom.id, basePosition: basePosition)
        }
        let atomByID = Dictionary(uniqueKeysWithValues: depictionMolecule.atoms.map { ($0.id, $0) })

        let degree = prepared.degreeByAtomID
        let aromaticRings = prepared.aromaticRings
        let aromaticBondIDs = prepared.aromaticBondIDs
        var aromaticEdgeCenters: [MetalRenderEdgeKey: [CGPoint]] = [:]
        var conjugatedDoubleEdgeCenters: [MetalRenderEdgeKey: [CGPoint]] = [:]
        let baseViewportBondWidth = max(1.0, style.bondWidth * viewportStrokeScale)
        let layoutBondWidth = autoFitBondWidth(style: style,
                                               baseViewportBondWidth: baseViewportBondWidth,
                                               representativeBondLength: representativeBaseBondLength,
                                               atomCount: depictionMolecule.atomCount,
                                               bondCount: depictionMolecule.bondCount,
                                               visibleAtomLabelCount: visibleAtomLabelCount)
        let baseBondWidth = max(1.0, layoutBondWidth * zoomedStrokeScale)
        let stereoScale = style.stereoAttenuation.clamped(to: 0.4...1.2)
        let hashedStereoScale = style.hashedStereoAttenuation.clamped(to: 0.5...1.05)
        // Match CDK's default standard generator more closely: the wide end of stereo wedges
        // is derived directly from the bond stroke width rather than a label-relative heuristic.
        let cdkDefaultStereoHalfWidth = baseBondWidth * 3.0
        let stereoSolidWedgeHalfWidth = cdkDefaultStereoHalfWidth * stereoScale
        let stereoHashedWedgeHalfWidth = cdkDefaultStereoHalfWidth * stereoScale * hashedStereoScale
        let doubleBondSeparation = max(3.4, style.bondWidth * 2.25 * viewportStrokeScale) * zoomedStrokeScale
        let doubleBondHalfSeparation = doubleBondSeparation * 0.5
        let tripleBondOffset = max(3.1, style.bondWidth * 2.2 * viewportStrokeScale) * zoomedStrokeScale
        let backgroundBoxes = CDKMarkushRendering.rGroupBoxes(molecule: depictionMolecule,
                                                              positionsByAtomID: positionByAtomID,
                                                              style: style,
                                                              padding: max(1.4, style.bondWidth * 0.55 * viewportStrokeScale) * zoomedStrokeScale,
                                                              fontSize: renderedLabelFontSize,
                                                              bondWidth: baseBondWidth,
                                                              scale: viewportStrokeScale * zoomedStrokeScale)
            .map {
                CDKMetalDepictionScene.BackgroundBox(rect: $0.boxRect,
                                                     cornerRadius: $0.cornerRadius,
                                                     color: $0.fillColor,
                                                     opacity: $0.fillColor.alpha)
            }
        let linkAnnotations = CDKMarkushRendering.linkAnnotations(molecule: depictionMolecule,
                                                                  positionsByAtomID: positionByAtomID,
                                                                  fontSize: renderedLabelFontSize)
        let sgroupBracketAnnotations = CDKSgroupRendering.bracketAnnotations(molecule: depictionMolecule,
                                                                             positionsByAtomID: positionByAtomID,
                                                                             transformPoint: renderedPoint(for:),
                                                                             fontSize: renderedLabelFontSize,
                                                                             bondWidth: baseBondWidth)

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
            let inset = max(1.6, style.bondWidth * 0.95 * viewportStrokeScale) * zoomedStrokeScale
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
                appendSegment(left,
                              right,
                              width: max(step * 1.35, baseBondWidth * 0.95 * stereoScale),
                              opacity: 0.985,
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

            let maxHalfWidth = min(stereoHashedWedgeHalfWidth, len * 0.34)
            let tickCount = max(9, min(22, Int(len / max(4.8, baseBondWidth * 1.55))))
            for i in 1...tickCount {
                let t = CGFloat(i) / CGFloat(tickCount + 1)
                let cx = start.x + dx * t
                let cy = start.y + dy * t
                let halfTick = max(baseBondWidth * 0.12 * stereoScale * hashedStereoScale, t * maxHalfWidth)
                let left = CGPoint(x: cx + px * halfTick, y: cy + py * halfTick)
                let right = CGPoint(x: cx - px * halfTick, y: cy - py * halfTick)
                appendSegment(left,
                              right,
                              width: max(0.75, baseBondWidth * 0.52 * stereoScale * hashedStereoScale),
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

        func appendAttachmentPointGlyph(center: CGPoint,
                                        other: CGPoint,
                                        color: CDKRenderColor,
                                        into segments: inout [CDKMetalDepictionScene.LineSegment]) {
            for segment in CDKMarkushRendering.attachmentPointSegments(center: center,
                                                                       other: other,
                                                                       style: style,
                                                                       scale: viewportStrokeScale * zoomedStrokeScale) {
                appendSegment(segment.from,
                              segment.to,
                              width: max(1.0, baseBondWidth * 0.92),
                              opacity: 0.95,
                              color: color,
                              into: &segments)
            }
        }

        func appendDashedSegment(_ from: CGPoint,
                                 _ to: CGPoint,
                                 width: CGFloat,
                                 opacity: CGFloat,
                                 color: CDKRenderColor,
                                 into segments: inout [CDKMetalDepictionScene.LineSegment]) {
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

        func appendGlowSegment(_ from: CGPoint,
                               _ to: CGPoint,
                               width: CGFloat,
                               into segments: inout [CDKMetalDepictionScene.LineSegment]) {
            let glowWidth = max(width * 2.35, width + max(1.6, style.bondWidth * 1.8))
            appendSegment(from,
                          to,
                          width: glowWidth,
                          opacity: 0.34,
                          color: .outerGlowHighlight,
                          into: &segments)
        }

        func appendQuerySecondaryLine(from p1: CGPoint,
                                      to p2: CGPoint,
                                      px: CGFloat,
                                      py: CGFloat,
                                      centers: [CGPoint],
                                      paritySeed: Int,
                                      color: CDKRenderColor,
                                      dashed: Bool,
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
            if dashed {
                appendDashedSegment(a,
                                    b,
                                    width: max(1.0, baseBondWidth * 0.82),
                                    opacity: 0.92,
                                    color: color,
                                    into: &segments)
            } else {
                appendSegment(a,
                              b,
                              width: max(1.0, baseBondWidth * 0.82),
                              opacity: 0.92,
                              color: color,
                              into: &segments)
            }
        }

        var bondSegments: [CDKMetalDepictionScene.LineSegment] = []
        bondSegments.reserveCapacity(max(16, depictionMolecule.bonds.count * 3))

        for bond in depictionMolecule.bonds {
            guard let p1 = positionByAtomID[bond.a1], let p2 = positionByAtomID[bond.a2] else { continue }
            let glowHighlighted = usesGlowHighlights(style: style) && highlightedBondIDs.contains(bond.id)
            let bondColor = CDKRenderingStyleResolver.bondColor(for: bond,
                                                                molecule: depictionMolecule,
                                                                style: style,
                                                                highlightedBondIDs: highlightedBondIDs)

            let dx = p2.x - p1.x
            let dy = p2.y - p1.y
            let len = max(0.0001, hypot(dx, dy))
            let px = -dy / len
            let py = dx / len
            let edgeKey = MetalRenderEdgeKey(bond.a1, bond.a2)
            let aromaticCenters = aromaticEdgeCenters[edgeKey] ?? []
            let conjugatedCenters = conjugatedDoubleEdgeCenters[edgeKey] ?? []
            let aromaticStyled = bond.order == .aromatic
                || (style.aromaticDisplayMode == .circle && aromaticBondIDs.contains(bond.id))

            func appendQueryBondSegments() {
                let queryColor = bondColor
                switch bond.queryType {
                case .any:
                    if glowHighlighted {
                        appendGlowSegment(p1, p2, width: baseBondWidth, into: &bondSegments)
                    }
                    appendDashedSegment(p1,
                                        p2,
                                        width: max(1.0, baseBondWidth * 0.92),
                                        opacity: 0.94,
                                        color: queryColor,
                                        into: &bondSegments)
                case .singleOrDouble:
                    if glowHighlighted {
                        appendGlowSegment(p1, p2, width: baseBondWidth, into: &bondSegments)
                    }
                    appendSegment(p1,
                                  p2,
                                  width: baseBondWidth,
                                  opacity: 0.95,
                                  color: queryColor,
                                  into: &bondSegments)
                    appendQuerySecondaryLine(from: p1,
                                             to: p2,
                                             px: px,
                                             py: py,
                                             centers: conjugatedCenters.isEmpty ? aromaticCenters : conjugatedCenters,
                                             paritySeed: bond.id,
                                             color: queryColor,
                                             dashed: true,
                                             into: &bondSegments)
                case .singleOrAromatic:
                    if glowHighlighted {
                        appendGlowSegment(p1, p2, width: max(1.0, baseBondWidth * 0.92), into: &bondSegments)
                    }
                    appendSegment(p1,
                                  p2,
                                  width: max(1.0, baseBondWidth * 0.92),
                                  opacity: 0.90,
                                  color: queryColor,
                                  into: &bondSegments)
                    appendQuerySecondaryLine(from: p1,
                                             to: p2,
                                             px: px,
                                             py: py,
                                             centers: aromaticCenters,
                                             paritySeed: bond.id,
                                             color: queryColor,
                                             dashed: true,
                                             into: &bondSegments)
                case .doubleOrAromatic:
                    if glowHighlighted {
                        appendGlowSegment(p1, p2, width: baseBondWidth, into: &bondSegments)
                    }
                    appendSegment(p1.offsetBy(dx: px * doubleBondHalfSeparation, dy: py * doubleBondHalfSeparation),
                                  p2.offsetBy(dx: px * doubleBondHalfSeparation, dy: py * doubleBondHalfSeparation),
                                  width: baseBondWidth,
                                  opacity: 0.95,
                                  color: queryColor,
                                  into: &bondSegments)
                    appendQuerySecondaryLine(from: p1,
                                             to: p2,
                                             px: px,
                                             py: py,
                                             centers: aromaticCenters,
                                             paritySeed: bond.id,
                                             color: queryColor,
                                             dashed: true,
                                             into: &bondSegments)
                case .none:
                    break
                }
            }

            if bond.queryType != nil {
                appendQueryBondSegments()
                if atomByID[bond.a1]?.attachmentPoint != nil {
                    appendAttachmentPointGlyph(center: p1, other: p2, color: bondColor, into: &bondSegments)
                }
                if atomByID[bond.a2]?.attachmentPoint != nil {
                    appendAttachmentPointGlyph(center: p2, other: p1, color: bondColor, into: &bondSegments)
                }
                continue
            }

            if aromaticStyled {
                if glowHighlighted {
                    appendGlowSegment(p1, p2, width: max(1.1, baseBondWidth * 0.88), into: &bondSegments)
                }
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
                if atomByID[bond.a1]?.attachmentPoint != nil {
                    appendAttachmentPointGlyph(center: p1, other: p2, color: bondColor, into: &bondSegments)
                }
                if atomByID[bond.a2]?.attachmentPoint != nil {
                    appendAttachmentPointGlyph(center: p2, other: p1, color: bondColor, into: &bondSegments)
                }
                continue
            }

            switch bond.order {
            case .single:
                switch bond.stereo {
                case .up:
                    if glowHighlighted {
                        appendGlowSegment(p1, p2, width: baseBondWidth, into: &bondSegments)
                    }
                    appendSolidWedge(from: p1, to: p2, color: bondColor, into: &bondSegments)
                case .upReversed:
                    if glowHighlighted {
                        appendGlowSegment(p1, p2, width: baseBondWidth, into: &bondSegments)
                    }
                    appendSolidWedge(from: p2, to: p1, color: bondColor, into: &bondSegments)
                case .down:
                    if glowHighlighted {
                        appendGlowSegment(p1, p2, width: baseBondWidth, into: &bondSegments)
                    }
                    appendHashedWedge(from: p1, to: p2, color: bondColor, into: &bondSegments)
                case .downReversed:
                    if glowHighlighted {
                        appendGlowSegment(p1, p2, width: baseBondWidth, into: &bondSegments)
                    }
                    appendHashedWedge(from: p2, to: p1, color: bondColor, into: &bondSegments)
                case .either:
                    if glowHighlighted {
                        appendGlowSegment(p1, p2, width: max(1.0, baseBondWidth * 0.72), into: &bondSegments)
                    }
                    appendSegment(p1, p2, width: max(1.0, baseBondWidth * 0.72), opacity: 0.78, color: bondColor, into: &bondSegments)
                case .none:
                    if glowHighlighted {
                        appendGlowSegment(p1, p2, width: baseBondWidth, into: &bondSegments)
                    }
                    appendSegment(p1, p2, width: baseBondWidth, opacity: 0.95, color: bondColor, into: &bondSegments)
                }
            case .double:
                if glowHighlighted {
                    appendGlowSegment(p1, p2, width: baseBondWidth, into: &bondSegments)
                }
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
                if glowHighlighted {
                    appendGlowSegment(p1, p2, width: baseBondWidth, into: &bondSegments)
                }
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
                if glowHighlighted {
                    appendGlowSegment(p1, p2, width: max(1.1, baseBondWidth * 0.88), into: &bondSegments)
                }
                appendSegment(p1,
                              p2,
                              width: max(1.1, baseBondWidth * 0.88),
                              opacity: 0.80,
                              color: bondColor,
                              into: &bondSegments)
            }

            if atomByID[bond.a1]?.attachmentPoint != nil {
                appendAttachmentPointGlyph(center: p1, other: p2, color: bondColor, into: &bondSegments)
            }
            if atomByID[bond.a2]?.attachmentPoint != nil {
                appendAttachmentPointGlyph(center: p2, other: p1, color: bondColor, into: &bondSegments)
            }
        }

        for annotation in sgroupBracketAnnotations {
            for segment in annotation.segments {
                appendSegment(segment.from,
                              segment.to,
                              width: max(0.95, baseBondWidth * 0.84),
                              opacity: 0.92,
                              color: .ink,
                              into: &bondSegments)
            }
        }

        if style.aromaticDisplayMode == .circle {
            for ring in aromaticRings where ring.count >= 5 {
                let ringPoints = ring.compactMap { positionByAtomID[$0] }
                guard ringPoints.count == ring.count else { continue }
                let ringColor = CDKRenderingStyleResolver.aromaticRingColor(atomIDs: ring,
                                                                            molecule: depictionMolecule,
                                                                            style: style,
                                                                            highlightedBondIDs: highlightedBondIDs)
                let ringIsHighlighted = usesGlowHighlights(style: style) && (0..<ring.count).contains { index in
                    let atom1 = ring[index]
                    let atom2 = ring[(index + 1) % ring.count]
                    guard let bond = depictionMolecule.bond(between: atom1, and: atom2) else { return false }
                    return highlightedBondIDs.contains(bond.id)
                }
                if ringIsHighlighted {
                    appendAromaticCircle(points: ringPoints,
                                         color: .outerGlowHighlight,
                                         width: max(1.0, baseBondWidth * 1.95),
                                         opacity: 0.26,
                                         into: &bondSegments)
                }
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

        for atom in depictionMolecule.atoms {
            let atomDegree = degree[atom.id] ?? 0
            guard CDKLabelText.shouldDrawLabel(atom: atom,
                                               degree: atomDegree,
                                               style: style,
                                               includeAromaticCarbonLabelsWhenCarbonsHidden: includeAromaticCarbonLabelsWhenCarbonsHidden,
                                               includeTerminalCarbonLabelsWhenCarbonsHidden: includeTerminalCarbonLabelsWhenCarbonsHidden,
                                               molecule: depictionMolecule,
                                               highlightedAtomIDs: highlightedAtomIDs,
                                               highlightedBondIDs: highlightedBondIDs),
                  let anchor = basePositionByAtomID[atom.id] else { continue }

            let implicitH = style.showImplicitHydrogens ? depictionMolecule.implicitHydrogenCount(for: atom.id) : 0
            let text = CDKLabelText.build(atom: atom, style: style, implicitHydrogenCount: implicitH)

            let neighborPositions = depictionMolecule.neighbors(of: atom.id).compactMap { basePositionByAtomID[$0] }
            let estimated = estimateLabelSize(text: text, fontSize: labelLayoutFontSize)
            let atomColor = CDKRenderingStyleResolver.atomColor(for: atom,
                                                                style: style,
                                                                highlightedAtomIDs: highlightedAtomIDs)
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

        let labelClipPadding = max(1.4, style.bondWidth * 0.55 * viewportStrokeScale) * zoomedStrokeScale
        let labelClipObstacles: [CDKLabelObstacle] = pendingLabels.compactMap { label in
            guard let center = resolvedCentersByAtomID[label.atomID] else { return nil }
            let renderedCenter = renderedPosition(for: label.atomID, basePosition: center)
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
        labels.reserveCapacity((pendingLabels.count * 2) + (linkAnnotations.count * 4) + (sgroupBracketAnnotations.count * 2) + backgroundBoxes.count)
        for label in pendingLabels {
            let resolvedBaseCenter = resolvedCentersByAtomID[label.atomID] ?? label.anchor
            let renderedCenter = renderedPosition(for: label.atomID, basePosition: resolvedBaseCenter)
            let atomHighlighted = usesGlowHighlights(style: style) && highlightedAtomIDs.contains(label.atomID)
            if atomHighlighted {
                labels.append(CDKMetalDepictionScene.AtomLabel(id: -(label.atomID + 1_000_000),
                                                               text: label.text,
                                                               position: renderedCenter,
                                                               fontSize: renderedLabelFontSize,
                                                               aromatic: false,
                                                               color: CDKRenderColor.outerGlowHighlight.withAlpha(0.92),
                                                               italicized: CDKLabelText.usesItalicRGroupFont(text: label.text),
                                                               drawsBackground: false,
                                                               usesGlowOverlay: true,
                                                               suppressesMatchedBackground: !preservesBackgroundDuringGlow(style: style)))
            }
            labels.append(CDKMetalDepictionScene.AtomLabel(id: label.atomID,
                                                           text: label.text,
                                                           position: renderedCenter,
                                                           fontSize: renderedLabelFontSize,
                                                           aromatic: label.aromatic,
                                                           color: label.color,
                                                           italicized: CDKLabelText.usesItalicRGroupFont(text: label.text),
                                                           drawsBackground: depictionMolecule.atom(id: label.atomID)?.rGroupMembership == nil))
        }

        var nextAnnotationID = -1

        func appendAnnotationLabel(text: String,
                                   position: CGPoint,
                                   fontSize: CGFloat,
                                   drawsBackground: Bool = true) {
            labels.append(CDKMetalDepictionScene.AtomLabel(id: nextAnnotationID,
                                                           text: text,
                                                           position: position,
                                                           fontSize: fontSize,
                                                           aromatic: false,
                                                           color: .ink,
                                                           italicized: CDKLabelText.usesItalicRGroupFont(text: text),
                                                           drawsBackground: drawsBackground))
            nextAnnotationID -= 1
        }

        for annotation in linkAnnotations {
            appendAnnotationLabel(text: "(",
                                  position: annotation.leftBracketPosition,
                                  fontSize: renderedLabelFontSize * 0.86,
                                  drawsBackground: false)
            appendAnnotationLabel(text: ")",
                                  position: annotation.rightBracketPosition,
                                  fontSize: renderedLabelFontSize * 0.86,
                                  drawsBackground: false)
            if let subscriptText = annotation.subscriptText,
               let subscriptPosition = annotation.subscriptPosition {
                appendAnnotationLabel(text: subscriptText,
                                      position: subscriptPosition,
                                      fontSize: renderedLabelFontSize * 0.78,
                                      drawsBackground: false)
            }
            if let superscriptText = annotation.superscriptText,
               let superscriptPosition = annotation.superscriptPosition {
                appendAnnotationLabel(text: superscriptText,
                                      position: superscriptPosition,
                                      fontSize: renderedLabelFontSize * 0.72,
                                      drawsBackground: false)
            }
        }

        for annotation in sgroupBracketAnnotations {
            if let subscriptText = annotation.subscriptText,
               let subscriptPosition = annotation.subscriptPosition {
                appendAnnotationLabel(text: subscriptText,
                                      position: subscriptPosition,
                                      fontSize: renderedLabelFontSize * 0.80,
                                      drawsBackground: false)
            }
            if let superscriptText = annotation.superscriptText,
               let superscriptPosition = annotation.superscriptPosition {
                appendAnnotationLabel(text: superscriptText,
                                      position: superscriptPosition,
                                      fontSize: renderedLabelFontSize * 0.74,
                                      drawsBackground: false)
            }
        }

        for annotation in CDKMarkushRendering.rGroupBoxes(molecule: depictionMolecule,
                                                          positionsByAtomID: positionByAtomID,
                                                          style: style,
                                                          padding: labelClipPadding,
                                                          fontSize: renderedLabelFontSize,
                                                          bondWidth: baseBondWidth,
                                                          scale: viewportStrokeScale * zoomedStrokeScale) {
            appendAnnotationLabel(text: annotation.label,
                                  position: annotation.labelPosition,
                                  fontSize: renderedLabelFontSize,
                                  drawsBackground: false)
        }

        return CDKMetalDepictionScene(gridSegments: gridSegments,
                                      backgroundBoxes: backgroundBoxes,
                                      bondSegments: bondSegments,
                                      labels: labels)
    }
}

private extension Comparable {
    func clamped(to range: ClosedRange<Self>) -> Self {
        min(max(self, range.lowerBound), range.upperBound)
    }
}
