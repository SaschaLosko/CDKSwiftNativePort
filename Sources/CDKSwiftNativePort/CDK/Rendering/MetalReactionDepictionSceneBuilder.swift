#if canImport(CoreGraphics)
import CoreGraphics
#endif
import Foundation

public struct CDKReactionParticipantSelection: Equatable, Hashable, Codable, Sendable {
    public let role: CDKReactionRole
    public let index: Int

    public init(role: CDKReactionRole, index: Int) {
        self.role = role
        self.index = index
    }
}

/// Metal scene builder for CDK-like 2D reaction depictions.
public enum CDKMetalReactionDepictionSceneBuilder {
    public static func build(reaction: CDKReaction,
                             style: RenderStyle,
                             canvasRect: CGRect,
                             zoom: CGFloat,
                             pan: CGSize,
                             rotationDegrees: CGFloat = 0,
                             highlightedParticipant: CDKReactionParticipantSelection? = nil,
                             withOuterGlowHighlight: Bool = true) -> CDKMetalDepictionScene {
        guard canvasRect.width > 1, canvasRect.height > 1 else {
            return CDKMetalDepictionScene(gridSegments: [], bondSegments: [], labels: [])
        }

        let gridSegments: [CDKMetalDepictionScene.LineSegment] = []
        let pad = style.padding
        let available = CGRect(x: canvasRect.minX + pad,
                               y: canvasRect.minY + pad,
                               width: max(1, canvasRect.width - 2 * pad),
                               height: max(1, canvasRect.height - 2 * pad))

        var reactants = reaction.reactants.map { CDKReactionParticipantLayoutRefiner.recomputeIfLowQuality($0) }
        let products = reaction.products.map { CDKReactionParticipantLayoutRefiner.recomputeIfLowQuality($0) }
        var agents = reaction.agents.map { CDKReactionParticipantLayoutRefiner.recomputeIfLowQuality($0) }

        let productMapTemplate = reactionMapTemplate(from: products)
        if !productMapTemplate.isEmpty {
            reactants = reactants.map {
                CDKReactionParticipantLayoutRefiner.reorientForMapTemplate($0, mapTemplate: productMapTemplate)
            }
            agents = agents.map {
                CDKReactionParticipantLayoutRefiner.reorientForMapTemplate($0, mapTemplate: productMapTemplate)
            }
        }

        if reactants.isEmpty && products.isEmpty && agents.isEmpty {
            return CDKMetalDepictionScene(gridSegments: gridSegments, bondSegments: [], labels: [])
        }

        let activeHighlight = validatedHighlight(highlightedParticipant,
                                                 reactantCount: reactants.count,
                                                 agentCount: agents.count,
                                                 productCount: products.count)

        // CDK-like reaction depictions keep labels comparatively compact, even if the
        // single-molecule depiction style uses a very large default atom font.
        let participantMolecules = reactants + agents + products
        let reactionStyle = compactReactionStyle(base: style,
                                                 participants: participantMolecules,
                                                 in: available)

        let mainRegion = available
        let splitRegions = reactionSideRegions(in: mainRegion,
                                               reactants: reactants,
                                               products: products,
                                               preferredArrowGap: preferredArrowGap(for: agents,
                                                                                   in: mainRegion,
                                                                                   style: reactionStyle))
        let arrowGap = splitRegions.arrowGap
        let leftRegion = splitRegions.leftRegion
        let rightRegion = splitRegions.rightRegion
        let agentRegion = agentRegion(in: mainRegion,
                                      splitRegions: splitRegions,
                                      agentCount: agents.count)

        let reactantLayout = layoutBoxes(for: reactants, in: leftRegion, style: reactionStyle)
        let productLayout = layoutBoxes(for: products, in: rightRegion, style: reactionStyle)
        let agentLayout = layoutBoxes(for: agents, in: agentRegion, style: reactionStyle)
        let participantScale = sharedParticipantScale(for: [reactantLayout, productLayout, agentLayout],
                                                       style: reactionStyle)

        var baseBondSegments: [CDKMetalDepictionScene.LineSegment] = []
        baseBondSegments.reserveCapacity((reactants.count + products.count + agents.count) * 32)
        var baseBackgroundBoxes: [CDKMetalDepictionScene.BackgroundBox] = []
        baseBackgroundBoxes.reserveCapacity((reactants.count + products.count + agents.count) * 2)
        var baseLabels: [CDKMetalDepictionScene.AtomLabel] = []
        baseLabels.reserveCapacity((reactants.count + products.count + agents.count) * 24)
        var nextLabelID = 1_000_000

        func shouldGlow(role: CDKReactionRole, index: Int) -> Bool {
            guard withOuterGlowHighlight, let activeHighlight else { return false }
            return activeHighlight.role == role && activeHighlight.index == index
        }

        func mergeMolecule(_ molecule: Molecule, into rect: CGRect, highlighted: Bool) {
            let renderRect = participantCanvasRect(for: molecule,
                                                   in: rect,
                                                   style: reactionStyle,
                                                   sharedScale: participantScale)
            let scene = CDKMetalDepictionSceneBuilder.build(molecule: molecule,
                                                            style: reactionStyle,
                                                            canvasRect: renderRect,
                                                            zoom: 1.0,
                                                            pan: .zero,
                                                            rotationDegrees: 0)
            let glowColor = CDKRenderColor.outerGlowHighlight.withAlpha(0.92)
            for box in scene.backgroundBoxes {
                baseBackgroundBoxes.append(box)
            }
            for segment in scene.bondSegments {
                if highlighted {
                    let glowWidth = max(segment.width * 2.45,
                                        segment.width + max(1.8, reactionStyle.bondWidth * 1.95))
                    let glowOpacity = (0.16 + (segment.opacity * 0.26)).clamped(to: 0.16...0.44)
                    baseBondSegments.append(
                        CDKMetalDepictionScene.LineSegment(from: segment.from,
                                                           to: segment.to,
                                                           width: glowWidth,
                                                           opacity: glowOpacity,
                                                           color: glowColor)
                    )
                }
                baseBondSegments.append(
                    CDKMetalDepictionScene.LineSegment(from: segment.from,
                                                       to: segment.to,
                                                       width: segment.width,
                                                       opacity: segment.opacity,
                                                       color: segment.color)
                )
            }
            for label in scene.labels {
                if highlighted {
                    nextLabelID += 1
                    baseLabels.append(
                        CDKMetalDepictionScene.AtomLabel(id: nextLabelID,
                                                         text: label.text,
                                                         position: label.position,
                                                         fontSize: max(7.2, label.fontSize),
                                                         aromatic: false,
                                                         color: CDKRenderColor.outerGlowHighlight.withAlpha(0.92),
                                                         italicized: label.italicized,
                                                         drawsBackground: false,
                                                         usesGlowOverlay: true,
                                                         suppressesMatchedBackground: true)
                    )
                }
                nextLabelID += 1
                baseLabels.append(
                    CDKMetalDepictionScene.AtomLabel(id: nextLabelID,
                                                     text: label.text,
                                                     position: label.position,
                                                     fontSize: max(7.2, label.fontSize),
                                                     aromatic: label.aromatic,
                                                     color: label.color,
                                                     italicized: label.italicized,
                                                     drawsBackground: label.drawsBackground,
                                                     usesGlowOverlay: label.usesGlowOverlay,
                                                     suppressesMatchedBackground: label.suppressesMatchedBackground)
                )
            }
        }

        for (index, item) in reactantLayout.items.enumerated() {
            mergeMolecule(item.molecule,
                          into: item.rect,
                          highlighted: shouldGlow(role: .reactant, index: index))
        }
        for (index, item) in productLayout.items.enumerated() {
            mergeMolecule(item.molecule,
                          into: item.rect,
                          highlighted: shouldGlow(role: .product, index: index))
        }
        for (index, item) in agentLayout.items.enumerated() {
            mergeMolecule(item.molecule,
                          into: item.rect,
                          highlighted: shouldGlow(role: .agent, index: index))
        }

        let symbolFontSize = min(12, max(9, reactionStyle.fontSize * 0.66))
        for point in reactantLayout.plusPositions + productLayout.plusPositions + agentLayout.plusPositions {
            nextLabelID += 1
            baseLabels.append(
                CDKMetalDepictionScene.AtomLabel(id: nextLabelID,
                                                 text: "+",
                                                 position: point,
                                                 fontSize: symbolFontSize,
                                                 aromatic: false,
                                                 color: .ink)
            )
        }

        let arrowY = mainRegion.midY
        var arrowStart = CGPoint(x: leftRegion.maxX + (arrowGap * 0.14), y: arrowY)
        var arrowEnd = CGPoint(x: rightRegion.minX - (arrowGap * 0.14), y: arrowY)
        if arrowEnd.x <= arrowStart.x {
            let halfWidth = min(60, mainRegion.width * 0.12)
            arrowStart = CGPoint(x: mainRegion.midX - halfWidth, y: arrowY)
            arrowEnd = CGPoint(x: mainRegion.midX + halfWidth, y: arrowY)
        }

        let arrowWidth = max(1.0, reactionStyle.bondWidth * 0.94)
        baseBondSegments.append(
            CDKMetalDepictionScene.LineSegment(from: arrowStart,
                                               to: arrowEnd,
                                               width: arrowWidth,
                                               opacity: 0.95,
                                               color: .ink)
        )

        let headLength = max(11, reactionStyle.fontSize * 0.54)
        let headHalf = headLength * 0.36
        let headBack = CGPoint(x: arrowEnd.x - headLength, y: arrowEnd.y)
        baseBondSegments.append(
            CDKMetalDepictionScene.LineSegment(from: CGPoint(x: headBack.x, y: headBack.y - headHalf),
                                               to: arrowEnd,
                                               width: arrowWidth,
                                               opacity: 0.95,
                                               color: .ink)
        )
        baseBondSegments.append(
            CDKMetalDepictionScene.LineSegment(from: CGPoint(x: headBack.x, y: headBack.y + headHalf),
                                               to: arrowEnd,
                                               width: arrowWidth,
                                               opacity: 0.95,
                                               color: .ink)
        )

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

        // Keep line and label sizing coupled to reaction viewport auto-fit so resizing
        // the window behaves similarly to explicit zooming.
        let referenceMainRegion = CGSize(width: 760, height: 320)
        let reactionFitScaleX = mainRegion.width / max(1, referenceMainRegion.width)
        let reactionFitScaleY = mainRegion.height / max(1, referenceMainRegion.height)
        // Keep an upper cap, but no lower clamp: compact windows should continue
        // scaling line and label sizing down with geometry.
        let reactionAutoFitScale = min(1.85, min(reactionFitScaleX, reactionFitScaleY))

        let transformedBonds = baseBondSegments.map { segment in
            CDKMetalDepictionScene.LineSegment(from: applyViewportTransform(segment.from),
                                               to: applyViewportTransform(segment.to),
                                               width: max(0.75, segment.width * zoom * reactionAutoFitScale),
                                               opacity: segment.opacity,
                                               color: segment.color)
        }
        let transformedBoxes = baseBackgroundBoxes.map { box in
            let minPoint = applyViewportTransform(CGPoint(x: box.rect.minX, y: box.rect.minY))
            let maxPoint = applyViewportTransform(CGPoint(x: box.rect.maxX, y: box.rect.maxY))
            let rect = CGRect(x: min(minPoint.x, maxPoint.x),
                              y: min(minPoint.y, maxPoint.y),
                              width: abs(maxPoint.x - minPoint.x),
                              height: abs(maxPoint.y - minPoint.y))
            return CDKMetalDepictionScene.BackgroundBox(rect: rect,
                                                        cornerRadius: max(2.0, box.cornerRadius * zoom * reactionAutoFitScale),
                                                        color: box.color,
                                                        opacity: box.opacity)
        }
        let transformedLabels = baseLabels.map { label in
            CDKMetalDepictionScene.AtomLabel(id: label.id,
                                             text: label.text,
                                             position: applyViewportTransform(label.position),
                                             fontSize: (label.fontSize * zoom * reactionAutoFitScale).clamped(to: 7.0...16.0),
                                             aromatic: label.aromatic,
                                             color: label.color,
                                             italicized: label.italicized,
                                             drawsBackground: label.drawsBackground,
                                             usesGlowOverlay: label.usesGlowOverlay,
                                             suppressesMatchedBackground: label.suppressesMatchedBackground)
        }

        return CDKMetalDepictionScene(gridSegments: gridSegments,
                                      backgroundBoxes: transformedBoxes,
                                      bondSegments: transformedBonds,
                                      labels: transformedLabels)
    }

    public static func participant(at point: CGPoint,
                                   in reaction: CDKReaction,
                                   style: RenderStyle,
                                   canvasRect: CGRect,
                                   zoom: CGFloat,
                                   pan: CGSize,
                                   rotationDegrees: CGFloat = 0,
                                   hitSlop: CGFloat = 12) -> CDKReactionParticipantSelection? {
        guard canvasRect.width > 1, canvasRect.height > 1 else {
            return nil
        }

        var reactants = reaction.reactants.map { CDKReactionParticipantLayoutRefiner.recomputeIfLowQuality($0) }
        let products = reaction.products.map { CDKReactionParticipantLayoutRefiner.recomputeIfLowQuality($0) }
        var agents = reaction.agents.map { CDKReactionParticipantLayoutRefiner.recomputeIfLowQuality($0) }

        let productMapTemplate = reactionMapTemplate(from: products)
        if !productMapTemplate.isEmpty {
            reactants = reactants.map {
                CDKReactionParticipantLayoutRefiner.reorientForMapTemplate($0, mapTemplate: productMapTemplate)
            }
            agents = agents.map {
                CDKReactionParticipantLayoutRefiner.reorientForMapTemplate($0, mapTemplate: productMapTemplate)
            }
        }

        if reactants.isEmpty && products.isEmpty && agents.isEmpty {
            return nil
        }

        let pad = style.padding
        let available = CGRect(x: canvasRect.minX + pad,
                               y: canvasRect.minY + pad,
                               width: max(1, canvasRect.width - 2 * pad),
                               height: max(1, canvasRect.height - 2 * pad))
        let reactionStyle = compactReactionStyle(base: style,
                                                 participants: reactants + agents + products,
                                                 in: available)

        let mainRegion = available
        let splitRegions = reactionSideRegions(in: mainRegion,
                                               reactants: reactants,
                                               products: products,
                                               preferredArrowGap: preferredArrowGap(for: agents,
                                                                                   in: mainRegion,
                                                                                   style: reactionStyle))
        let leftRegion = splitRegions.leftRegion
        let rightRegion = splitRegions.rightRegion
        let agentRegion = agentRegion(in: mainRegion,
                                      splitRegions: splitRegions,
                                      agentCount: agents.count)

        let reactantLayout = layoutBoxes(for: reactants, in: leftRegion, style: reactionStyle)
        let productLayout = layoutBoxes(for: products, in: rightRegion, style: reactionStyle)
        let agentLayout = layoutBoxes(for: agents, in: agentRegion, style: reactionStyle)

        let untransformedPoint = inverseViewportTransform(point,
                                                          canvasRect: canvasRect,
                                                          zoom: zoom,
                                                          pan: pan,
                                                          rotationDegrees: rotationDegrees)
        let baseSlop = max(2, hitSlop / max(0.25, abs(zoom)))

        if let selection = hitTestParticipant(in: reactantLayout.items,
                                              role: .reactant,
                                              point: untransformedPoint,
                                              slop: baseSlop) {
            return selection
        }
        if let selection = hitTestParticipant(in: agentLayout.items,
                                              role: .agent,
                                              point: untransformedPoint,
                                              slop: baseSlop) {
            return selection
        }
        if let selection = hitTestParticipant(in: productLayout.items,
                                              role: .product,
                                              point: untransformedPoint,
                                              slop: baseSlop) {
            return selection
        }
        return nil
    }

    private struct MoleculeBoxLayout {
        struct Item {
            let molecule: Molecule
            let rect: CGRect
        }

        let items: [Item]
        let plusPositions: [CGPoint]
    }

    private struct ReactionSideRegions {
        let leftRegion: CGRect
        let rightRegion: CGRect
        let arrowGap: CGFloat
    }

    private static func sharedParticipantScale(for layouts: [MoleculeBoxLayout],
                                               style: RenderStyle) -> CGFloat? {
        let candidates = layouts.flatMap(\.items).compactMap { item in
            participantScaleCandidate(for: item.molecule,
                                      in: item.rect,
                                      style: style)
        }
        guard !candidates.isEmpty else {
            return nil
        }
        return candidates.min()
    }

    private static func participantScaleCandidate(for molecule: Molecule,
                                                  in rect: CGRect,
                                                  style: RenderStyle) -> CGFloat? {
        guard let box = molecule.boundingBox() else {
            return nil
        }

        let xPadding = min(style.padding, rect.width * 0.25)
        let yPadding = min(style.padding, rect.height * 0.25)
        let availableWidth = max(1, rect.width - (2 * xPadding))
        let availableHeight = max(1, rect.height - (2 * yPadding))

        var candidates: [CGFloat] = []
        if box.width > 0.0001 {
            candidates.append(availableWidth / box.width)
        }
        if box.height > 0.0001 {
            candidates.append(availableHeight / box.height)
        }

        return candidates.min()
    }

    private static func participantCanvasRect(for molecule: Molecule,
                                              in rect: CGRect,
                                              style: RenderStyle,
                                              sharedScale: CGFloat?) -> CGRect {
        guard let sharedScale,
              sharedScale.isFinite,
              sharedScale > 0,
              let box = molecule.boundingBox() else {
            return rect
        }

        let minimumCoordinateExtent: CGFloat = 1.2
        let contentWidth = max(box.width, minimumCoordinateExtent) * sharedScale
        let contentHeight = max(box.height, minimumCoordinateExtent) * sharedScale
        let width = max(1, contentWidth + (2 * style.padding))
        let height = max(1, contentHeight + (2 * style.padding))

        return CGRect(x: rect.midX - (width * 0.5),
                      y: rect.midY - (height * 0.5),
                      width: width,
                      height: height)
    }

    private static func layoutBoxes(for molecules: [Molecule], in region: CGRect, style: RenderStyle) -> MoleculeBoxLayout {
        guard !molecules.isEmpty else {
            return MoleculeBoxLayout(items: [], plusPositions: [])
        }

        let count = molecules.count
        let plusCount = max(0, count - 1)
        let plusGlyphWidth = max(14, style.fontSize * 0.58)
        let interItemSpacing = max(8, style.fontSize * 0.30)
        let plusSpace = CGFloat(plusCount) * plusGlyphWidth
        let spacing = CGFloat(plusCount) * interItemSpacing
        let availableWidth = max(1, region.width - plusSpace - spacing)
        let itemHeight = max(1, region.height * 0.90)
        let itemY = region.midY - (itemHeight * 0.5)

        let rawWeights = molecules.map { preferredWidthWeight(for: $0) }
        if shouldWrapParticipants(molecules, in: region, style: style) {
            return wrappedLayoutBoxes(for: molecules,
                                      weights: rawWeights,
                                      in: region,
                                      style: style)
        }
        let weightSum = max(0.0001, rawWeights.reduce(0, +))
        let minimumWidth = min(max(38, availableWidth * 0.08), availableWidth / CGFloat(count))
        let remainderWidth = max(0, availableWidth - (minimumWidth * CGFloat(count)))

        var items: [MoleculeBoxLayout.Item] = []
        items.reserveCapacity(count)
        var plusPositions: [CGPoint] = []
        plusPositions.reserveCapacity(plusCount)

        var x = region.minX
        for idx in molecules.indices {
            let weightedShare = rawWeights[idx] / weightSum
            let proposedWidth = minimumWidth + (remainderWidth * weightedShare)
            let width = min(proposedWidth, region.maxX - x)
            let rect = CGRect(x: x,
                              y: itemY,
                              width: max(1, width),
                              height: itemHeight)
            items.append(MoleculeBoxLayout.Item(molecule: molecules[idx], rect: rect))
            x = rect.maxX
            if idx < molecules.count - 1 {
                let plusX = x + (plusGlyphWidth * 0.5)
                plusPositions.append(CGPoint(x: plusX, y: region.midY))
                x += plusGlyphWidth + interItemSpacing
            }
        }

        return MoleculeBoxLayout(items: items, plusPositions: plusPositions)
    }

    private static func wrappedLayoutBoxes(for molecules: [Molecule],
                                           weights: [CGFloat],
                                           in region: CGRect,
                                           style: RenderStyle) -> MoleculeBoxLayout {
        let rowAssignments = balancedRowAssignments(for: weights, maxRows: min(2, molecules.count))
        let rowCount = max(1, rowAssignments.count)
        let rowGap = max(14, style.fontSize * 0.72)
        let totalGap = CGFloat(max(0, rowCount - 1)) * rowGap
        let rowHeight = max(1, (region.height - totalGap) / CGFloat(rowCount))
        let innerHeight = max(1, rowHeight * 0.88)
        let totalInnerHeight = (CGFloat(rowCount) * innerHeight) + totalGap
        var y = region.midY - (totalInnerHeight * 0.5)

        var items: [MoleculeBoxLayout.Item] = []
        items.reserveCapacity(molecules.count)

        for row in rowAssignments {
            let rowIndices = Array(row)
            let rowSpacing = max(14, style.fontSize * 0.48)
            let dampedWeights = rowIndices.map { sqrt(max(0.25, weights[$0])) }
            let rowWeightSum = max(0.0001, dampedWeights.reduce(CGFloat.zero, +))
            let baseAvailableRowWidth = max(1, region.width - (CGFloat(max(0, rowIndices.count - 1)) * rowSpacing))
            let availableRowWidth: CGFloat
            if rowIndices.count == 1 {
                availableRowWidth = min(baseAvailableRowWidth, region.width * 0.74)
            } else {
                availableRowWidth = baseAvailableRowWidth
            }
            let minimumWidth = min(max(54, availableRowWidth * 0.16),
                                   availableRowWidth / CGFloat(max(1, rowIndices.count)))
            let remainderWidth = max(0, availableRowWidth - (minimumWidth * CGFloat(rowIndices.count)))
            var widths = rowIndices.enumerated().map { offset, idx in
                minimumWidth + (remainderWidth * (dampedWeights[offset] / rowWeightSum))
            }
            if rowIndices.count == 1 {
                widths[0] = min(widths[0], region.width * 0.62)
            } else {
                let maxWidth = region.width * 0.52
                let cappedWidths = widths.map { min($0, maxWidth) }
                let freedWidth = widths.reduce(CGFloat.zero, +) - cappedWidths.reduce(CGFloat.zero, +)
                widths = cappedWidths
                if freedWidth > 0.5 {
                    let uncappedSlots = widths.indices.filter { widths[$0] < maxWidth - 0.5 }
                    if !uncappedSlots.isEmpty {
                        let bonus = freedWidth / CGFloat(uncappedSlots.count)
                        for index in uncappedSlots {
                            widths[index] = min(maxWidth, widths[index] + bonus)
                        }
                    }
                }
            }
            let totalRowWidth = widths.reduce(0, +) + (CGFloat(max(0, rowIndices.count - 1)) * rowSpacing)
            var x = region.midX - (totalRowWidth * 0.5)

            for (offset, idx) in rowIndices.enumerated() {
                let width = widths[offset]
                let rect = CGRect(x: x,
                                  y: y,
                                  width: max(1, width),
                                  height: innerHeight)
                items.append(MoleculeBoxLayout.Item(molecule: molecules[idx], rect: rect))
                x = rect.maxX + rowSpacing
            }

            y += innerHeight + rowGap
        }

        return MoleculeBoxLayout(items: items, plusPositions: [])
    }

    private static func shouldWrapParticipants(_ molecules: [Molecule],
                                               in region: CGRect,
                                               style: RenderStyle) -> Bool {
        guard molecules.count >= 3 else { return false }
        return region.height >= max(132, style.fontSize * 8.5)
    }

    private static func balancedRowAssignments(for weights: [CGFloat],
                                               maxRows: Int) -> [[Int]] {
        let count = weights.count
        guard count > 0 else { return [] }
        guard maxRows > 1, count > 1 else {
            return [Array(0..<count)]
        }

        if maxRows == 2 {
            var prefixSums: [CGFloat] = Array(repeating: 0, count: count + 1)
            for index in weights.indices {
                prefixSums[index + 1] = prefixSums[index] + weights[index]
            }

            let total = prefixSums[count]
            var bestBreak = max(1, min(count - 1, (count + 1) / 2))
            var bestScore = CGFloat.greatestFiniteMagnitude

            for breakIndex in 1..<count {
                let first = prefixSums[breakIndex]
                let second = total - first
                let rowImbalance = abs(first - second)
                let countImbalance = abs(CGFloat(breakIndex) - CGFloat(count - breakIndex))
                let singletonTopPenalty: CGFloat = (breakIndex == 1 && count > 2) ? max(0.45, total * 0.16) : 0
                let score = max(first, second) + (rowImbalance * 0.24) + (countImbalance * 0.30) + singletonTopPenalty
                if score < bestScore {
                    bestScore = score
                    bestBreak = breakIndex
                }
            }

            return [Array(0..<bestBreak), Array(bestBreak..<count)]
        }

        return [Array(0..<count)]
    }

    private static func agentRegion(in mainRegion: CGRect,
                                    splitRegions: ReactionSideRegions,
                                    agentCount: Int) -> CGRect {
        guard agentCount > 0 else {
            return .zero
        }

        let arrowY = mainRegion.midY
        let sectionSpacing: CGFloat = 12
        let preferredHeight = min(max(52, mainRegion.height * 0.18), 116)
        let top = max(mainRegion.minY, arrowY - preferredHeight - sectionSpacing)
        let height = max(1, min(preferredHeight, arrowY - sectionSpacing - top))
        return CGRect(x: splitRegions.leftRegion.maxX,
                      y: top,
                      width: max(1, splitRegions.arrowGap),
                      height: height)
    }

    private static func preferredArrowGap(for agents: [Molecule],
                                          in mainRegion: CGRect,
                                          style: RenderStyle) -> CGFloat? {
        guard !agents.isEmpty else {
            return nil
        }

        let plusCount = max(0, agents.count - 1)
        let plusGlyphWidth = max(14, style.fontSize * 0.58)
        let interItemSpacing = max(8, style.fontSize * 0.30)
        let separatorSpace = CGFloat(plusCount) * (plusGlyphWidth + interItemSpacing)
        let agentDemand = max(CGFloat(agents.count), layoutDemandWeight(for: agents))
        let unitWidth = max(42, min(92, mainRegion.width * 0.06))
        return separatorSpace + (agentDemand * unitWidth)
    }

    private static func reactionSideRegions(in mainRegion: CGRect,
                                            reactants: [Molecule],
                                            products: [Molecule],
                                            preferredArrowGap: CGFloat? = nil) -> ReactionSideRegions {
        let baseArrowGap = max(56, min(120, mainRegion.width * 0.16))
        let maximumArrowGap = max(1, mainRegion.width * 0.46)
        let arrowGap = min(max(baseArrowGap, preferredArrowGap ?? baseArrowGap), maximumArrowGap)
        let sideWidth = max(1, mainRegion.width - arrowGap)

        let reactantDemand = layoutDemandWeight(for: reactants)
        let productDemand = layoutDemandWeight(for: products)
        let leftShareRaw: CGFloat
        if reactants.isEmpty && !products.isEmpty {
            leftShareRaw = 0.28
        } else if products.isEmpty && !reactants.isEmpty {
            leftShareRaw = 0.72
        } else {
            leftShareRaw = reactantDemand / max(0.0001, reactantDemand + productDemand)
        }

        let wrapAware = reactants.count >= 3 || products.count >= 3
        let leftShare = leftShareRaw.clamped(to: wrapAware ? 0.36...0.64 : 0.33...0.67)
        let leftWidth = sideWidth * leftShare
        let rightWidth = sideWidth - leftWidth

        let leftRegion = CGRect(x: mainRegion.minX,
                                y: mainRegion.minY,
                                width: max(1, leftWidth),
                                height: mainRegion.height)
        let rightRegion = CGRect(x: leftRegion.maxX + arrowGap,
                                 y: mainRegion.minY,
                                 width: max(1, rightWidth),
                                 height: mainRegion.height)
        return ReactionSideRegions(leftRegion: leftRegion,
                                   rightRegion: rightRegion,
                                   arrowGap: arrowGap)
    }

    private static func inverseViewportTransform(_ point: CGPoint,
                                                 canvasRect: CGRect,
                                                 zoom: CGFloat,
                                                 pan: CGSize,
                                                 rotationDegrees: CGFloat) -> CGPoint {
        let viewportCenter = CGPoint(x: canvasRect.midX, y: canvasRect.midY)
        let safeZoom = max(0.0001, abs(zoom))

        let dePanned = CGPoint(x: point.x - pan.width, y: point.y - pan.height)
        let scaled = CGPoint(x: ((dePanned.x - viewportCenter.x) / safeZoom) + viewportCenter.x,
                             y: ((dePanned.y - viewportCenter.y) / safeZoom) + viewportCenter.y)

        let normalizedRotation = rotationDegrees.truncatingRemainder(dividingBy: 360)
        let radians = normalizedRotation * (.pi / 180)
        let c = cos(radians)
        let s = sin(radians)
        let rx = scaled.x - viewportCenter.x
        let ry = scaled.y - viewportCenter.y
        return CGPoint(x: (rx * c) + (ry * s) + viewportCenter.x,
                       y: (-rx * s) + (ry * c) + viewportCenter.y)
    }

    private static func hitTestParticipant(in items: [MoleculeBoxLayout.Item],
                                           role: CDKReactionRole,
                                           point: CGPoint,
                                           slop: CGFloat) -> CDKReactionParticipantSelection? {
        var best: (index: Int, distance: CGFloat)?
        for (index, item) in items.enumerated() {
            let expanded = item.rect.insetBy(dx: -slop, dy: -slop)
            guard expanded.contains(point) else { continue }
            let distance = distance(from: point, to: item.rect)
            if let current = best {
                if distance < current.distance {
                    best = (index, distance)
                }
            } else {
                best = (index, distance)
            }
        }

        guard let match = best else { return nil }
        return CDKReactionParticipantSelection(role: role, index: match.index)
    }

    private static func distance(from point: CGPoint, to rect: CGRect) -> CGFloat {
        let dx: CGFloat
        if point.x < rect.minX {
            dx = rect.minX - point.x
        } else if point.x > rect.maxX {
            dx = point.x - rect.maxX
        } else {
            dx = 0
        }

        let dy: CGFloat
        if point.y < rect.minY {
            dy = rect.minY - point.y
        } else if point.y > rect.maxY {
            dy = point.y - rect.maxY
        } else {
            dy = 0
        }

        return hypot(dx, dy)
    }

    private static func compactReactionStyle(base style: RenderStyle,
                                             participants: [Molecule],
                                             in available: CGRect) -> RenderStyle {
        var compact = style
        let participantCount = max(1, participants.count)
        let totalAtoms = max(1, participants.reduce(0) { $0 + $1.atomCount })
        let avgAtoms = CGFloat(totalAtoms) / CGFloat(participantCount)
        let areaPerParticipant = max(1, (available.width * available.height) / CGFloat(participantCount))

        // More participants + larger molecules => proportionally smaller labels.
        let countScale = (1.0 / sqrt(CGFloat(participantCount) / 1.4)).clamped(to: 0.44...1.0)
        let complexityScale = (14.0 / max(8, avgAtoms)).clamped(to: 0.46...1.0)
        let areaScale = (sqrt(areaPerParticipant / 180_000.0)).clamped(to: 0.52...1.0)
        let compactScale = min(countScale, complexityScale, areaScale) * 0.90

        compact.fontSize = (style.fontSize * compactScale).clamped(to: 7.6...11.6)
        compact.bondWidth = (style.bondWidth * (0.72 + (0.12 * compactScale))).clamped(to: 0.92...style.bondWidth)
        compact.padding = max(6, style.padding * 0.42)
        // Reaction depictions are denser than single-molecule scenes.
        // Apply a reaction-only stereo attenuation to keep wedges/hashes CDK-like.
        compact.stereoAttenuation = (0.62 + (0.14 * compactScale)).clamped(to: 0.56...0.78)
        // Additional micro-attenuation for hashed wedges only in reaction mode.
        compact.hashedStereoAttenuation = 0.82
        return compact
    }

    private static func preferredWidthWeight(for molecule: Molecule) -> CGFloat {
        let box = molecule.boundingBox()
        let atomFactor = sqrt(CGFloat(max(4, molecule.atomCount)) / 6.0).clamped(to: 0.8...3.4)
        let heteroCount = molecule.atoms.reduce(0) { partial, atom in
            let element = atom.element.uppercased()
            if element == "C" || element == "H" {
                return partial
            }
            return partial + 1
        }
        let labelPressure = (1.0 + (CGFloat(heteroCount) * 0.065)).clamped(to: 1.0...2.8)
        guard let box else {
            return max(0.5, atomFactor * labelPressure)
        }

        let aspect = (box.width / max(0.1, box.height)).clamped(to: 0.42...3.2)
        let extent = sqrt(max(1, box.width * box.height)).clamped(to: 1.0...10.0)
        let extentFactor = (extent / 2.2).clamped(to: 0.7...2.4)
        return max(0.5, aspect * atomFactor * labelPressure * extentFactor)
    }

    private static func layoutDemandWeight(for molecules: [Molecule]) -> CGFloat {
        guard !molecules.isEmpty else { return 1 }
        let weights = molecules.map { preferredWidthWeight(for: $0) }
        let total = max(1, weights.reduce(CGFloat.zero, +))
        guard molecules.count >= 3 else { return total }

        let rows = balancedRowAssignments(for: weights, maxRows: 2)
        let rowMax = rows.reduce(CGFloat.zero) { partial, row in
            max(partial, row.reduce(CGFloat.zero) { subtotal, idx in
                subtotal + weights[idx]
            })
        }
        return max(1, rowMax)
    }

    private static func reactionMapTemplate(from molecules: [Molecule]) -> [Int: CGPoint] {
        var pointsByMap: [Int: [CGPoint]] = [:]
        for molecule in molecules {
            guard let box = molecule.boundingBox() else { continue }
            let width = max(0.0001, box.width)
            let height = max(0.0001, box.height)
            for atom in molecule.atoms {
                guard let mapNo = atom.atomMapNumber, mapNo > 0 else { continue }
                let normalized = CGPoint(x: (atom.position.x - box.minX) / width,
                                         y: (atom.position.y - box.minY) / height)
                pointsByMap[mapNo, default: []].append(normalized)
            }
        }

        var out: [Int: CGPoint] = [:]
        out.reserveCapacity(pointsByMap.count)
        for (mapNo, values) in pointsByMap where !values.isEmpty {
            let sx = values.reduce(CGFloat.zero) { partial, point in
                partial + point.x
            }
            let sy = values.reduce(CGFloat.zero) { partial, point in
                partial + point.y
            }
            out[mapNo] = CGPoint(x: sx / CGFloat(values.count),
                                 y: sy / CGFloat(values.count))
        }
        return out
    }

    private static func validatedHighlight(_ selection: CDKReactionParticipantSelection?,
                                           reactantCount: Int,
                                           agentCount: Int,
                                           productCount: Int) -> CDKReactionParticipantSelection? {
        guard let selection, selection.index >= 0 else { return nil }
        switch selection.role {
        case .reactant where selection.index < reactantCount:
            return selection
        case .agent where selection.index < agentCount:
            return selection
        case .product where selection.index < productCount:
            return selection
        default:
            return nil
        }
    }
}

enum CDKReactionParticipantLayoutRefiner {
    private static let lowQualityThreshold: CGFloat = 60
    private static let requiredImprovementFactor: CGFloat = 0.9
    private static let requiredAbsoluteImprovement: CGFloat = 10
    private static let cache = ParticipantLayoutCache(maxEntries: 256)

    static func recomputeIfLowQuality(_ molecule: Molecule) -> Molecule {
        if let cached = cache.value(for: molecule) {
            return cached
        }

        guard shouldRecompute2DLayout(for: molecule) else {
            let oriented = orientForReactionPresentation(molecule)
            cache.store(oriented, for: molecule)
            return oriented
        }

        let originalPenalty = qualityPenalty(for: molecule)
        let regenerated = Depiction2DGenerator.generate(for: molecule)
        let regeneratedPenalty = qualityPenalty(for: regenerated)

        let keepRegenerated =
            regeneratedPenalty <= (originalPenalty * requiredImprovementFactor) ||
            regeneratedPenalty <= (originalPenalty - requiredAbsoluteImprovement)

        let coarse = keepRegenerated ? regenerated : molecule
        let result = orientForReactionPresentation(coarse)
        cache.store(result, for: molecule)
        return result
    }

    static func shouldRecompute2DLayout(for molecule: Molecule) -> Bool {
        guard molecule.atomCount > 1, molecule.bondCount > 0 else {
            return false
        }
        return qualityPenalty(for: molecule) >= lowQualityThreshold
    }

    static func reorientForMapTemplate(_ molecule: Molecule, mapTemplate: [Int: CGPoint]) -> Molecule {
        guard !mapTemplate.isEmpty else { return molecule }
        return orientForReactionPresentation(molecule, mapTemplate: mapTemplate)
    }

    // Backward-compatible shim for tests/callers still providing x-only templates.
    static func reorientForMapTemplate(_ molecule: Molecule, mapXTemplate: [Int: CGFloat]) -> Molecule {
        let mapTemplate = Dictionary(uniqueKeysWithValues: mapXTemplate.map { key, value in
            (key, CGPoint(x: value, y: 0.5))
        })
        return reorientForMapTemplate(molecule, mapTemplate: mapTemplate)
    }

    static func qualityPenalty(for molecule: Molecule) -> CGFloat {
        guard molecule.atomCount > 1 else { return 0 }

        let positionsByAtom = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0.position) })
        let bondLengths = molecule.bonds.compactMap { bond -> CGFloat? in
            guard let p1 = positionsByAtom[bond.a1], let p2 = positionsByAtom[bond.a2] else { return nil }
            let d = p1.distance(to: p2)
            return d.isFinite && d > 0.000_01 ? d : nil
        }

        guard !bondLengths.isEmpty else {
            return 1_000
        }

        let meanBondLength = bondLengths.reduce(0, +) / CGFloat(bondLengths.count)
        if meanBondLength <= 0.000_01 || !meanBondLength.isFinite {
            return 1_000
        }

        var penalty: CGFloat = 0

        let variance = bondLengths.reduce(0) { partial, value in
            let delta = value - meanBondLength
            return partial + (delta * delta)
        } / CGFloat(bondLengths.count)
        let stdDev = sqrt(max(0, variance))
        let coefficientOfVariation = stdDev / max(0.000_1, meanBondLength)
        penalty += max(0, coefficientOfVariation - 0.22) * 300

        if let minBond = bondLengths.min(), minBond < meanBondLength * 0.35 {
            penalty += (1 - (minBond / (meanBondLength * 0.35))) * 140
        }
        if let maxBond = bondLengths.max(), maxBond > meanBondLength * 2.8 {
            penalty += ((maxBond / (meanBondLength * 2.8)) - 1) * 120
        }

        if let box = molecule.boundingBox() {
            let minSide = min(box.width, box.height)
            if minSide < meanBondLength * 0.5 {
                penalty += ((meanBondLength * 0.5) - minSide) / max(0.001, meanBondLength) * 120
            }

            let aspect = max(box.width, box.height) / max(0.0001, minSide)
            if aspect > 8 {
                penalty += (aspect - 8) * 20
            }
        } else {
            penalty += 200
        }

        let overlapTolerance = max(0.05, meanBondLength * 0.06)
        var uniquePoints = Set<QuantizedPoint>()
        for atom in molecule.atoms {
            uniquePoints.insert(QuantizedPoint(point: atom.position, tolerance: overlapTolerance))
        }
        let duplicateCount = molecule.atomCount - uniquePoints.count
        if duplicateCount > 0 {
            penalty += CGFloat(duplicateCount) * 60
        }

        var bondedPairs = Set<BondedPair>()
        bondedPairs.reserveCapacity(molecule.bondCount)
        for bond in molecule.bonds {
            bondedPairs.insert(BondedPair(a: bond.a1, b: bond.a2))
        }

        let atoms = molecule.atoms
        let severeOverlapDistance = meanBondLength * 0.12
        let closeContactDistance = meanBondLength * 0.24

        for i in atoms.indices {
            let a = atoms[i]
            for j in atoms.indices where j > i {
                let b = atoms[j]
                if bondedPairs.contains(BondedPair(a: a.id, b: b.id)) {
                    continue
                }

                let distance = a.position.distance(to: b.position)
                if distance < severeOverlapDistance {
                    penalty += 120
                } else if distance < closeContactDistance {
                    penalty += ((closeContactDistance - distance) / max(0.001, closeContactDistance)) * 40
                }
            }
        }

        return penalty
    }

    private static func orientForReactionPresentation(_ molecule: Molecule,
                                                      mapTemplate: [Int: CGPoint]? = nil) -> Molecule {
        guard molecule.atomCount > 1, let baseBox = molecule.boundingBox() else {
            return molecule
        }

        let mappedPairs: [(map: Int, id: Int)] = molecule.atoms.compactMap { atom in
            guard let map = atom.atomMapNumber, map > 0 else { return nil }
            return (map: map, id: atom.id)
        }
        let hasUsefulMaps = mappedPairs.count >= 3
        let hasTemplate = !(mapTemplate?.isEmpty ?? true)
        let center = CGPoint(x: baseBox.midX, y: baseBox.midY)

        struct OrientationCandidate {
            let rotationRadians: CGFloat
            let mirrorX: Bool
        }

        var candidates: [OrientationCandidate] = [
            OrientationCandidate(rotationRadians: 0, mirrorX: false),
            OrientationCandidate(rotationRadians: .pi / 2, mirrorX: false),
            OrientationCandidate(rotationRadians: -.pi / 2, mirrorX: false),
            OrientationCandidate(rotationRadians: .pi, mirrorX: false)
        ]
        if hasUsefulMaps || hasTemplate {
            candidates.append(contentsOf: [
                OrientationCandidate(rotationRadians: 0, mirrorX: true),
                OrientationCandidate(rotationRadians: .pi / 2, mirrorX: true),
                OrientationCandidate(rotationRadians: -.pi / 2, mirrorX: true),
                OrientationCandidate(rotationRadians: .pi, mirrorX: true)
            ])
        }

        if let axisAngle = principalAxisAngle(of: molecule) {
            // CDK-like readability: additionally evaluate orientations that make the
            // dominant participant axis horizontal, including ordinary unmapped reactions.
            let align = -axisAngle
            let opposite = align + .pi
            candidates.append(OrientationCandidate(rotationRadians: align, mirrorX: false))
            candidates.append(OrientationCandidate(rotationRadians: opposite, mirrorX: false))
            if hasUsefulMaps || hasTemplate {
                candidates.append(OrientationCandidate(rotationRadians: align, mirrorX: true))
                candidates.append(OrientationCandidate(rotationRadians: opposite, mirrorX: true))
            }
        }

        var best = molecule
        var bestScore = scoreOrientation(molecule,
                                         mappedPairs: mappedPairs,
                                         hasUsefulMaps: hasUsefulMaps,
                                         mapTemplate: mapTemplate)

        for candidate in candidates {
            let transformed = applyOrientationTransform(molecule,
                                                        around: center,
                                                        rotationRadians: candidate.rotationRadians,
                                                        mirrorX: candidate.mirrorX)
            let score = scoreOrientation(transformed,
                                         mappedPairs: mappedPairs,
                                         hasUsefulMaps: hasUsefulMaps,
                                         mapTemplate: mapTemplate)
            if score > bestScore {
                bestScore = score
                best = transformed
            }
        }

        return best
    }

    private static func applyOrientationTransform(_ molecule: Molecule,
                                                  around center: CGPoint,
                                                  rotationRadians: CGFloat,
                                                  mirrorX: Bool) -> Molecule {
        let c = cos(rotationRadians)
        let s = sin(rotationRadians)

        var out = molecule
        out.atoms = molecule.atoms.map { atom in
            var updated = atom
            var x = atom.position.x - center.x
            let y = atom.position.y - center.y

            if mirrorX {
                x = -x
            }

            let rx = (x * c) - (y * s)
            let ry = (x * s) + (y * c)
            updated.position = CGPoint(x: rx + center.x, y: ry + center.y)
            return updated
        }
        return out
    }

    private static func scoreOrientation(_ molecule: Molecule,
                                         mappedPairs: [(map: Int, id: Int)],
                                         hasUsefulMaps: Bool,
                                         mapTemplate: [Int: CGPoint]?) -> CGFloat {
        guard let box = molecule.boundingBox() else {
            return -.greatestFiniteMagnitude
        }

        let width = max(0.0001, box.width)
        let height = max(0.0001, box.height)
        let horizontalPreference = (width / height).clamped(to: 0.35...6.5)

        // Prefer compact participants after orientation transforms.
        let penalty = qualityPenalty(for: molecule)
        var score = (horizontalPreference * 7.5) - (penalty * 0.018)

        if hasUsefulMaps {
            var maps: [CGFloat] = []
            var xs: [CGFloat] = []
            for pair in mappedPairs {
                guard let atom = molecule.atom(id: pair.id) else { continue }
                maps.append(CGFloat(pair.map))
                xs.append(atom.position.x)
            }
            if let correlation = pearsonCorrelation(xs: xs, ys: maps) {
                // Positive map-number progression left->right is easier to read in reaction depictions.
                score += (correlation * 12.0)
            }
            if let templateScore = mapTemplateAlignmentScore(molecule,
                                                             mappedPairs: mappedPairs,
                                                             mapTemplate: mapTemplate) {
                // Prefer reactant/agent orientation that aligns with product-side mapped atom progression.
                // X carries stronger weight for reaction left-to-right readability; Y helps
                // avoid mirrored/rotated placements that look geometrically "off".
                score += templateScore
            }
        }

        return score
    }

    private static func pearsonCorrelation(xs: [CGFloat], ys: [CGFloat]) -> CGFloat? {
        guard xs.count == ys.count, xs.count >= 3 else { return nil }
        let n = CGFloat(xs.count)
        let meanX = xs.reduce(0, +) / n
        let meanY = ys.reduce(0, +) / n

        var covariance: CGFloat = 0
        var varianceX: CGFloat = 0
        var varianceY: CGFloat = 0
        for idx in xs.indices {
            let dx = xs[idx] - meanX
            let dy = ys[idx] - meanY
            covariance += dx * dy
            varianceX += dx * dx
            varianceY += dy * dy
        }

        let denom = sqrt(max(0.000_000_1, varianceX * varianceY))
        guard denom.isFinite, denom > 0 else { return nil }
        let r = covariance / denom
        return r.isFinite ? r.clamped(to: -1...1) : nil
    }

    private static func principalAxisAngle(of molecule: Molecule) -> CGFloat? {
        guard molecule.atomCount >= 2 else { return nil }
        let xs = molecule.atoms.map { $0.position.x }
        let ys = molecule.atoms.map { $0.position.y }
        let n = CGFloat(molecule.atomCount)
        let meanX = xs.reduce(0, +) / n
        let meanY = ys.reduce(0, +) / n

        var sxx: CGFloat = 0
        var syy: CGFloat = 0
        var sxy: CGFloat = 0
        for atom in molecule.atoms {
            let dx = atom.position.x - meanX
            let dy = atom.position.y - meanY
            sxx += dx * dx
            syy += dy * dy
            sxy += dx * dy
        }

        if abs(sxy) < 0.000_001 && abs(sxx - syy) < 0.000_001 {
            return nil
        }
        return 0.5 * atan2(2 * sxy, sxx - syy)
    }

    private static func mapTemplateAlignmentScore(_ molecule: Molecule,
                                                  mappedPairs: [(map: Int, id: Int)],
                                                  mapTemplate: [Int: CGPoint]?) -> CGFloat? {
        guard let mapTemplate, !mapTemplate.isEmpty else { return nil }

        var xs: [CGFloat] = []
        var ys: [CGFloat] = []
        var xTemplate: [CGFloat] = []
        var yTemplate: [CGFloat] = []
        xs.reserveCapacity(mappedPairs.count)
        ys.reserveCapacity(mappedPairs.count)
        xTemplate.reserveCapacity(mappedPairs.count)
        yTemplate.reserveCapacity(mappedPairs.count)

        guard let box = molecule.boundingBox() else { return nil }
        let width = max(0.0001, box.width)
        let height = max(0.0001, box.height)

        for pair in mappedPairs {
            guard let template = mapTemplate[pair.map], let atom = molecule.atom(id: pair.id) else { continue }
            xs.append((atom.position.x - box.minX) / width)
            ys.append((atom.position.y - box.minY) / height)
            xTemplate.append(template.x)
            yTemplate.append(template.y)
        }

        guard xs.count >= 3 else { return nil }
        let xCorr = pearsonCorrelation(xs: xs, ys: xTemplate) ?? 0
        let yCorr = pearsonCorrelation(xs: ys, ys: yTemplate) ?? 0
        return (xCorr * 16.0) + (yCorr * 8.5)
    }
}

private final class ParticipantLayoutCache: @unchecked Sendable {
    private let lock = NSLock()
    private let maxEntries: Int
    private var storage: [Molecule: Molecule] = [:]
    private var insertionOrder: [Molecule] = []

    init(maxEntries: Int) {
        self.maxEntries = max(8, maxEntries)
    }

    func value(for molecule: Molecule) -> Molecule? {
        lock.lock()
        defer { lock.unlock() }
        return storage[molecule]
    }

    func store(_ value: Molecule, for key: Molecule) {
        lock.lock()
        defer { lock.unlock() }

        if storage[key] == nil {
            insertionOrder.append(key)
        }
        storage[key] = value

        while insertionOrder.count > maxEntries {
            let removed = insertionOrder.removeFirst()
            storage.removeValue(forKey: removed)
        }
    }
}

private struct QuantizedPoint: Hashable {
    let x: Int
    let y: Int

    init(point: CGPoint, tolerance: CGFloat) {
        let t = max(0.000_1, tolerance)
        self.x = Int((point.x / t).rounded())
        self.y = Int((point.y / t).rounded())
    }
}

private struct BondedPair: Hashable {
    let a: Int
    let b: Int

    init(a: Int, b: Int) {
        if a <= b {
            self.a = a
            self.b = b
        } else {
            self.a = b
            self.b = a
        }
    }
}

private extension Comparable {
    func clamped(to range: ClosedRange<Self>) -> Self {
        min(max(self, range.lowerBound), range.upperBound)
    }
}
