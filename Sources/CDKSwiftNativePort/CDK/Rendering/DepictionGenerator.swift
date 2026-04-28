#if canImport(CoreGraphics)
import CoreGraphics
#endif
import Foundation

/// CDK-style 2D depiction export to SVG.
///
/// This is a native Swift port of the depiction/export stage used in this project:
/// layout is assumed to be generated already (StructureDiagramGenerator), then rendered
/// to vector primitives with CDK-like bond and label conventions.
public enum CDKDepictionGenerator {
    public static func toSVG(molecule: Molecule,
                             style: RenderStyle = RenderStyle(),
                             canvasSize: CGSize = CGSize(width: 900, height: 620),
                             includeBackground: Bool = true) -> String {
        let width = max(1, canvasSize.width)
        let height = max(1, canvasSize.height)

        let depictionMolecule = CDKDepictionPreprocessor.prepareForRendering(molecule: molecule, style: style)
        let highlightedAtomIDs = style.highlightStyle == .none ? Set<Int>() : Set(depictionMolecule.highlightedAtomIDs)
        let highlightedBondIDs = style.highlightStyle == .none ? Set<Int>() : Set(depictionMolecule.highlightedBondIDs)
        guard let box = depictionMolecule.boundingBox() else {
            return """
            <svg xmlns="http://www.w3.org/2000/svg" width="\(Int(width))" height="\(Int(height))" viewBox="0 0 \(Int(width)) \(Int(height))"></svg>
            """
        }

        let rect = CGRect(origin: .zero, size: CGSize(width: width, height: height))
        let pad = style.padding
        let available = CGRect(x: rect.minX + pad,
                               y: rect.minY + pad,
                               width: max(1, rect.width - 2 * pad),
                               height: max(1, rect.height - 2 * pad))
        let scaleX = available.width / max(0.0001, box.width)
        let scaleY = available.height / max(0.0001, box.height)
        let scale = min(scaleX, scaleY)

        let center = CGPoint(x: available.midX, y: available.midY)
        let transform = CGAffineTransform.identity
            .translatedBy(x: center.x, y: center.y)
            .scaledBy(x: scale, y: -scale)
            .translatedBy(x: -box.midX, y: -box.midY)

        let positionsByAtomID = Dictionary(uniqueKeysWithValues: depictionMolecule.atoms.map {
            ($0.id, $0.position.applying(transform))
        })
        let atomByID = Dictionary(uniqueKeysWithValues: depictionMolecule.atoms.map { ($0.id, $0) })

        func p(_ atomID: Int) -> CGPoint? {
            positionsByAtomID[atomID]
        }

        var degree: [Int: Int] = [:]
        for b in depictionMolecule.bonds {
            degree[b.a1, default: 0] += 1
            degree[b.a2, default: 0] += 1
        }

        let aromaticRings = depictionMolecule.aromaticDisplayRings()
        let aromaticBondIDs = depictionMolecule.aromaticDisplayBondIDs()
        var aromaticEdgeCenters: [RenderEdgeKey: [CGPoint]] = [:]
        var conjugatedDoubleEdgeCenters: [RenderEdgeKey: [CGPoint]] = [:]

        for ring in aromaticRings where ring.count >= 3 {
            let ringPoints = ring.compactMap(p)
            guard ringPoints.count == ring.count else { continue }
            let c = ringPoints.reduce(CGPoint.zero) { acc, q in
                CGPoint(x: acc.x + q.x / CGFloat(ringPoints.count),
                        y: acc.y + q.y / CGFloat(ringPoints.count))
            }
            for i in 0..<ring.count {
                let a = ring[i]
                let b = ring[(i + 1) % ring.count]
                guard let bond = depictionMolecule.bond(between: a, and: b), aromaticBondIDs.contains(bond.id) else { continue }
                aromaticEdgeCenters[RenderEdgeKey(a, b), default: []].append(c)
            }
        }

        let conjugatedRings = depictionMolecule.simpleCycles(maxSize: 10)
        for ring in conjugatedRings where ring.count >= 5 {
            let ringBonds: [Bond] = (0..<ring.count).compactMap { i in
                depictionMolecule.bond(between: ring[i], and: ring[(i + 1) % ring.count])
            }
            guard ringBonds.count == ring.count else { continue }

            let piCount = ringBonds.reduce(0) { acc, b in
                acc + ((b.order == .double || b.order == .aromatic) ? 1 : 0)
            }
            let hasConjugation = piCount >= 2 && ringBonds.allSatisfy { $0.order != .triple }
            guard hasConjugation else { continue }

            let ringPoints = ring.compactMap(p)
            guard ringPoints.count == ring.count else { continue }
            let c = ringPoints.reduce(CGPoint.zero) { acc, q in
                CGPoint(x: acc.x + q.x / CGFloat(ringPoints.count),
                        y: acc.y + q.y / CGFloat(ringPoints.count))
            }

            for i in 0..<ring.count {
                let a = ring[i]
                let b = ring[(i + 1) % ring.count]
                guard let bond = depictionMolecule.bond(between: a, and: b) else { continue }
                if bond.order == .double || aromaticBondIDs.contains(bond.id) {
                    conjugatedDoubleEdgeCenters[RenderEdgeKey(a, b), default: []].append(c)
                }
            }
        }

        let labelClipPadding = max(1.4, style.bondWidth * 0.55)
        let labelClipObstacles: [CDKLabelObstacle] = depictionMolecule.atoms.compactMap { atom in
            let atomDegree = degree[atom.id] ?? 0
            guard CDKLabelText.shouldDrawLabel(atom: atom,
                                               degree: atomDegree,
                                               style: style,
                                               molecule: depictionMolecule,
                                               highlightedAtomIDs: highlightedAtomIDs,
                                               highlightedBondIDs: highlightedBondIDs),
                  let pos = p(atom.id) else {
                return nil
            }
            let implicitH = style.showImplicitHydrogens ? depictionMolecule.implicitHydrogenCount(for: atom.id) : 0
            let text = CDKLabelText.build(atom: atom, style: style, implicitHydrogenCount: implicitH)
            let centerOffset = CDKLabelText.centerOffset(atom: atom,
                                                         style: style,
                                                         implicitHydrogenCount: implicitH,
                                                         fontSize: style.fontSize)
            let labelCenter = pos.offsetBy(dx: centerOffset.dx, dy: centerOffset.dy)
            return CDKLabelClipping.makeGlyphObstacle(text: text,
                                                      center: labelCenter,
                                                      fontSize: style.fontSize,
                                                      padding: labelClipPadding)
        }

        let markushBoxes = CDKMarkushRendering.rGroupBoxes(molecule: depictionMolecule,
                                                           positionsByAtomID: positionsByAtomID,
                                                           style: style,
                                                           padding: labelClipPadding,
                                                           fontSize: style.fontSize,
                                                           bondWidth: style.bondWidth)
        let linkAnnotations = CDKMarkushRendering.linkAnnotations(molecule: depictionMolecule,
                                                                  positionsByAtomID: positionsByAtomID,
                                                                  fontSize: style.fontSize)
        let sgroupBracketAnnotations = CDKSgroupRendering.bracketAnnotations(molecule: depictionMolecule,
                                                                             positionsByAtomID: positionsByAtomID,
                                                                             transformPoint: { $0.applying(transform) },
                                                                             fontSize: style.fontSize,
                                                                             bondWidth: style.bondWidth)

        var lines: [String] = []
        lines.reserveCapacity(max(32, depictionMolecule.bonds.count * 4 + depictionMolecule.atoms.count))
        lines.append("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"\(Int(width))\" height=\"\(Int(height))\" viewBox=\"0 0 \(Int(width)) \(Int(height))\">")
        if includeBackground {
            lines.append("<rect x=\"0\" y=\"0\" width=\"\(Int(width))\" height=\"\(Int(height))\" fill=\"#ffffff\"/>")
        }

        for annotation in markushBoxes {
            lines.append(
                "<rect x=\"\(fmt(annotation.boxRect.minX))\" y=\"\(fmt(annotation.boxRect.minY))\" width=\"\(fmt(annotation.boxRect.width))\" height=\"\(fmt(annotation.boxRect.height))\" rx=\"\(fmt(annotation.cornerRadius))\" ry=\"\(fmt(annotation.cornerRadius))\" fill=\"\(annotation.fillColor.svgHexRGB())\" fill-opacity=\"\(fmt(annotation.fillColor.alpha))\"/>"
            )
        }

        @inline(__always)
        func addLine(_ a: CGPoint,
                     _ b: CGPoint,
                     width: CGFloat = style.bondWidth,
                     opacity: CGFloat = 0.95,
                     color: String = CDKRenderColor.ink.svgHexRGB(),
                     dashArray: String? = nil) {
            guard let (start, end) = CDKLabelClipping.clipSegmentEndpoints(a,
                                                                            b,
                                                                            labelObstacles: labelClipObstacles,
                                                                            padding: labelClipPadding) else {
                return
            }
            let dashAttribute = dashArray.map { " stroke-dasharray=\"\($0)\"" } ?? ""
            lines.append(
                "<line x1=\"\(fmt(start.x))\" y1=\"\(fmt(start.y))\" x2=\"\(fmt(end.x))\" y2=\"\(fmt(end.y))\" stroke=\"\(color)\" stroke-width=\"\(fmt(width))\" stroke-opacity=\"\(fmt(opacity))\" stroke-linecap=\"round\"\(dashAttribute)/>"
            )
        }

        @inline(__always)
        func addGlowLine(_ a: CGPoint,
                         _ b: CGPoint,
                         width: CGFloat) {
            addLine(a,
                    b,
                    width: max(width * 2.35, width + max(1.6, style.bondWidth * 1.8)),
                    opacity: 0.34,
                    color: CDKRenderColor.outerGlowHighlight.svgHexRGB())
        }

        @inline(__always)
        func addPolygon(points: [CGPoint], opacity: CGFloat = 0.96, color: String = CDKRenderColor.ink.svgHexRGB()) {
            let pointsText = points.map { "\(fmt($0.x)),\(fmt($0.y))" }.joined(separator: " ")
            lines.append(
                "<polygon points=\"\(pointsText)\" fill=\"\(color)\" fill-opacity=\"\(fmt(opacity))\"/>"
            )
        }

        func normalize(_ vec: CGVector) -> CGVector? {
            let l = hypot(vec.dx, vec.dy)
            guard l > 0.0001 else { return nil }
            return CGVector(dx: vec.dx / l, dy: vec.dy / l)
        }

        for bond in depictionMolecule.bonds {
            guard let p1 = p(bond.a1), let p2 = p(bond.a2) else { continue }

            let v = CGVector(dx: p2.x - p1.x, dy: p2.y - p1.y)
            let len = max(0.0001, hypot(v.dx, v.dy))
            let px = -v.dy / len
            let py = v.dx / len
            let edgeKey = RenderEdgeKey(bond.a1, bond.a2)
            let aromaticCenters = aromaticEdgeCenters[edgeKey] ?? []
            let conjugatedCenters = conjugatedDoubleEdgeCenters[edgeKey] ?? []
            let aromaticStyled = bond.order == .aromatic
                || (style.aromaticDisplayMode == .circle && aromaticBondIDs.contains(bond.id))

            let baseColor = CDKRenderingStyleResolver
                .bondColor(for: bond,
                           molecule: depictionMolecule,
                           style: style,
                           highlightedBondIDs: highlightedBondIDs)
                .svgHexRGB()
            let doubleBondSeparation = max(3.4, style.bondWidth * 2.25)
            let doubleBondHalfSeparation = doubleBondSeparation * 0.5
            let tripleBondOffset = max(3.1, style.bondWidth * 2.2)

            func drawAromaticInnerLine(center: CGPoint) {
                let mid = CGPoint(x: (p1.x + p2.x) * 0.5, y: (p1.y + p2.y) * 0.5)
                guard let inward = normalize(CGVector(dx: center.x - mid.x, dy: center.y - mid.y)) else { return }
                let trim = CGFloat(0.16)
                let inset = max(1.6, style.bondWidth * 0.95)
                let a = CGPoint(x: p1.x + (p2.x - p1.x) * trim + inward.dx * inset,
                                y: p1.y + (p2.y - p1.y) * trim + inward.dy * inset)
                let b = CGPoint(x: p2.x + (p1.x - p2.x) * trim + inward.dx * inset,
                                y: p2.y + (p1.y - p2.y) * trim + inward.dy * inset)
                addLine(a, b, width: max(1.0, style.bondWidth * 0.8), opacity: 0.76, color: baseColor)
            }

            func drawConjugatedDoubleInnerLine(centers: [CGPoint]) {
                let mid = CGPoint(x: (p1.x + p2.x) * 0.5, y: (p1.y + p2.y) * 0.5)
                let sum = centers.reduce(CGVector.zero) { acc, c in
                    CGVector(dx: acc.dx + (c.x - mid.x), dy: acc.dy + (c.y - mid.y))
                }
                let preferred = normalize(sum)
                let perp = CGVector(dx: px, dy: py)
                let sign: CGFloat
                if let preferred {
                    sign = (preferred.dx * perp.dx + preferred.dy * perp.dy) >= 0 ? 1 : -1
                } else {
                    sign = (bond.id % 2 == 0) ? 1 : -1
                }

                let trim = CGFloat(0.15)
                let inset = doubleBondSeparation
                let a = CGPoint(
                    x: p1.x + (p2.x - p1.x) * trim + perp.dx * inset * sign,
                    y: p1.y + (p2.y - p1.y) * trim + perp.dy * inset * sign
                )
                let b = CGPoint(
                    x: p2.x + (p1.x - p2.x) * trim + perp.dx * inset * sign,
                    y: p2.y + (p1.y - p2.y) * trim + perp.dy * inset * sign
                )
                addLine(a, b, width: max(1.05, style.bondWidth * 0.85), opacity: 0.88, color: baseColor)
            }

            func drawQuerySecondaryLine(centers: [CGPoint], dashed: Bool) {
                let mid = CGPoint(x: (p1.x + p2.x) * 0.5, y: (p1.y + p2.y) * 0.5)
                let sum = centers.reduce(CGVector.zero) { acc, c in
                    CGVector(dx: acc.dx + (c.x - mid.x), dy: acc.dy + (c.y - mid.y))
                }
                let preferred = normalize(sum)
                let perp = CGVector(dx: px, dy: py)
                let sign: CGFloat
                if let preferred {
                    sign = (preferred.dx * perp.dx + preferred.dy * perp.dy) >= 0 ? 1 : -1
                } else {
                    sign = (bond.id % 2 == 0) ? 1 : -1
                }

                let trim = CGFloat(0.15)
                let inset = doubleBondSeparation
                let a = CGPoint(
                    x: p1.x + (p2.x - p1.x) * trim + perp.dx * inset * sign,
                    y: p1.y + (p2.y - p1.y) * trim + perp.dy * inset * sign
                )
                let b = CGPoint(
                    x: p2.x + (p1.x - p2.x) * trim + perp.dx * inset * sign,
                    y: p2.y + (p1.y - p2.y) * trim + perp.dy * inset * sign
                )
                addLine(a,
                        b,
                        width: max(1.0, style.bondWidth * 0.82),
                        opacity: 0.92,
                        color: baseColor,
                        dashArray: dashed ? "\(fmt(max(style.bondWidth * 3.6, 8.0))) \(fmt(max(style.bondWidth * 1.8, 4.8)))" : nil)
            }

            func drawAttachmentPointGlyph(center: CGPoint, other: CGPoint) {
                for segment in CDKMarkushRendering.attachmentPointSegments(center: center,
                                                                           other: other,
                                                                           style: style) {
                    addLine(segment.from,
                            segment.to,
                            width: max(1.0, style.bondWidth * 0.92),
                            opacity: 0.95,
                            color: baseColor)
                }
            }

            func drawSolidWedge(from a: CGPoint, to b: CGPoint) {
                guard let (start, end) = CDKLabelClipping.clipSegmentEndpoints(a,
                                                                                b,
                                                                                labelObstacles: labelClipObstacles,
                                                                                padding: labelClipPadding) else {
                    return
                }
                let w = max(7.2, style.bondWidth * 3.9)
                let edge1 = start.offsetBy(dx: px * w, dy: py * w)
                let edge2 = start.offsetBy(dx: -px * w, dy: -py * w)
                addPolygon(points: [end, edge1, edge2], opacity: 0.96, color: baseColor)
            }

            func drawHashedWedge(from a: CGPoint, to b: CGPoint) {
                guard let (start, end) = CDKLabelClipping.clipSegmentEndpoints(a,
                                                                                b,
                                                                                labelObstacles: labelClipObstacles,
                                                                                padding: labelClipPadding) else {
                    return
                }
                let segments = 9
                let w = max(7.6, style.bondWidth * 4.15)
                for i in 1...segments {
                    let t = CGFloat(i) / CGFloat(segments + 1)
                    let c = CGPoint(x: start.x + (end.x - start.x) * t, y: start.y + (end.y - start.y) * t)
                    let hw = w * (1 - t)
                    let l = c.offsetBy(dx: px * hw, dy: py * hw)
                    let r = c.offsetBy(dx: -px * hw, dy: -py * hw)
                    addLine(l, r, width: max(1.15, style.bondWidth * 0.86), opacity: 0.95, color: baseColor)
                }
            }

            let glowHighlighted = usesGlowHighlights(style: style) && highlightedBondIDs.contains(bond.id)
            if aromaticStyled {
                if glowHighlighted {
                    addGlowLine(p1, p2, width: max(1.2, style.bondWidth * 0.9))
                }
                addLine(p1, p2, width: max(1.2, style.bondWidth * 0.9), opacity: 0.82, color: baseColor)
                if style.aromaticDisplayMode == .innerLine,
                   aromaticCenters.count == 1,
                   let center = aromaticCenters.first {
                    drawAromaticInnerLine(center: center)
                }
                continue
            }

            if bond.queryType != nil {
                switch bond.queryType {
                case .any:
                    if glowHighlighted {
                        addGlowLine(p1, p2, width: max(1.0, style.bondWidth * 0.92))
                    }
                    addLine(p1,
                            p2,
                            width: max(1.0, style.bondWidth * 0.92),
                            opacity: 0.94,
                            color: baseColor,
                            dashArray: "\(fmt(max(style.bondWidth * 3.6, 8.0))) \(fmt(max(style.bondWidth * 1.8, 4.8)))")
                case .singleOrDouble:
                    if glowHighlighted {
                        addGlowLine(p1, p2, width: style.bondWidth)
                    }
                    addLine(p1, p2, width: style.bondWidth, opacity: 0.95, color: baseColor)
                    drawQuerySecondaryLine(centers: conjugatedCenters.isEmpty ? aromaticCenters : conjugatedCenters, dashed: true)
                case .singleOrAromatic:
                    if glowHighlighted {
                        addGlowLine(p1, p2, width: max(1.0, style.bondWidth * 0.92))
                    }
                    addLine(p1, p2, width: max(1.0, style.bondWidth * 0.92), opacity: 0.90, color: baseColor)
                    drawQuerySecondaryLine(centers: aromaticCenters, dashed: true)
                case .doubleOrAromatic:
                    if glowHighlighted {
                        addGlowLine(p1, p2, width: style.bondWidth)
                    }
                    addLine(p1.offsetBy(dx: px * doubleBondHalfSeparation, dy: py * doubleBondHalfSeparation),
                            p2.offsetBy(dx: px * doubleBondHalfSeparation, dy: py * doubleBondHalfSeparation),
                            width: style.bondWidth,
                            opacity: 0.95,
                            color: baseColor)
                    drawQuerySecondaryLine(centers: aromaticCenters, dashed: true)
                case .none:
                    break
                }
                if atomByID[bond.a1]?.attachmentPoint != nil {
                    drawAttachmentPointGlyph(center: p1, other: p2)
                }
                if atomByID[bond.a2]?.attachmentPoint != nil {
                    drawAttachmentPointGlyph(center: p2, other: p1)
                }
                continue
            }

            if bond.order == .single {
                switch bond.stereo {
                case .up:
                    if glowHighlighted {
                        addGlowLine(p1, p2, width: style.bondWidth)
                    }
                    drawSolidWedge(from: p1, to: p2)
                    continue
                case .down:
                    if glowHighlighted {
                        addGlowLine(p1, p2, width: style.bondWidth)
                    }
                    drawHashedWedge(from: p1, to: p2)
                    continue
                case .upReversed:
                    if glowHighlighted {
                        addGlowLine(p1, p2, width: style.bondWidth)
                    }
                    drawSolidWedge(from: p2, to: p1)
                    continue
                case .downReversed:
                    if glowHighlighted {
                        addGlowLine(p1, p2, width: style.bondWidth)
                    }
                    drawHashedWedge(from: p2, to: p1)
                    continue
                case .either:
                    if glowHighlighted {
                        addGlowLine(p1, p2, width: max(1.0, style.bondWidth * 0.75))
                    }
                    addLine(p1, p2, width: max(1.0, style.bondWidth * 0.75), opacity: 0.78, color: baseColor)
                    continue
                case .none:
                    break
                }
            }

            switch bond.order {
            case .single:
                if glowHighlighted {
                    addGlowLine(p1, p2, width: style.bondWidth)
                }
                addLine(p1, p2, width: style.bondWidth, opacity: 0.95, color: baseColor)
            case .double:
                if glowHighlighted {
                    addGlowLine(p1, p2, width: style.bondWidth)
                }
                if !conjugatedCenters.isEmpty {
                    addLine(p1, p2, width: style.bondWidth, opacity: 0.95, color: baseColor)
                    drawConjugatedDoubleInnerLine(centers: conjugatedCenters)
                } else {
                    addLine(p1.offsetBy(dx: px * doubleBondHalfSeparation, dy: py * doubleBondHalfSeparation),
                            p2.offsetBy(dx: px * doubleBondHalfSeparation, dy: py * doubleBondHalfSeparation),
                            width: style.bondWidth,
                            opacity: 0.95,
                            color: baseColor)
                    addLine(p1.offsetBy(dx: -px * doubleBondHalfSeparation, dy: -py * doubleBondHalfSeparation),
                            p2.offsetBy(dx: -px * doubleBondHalfSeparation, dy: -py * doubleBondHalfSeparation),
                            width: style.bondWidth,
                            opacity: 0.95,
                            color: baseColor)
                }
            case .triple:
                if glowHighlighted {
                    addGlowLine(p1, p2, width: style.bondWidth)
                }
                addLine(p1, p2, width: style.bondWidth, opacity: 0.95, color: baseColor)
                addLine(p1.offsetBy(dx: px * tripleBondOffset, dy: py * tripleBondOffset),
                        p2.offsetBy(dx: px * tripleBondOffset, dy: py * tripleBondOffset),
                        width: style.bondWidth,
                        opacity: 0.95,
                        color: baseColor)
                addLine(p1.offsetBy(dx: -px * tripleBondOffset, dy: -py * tripleBondOffset),
                        p2.offsetBy(dx: -px * tripleBondOffset, dy: -py * tripleBondOffset),
                        width: style.bondWidth,
                        opacity: 0.95,
                        color: baseColor)
            case .aromatic:
                if glowHighlighted {
                    addGlowLine(p1, p2, width: max(1.2, style.bondWidth * 0.85))
                }
                addLine(p1, p2, width: max(1.2, style.bondWidth * 0.85), opacity: 0.8, color: baseColor)
            }

            if atomByID[bond.a1]?.attachmentPoint != nil {
                drawAttachmentPointGlyph(center: p1, other: p2)
            }
            if atomByID[bond.a2]?.attachmentPoint != nil {
                drawAttachmentPointGlyph(center: p2, other: p1)
            }
        }

        for annotation in sgroupBracketAnnotations {
            for segment in annotation.segments {
                addLine(segment.from,
                        segment.to,
                        width: max(0.95, style.bondWidth * 0.84),
                        opacity: 0.92,
                        color: CDKRenderColor.ink.svgHexRGB())
            }
        }

        if style.aromaticDisplayMode == .circle {
            for ring in aromaticRings where ring.count >= 5 {
                let ringPoints = ring.compactMap(p)
                guard ringPoints.count == ring.count else { continue }
                let center = ringPoints.reduce(CGPoint.zero) { acc, q in
                    CGPoint(x: acc.x + q.x / CGFloat(ringPoints.count),
                            y: acc.y + q.y / CGFloat(ringPoints.count))
                }
                let minRadius = ringPoints.map { hypot($0.x - center.x, $0.y - center.y) }.min() ?? 0
                guard minRadius > 1 else { continue }
                let radius = minRadius * 0.53
                let circleColor = CDKRenderingStyleResolver
                    .aromaticRingColor(atomIDs: ring,
                                       molecule: depictionMolecule,
                                       style: style,
                                       highlightedBondIDs: highlightedBondIDs)
                    .svgHexRGB()
                let ringIsHighlighted = usesGlowHighlights(style: style) && (0..<ring.count).contains { index in
                    let atom1 = ring[index]
                    let atom2 = ring[(index + 1) % ring.count]
                    guard let bond = depictionMolecule.bond(between: atom1, and: atom2) else { return false }
                    return highlightedBondIDs.contains(bond.id)
                }
                if ringIsHighlighted {
                    lines.append(
                        "<circle cx=\"\(fmt(center.x))\" cy=\"\(fmt(center.y))\" r=\"\(fmt(radius))\" fill=\"none\" stroke=\"\(CDKRenderColor.outerGlowHighlight.svgHexRGB())\" stroke-width=\"\(fmt(max(1.0, style.bondWidth * 1.95)))\" stroke-opacity=\"0.26\"/>"
                    )
                }
                lines.append(
                    "<circle cx=\"\(fmt(center.x))\" cy=\"\(fmt(center.y))\" r=\"\(fmt(radius))\" fill=\"none\" stroke=\"\(circleColor)\" stroke-width=\"\(fmt(max(1.0, style.bondWidth * 0.84)))\" stroke-opacity=\"0.82\"/>"
                )
            }
        }

        for atom in depictionMolecule.atoms {
            guard CDKLabelText.shouldDrawLabel(atom: atom,
                                               degree: degree[atom.id] ?? 0,
                                               style: style,
                                               molecule: depictionMolecule,
                                               highlightedAtomIDs: highlightedAtomIDs,
                                               highlightedBondIDs: highlightedBondIDs) else { continue }
            guard let pos = p(atom.id) else { continue }

            let implicitH = style.showImplicitHydrogens ? depictionMolecule.implicitHydrogenCount(for: atom.id) : 0
            let label = CDKLabelText.build(atom: atom, style: style, implicitHydrogenCount: implicitH)
            let centerOffset = CDKLabelText.centerOffset(atom: atom,
                                                         style: style,
                                                         implicitHydrogenCount: implicitH,
                                                         fontSize: style.fontSize)
            let labelCenter = pos.offsetBy(dx: centerOffset.dx, dy: centerOffset.dy)
            let atomColor = CDKRenderingStyleResolver
                .atomColor(for: atom,
                           style: style,
                           highlightedAtomIDs: highlightedAtomIDs)
                .svgHexRGB()
            let atomHighlighted = usesGlowHighlights(style: style) && highlightedAtomIDs.contains(atom.id)
            if atomHighlighted {
                lines.append(
                    svgGlowText(label,
                                position: labelCenter,
                                fontSize: style.fontSize,
                                fill: CDKRenderColor.outerGlowHighlight.svgHexRGB())
                )
            }
            lines.append(
                svgText(label,
                        position: labelCenter,
                        fontSize: style.fontSize,
                        fill: atomColor,
                        outlineWhite: atomHighlighted && preservesBackgroundDuringGlow(style: style))
            )
        }

        for annotation in linkAnnotations {
            lines.append(annotationText("(", position: annotation.leftBracketPosition, fontSize: style.fontSize * 0.86))
            lines.append(annotationText(")", position: annotation.rightBracketPosition, fontSize: style.fontSize * 0.86))
            if let subscriptText = annotation.subscriptText,
               let subscriptPosition = annotation.subscriptPosition {
                lines.append(annotationText(subscriptText, position: subscriptPosition, fontSize: style.fontSize * 0.78))
            }
            if let superscriptText = annotation.superscriptText,
               let superscriptPosition = annotation.superscriptPosition {
                lines.append(annotationText(superscriptText, position: superscriptPosition, fontSize: style.fontSize * 0.72))
            }
        }

        for annotation in sgroupBracketAnnotations {
            if let subscriptText = annotation.subscriptText,
               let subscriptPosition = annotation.subscriptPosition {
                lines.append(annotationText(subscriptText, position: subscriptPosition, fontSize: style.fontSize * 0.80))
            }
            if let superscriptText = annotation.superscriptText,
               let superscriptPosition = annotation.superscriptPosition {
                lines.append(annotationText(superscriptText, position: superscriptPosition, fontSize: style.fontSize * 0.74))
            }
        }

        for annotation in markushBoxes {
            lines.append(annotationText(annotation.label, position: annotation.labelPosition, fontSize: style.fontSize))
        }

        lines.append("</svg>")
        return lines.joined(separator: "\n")
    }

    @inline(__always)
    private static func fmt(_ x: CGFloat) -> String {
        String(format: "%.2f", Double(x))
    }

    private static func svgEsc(_ s: String) -> String {
        s.replacingOccurrences(of: "&", with: "&amp;")
            .replacingOccurrences(of: "<", with: "&lt;")
            .replacingOccurrences(of: ">", with: "&gt;")
            .replacingOccurrences(of: "\"", with: "&quot;")
            .replacingOccurrences(of: "'", with: "&apos;")
    }

    private static func annotationText(_ text: String, position: CGPoint, fontSize: CGFloat) -> String {
        svgText(text,
                position: position,
                fontSize: fontSize,
                fill: CDKRenderColor.ink.svgHexRGB())
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

    private static func svgGlowText(_ text: String,
                                    position: CGPoint,
                                    fontSize: CGFloat,
                                    fill: String) -> String {
        let italicized = CDKLabelText.usesItalicRGroupFont(text: text)
        let fontFamily = italicized
            ? "Times New Roman, Georgia, serif"
            : "Helvetica, Arial, sans-serif"
        let fontStyle = italicized ? " font-style=\"italic\"" : ""
        let glowStrokeWidth = max(1.2, fontSize * 0.16)
        return "<text x=\"\(fmt(position.x))\" y=\"\(fmt(position.y))\" text-anchor=\"middle\" dominant-baseline=\"middle\" font-family=\"\(fontFamily)\" font-size=\"\(fmt(fontSize))\" font-weight=\"600\"\(fontStyle) fill=\"\(fill)\" fill-opacity=\"0.78\" stroke=\"\(fill)\" stroke-opacity=\"0.58\" stroke-width=\"\(fmt(glowStrokeWidth))\" paint-order=\"stroke fill\">\(svgEsc(text))</text>"
    }

    private static func svgText(_ text: String,
                                position: CGPoint,
                                fontSize: CGFloat,
                                fill: String,
                                outlineWhite: Bool = false) -> String {
        let italicized = CDKLabelText.usesItalicRGroupFont(text: text)
        let fontFamily = italicized
            ? "Times New Roman, Georgia, serif"
            : "Helvetica, Arial, sans-serif"
        let fontStyle = italicized ? " font-style=\"italic\"" : ""
        let outlineAttributes: String
        if outlineWhite {
            outlineAttributes = " stroke=\"#FFFFFF\" stroke-width=\"\(fmt(max(0.9, fontSize * 0.12)))\" paint-order=\"stroke fill\""
        } else {
            outlineAttributes = ""
        }
        return "<text x=\"\(fmt(position.x))\" y=\"\(fmt(position.y))\" text-anchor=\"middle\" dominant-baseline=\"middle\" font-family=\"\(fontFamily)\" font-size=\"\(fmt(fontSize))\" font-weight=\"600\"\(fontStyle) fill=\"\(fill)\"\(outlineAttributes)>\(svgEsc(text))</text>"
    }

    private struct RenderEdgeKey: Hashable {
        let a: Int
        let b: Int

        init(_ u: Int, _ v: Int) {
            a = min(u, v)
            b = max(u, v)
        }
    }
}
