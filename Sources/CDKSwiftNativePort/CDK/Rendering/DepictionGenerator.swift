import CoreGraphics
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

        func p(_ atomID: Int) -> CGPoint? {
            depictionMolecule.atoms.first(where: { $0.id == atomID }).map { $0.position.applying(transform) }
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
            guard CDKLabelText.shouldDrawLabel(atom: atom, degree: atomDegree, style: style),
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

        var lines: [String] = []
        lines.reserveCapacity(max(32, depictionMolecule.bonds.count * 4 + depictionMolecule.atoms.count))
        lines.append("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"\(Int(width))\" height=\"\(Int(height))\" viewBox=\"0 0 \(Int(width)) \(Int(height))\">")
        if includeBackground {
            lines.append("<rect x=\"0\" y=\"0\" width=\"\(Int(width))\" height=\"\(Int(height))\" fill=\"#ffffff\"/>")
        }

        @inline(__always)
        func addLine(_ a: CGPoint,
                     _ b: CGPoint,
                     width: CGFloat = style.bondWidth,
                     opacity: CGFloat = 0.95,
                     color: String = CDKRenderColor.ink.svgHexRGB()) {
            guard let (start, end) = CDKLabelClipping.clipSegmentEndpoints(a,
                                                                            b,
                                                                            labelObstacles: labelClipObstacles,
                                                                            padding: labelClipPadding) else {
                return
            }
            lines.append(
                "<line x1=\"\(fmt(start.x))\" y1=\"\(fmt(start.y))\" x2=\"\(fmt(end.x))\" y2=\"\(fmt(end.y))\" stroke=\"\(color)\" stroke-width=\"\(fmt(width))\" stroke-opacity=\"\(fmt(opacity))\" stroke-linecap=\"round\"/>"
            )
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

            let baseColor = CDKRenderingStyleResolver
                .bondColor(for: bond, molecule: depictionMolecule, style: style)
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

            func drawSolidWedge(from a: CGPoint, to b: CGPoint) {
                guard let (start, end) = CDKLabelClipping.clipSegmentEndpoints(a,
                                                                                b,
                                                                                labelObstacles: labelClipObstacles,
                                                                                padding: labelClipPadding) else {
                    return
                }
                let w = max(7.2, style.bondWidth * 3.9)
                let edge1 = end.offsetBy(dx: px * w, dy: py * w)
                let edge2 = end.offsetBy(dx: -px * w, dy: -py * w)
                addPolygon(points: [start, edge1, edge2], opacity: 0.96, color: baseColor)
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
                    let hw = w * t
                    let l = c.offsetBy(dx: px * hw, dy: py * hw)
                    let r = c.offsetBy(dx: -px * hw, dy: -py * hw)
                    addLine(l, r, width: max(1.15, style.bondWidth * 0.86), opacity: 0.95, color: baseColor)
                }
            }

            if bond.order == .aromatic {
                addLine(p1, p2, width: max(1.2, style.bondWidth * 0.9), opacity: 0.82, color: baseColor)
                if style.aromaticDisplayMode == .innerLine,
                   aromaticCenters.count == 1,
                   let center = aromaticCenters.first {
                    drawAromaticInnerLine(center: center)
                }
                continue
            }

            if bond.order == .single {
                switch bond.stereo {
                case .up:
                    drawSolidWedge(from: p1, to: p2)
                    continue
                case .down:
                    drawHashedWedge(from: p1, to: p2)
                    continue
                case .upReversed:
                    drawSolidWedge(from: p2, to: p1)
                    continue
                case .downReversed:
                    drawHashedWedge(from: p2, to: p1)
                    continue
                case .either:
                    addLine(p1, p2, width: max(1.0, style.bondWidth * 0.75), opacity: 0.78, color: baseColor)
                    continue
                case .none:
                    break
                }
            }

            switch bond.order {
            case .single:
                addLine(p1, p2, width: style.bondWidth, opacity: 0.95, color: baseColor)
            case .double:
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
                addLine(p1, p2, width: max(1.2, style.bondWidth * 0.85), opacity: 0.8, color: baseColor)
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
                    .aromaticRingColor(atomIDs: ring, molecule: depictionMolecule, style: style)
                    .svgHexRGB()
                lines.append(
                    "<circle cx=\"\(fmt(center.x))\" cy=\"\(fmt(center.y))\" r=\"\(fmt(radius))\" fill=\"none\" stroke=\"\(circleColor)\" stroke-width=\"\(fmt(max(1.0, style.bondWidth * 0.84)))\" stroke-opacity=\"0.82\"/>"
                )
            }
        }

        for atom in depictionMolecule.atoms {
            guard CDKLabelText.shouldDrawLabel(atom: atom, degree: degree[atom.id] ?? 0, style: style) else { continue }
            guard let pos = p(atom.id) else { continue }

            let implicitH = style.showImplicitHydrogens ? depictionMolecule.implicitHydrogenCount(for: atom.id) : 0
            let label = CDKLabelText.build(atom: atom, style: style, implicitHydrogenCount: implicitH)
            let centerOffset = CDKLabelText.centerOffset(atom: atom,
                                                         style: style,
                                                         implicitHydrogenCount: implicitH,
                                                         fontSize: style.fontSize)
            let labelCenter = pos.offsetBy(dx: centerOffset.dx, dy: centerOffset.dy)
            let atomColor = CDKRenderingStyleResolver
                .atomColor(for: atom, style: style)
                .svgHexRGB()
            lines.append(
                "<text x=\"\(fmt(labelCenter.x))\" y=\"\(fmt(labelCenter.y))\" text-anchor=\"middle\" dominant-baseline=\"middle\" font-family=\"Helvetica, Arial, sans-serif\" font-size=\"\(fmt(style.fontSize))\" font-weight=\"700\" fill=\"\(atomColor)\">\(svgEsc(label))</text>"
            )
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

    private struct RenderEdgeKey: Hashable {
        let a: Int
        let b: Int

        init(_ u: Int, _ v: Int) {
            a = min(u, v)
            b = max(u, v)
        }
    }
}
