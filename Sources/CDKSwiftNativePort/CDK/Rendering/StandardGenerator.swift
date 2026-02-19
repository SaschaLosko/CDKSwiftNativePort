import SwiftUI

public struct RenderStyle: Hashable {
    public var showCarbons: Bool = false
    public var showImplicitHydrogens: Bool = true
    public var showExplicitHydrogens: Bool = false
    public var showAtomIDs: Bool = false
    public var atomColoringMode: CDKAtomColoringMode = .monochrome
    public var colorBondsByAtom: Bool = false
    public var aromaticDisplayMode: CDKAromaticDisplayMode = .innerLine
    public var bondWidth: CGFloat = 2.35
    public var fontSize: CGFloat = 28.0
    public var padding: CGFloat = 24.0

    public init(showCarbons: Bool = false,
                showImplicitHydrogens: Bool = true,
                showExplicitHydrogens: Bool = false,
                showAtomIDs: Bool = false,
                atomColoringMode: CDKAtomColoringMode = .monochrome,
                colorBondsByAtom: Bool = false,
                aromaticDisplayMode: CDKAromaticDisplayMode = .innerLine,
                bondWidth: CGFloat = 2.35,
                fontSize: CGFloat = 28.0,
                padding: CGFloat = 24.0) {
        self.showCarbons = showCarbons
        self.showImplicitHydrogens = showImplicitHydrogens
        self.showExplicitHydrogens = showExplicitHydrogens
        self.showAtomIDs = showAtomIDs
        self.atomColoringMode = atomColoringMode
        self.colorBondsByAtom = colorBondsByAtom
        self.aromaticDisplayMode = aromaticDisplayMode
        self.bondWidth = bondWidth
        self.fontSize = fontSize
        self.padding = padding
    }
}

public enum CDKStandardGenerator {

    private struct RenderEdgeKey: Hashable {
        let a: Int
        let b: Int
        init(_ u: Int, _ v: Int) {
            a = min(u, v)
            b = max(u, v)
        }
    }

    public static func draw(molecule: Molecule, in rect: CGRect, style: RenderStyle, context: inout GraphicsContext) {
        let depictionMolecule = CDKDepictionPreprocessor.prepareForRendering(molecule: molecule, style: style)
        guard let box = depictionMolecule.boundingBox() else { return }

        let pad = style.padding
        let available = CGRect(x: rect.minX + pad, y: rect.minY + pad, width: max(1, rect.width - 2*pad), height: max(1, rect.height - 2*pad))

        // Determine transform: center molecule in available rect, fit scale, and flip Y (chem coords usually Y-up).
        let scaleX = available.width / max(0.0001, box.width)
        let scaleY = available.height / max(0.0001, box.height)
        let scale = min(scaleX, scaleY)

        let center = CGPoint(x: available.midX, y: available.midY)

        // Order matters: with CoreGraphics concatenation semantics, this sequence
        // maps molecule coordinates into view coordinates correctly.
        let t = CGAffineTransform.identity
            .translatedBy(x: center.x, y: center.y)   // Move origin to canvas center
            .scaledBy(x: scale, y: -scale)            // Scale + flip Y to screen coords
            .translatedBy(x: -box.midX, y: -box.midY) // Center molecule at origin

        func p(_ atomID: Int) -> CGPoint? {
            depictionMolecule.atoms.first(where: { $0.id == atomID }).map { $0.position.applying(t) }
        }

        // Precompute degrees
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
            let center = ringPoints.reduce(CGPoint.zero) { acc, p in
                CGPoint(x: acc.x + p.x / CGFloat(ringPoints.count),
                        y: acc.y + p.y / CGFloat(ringPoints.count))
            }
            for i in 0..<ring.count {
                let a = ring[i]
                let b = ring[(i + 1) % ring.count]
                guard let bond = depictionMolecule.bond(between: a, and: b), aromaticBondIDs.contains(bond.id) else { continue }
                aromaticEdgeCenters[RenderEdgeKey(a, b), default: []].append(center)
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
            let center = ringPoints.reduce(CGPoint.zero) { acc, q in
                CGPoint(x: acc.x + q.x / CGFloat(ringPoints.count),
                        y: acc.y + q.y / CGFloat(ringPoints.count))
            }
            for i in 0..<ring.count {
                let a = ring[i]
                let b = ring[(i + 1) % ring.count]
                guard let bond = depictionMolecule.bond(between: a, and: b) else { continue }
                if bond.order == .double || aromaticBondIDs.contains(bond.id) {
                    conjugatedDoubleEdgeCenters[RenderEdgeKey(a, b), default: []].append(center)
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

        // Bonds
        for bond in depictionMolecule.bonds {
            guard let p1 = p(bond.a1), let p2 = p(bond.a2) else { continue }
            let edgeKey = RenderEdgeKey(bond.a1, bond.a2)
            let aromaticCenters = aromaticEdgeCenters[edgeKey] ?? []
            let conjugatedCenters = conjugatedDoubleEdgeCenters[edgeKey] ?? []
            let bondColor = CDKRenderingStyleResolver.bondColor(for: bond, molecule: depictionMolecule, style: style)
            drawBond(from: p1,
                     to: p2,
                     bond: bond,
                     bondColor: bondColor,
                     aromaticStyled: bond.order == .aromatic,
                     aromaticRingCenters: aromaticCenters,
                     conjugatedRingCenters: conjugatedCenters,
                     labelClipObstacles: labelClipObstacles,
                     labelClipPadding: labelClipPadding,
                     style: style,
                     context: &context)
        }

        if style.aromaticDisplayMode == .circle {
            for ring in aromaticRings where ring.count >= 5 {
                let ringPoints = ring.compactMap(p)
                guard ringPoints.count == ring.count else { continue }
                let ringColor = CDKRenderingStyleResolver.aromaticRingColor(atomIDs: ring, molecule: depictionMolecule, style: style)
                drawAromaticCircle(points: ringPoints,
                                   color: ringColor,
                                   bondWidth: style.bondWidth,
                                   context: &context)
            }
        }

        // Atoms (labels)
        for atom in depictionMolecule.atoms {
            guard CDKLabelText.shouldDrawLabel(atom: atom, degree: degree[atom.id] ?? 0, style: style) else { continue }
            let pos = atom.position.applying(t)
            let implicitH = style.showImplicitHydrogens ? depictionMolecule.implicitHydrogenCount(for: atom.id) : 0
            let atomColor = CDKRenderingStyleResolver.atomColor(for: atom, style: style)
            drawAtomLabel(atom: atom,
                          at: pos,
                          color: atomColor,
                          style: style,
                          context: &context,
                          degree: degree[atom.id] ?? 0,
                          implicitHydrogenCount: implicitH)
        }
    }

    private static func drawBond(from p1: CGPoint,
                                 to p2: CGPoint,
                                 bond: Bond,
                                 bondColor: CDKRenderColor,
                                 aromaticStyled: Bool,
                                 aromaticRingCenters: [CGPoint],
                                 conjugatedRingCenters: [CGPoint],
                                 labelClipObstacles: [CDKLabelObstacle],
                                 labelClipPadding: CGFloat,
                                 style: RenderStyle,
                                 context: inout GraphicsContext) {
        let v = CGVector(dx: p2.x - p1.x, dy: p2.y - p1.y)
        let len = max(0.0001, hypot(v.dx, v.dy))
        let ny = v.dy / len

        // Perpendicular for multiple lines
        let px = -ny
        let py = v.dx / len

        let doubleBondSeparation = max(3.4, style.bondWidth * 2.25)
        let doubleBondHalfSeparation = doubleBondSeparation * 0.5
        let tripleBondOffset = max(3.1, style.bondWidth * 2.2)

        func normalize(_ vec: CGVector) -> CGVector? {
            let l = hypot(vec.dx, vec.dy)
            guard l > 0.0001 else { return nil }
            return CGVector(dx: vec.dx / l, dy: vec.dy / l)
        }

        func strokeLine(_ a: CGPoint, _ b: CGPoint, width: CGFloat = style.bondWidth, opacity: CGFloat = 0.95) {
            guard let (start, end) = CDKLabelClipping.clipSegmentEndpoints(a,
                                                                            b,
                                                                            labelObstacles: labelClipObstacles,
                                                                            padding: labelClipPadding) else {
                return
            }
            var path = Path()
            path.move(to: start)
            path.addLine(to: end)
            context.stroke(path,
                           with: .color(swiftUIColor(from: bondColor.withAlpha(opacity))),
                           lineWidth: width)
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
            var path = Path()
            path.move(to: start)
            path.addLine(to: edge1)
            path.addLine(to: edge2)
            path.closeSubpath()
            context.fill(path, with: .color(swiftUIColor(from: bondColor.withAlpha(0.96))))
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
                strokeLine(l, r, width: max(1.15, style.bondWidth * 0.86))
            }
        }

        func drawAromaticInnerLine(center: CGPoint) {
            let mid = CGPoint(x: (p1.x + p2.x) * 0.5, y: (p1.y + p2.y) * 0.5)
            guard let inward = normalize(CGVector(dx: center.x - mid.x, dy: center.y - mid.y)) else { return }
            let trim = CGFloat(0.16)
            let inset = max(1.6, style.bondWidth * 0.95)
            let a = CGPoint(x: p1.x + (p2.x - p1.x) * trim + inward.dx * inset,
                            y: p1.y + (p2.y - p1.y) * trim + inward.dy * inset)
            let b = CGPoint(x: p2.x + (p1.x - p2.x) * trim + inward.dx * inset,
                            y: p2.y + (p1.y - p2.y) * trim + inward.dy * inset)
            strokeLine(a, b, width: max(1.0, style.bondWidth * 0.8), opacity: 0.76)
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
            strokeLine(a, b, width: max(1.05, style.bondWidth * 0.85), opacity: 0.88)
        }

        if aromaticStyled {
            strokeLine(p1, p2, width: max(1.2, style.bondWidth * 0.9), opacity: 0.82)
            if style.aromaticDisplayMode == .innerLine {
                // CDK-like aromatic rendering: inner line towards ring center.
                // For fused/shared aromatic edges, avoid overdraw by skipping inner line.
                if aromaticRingCenters.count == 1, let center = aromaticRingCenters.first {
                    drawAromaticInnerLine(center: center)
                }
            }
            return
        }

        if bond.order == .single {
            switch bond.stereo {
            case .up:
                drawSolidWedge(from: p1, to: p2)
                return
            case .down:
                drawHashedWedge(from: p1, to: p2)
                return
            case .upReversed:
                drawSolidWedge(from: p2, to: p1)
                return
            case .downReversed:
                drawHashedWedge(from: p2, to: p1)
                return
            case .either:
                strokeLine(p1, p2, width: max(1.0, style.bondWidth * 0.75), opacity: 0.78)
                return
            case .none:
                break
            }
        }

        switch bond.order {
        case .single:
            strokeLine(p1, p2)
        case .double:
            if !conjugatedRingCenters.isEmpty {
                // CDK-like style: primary bond line + one inward secondary line.
                strokeLine(p1, p2)
                drawConjugatedDoubleInnerLine(centers: conjugatedRingCenters)
            } else {
                strokeLine(p1.offsetBy(dx: px * doubleBondHalfSeparation, dy: py * doubleBondHalfSeparation),
                           p2.offsetBy(dx: px * doubleBondHalfSeparation, dy: py * doubleBondHalfSeparation))
                strokeLine(p1.offsetBy(dx: -px * doubleBondHalfSeparation, dy: -py * doubleBondHalfSeparation),
                           p2.offsetBy(dx: -px * doubleBondHalfSeparation, dy: -py * doubleBondHalfSeparation))
            }
        case .triple:
            strokeLine(p1, p2)
            strokeLine(p1.offsetBy(dx: px * tripleBondOffset, dy: py * tripleBondOffset),
                       p2.offsetBy(dx: px * tripleBondOffset, dy: py * tripleBondOffset))
            strokeLine(p1.offsetBy(dx: -px * tripleBondOffset, dy: -py * tripleBondOffset),
                       p2.offsetBy(dx: -px * tripleBondOffset, dy: -py * tripleBondOffset))
        case .aromatic:
            strokeLine(p1, p2, width: max(1.2, style.bondWidth * 0.85), opacity: 0.8)
        }
    }

    private static func drawAtomLabel(atom: Atom,
                                      at pos: CGPoint,
                                      color: CDKRenderColor,
                                      style: RenderStyle,
                                      context: inout GraphicsContext,
                                      degree: Int,
                                      implicitHydrogenCount: Int) {
        let label = CDKLabelText.build(atom: atom, style: style, implicitHydrogenCount: implicitHydrogenCount)
        let centerOffset = CDKLabelText.centerOffset(atom: atom,
                                                     style: style,
                                                     implicitHydrogenCount: implicitHydrogenCount,
                                                     fontSize: style.fontSize)
        let labelCenter = pos.offsetBy(dx: centerOffset.dx, dy: centerOffset.dy)

        let text = Text(label)
            .font(.system(size: style.fontSize, weight: .semibold, design: .rounded))
            .foregroundStyle(swiftUIColor(from: color.withAlpha(atom.aromatic ? 0.90 : 1.0)))

        let resolved = context.resolve(text)
        let size = resolved.measure(in: CGSize(width: 200, height: 200))

        let box = CGRect(
            x: labelCenter.x - size.width/2 - 4,
            y: labelCenter.y - size.height/2 - 2,
            width: size.width + 8,
            height: size.height + 4
        )

        // A subtle background behind labels to keep them legible over bonds.
        let labelBackground = Color(.sRGB, red: 1.0, green: 1.0, blue: 1.0, opacity: 0.82)
        context.fill(Path(roundedRect: box, cornerRadius: 6), with: .color(labelBackground))

        context.draw(resolved, at: CGPoint(x: labelCenter.x, y: labelCenter.y), anchor: .center)

        // Tiny dot for carbon that’s shown only because it’s terminal.
        if atom.element == "C" && !style.showCarbons && atom.charge == 0 && degree <= 1 {
            let r: CGFloat = 2.2
            context.fill(Path(ellipseIn: CGRect(x: pos.x - r, y: pos.y - r, width: 2*r, height: 2*r)),
                         with: .color(swiftUIColor(from: color.withAlpha(0.55))))
        }
    }

    private static func drawAromaticCircle(points: [CGPoint],
                                           color: CDKRenderColor,
                                           bondWidth: CGFloat,
                                           context: inout GraphicsContext) {
        guard points.count >= 5 else { return }
        let center = points.reduce(CGPoint.zero) { acc, p in
            CGPoint(x: acc.x + p.x / CGFloat(points.count),
                    y: acc.y + p.y / CGFloat(points.count))
        }
        let minRadius = points.map { hypot($0.x - center.x, $0.y - center.y) }.min() ?? 0
        guard minRadius > 1 else { return }
        let radius = minRadius * 0.53
        let rect = CGRect(x: center.x - radius, y: center.y - radius, width: radius * 2, height: radius * 2)
        context.stroke(Path(ellipseIn: rect),
                       with: .color(swiftUIColor(from: color.withAlpha(0.82))),
                       lineWidth: max(1.0, bondWidth * 0.84))
    }

    private static func swiftUIColor(from color: CDKRenderColor) -> Color {
        Color(.sRGB,
              red: Double(color.red),
              green: Double(color.green),
              blue: Double(color.blue),
              opacity: Double(color.alpha))
    }
}
