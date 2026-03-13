import CoreGraphics
import Foundation

struct CDKMarkushBoxAnnotation: Hashable {
    let label: String
    let boxRect: CGRect
    let labelPosition: CGPoint
    let cornerRadius: CGFloat
    let fillColor: CDKRenderColor
}

struct CDKMarkushLinkAnnotation: Hashable {
    let leftBracketPosition: CGPoint
    let rightBracketPosition: CGPoint
    let subscriptText: String?
    let subscriptPosition: CGPoint?
    let superscriptText: String?
    let superscriptPosition: CGPoint?
}

struct CDKMarkushWaveSegment: Hashable {
    let from: CGPoint
    let to: CGPoint
}

enum CDKMarkushRendering {
    private static let rGroupBoxColor = CDKRenderColor(red: 0.86,
                                                       green: 0.89,
                                                       blue: 0.97,
                                                       alpha: 0.90)

    static func rGroupBoxes(molecule: Molecule,
                            positionsByAtomID: [Int: CGPoint],
                            style: RenderStyle,
                            padding: CGFloat,
                            fontSize: CGFloat,
                            bondWidth: CGFloat,
                            scale: CGFloat = 1.0) -> [CDKMarkushBoxAnnotation] {
        var boundsByLabel: [String: CGRect] = [:]
        let atomByID = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0) })
        var degreeByAtomID: [Int: Int] = [:]
        degreeByAtomID.reserveCapacity(molecule.atoms.count)

        for bond in molecule.bonds {
            degreeByAtomID[bond.a1, default: 0] += 1
            degreeByAtomID[bond.a2, default: 0] += 1
        }

        func union(_ rect: CGRect, with label: String) {
            boundsByLabel[label] = boundsByLabel[label]?.union(rect) ?? rect
        }

        for atom in molecule.atoms {
            guard let label = atom.rGroupMembership,
                  let position = positionsByAtomID[atom.id] else {
                continue
            }

            let degree = degreeByAtomID[atom.id] ?? 0
            if atom.attachmentPoint == nil,
               CDKLabelText.shouldDrawLabel(atom: atom, degree: degree, style: style) {
                let implicitH = style.showImplicitHydrogens ? molecule.implicitHydrogenCount(for: atom.id) : 0
                let text = CDKLabelText.build(atom: atom, style: style, implicitHydrogenCount: implicitH)
                let centerOffset = CDKLabelText.centerOffset(atom: atom,
                                                             style: style,
                                                             implicitHydrogenCount: implicitH,
                                                             fontSize: fontSize)
                let center = position.offsetBy(dx: centerOffset.dx, dy: centerOffset.dy)
                let estimatedTextSize = CDKLabelClipping.estimateLabelSize(text: text, fontSize: fontSize)
                union(CDKLabelClipping.makeLabelRect(center: center, estimatedTextSize: estimatedTextSize),
                      with: label)
            } else {
                let anchorPad = max(padding * 0.5, bondWidth * 0.6)
                union(CGRect(x: position.x, y: position.y, width: 0, height: 0)
                    .insetBy(dx: -anchorPad, dy: -anchorPad),
                      with: label)
            }
        }

        for bond in molecule.bonds {
            guard let atom1 = atomByID[bond.a1],
                  let atom2 = atomByID[bond.a2],
                  let p1 = positionsByAtomID[bond.a1],
                  let p2 = positionsByAtomID[bond.a2] else {
                continue
            }

            let labels = [atom1.rGroupMembership, atom2.rGroupMembership].compactMap { $0 }
            let sharedLabel = labels.first { label in
                labels.filter { $0 == label }.count == labels.count
            }

            if let label = sharedLabel {
                var expansion = max(padding * 0.55, bondWidth * 0.75)
                switch bond.order {
                case .single, .aromatic:
                    break
                case .double:
                    expansion += max(3.4, bondWidth * 2.25) * 0.5
                case .triple:
                    expansion += max(3.1, bondWidth * 2.2)
                }

                let bondRect = CGRect(x: min(p1.x, p2.x),
                                      y: min(p1.y, p2.y),
                                      width: abs(p1.x - p2.x),
                                      height: abs(p1.y - p2.y))
                    .insetBy(dx: -expansion, dy: -expansion)
                union(bondRect, with: label)

                if atom1.attachmentPoint != nil {
                    let wavePad = max(padding * 0.55, bondWidth * 0.72)
                    for segment in attachmentPointSegments(center: p1,
                                                           other: p2,
                                                           style: style,
                                                           scale: scale) {
                        let segmentRect = CGRect(x: min(segment.from.x, segment.to.x),
                                                 y: min(segment.from.y, segment.to.y),
                                                 width: abs(segment.from.x - segment.to.x),
                                                 height: abs(segment.from.y - segment.to.y))
                            .insetBy(dx: -wavePad, dy: -wavePad)
                        union(segmentRect, with: label)
                    }
                }

                if atom2.attachmentPoint != nil {
                    let wavePad = max(padding * 0.55, bondWidth * 0.72)
                    for segment in attachmentPointSegments(center: p2,
                                                           other: p1,
                                                           style: style,
                                                           scale: scale) {
                        let segmentRect = CGRect(x: min(segment.from.x, segment.to.x),
                                                 y: min(segment.from.y, segment.to.y),
                                                 width: abs(segment.from.x - segment.to.x),
                                                 height: abs(segment.from.y - segment.to.y))
                            .insetBy(dx: -wavePad, dy: -wavePad)
                        union(segmentRect, with: label)
                    }
                }
            }
        }

        return boundsByLabel.keys.sorted(by: compareRGroupLabels).compactMap { label in
            guard let rawBounds = boundsByLabel[label] else { return nil }
            let insetX = max(padding * 2.6, fontSize * 0.76, bondWidth * 2.6)
            let insetY = max(padding * 2.2, fontSize * 0.58, bondWidth * 2.2)
            let rect = rawBounds.insetBy(dx: -insetX, dy: -insetY)
            let labelPosition = CGPoint(x: rect.minX - max(fontSize * 1.85, insetX * 2.0),
                                        y: rect.midY)
            return CDKMarkushBoxAnnotation(label: "\(CDKLabelText.displayText(label)) is",
                                           boxRect: rect,
                                           labelPosition: labelPosition,
                                           cornerRadius: max(10, insetY * 1.8),
                                           fillColor: rGroupBoxColor)
        }
    }

    static func linkAnnotations(molecule: Molecule,
                                positionsByAtomID: [Int: CGPoint],
                                fontSize: CGFloat) -> [CDKMarkushLinkAnnotation] {
        molecule.sgroups.compactMap { sgroup in
            guard sgroup.kind == .structureRepeatUnit else { return nil }

            let points = sgroup.atomIDs.compactMap { positionsByAtomID[$0] }
            guard !points.isEmpty else { return nil }

            let bounds = points.reduce(into: CGRect.null) { partial, point in
                partial = partial.union(CGRect(x: point.x, y: point.y, width: 0, height: 0))
            }

            let bracketOffset = max(fontSize * 0.58, 9)
            let left = CGPoint(x: bounds.minX - bracketOffset,
                               y: bounds.midY - (fontSize * 0.08))
            let right = CGPoint(x: bounds.maxX + bracketOffset,
                                y: bounds.midY - (fontSize * 0.08))

            let subscriptPosition = sgroup.subscriptText.map { _ in
                CGPoint(x: right.x + max(fontSize * 0.18, 2),
                        y: bounds.maxY + max(fontSize * 0.68, 10))
            }

            let shouldShowSuperscript = (sgroup.superscriptText?.isEmpty == false)
                && sgroup.superscriptText != "ht"
                && sgroup.atomIDs.count > 1
            let superscriptPosition = shouldShowSuperscript ? CGPoint(
                x: right.x + max(fontSize * 0.18, 2),
                y: bounds.minY - max(fontSize * 0.62, 9)
            ) : nil

            return CDKMarkushLinkAnnotation(leftBracketPosition: left,
                                            rightBracketPosition: right,
                                            subscriptText: sgroup.subscriptText,
                                            subscriptPosition: subscriptPosition,
                                            superscriptText: shouldShowSuperscript ? sgroup.superscriptText : nil,
                                            superscriptPosition: superscriptPosition)
        }
    }

    static func attachmentPointSegments(molecule: Molecule,
                                        positionsByAtomID: [Int: CGPoint],
                                        style: RenderStyle,
                                        scale: CGFloat = 1.0) -> [CDKMarkushWaveSegment] {
        let atomByID = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0) })
        var segments: [CDKMarkushWaveSegment] = []

        for bond in molecule.bonds {
            if let atom = atomByID[bond.a1],
               atom.attachmentPoint != nil,
               let center = positionsByAtomID[bond.a1],
               let other = positionsByAtomID[bond.a2] {
                segments.append(contentsOf: attachmentPointSegments(center: center,
                                                                   other: other,
                                                                   style: style,
                                                                   scale: scale))
            }
            if let atom = atomByID[bond.a2],
               atom.attachmentPoint != nil,
               let center = positionsByAtomID[bond.a2],
               let other = positionsByAtomID[bond.a1] {
                segments.append(contentsOf: attachmentPointSegments(center: center,
                                                                   other: other,
                                                                   style: style,
                                                                   scale: scale))
            }
        }

        return segments
    }

    static func attachmentPointSegments(center: CGPoint,
                                        other: CGPoint,
                                        style: RenderStyle,
                                        scale: CGFloat = 1.0) -> [CDKMarkushWaveSegment] {
        let dx = other.x - center.x
        let dy = other.y - center.y
        let length = hypot(dx, dy)
        guard length > 0.0001 else { return [] }

        let ux = dx / length
        let uy = dy / length
        let px = -uy
        let py = ux

        let span = min(max(length * 0.72, style.fontSize * 0.92 * scale), length)
        let amplitude = min(max(style.bondWidth * 1.7 * scale, 2.4), length * 0.18)
        let halfSpan = span * 0.5
        let start = center.offsetBy(dx: px * halfSpan, dy: py * halfSpan)

        let pointCount = 13
        var points: [CGPoint] = []
        points.reserveCapacity(pointCount)

        for index in 0..<pointCount {
            let t = CGFloat(index) / CGFloat(pointCount - 1)
            let base = start.offsetBy(dx: -px * span * t, dy: -py * span * t)
            let phase = (t - 0.5) * .pi * 4
            let offset = sin(phase) * amplitude
            points.append(base.offsetBy(dx: ux * offset, dy: uy * offset))
        }

        return zip(points, points.dropFirst()).map {
            CDKMarkushWaveSegment(from: $0.0, to: $0.1)
        }
    }

    private static func compareRGroupLabels(_ lhs: String, _ rhs: String) -> Bool {
        let leftNumber = Int(lhs.drop(while: { $0.isLetter }))
        let rightNumber = Int(rhs.drop(while: { $0.isLetter }))
        if let leftNumber, let rightNumber, leftNumber != rightNumber {
            return leftNumber < rightNumber
        }
        return lhs < rhs
    }
}
