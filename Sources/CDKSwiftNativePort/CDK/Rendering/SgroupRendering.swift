import CoreGraphics
import Foundation

struct CDKBracketSegment: Hashable {
    let from: CGPoint
    let to: CGPoint
}

struct CDKSgroupBracketAnnotation: Hashable {
    let segments: [CDKBracketSegment]
    let subscriptText: String?
    let subscriptPosition: CGPoint?
    let superscriptText: String?
    let superscriptPosition: CGPoint?
}

enum CDKSgroupRendering {
    static func bracketAnnotations(molecule: Molecule,
                                   positionsByAtomID: [Int: CGPoint],
                                   transformPoint: ((CGPoint) -> CGPoint)? = nil,
                                   fontSize: CGFloat,
                                   bondWidth: CGFloat) -> [CDKSgroupBracketAnnotation] {
        molecule.sgroups.compactMap { sgroup in
            guard shouldRenderBracket(for: sgroup) else { return nil }
            guard let geometry = bracketGeometry(for: sgroup,
                                                molecule: molecule,
                                                positionsByAtomID: positionsByAtomID,
                                                transformPoint: transformPoint,
                                                fontSize: fontSize,
                                                bondWidth: bondWidth) else {
                return nil
            }
            let descriptor = displayDescriptor(for: sgroup, in: molecule)
            if descriptor.subscriptText == nil,
               descriptor.superscriptText == nil,
               geometry.segments.isEmpty {
                return nil
            }
            return CDKSgroupBracketAnnotation(segments: geometry.segments,
                                              subscriptText: descriptor.subscriptText,
                                              subscriptPosition: geometry.subscriptPosition,
                                              superscriptText: descriptor.superscriptText,
                                              superscriptPosition: geometry.superscriptPosition)
        }
    }

    private struct AnnotationGeometry {
        let segments: [CDKBracketSegment]
        let subscriptPosition: CGPoint?
        let superscriptPosition: CGPoint?
    }

    private struct DisplayDescriptor {
        let subscriptText: String?
        let superscriptText: String?
    }

    private static func shouldRenderBracket(for sgroup: MoleculeSgroup) -> Bool {
        if sgroup.kind == .structureRepeatUnit {
            // CXSMILES link nodes are represented as round-bracket repeat-unit
            // Sgroups and are rendered by the dedicated Markush/link-node path.
            // Keep square-bracket SRUs in the generic Sgroup renderer.
            return !sgroup.roundBrackets
        }
        let keyword = normalizedKeyword(sgroup.keyword)
        switch keyword {
        case "SUP", "DAT", "M":
            return false
        default:
            return !sgroup.atomIDs.isEmpty
        }
    }

    private static func bracketGeometry(for sgroup: MoleculeSgroup,
                                        molecule: Molecule,
                                        positionsByAtomID: [Int: CGPoint],
                                        transformPoint: ((CGPoint) -> CGPoint)?,
                                        fontSize: CGFloat,
                                        bondWidth: CGFloat) -> AnnotationGeometry? {
        let points = sgroup.atomIDs.compactMap { positionsByAtomID[$0] }
        guard !points.isEmpty else { return nil }
        let groupBounds = bounds(of: points)

        let bracketSegments: [CDKBracketSegment]
        let rightMid: CGPoint

        if sgroup.brackets.count >= 2, let transformPoint {
            let transformedBrackets = sgroup.brackets.map { bracket in
                MoleculeSgroupBracket(firstPoint: transformPoint(bracket.firstPoint),
                                      secondPoint: transformPoint(bracket.secondPoint))
            }
            bracketSegments = transformedBrackets.flatMap { bracket in
                bracketSegmentsForBracket(bracket,
                                          groupCenter: CGPoint(x: groupBounds.midX, y: groupBounds.midY),
                                          fontSize: fontSize,
                                          bondWidth: bondWidth)
            }
            let ordered = transformedBrackets.sorted { lhs, rhs in
                ((lhs.firstPoint.x + lhs.secondPoint.x) * 0.5) < ((rhs.firstPoint.x + rhs.secondPoint.x) * 0.5)
            }
            guard let last = ordered.last else { return nil }
            rightMid = midpoint(last.firstPoint, last.secondPoint)
        } else {
            let inset = max(fontSize * 0.22, bondWidth * 1.8, 5.5)
            let halfHeight = max(groupBounds.height * 0.52, fontSize * 1.1)
            let topY = groupBounds.midY - halfHeight
            let bottomY = groupBounds.midY + halfHeight
            let leftX = groupBounds.minX - max(fontSize * 0.52, bondWidth * 2.2, 9)
            let rightX = groupBounds.maxX + max(fontSize * 0.52, bondWidth * 2.2, 9)
            let leftTop = CGPoint(x: leftX, y: topY)
            let leftBottom = CGPoint(x: leftX, y: bottomY)
            let rightTop = CGPoint(x: rightX, y: topY)
            let rightBottom = CGPoint(x: rightX, y: bottomY)
            bracketSegments = [
                CDKBracketSegment(from: leftTop, to: leftBottom),
                CDKBracketSegment(from: leftTop, to: CGPoint(x: leftTop.x + inset, y: leftTop.y)),
                CDKBracketSegment(from: leftBottom, to: CGPoint(x: leftBottom.x + inset, y: leftBottom.y)),
                CDKBracketSegment(from: rightTop, to: rightBottom),
                CDKBracketSegment(from: rightTop, to: CGPoint(x: rightTop.x - inset, y: rightTop.y)),
                CDKBracketSegment(from: rightBottom, to: CGPoint(x: rightBottom.x - inset, y: rightBottom.y))
            ]
            rightMid = midpoint(rightTop, rightBottom)
        }

        let subscriptPosition = CGPoint(x: rightMid.x + max(fontSize * 0.38, 6),
                                        y: groupBounds.maxY + max(fontSize * 0.70, 10))
        let superscriptPosition = CGPoint(x: rightMid.x + max(fontSize * 0.30, 5),
                                          y: groupBounds.minY - max(fontSize * 0.68, 10))

        return AnnotationGeometry(segments: bracketSegments,
                                  subscriptPosition: subscriptPosition,
                                  superscriptPosition: superscriptPosition)
    }

    private static func displayDescriptor(for sgroup: MoleculeSgroup, in molecule: Molecule) -> DisplayDescriptor {
        let keyword = normalizedKeyword(sgroup.keyword)
        let subtype = sgroup.subtype?.trimmingCharacters(in: .whitespacesAndNewlines).uppercased() ?? ""
        let defaultSubscript = sgroup.subscriptText?.trimmingCharacters(in: .whitespacesAndNewlines)
        let connectivity = sgroup.connectivity?.trimmingCharacters(in: .whitespacesAndNewlines).lowercased()

        switch keyword {
        case "SRU", "N":
            let subscriptLabel = normalizedText(defaultSubscript) ?? "n"
            let superscript: String?
            if sgroup.atomIDs.count == 1 || connectivity == "ht" {
                superscript = nil
            } else {
                superscript = normalizedText(connectivity ?? "eu")
            }
            return DisplayDescriptor(subscriptText: subscriptLabel, superscriptText: superscript)
        case "COP", "CO":
            let copolymerLabel: String
            switch subtype {
            case "RAN":
                copolymerLabel = "ran"
            case "BLK":
                copolymerLabel = "blk"
            case "ALT":
                copolymerLabel = "alt"
            default:
                copolymerLabel = "co"
            }
            return DisplayDescriptor(subscriptText: copolymerLabel, superscriptText: nil)
        case "CRO":
            return DisplayDescriptor(subscriptText: "xl", superscriptText: nil)
        case "MON":
            return DisplayDescriptor(subscriptText: "mon", superscriptText: nil)
        case "MER":
            return DisplayDescriptor(subscriptText: "mer", superscriptText: nil)
        case "GRA":
            return DisplayDescriptor(subscriptText: "grf", superscriptText: nil)
        case "MOD":
            return DisplayDescriptor(subscriptText: "mod", superscriptText: nil)
        case "MIX":
            return DisplayDescriptor(subscriptText: "mix", superscriptText: childDataSuperscript(for: sgroup, in: molecule))
        case "FOR", "F":
            return DisplayDescriptor(subscriptText: "f", superscriptText: childDataSuperscript(for: sgroup, in: molecule))
        case "COM", "C":
            let componentLabel: String
            if let componentNumber = sgroup.componentNumber, componentNumber > 0 {
                componentLabel = "c\(componentNumber)"
            } else {
                componentLabel = "c"
            }
            return DisplayDescriptor(subscriptText: componentLabel, superscriptText: childDataSuperscript(for: sgroup, in: molecule))
        case "MUL":
            return DisplayDescriptor(subscriptText: normalizedText(defaultSubscript),
                                     superscriptText: nil)
        case "ANY":
            return DisplayDescriptor(subscriptText: "any", superscriptText: nil)
        case "GEN":
            return DisplayDescriptor(subscriptText: normalizedText(defaultSubscript) ?? "gen",
                                     superscriptText: nil)
        default:
            return DisplayDescriptor(subscriptText: normalizedText(defaultSubscript),
                                     superscriptText: normalizedText(connectivity))
        }
    }

    private static func childDataSuperscript(for sgroup: MoleculeSgroup, in molecule: Molecule) -> String? {
        for childIndex in sgroup.childGroupIndices {
            guard molecule.sgroups.indices.contains(childIndex) else { continue }
            let child = molecule.sgroups[childIndex]
            guard child.kind == .data else { continue }
            if let dataValue = normalizedText(child.dataValue) {
                if let unit = normalizedText(child.dataUnit) {
                    return dataValue + unit
                }
                return dataValue
            }
        }
        return nil
    }

    private static func bracketSegmentsForBracket(_ bracket: MoleculeSgroupBracket,
                                                  groupCenter: CGPoint,
                                                  fontSize: CGFloat,
                                                  bondWidth: CGFloat) -> [CDKBracketSegment] {
        let start = bracket.firstPoint
        let end = bracket.secondPoint
        let midpoint = midpoint(start, end)
        let line = CGVector(dx: end.x - start.x, dy: end.y - start.y)
        let lineLength = max(0.0001, hypot(line.dx, line.dy))
        let normalA = CGVector(dx: -line.dy / lineLength, dy: line.dx / lineLength)
        let normalB = CGVector(dx: line.dy / lineLength, dy: -line.dx / lineLength)
        let toCenter = CGVector(dx: groupCenter.x - midpoint.x, dy: groupCenter.y - midpoint.y)
        let inward = ((toCenter.dx * normalA.dx) + (toCenter.dy * normalA.dy)) >= 0 ? normalA : normalB
        let depth = max(fontSize * 0.34, bondWidth * 2.2, 5.5)
        let firstInset = CGPoint(x: start.x + inward.dx * depth, y: start.y + inward.dy * depth)
        let secondInset = CGPoint(x: end.x + inward.dx * depth, y: end.y + inward.dy * depth)
        return [
            CDKBracketSegment(from: start, to: end),
            CDKBracketSegment(from: start, to: firstInset),
            CDKBracketSegment(from: end, to: secondInset)
        ]
    }

    private static func bounds(of points: [CGPoint]) -> CGRect {
        points.reduce(into: CGRect.null) { bounds, point in
            bounds = bounds.union(CGRect(x: point.x, y: point.y, width: 0, height: 0))
        }
    }

    private static func midpoint(_ a: CGPoint, _ b: CGPoint) -> CGPoint {
        CGPoint(x: (a.x + b.x) * 0.5, y: (a.y + b.y) * 0.5)
    }

    private static func normalizedKeyword(_ keyword: String?) -> String {
        keyword?.trimmingCharacters(in: .whitespacesAndNewlines).uppercased() ?? ""
    }

    private static func normalizedText(_ text: String?) -> String? {
        guard let text else { return nil }
        let trimmed = text.trimmingCharacters(in: .whitespacesAndNewlines)
        return trimmed.isEmpty ? nil : trimmed.lowercased()
    }
}
