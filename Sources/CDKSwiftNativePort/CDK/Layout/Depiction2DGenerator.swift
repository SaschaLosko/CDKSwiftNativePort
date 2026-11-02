import Foundation

#if canImport(CoreGraphics)
    import CoreGraphics
#endif

public enum Depiction2DGenerator {
    public static func generate(for molecule: Molecule) -> Molecule {
        guard let collapsed = CDKHydrogenLayoutCollapse.collapse(in: molecule) else {
            return CDKStructureDiagramGenerator.apply(to: molecule)
        }

        let laidOutScaffold = CDKStructureDiagramGenerator.apply(to: collapsed.scaffold)
        return CDKHydrogenLayoutCollapse.restore(
            original: molecule,
            from: laidOutScaffold,
            hydrogenAnchorByID: collapsed.hydrogenAnchorByID)
    }
}

private struct CDKHydrogenLayoutCollapseResult {
    let scaffold: Molecule
    let hydrogenAnchorByID: [Int: Int]
}

private enum CDKHydrogenLayoutCollapse {
    static func collapse(in molecule: Molecule) -> CDKHydrogenLayoutCollapseResult? {
        let hydrogenAnchorByID = Dictionary(
            uniqueKeysWithValues: molecule.atoms.compactMap { atom -> (Int, Int)? in
                guard isSuppressibleHydrogen(atomID: atom.id, in: molecule),
                    let bond = molecule.bonds(forAtom: atom.id).first
                else {
                    return nil
                }
                let anchorID = (bond.a1 == atom.id) ? bond.a2 : bond.a1
                return (atom.id, anchorID)
            })

        guard !hydrogenAnchorByID.isEmpty else { return nil }

        let suppressedHydrogenIDs = Set(hydrogenAnchorByID.keys)
        var removedHydrogenCountByNeighbor: [Int: Int] = [:]
        for anchorID in hydrogenAnchorByID.values {
            removedHydrogenCountByNeighbor[anchorID, default: 0] += 1
        }

        var scaffold = molecule
        scaffold.atoms = molecule.atoms.compactMap { atom in
            guard !suppressedHydrogenIDs.contains(atom.id) else { return nil }
            var updated = atom
            if let removed = removedHydrogenCountByNeighbor[atom.id], removed > 0 {
                updated.explicitHydrogenCount = max(0, updated.explicitHydrogenCount ?? 0) + removed
            }
            return updated
        }
        scaffold.bonds = molecule.bonds.filter { bond in
            !suppressedHydrogenIDs.contains(bond.a1) && !suppressedHydrogenIDs.contains(bond.a2)
        }

        return CDKHydrogenLayoutCollapseResult(
            scaffold: scaffold,
            hydrogenAnchorByID: hydrogenAnchorByID)
    }

    static func restore(
        original: Molecule,
        from scaffold: Molecule,
        hydrogenAnchorByID: [Int: Int]
    ) -> Molecule {
        var restored = original
        let scaffoldPositionByAtomID = Dictionary(
            uniqueKeysWithValues: scaffold.atoms.map { ($0.id, $0.position) })

        for index in restored.atoms.indices {
            if let position = scaffoldPositionByAtomID[restored.atoms[index].id] {
                restored.atoms[index].position = position
            }
        }

        placeSuppressedHydrogens(in: &restored, hydrogenAnchorByID: hydrogenAnchorByID)
        return restored
    }

    private static func placeSuppressedHydrogens(
        in molecule: inout Molecule,
        hydrogenAnchorByID: [Int: Int]
    ) {
        let groupedByAnchor = Dictionary(
            grouping: hydrogenAnchorByID.keys, by: { hydrogenAnchorByID[$0] ?? -1 })

        for (anchorID, hydrogenIDs) in groupedByAnchor {
            guard let anchorIndex = molecule.atoms.firstIndex(where: { $0.id == anchorID }) else {
                continue
            }
            let anchor = molecule.atoms[anchorIndex]
            let anchorPoint = anchor.position

            let placedNeighborIDs = molecule.neighbors(of: anchorID).filter { neighborID in
                !hydrogenIDs.contains(neighborID)
            }
            let placedNeighborPoints = placedNeighborIDs.compactMap { molecule.atom(id: $0)?.position }
            let heavyNeighborPoints = placedNeighborIDs.compactMap { neighborID -> CGPoint? in
                guard let atom = molecule.atom(id: neighborID),
                    !isHydrogen(atom)
                else {
                    return nil
                }
                return atom.position
            }

            let baseAngle = preferredHydrogenBaseAngle(
                anchor: anchorPoint,
                neighborPoints: heavyNeighborPoints.isEmpty ? placedNeighborPoints : heavyNeighborPoints)
            let step = CGFloat.pi / 3.0
            let centeredOffsets = hydrogenOffsets(count: hydrogenIDs.count, step: step)
            let referenceNeighbors =
                heavyNeighborPoints.isEmpty ? placedNeighborPoints : heavyNeighborPoints
            let meanNeighborDistance =
                referenceNeighbors.isEmpty
                ? CGFloat(1.4)
                : referenceNeighbors.reduce(CGFloat.zero) { partial, point in
                    partial + anchorPoint.distance(to: point)
                } / CGFloat(referenceNeighbors.count)
            let hydrogenBondLength = max(0.72, min(1.15, meanNeighborDistance * 0.78))

            for (offset, hydrogenID) in zip(centeredOffsets, hydrogenIDs.sorted()) {
                guard let hydrogenIndex = molecule.atoms.firstIndex(where: { $0.id == hydrogenID }) else {
                    continue
                }
                let angle = baseAngle + offset
                molecule.atoms[hydrogenIndex].position = CGPoint(
                    x: anchorPoint.x + cos(angle) * hydrogenBondLength,
                    y: anchorPoint.y + sin(angle) * hydrogenBondLength)
            }
        }
    }

    private static func preferredHydrogenBaseAngle(
        anchor: CGPoint,
        neighborPoints: [CGPoint]
    ) -> CGFloat {
        guard !neighborPoints.isEmpty else { return 0 }

        let unitVectors = neighborPoints.compactMap { point -> CGVector? in
            let dx = point.x - anchor.x
            let dy = point.y - anchor.y
            let length = hypot(dx, dy)
            guard length > 0.0001 else { return nil }
            return CGVector(dx: dx / length, dy: dy / length)
        }

        guard !unitVectors.isEmpty else { return 0 }

        let summed = unitVectors.reduce(CGVector.zero) { partial, vector in
            CGVector(dx: partial.dx + vector.dx, dy: partial.dy + vector.dy)
        }
        let magnitude = hypot(summed.dx, summed.dy)
        if magnitude > 0.2 {
            return atan2(-summed.dy, -summed.dx)
        }

        let sortedAngles =
            unitVectors
            .map { atan2($0.dy, $0.dx) }
            .sorted()
        guard sortedAngles.count >= 2 else {
            return sortedAngles[0] + .pi
        }

        var widestGap: CGFloat = -.greatestFiniteMagnitude
        var bestMidpoint = sortedAngles[0] + .pi
        for index in sortedAngles.indices {
            let start = sortedAngles[index]
            let end =
                sortedAngles[(index + 1) % sortedAngles.count]
                + (index + 1 == sortedAngles.count ? (.pi * 2.0) : 0)
            let gap = end - start
            if gap > widestGap {
                widestGap = gap
                bestMidpoint = start + (gap * 0.5)
            }
        }

        return normalizedAngle(bestMidpoint)
    }

    private static func hydrogenOffsets(count: Int, step: CGFloat) -> [CGFloat] {
        guard count > 1 else { return [0] }
        let center = CGFloat(count - 1) * 0.5
        return (0..<count).map { index in
            (CGFloat(index) - center) * step
        }
    }

    private static func normalizedAngle(_ angle: CGFloat) -> CGFloat {
        var normalized = angle.truncatingRemainder(dividingBy: .pi * 2.0)
        if normalized <= -.pi {
            normalized += .pi * 2.0
        } else if normalized > .pi {
            normalized -= .pi * 2.0
        }
        return normalized
    }

    private static func isSuppressibleHydrogen(atomID: Int, in molecule: Molecule) -> Bool {
        guard let atom = molecule.atom(id: atomID), isHydrogen(atom) else { return false }
        guard atom.charge == 0,
            atom.isotopeMassNumber == nil,
            atom.radical == nil,
            atom.queryType == nil,
            atom.atomList == nil
        else {
            return false
        }

        let attachedBonds = molecule.bonds(forAtom: atom.id)
        guard attachedBonds.count == 1, let bond = attachedBonds.first else { return false }
        let neighborID = (bond.a1 == atom.id) ? bond.a2 : bond.a1
        guard bond.order == .single,
            bond.queryType == nil,
            bond.stereo == .none
                || isDirectionalDoubleBondHydrogen(
                    atomID: atomID,
                    neighborID: neighborID,
                    in: molecule)
        else {
            return false
        }
        guard let neighbor = molecule.atom(id: neighborID), !isHydrogen(neighbor) else { return false }
        return true
    }

    private static func isDirectionalDoubleBondHydrogen(
        atomID: Int,
        neighborID: Int,
        in molecule: Molecule
    ) -> Bool {
        let hydrogenNeighborCount = molecule.neighbors(of: neighborID).filter { candidateID in
            molecule.atom(id: candidateID).map(isHydrogen) == true
        }.count
        guard hydrogenNeighborCount >= 2 else { return false }
        return molecule.bonds(forAtom: neighborID).contains { bond in
            bond.order == .double && bond.a1 != atomID && bond.a2 != atomID
        }
    }

    private static func isHydrogen(_ atom: Atom) -> Bool {
        atom.element.trimmingCharacters(in: .whitespacesAndNewlines).uppercased() == "H"
    }
}
