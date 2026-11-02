import Foundation

#if canImport(CoreGraphics)
    import CoreGraphics
#endif

extension CDKStructureDiagramGenerator {
    static func arrangeCoordinateBondedSandwichRings(
        molecule: Molecule,
        rings: [[Int]],
        positions: inout [Int: CGPoint],
        bondLength: CGFloat
    ) {
        let coordinateBondsByReference = Dictionary(
            grouping: molecule.bonds.compactMap { bond -> Bond? in
                bond.coordinateBondReferenceAtomID == nil ? nil : bond
            }, by: { $0.coordinateBondReferenceAtomID ?? -1 })

        for (referenceAtomID, coordinateBonds) in coordinateBondsByReference where referenceAtomID >= 0 {
            let donorAtomIDs = Set(
                coordinateBonds.compactMap { bond -> Int? in
                    if bond.a1 == referenceAtomID { return bond.a2 }
                    if bond.a2 == referenceAtomID { return bond.a1 }
                    return nil
                })
            let attachedRings =
                rings
                .filter { ring in ring.contains(where: donorAtomIDs.contains) }
                .sorted { lhs, rhs in
                    if lhs.count != rhs.count { return lhs.count < rhs.count }
                    return lhs.lexicographicallyPrecedes(rhs)
                }
            guard attachedRings.count == 2,
                Set(attachedRings[0]).isDisjoint(with: attachedRings[1])
            else {
                continue
            }

            let metalPosition = positions[referenceAtomID] ?? .zero
            let verticalOffset = bondLength * 1.30
            for (index, ring) in attachedRings.enumerated() {
                guard let donorAtomID = ring.first(where: donorAtomIDs.contains) else { continue }
                let center = CGPoint(
                    x: metalPosition.x,
                    y: metalPosition.y + (index == 0 ? verticalOffset : -verticalOffset))
                placeSandwichRing(
                    ring,
                    donorAtomID: donorAtomID,
                    center: center,
                    pointsToward: metalPosition,
                    bondLength: bondLength,
                    positions: &positions)
            }
            positions[referenceAtomID] = metalPosition
        }
    }

    static func placeSandwichRing(
        _ ring: [Int],
        donorAtomID: Int,
        center: CGPoint,
        pointsToward target: CGPoint,
        bondLength: CGFloat,
        positions: inout [Int: CGPoint]
    ) {
        guard ring.count >= 3,
            let donorIndex = ring.firstIndex(of: donorAtomID)
        else {
            return
        }
        let count = CGFloat(ring.count)
        let radius = bondLength / (2 * sin(.pi / count))
        let targetAngle = atan2(target.y - center.y, target.x - center.x)
        let step = (2 * .pi) / count
        let base = targetAngle - (CGFloat(donorIndex) * step)

        for (index, atomID) in ring.enumerated() {
            let angle = base + (CGFloat(index) * step)
            positions[atomID] = CGPoint(
                x: center.x + cos(angle) * radius,
                y: center.y + sin(angle) * radius)
        }
    }
}
