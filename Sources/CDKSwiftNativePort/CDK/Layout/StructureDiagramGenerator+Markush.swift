import Foundation

#if canImport(CoreGraphics)
    import CoreGraphics
#endif

extension CDKStructureDiagramGenerator {
    static func compareComponents(
        _ lhs: Set<Int>,
        _ rhs: Set<Int>,
        molecule: Molecule,
        markush: Bool
    ) -> Bool {
        if markush {
            let leftLabel = componentRGroupMembership(lhs, molecule: molecule)
            let rightLabel = componentRGroupMembership(rhs, molecule: molecule)
            if leftLabel != rightLabel {
                if leftLabel == nil { return true }
                if rightLabel == nil { return false }
                return (leftLabel ?? "") < (rightLabel ?? "")
            }

            let leftGroupID = componentGroupID(lhs, molecule: molecule)
            let rightGroupID = componentGroupID(rhs, molecule: molecule)
            if leftGroupID != rightGroupID {
                if leftGroupID == nil { return true }
                if rightGroupID == nil { return false }
                return (leftGroupID ?? 0) < (rightGroupID ?? 0)
            }
        }

        return (lhs.min() ?? 0) < (rhs.min() ?? 0)
    }

    static func componentRGroupMembership(_ component: Set<Int>, molecule: Molecule)
        -> String?
    {
        component
            .compactMap { molecule.atom(id: $0)?.rGroupMembership }
            .sorted()
            .first
    }

    static func componentGroupID(_ component: Set<Int>, molecule: Molecule) -> Int? {
        component
            .compactMap { molecule.atom(id: $0)?.componentGroupID }
            .sorted()
            .first
    }

    static func arrangeMarkushComponents(
        components: [Set<Int>],
        molecule: Molecule,
        positions: inout [Int: CGPoint],
        spacing: CGFloat
    ) {
        guard !components.isEmpty else { return }

        let rows = buildMarkushRows(from: components, molecule: molecule)
        guard !rows.isEmpty else { return }

        var rowBounds: [CGRect] = []
        rowBounds.reserveCapacity(rows.count)

        for row in rows {
            layoutMarkushGrid(
                cells: row,
                positions: &positions,
                spacing: spacing)
            rowBounds.append(markushBounds(for: row, positions: positions) ?? .zero)
        }

        var yOffsets = Array(repeating: CGFloat(0), count: rows.count + 1)
        for row in 0..<rows.count {
            let rowHeight = spacing + rowBounds[row].height
            if rowHeight > yOffsets[row] {
                yOffsets[row] = rowHeight
            }
        }
        for row in stride(from: rows.count - 1, through: 0, by: -1) {
            yOffsets[row] += yOffsets[row + 1]
        }

        for row in 0..<rows.count {
            guard let currentBounds = markushBounds(for: rows[row], positions: positions) else {
                continue
            }
            let dest = CGPoint(x: 0, y: (yOffsets[row] + yOffsets[row + 1]) * 0.5)
            translateMarkushCells(
                rows[row],
                positions: &positions,
                dx: dest.x - currentBounds.midX,
                dy: dest.y - currentBounds.midY)
        }
    }

    static func buildMarkushRows(
        from components: [Set<Int>],
        molecule: Molecule
    ) -> [[[Set<Int>]]] {
        var rows: [[[Set<Int>]]] = []
        var currentRow: [[Set<Int>]] = []
        var currentCell: [Set<Int>] = []
        var currentLabel: String?
        var currentGroupID: Int?
        var started = false

        func flushCell() {
            guard !currentCell.isEmpty else { return }
            currentRow.append(currentCell)
            currentCell = []
            currentGroupID = nil
        }

        func flushRow() {
            flushCell()
            guard !currentRow.isEmpty else { return }
            rows.append(currentRow)
            currentRow = []
        }

        for component in components {
            let label = componentRGroupMembership(component, molecule: molecule)
            let groupID = componentGroupID(component, molecule: molecule)

            if started, label != currentLabel {
                flushRow()
            }

            if currentCell.isEmpty {
                currentCell = [component]
                currentGroupID = groupID
            } else if let groupID, groupID == currentGroupID {
                currentCell.append(component)
            } else {
                flushCell()
                currentCell = [component]
                currentGroupID = groupID
            }

            currentLabel = label
            started = true
        }

        flushRow()
        return rows
    }

    static func layoutMarkushGrid(
        cells: [[Set<Int>]],
        positions: inout [Int: CGPoint],
        spacing: CGFloat
    ) {
        guard !cells.isEmpty else { return }

        let limits = cells.map { markushBounds(for: [$0], positions: positions) ?? .zero }
        let fragmentCount = cells.count
        let rowCount = max(1, Int(floor(sqrt(Double(fragmentCount)))))
        let columnCount = Int(ceil(Double(fragmentCount) / Double(rowCount)))

        var xOffsets = Array(repeating: CGFloat(0), count: columnCount + 1)
        var yOffsets = Array(repeating: CGFloat(0), count: rowCount + 1)

        for index in 0..<fragmentCount {
            let column = 1 + index % columnCount
            let row = index / columnCount
            let bounds = limits[index]
            let width = spacing + bounds.width
            let height = spacing + bounds.height
            if width > xOffsets[column] {
                xOffsets[column] = width
            }
            if height > yOffsets[row] {
                yOffsets[row] = height
            }
        }

        for index in 1..<xOffsets.count {
            xOffsets[index] += xOffsets[index - 1]
        }
        for index in stride(from: yOffsets.count - 2, through: 0, by: -1) {
            yOffsets[index] += yOffsets[index + 1]
        }

        for index in 0..<fragmentCount {
            let row = index / columnCount
            let column = index % columnCount
            let bounds = limits[index]
            let dest = CGPoint(
                x: (xOffsets[column] + xOffsets[column + 1]) * 0.5,
                y: (yOffsets[row] + yOffsets[row + 1]) * 0.5)
            translateMarkushCells(
                [cells[index]],
                positions: &positions,
                dx: dest.x - bounds.midX,
                dy: dest.y - bounds.midY)
        }
    }

    static func markushBounds(
        for cells: [[Set<Int>]],
        positions: [Int: CGPoint]
    ) -> CGRect? {
        var result: CGRect?
        for cell in cells {
            for component in cell {
                guard let box = boundingBox(component: component, positions: positions) else { continue }
                result = result?.union(box) ?? box
            }
        }
        return result
    }

    static func translateMarkushCells(
        _ cells: [[Set<Int>]],
        positions: inout [Int: CGPoint],
        dx: CGFloat,
        dy: CGFloat
    ) {
        guard abs(dx) > 0.0001 || abs(dy) > 0.0001 else { return }
        for cell in cells {
            for component in cell {
                for atomID in component {
                    guard let point = positions[atomID] else { continue }
                    positions[atomID] = CGPoint(x: point.x + dx, y: point.y + dy)
                }
            }
        }
    }

    static func selectOrientation(
        component: Set<Int>,
        graph: CDKMolecularGraph,
        molecule: Molecule,
        positions: inout [Int: CGPoint],
        widthDiff: CGFloat,
        alignDiff: Int
    ) {
        let bonds = orientableBonds(
            in: component,
            graph: graph,
            molecule: molecule)
        guard !bonds.isEmpty,
            let bounds = boundingBox(component: component, positions: positions)
        else {
            return
        }

        let pivot = CGPoint(x: bounds.midX, y: bounds.midY)

        if bonds.count == 1,
            let p1 = positions[bonds[0].a1],
            let p2 = positions[bonds[0].a2]
        {
            let dx = p1.x - p2.x
            let dy = p1.y - p2.y
            rotate(
                component: component,
                around: pivot,
                by: atan2(dx, dy) + (.pi / 2),
                positions: &positions)
            return
        }

        var directionHistogram = Array(repeating: 0, count: 180)
        calcDirectionHistogram(
            bonds: bonds,
            positions: positions,
            counts: &directionHistogram,
            limit: 60)
        if let maxBucket = directionHistogram.indices.max(by: {
            directionHistogram[$0] < directionHistogram[$1]
        }),
            maxBucket != 0,
            Double(directionHistogram[maxBucket]) / Double(max(1, bonds.count)) > 0.5
        {
            rotate(
                component: component,
                around: pivot,
                by: CGFloat((60 - maxBucket)) * .pi / 180.0,
                positions: &positions)
        }

        guard let startingBounds = boundingBox(component: component, positions: positions) else {
            return
        }
        let baselineWidth = startingBounds.width
        var bestWidth = baselineWidth

        calcDirectionHistogram(
            bonds: bonds,
            positions: positions,
            counts: &directionHistogram,
            limit: 180)
        var maxAligned = directionHistogram[60] + directionHistogram[120]
        var bestCoordinates = Dictionary(
            uniqueKeysWithValues: component.compactMap { atomID in
                positions[atomID].map { (atomID, $0) }
            })

        let step = CGFloat.pi / 3
        var total: CGFloat = 0

        while total < (2 * .pi) {
            total += step
            rotate(
                component: component,
                around: pivot,
                by: step,
                positions: &positions)

            guard let currentBounds = boundingBox(component: component, positions: positions) else {
                continue
            }
            let width = currentBounds.width
            let widthDelta = abs(width - baselineWidth)

            if widthDelta >= widthDiff && width > bestWidth {
                bestWidth = width
                bestCoordinates = Dictionary(
                    uniqueKeysWithValues: component.compactMap { atomID in
                        positions[atomID].map { (atomID, $0) }
                    })
                continue
            }

            guard widthDelta <= widthDiff else { continue }

            calcDirectionHistogram(
                bonds: bonds,
                positions: positions,
                counts: &directionHistogram,
                limit: 180)
            let aligned = directionHistogram[60] + directionHistogram[120]
            let alignDelta = aligned - maxAligned

            if alignDelta > alignDiff || aligned == bonds.count || (alignDelta == 0 && width > bestWidth) {
                maxAligned = aligned
                bestWidth = width
                bestCoordinates = Dictionary(
                    uniqueKeysWithValues: component.compactMap { atomID in
                        positions[atomID].map { (atomID, $0) }
                    })
            }
        }

        for atomID in component {
            if let point = bestCoordinates[atomID] {
                positions[atomID] = point
            }
        }
    }

    static func orientMarkushRootComponent(
        component: Set<Int>,
        molecule: Molecule,
        positions: inout [Int: CGPoint]
    ) {
        guard componentRGroupMembership(component, molecule: molecule) == nil,
            let bounds = boundingBox(component: component, positions: positions)
        else {
            return
        }

        let repeatAtomID = molecule.sgroups
            .filter {
                $0.kind == .structureRepeatUnit && $0.atomIDs.count == 1
                    && component.contains($0.atomIDs[0])
            }
            .compactMap(\.atomIDs.first)
            .first

        guard let repeatAtomID,
            let repeatPosition = positions[repeatAtomID]
        else {
            return
        }

        let center = CGPoint(x: bounds.midX, y: bounds.midY)
        let dx = repeatPosition.x - center.x
        let dy = repeatPosition.y - center.y
        guard hypot(dx, dy) > 0.0001 else { return }

        rotate(
            component: component,
            around: center,
            by: .pi - atan2(dy, dx),
            positions: &positions)
    }

    static func orientMarkushAttachmentComponent(
        component: Set<Int>,
        graph: CDKMolecularGraph,
        molecule: Molecule,
        positions: inout [Int: CGPoint]
    ) -> Bool {
        guard componentRGroupMembership(component, molecule: molecule) != nil else {
            return false
        }

        let attachmentAtoms = component.compactMap { atomID in
            molecule.atom(id: atomID).flatMap { atom -> Atom? in
                atom.attachmentPoint == 1 ? atom : nil
            }
        }
        guard let attachment = attachmentAtoms.first,
            component.count > 2,
            let bounds = boundingBox(component: component, positions: positions),
            let attachmentNeighbors = positions[attachment.id].map({ _ in
                graph.neighbors(of: attachment.id).filter(component.contains)
            }),
            attachmentNeighbors.count == 1,
            let attachmentPosition = positions[attachment.id],
            let neighborPosition = positions[attachmentNeighbors[0]]
        else {
            return false
        }

        let pivot = CGPoint(x: bounds.midX, y: bounds.midY)
        rotate(
            component: component,
            around: pivot,
            by: -atan2(
                neighborPosition.y - attachmentPosition.y,
                neighborPosition.x - attachmentPosition.x),
            positions: &positions)

        guard let rotatedAttachmentPosition = positions[attachment.id] else {
            return false
        }

        var lowestBelow = CGFloat.zero
        var highestAbove = CGFloat.zero
        for atomID in component {
            guard let point = positions[atomID] else { continue }
            let delta = rotatedAttachmentPosition.y - point.y
            if delta > 0, delta > highestAbove {
                highestAbove = delta
            } else if delta < 0, delta < lowestBelow {
                lowestBelow = delta
            }
        }

        if abs(lowestBelow) < highestAbove {
            for atomID in component {
                guard let point = positions[atomID] else { continue }
                positions[atomID] = CGPoint(
                    x: point.x,
                    y: (2 * rotatedAttachmentPosition.y) - point.y)
            }
        }

        rotate(
            component: component,
            around: pivot,
            by: -.pi / 6,
            positions: &positions)
        return true
    }

    static func orientableBonds(
        in component: Set<Int>,
        graph: CDKMolecularGraph,
        molecule: Molecule
    ) -> [Bond] {
        let componentBonds = molecule.bonds.filter {
            component.contains($0.a1) && component.contains($0.a2)
        }
        var filtered: [Bond] = []
        filtered.reserveCapacity(componentBonds.count)

        for bond in componentBonds {
            let beginDegree = graph.neighbors(of: bond.a1).filter { component.contains($0) }.count
            let endDegree = graph.neighbors(of: bond.a2).filter { component.contains($0) }.count

            if (beginDegree == 1 && endDegree > 2) || (endDegree == 1 && beginDegree > 2) {
                continue
            }
            if beginDegree != 2 && endDegree != 2 {
                continue
            }
            filtered.append(bond)
        }

        return filtered.isEmpty ? componentBonds : filtered
    }

    static func calcDirectionHistogram(
        bonds: [Bond],
        positions: [Int: CGPoint],
        counts: inout [Int],
        limit: Int
    ) {
        guard limit > 0, limit <= counts.count else { return }
        for index in counts.indices {
            counts[index] = 0
        }

        for bond in bonds {
            guard let begin = positions[bond.a1], let end = positions[bond.a2] else { continue }
            var vector = CGVector(dx: end.x - begin.x, dy: end.y - begin.y)
            if vector.dx < 0 {
                vector.dx.negate()
                vector.dy.negate()
            }
            let angle = (.pi / 2) + atan2(vector.dy, vector.dx)
            let bucket = Int(lround(Double(angle * 180.0 / .pi))) % limit
            counts[bucket] += 1
        }
    }

    static func rotate(
        component: Set<Int>,
        around pivot: CGPoint,
        by angle: CGFloat,
        positions: inout [Int: CGPoint]
    ) {
        guard abs(angle) > 0.0001 else { return }
        let cosAngle = cos(angle)
        let sinAngle = sin(angle)
        for atomID in component {
            guard let point = positions[atomID] else { continue }
            let dx = point.x - pivot.x
            let dy = point.y - pivot.y
            positions[atomID] = CGPoint(
                x: pivot.x + (dx * cosAngle) - (dy * sinAngle),
                y: pivot.y + (dx * sinAngle) + (dy * cosAngle))
        }
    }
}
