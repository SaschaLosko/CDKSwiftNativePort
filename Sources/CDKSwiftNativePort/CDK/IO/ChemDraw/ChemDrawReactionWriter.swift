import Foundation

extension ChemDrawBuilder {
    mutating func reactionPage(_ reaction: CDKReaction) throws -> ChemDrawObject {
        guard !reaction.participants.isEmpty else { throw ChemError.emptyInput }
        guard reaction.participants.allSatisfy({ $0.stoichiometry == nil || $0.stoichiometry == 1 }),
            reaction.direction != .noGo
        else {
            throw ChemError.unsupported(
                "ChemDraw export does not yet support stoichiometric coefficients or no-go arrows.")
        }
        // Reaction-level CX annotations are not necessarily applied to participants.
        // Refuse them until their chemistry can be represented without information loss.
        if let state = reaction.cxState, state != CDKCxSmilesState() {
            throw ChemError.unsupported(
                "ChemDraw export does not yet support reaction-level CXSMILES annotations. Use CML.")
        }
        var page = try object(0x8001, "page")
        if let name = reaction.name { page.children.append(try text(name, x: 40, y: 25)) }
        // Place agents above the arrow and substrates/products below, without overlap.
        var agentIDs: [UInt32] = []
        var agentWidth = 0.0
        var agentHeight = 0.0
        for molecule in reaction.agents {
            let item = try fragment(molecule, x: 40 + agentWidth, y: 50)
            page.children.append(item.object)
            agentIDs.append(item.object.id)
            agentWidth += item.width + 45
            agentHeight = max(agentHeight, item.height)
        }
        let y = 100 + agentHeight
        var x = 40.0
        var reactantIDs: [UInt32] = []
        var productIDs: [UInt32] = []
        var plusIDs: [UInt32] = []
        var reactantMaps: [Int: UInt32] = [:]
        var productMaps: [Int: UInt32] = [:]
        for (index, molecule) in reaction.reactants.enumerated() {
            let item = try fragment(molecule, x: x, y: y)
            page.children.append(item.object)
            reactantIDs.append(item.object.id)
            try merge(item.maps, into: &reactantMaps)
            x += item.width + 25
            if index < reaction.reactants.count - 1 {
                let plus = try text("+", x: x, y: y + 15)
                page.children.append(plus)
                plusIDs.append(plus.id)
                x += 30
            }
        }
        let arrowStart = x
        x += 75
        let arrowEnd = x
        x += 30
        let agentShift = (arrowStart + arrowEnd) / 2 - (40 + max(0, agentWidth - 45) / 2)
        for index in page.children.indices where agentIDs.contains(page.children[index].id) {
            for nodeIndex in page.children[index].children.indices {
                for propertyIndex in page.children[index].children[nodeIndex].properties.indices {
                    let property = page.children[index].children[nodeIndex].properties[propertyIndex]
                    if property.tag == 0x0200 {
                        let values = property.text.split(separator: " ").compactMap { Double($0) }
                        page.children[index].children[nodeIndex].properties[propertyIndex] = try .coordinates(
                            0x0200, "p", [values[0] + agentShift, values[1]], binaryOrder: [1, 0])
                    }
                }
            }
        }
        for (index, molecule) in reaction.products.enumerated() {
            let item = try fragment(molecule, x: x, y: y)
            page.children.append(item.object)
            productIDs.append(item.object.id)
            try merge(item.maps, into: &productMaps)
            x += item.width + 25
            if index < reaction.products.count - 1 {
                let plus = try text("+", x: x, y: y + 15)
                page.children.append(plus)
                plusIDs.append(plus.id)
                x += 30
            }
        }
        var arrow = try object(0x8007, "graphic")
        let arrowType: (UInt16, String) =
            switch reaction.direction {
            case .forward, .backward: (2, "FullHead")
            case .bidirectional: (8, "Equilibrium")
            case .resonance: (4, "Resonance")
            case .retroSynthetic: (32, "RetroSynthetic")
            case .undirected: (0, "NoHead")
            case .noGo: (0, "NoHead")  // Rejected above.
            }
        let head = reaction.direction == .backward ? arrowStart : arrowEnd
        let tail = reaction.direction == .backward ? arrowEnd : arrowStart
        // Graphic BoundingBox encodes the directed head then tail, not sorted bounds.
        arrow.properties = [
            try .coordinates(
                0x0204, "BoundingBox", [head, y + 15, tail, y + 15], binaryOrder: [1, 0, 3, 2]),
            .integer(0x0A00, "GraphicType", UInt16(1), text: "Line"),
            .integer(0x0A02, "ArrowType", arrowType.0, text: arrowType.1),
        ]
        page.children.append(arrow)
        var scheme = try object(0x800D, "scheme")
        var step = try object(0x800E, "step")
        var mapping: [UInt32] = []
        for number in reactantMaps.keys.sorted() {
            if let product = productMaps[number] { mapping += [reactantMaps[number]!, product] }
        }
        step.properties = [
            .references(0x0C01, "ReactionStepReactants", reactantIDs),
            .references(0x0C02, "ReactionStepProducts", productIDs),
            .references(0x0C03, "ReactionStepPlusses", plusIDs),
            .references(0x0C04, "ReactionStepArrows", [arrow.id]),
            .references(0x0C05, "ReactionStepObjectsAboveArrow", agentIDs),
            .references(0x0C00, "ReactionStepAtomMap", mapping),
            .references(0x0C07, "ReactionStepAtomMapManual", mapping),
        ]
        scheme.children.append(step)
        page.children.append(scheme)
        return page
    }

    private func merge(_ source: [Int: UInt32], into destination: inout [Int: UInt32]) throws {
        for (number, id) in source {
            guard destination.updateValue(id, forKey: number) == nil else {
                throw ChemError.unsupported(
                    "ChemDraw export requires unique atom map numbers on each reaction side.")
            }
        }
    }
}
