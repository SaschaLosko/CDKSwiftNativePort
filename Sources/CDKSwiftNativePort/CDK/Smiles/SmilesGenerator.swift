import Foundation

/// Swift counterpart of CDK's `org.openscience.cdk.smiles.SmilesGenerator`.
public final class CDKSmilesGenerator {
    public let flavor: CDKSmiFlavor

    private enum RenderEvent {
        case atom(id: Int, parent: Int?)
        case text(String)
    }

    public init(flavor: CDKSmiFlavor = .cdkDefault) {
        self.flavor = flavor
    }

    public func create(_ molecule: Molecule) -> String {
        guard !molecule.atoms.isEmpty else { return "" }
        let core = createCoreSmiles(molecule)
        guard let cxSuffix = createCxLayer(for: molecule, atomOutputOrder: core.atomOutputOrder) else {
            return core.smiles
        }
        return core.smiles + cxSuffix
    }

    public func create(_ reaction: CDKReaction) -> String {
        let reactants = reaction.reactantParticipants.map { createCoreSmiles($0.molecule).smiles }
            .joined(separator: ".")
        let agents = reaction.agentParticipants.map { createCoreSmiles($0.molecule).smiles }.joined(
            separator: ".")
        let products = reaction.productParticipants.map { createCoreSmiles($0.molecule).smiles }.joined(
            separator: ".")

        let base = reactants + ">" + agents + ">" + products
        guard wantsCxOutput, let cxState = reaction.cxState else { return base }
        guard let cxSuffix = createCxLayer(from: cxState, atomOutputOrder: []) else {
            return base
        }
        return base + cxSuffix
    }

    private func createCoreSmiles(_ molecule: Molecule) -> (smiles: String, atomOutputOrder: [Int]) {
        let atomByID = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0) })
        let adjacency = buildAdjacency(molecule)
        let bondByEdge = buildBondLookup(molecule)
        let aromaticBondIDs =
            flavor.contains(.useAromaticSymbols) ? molecule.aromaticDisplayBondIDs() : Set<Int>()

        let components = connectedComponents(molecule, adjacency: adjacency)

        var componentStrings: [String] = []
        componentStrings.reserveCapacity(components.count)
        var atomOutputOrder: [Int] = []

        for component in components {
            let componentSet = Set(component)
            guard let root = component.min() else { continue }

            let (treeAdjacency, treeEdges) = buildTree(
                component: componentSet,
                root: root,
                adjacency: adjacency)

            let componentEdges = Set(
                bondByEdge.keys.filter { componentSet.contains($0.a) && componentSet.contains($0.b) }
            )
            let ringEdges = componentEdges.subtracting(treeEdges)
            let ringNumbers = assignRingNumbers(to: ringEdges)

            var ringClosuresByAtom: [Int: [RingClosure]] = [:]
            for edge in ringEdges {
                guard let bond = bondByEdge[edge], let ringNumber = ringNumbers[edge] else { continue }
                let symbol = bondToken(bond, atomByID: atomByID, aromaticBondIDs: aromaticBondIDs)
                ringClosuresByAtom[edge.a, default: []].append(
                    RingClosure(number: ringNumber, bondSymbol: symbol, partner: edge.b)
                )
                ringClosuresByAtom[edge.b, default: []].append(
                    RingClosure(number: ringNumber, bondSymbol: symbol, partner: edge.a)
                )
            }

            for atomID in ringClosuresByAtom.keys {
                ringClosuresByAtom[atomID]?.sort { lhs, rhs in
                    if lhs.number != rhs.number { return lhs.number < rhs.number }
                    return lhs.partner < rhs.partner
                }
            }

            var visitedTreeAtoms: Set<Int> = []
            let rendered = renderTreeAtom(
                root,
                parent: nil,
                treeAdjacency: treeAdjacency,
                ringClosuresByAtom: ringClosuresByAtom,
                bondByEdge: bondByEdge,
                atomByID: atomByID,
                aromaticBondIDs: aromaticBondIDs,
                atomOutputOrder: &atomOutputOrder,
                visited: &visitedTreeAtoms)
            componentStrings.append(rendered)
        }

        return (componentStrings.joined(separator: "."), atomOutputOrder)
    }

    private var wantsCxOutput: Bool {
        flavor.contains(.cxsmiles) || !flavor.intersection(.cxAll).isEmpty
    }

    private func createCxLayer(for molecule: Molecule, atomOutputOrder: [Int]) -> String? {
        guard wantsCxOutput else { return nil }

        let atomIDs = atomOutputOrder.isEmpty ? molecule.atoms.map(\.id) : atomOutputOrder
        let atomIndexByID = Dictionary(
            uniqueKeysWithValues: atomIDs.enumerated().map { ($0.element, $0.offset) })
        let atomByID = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0) })
        let orderedBonds = molecule.bonds.sorted { $0.id < $1.id }
        let bondIndexByID = Dictionary(
            uniqueKeysWithValues: orderedBonds.enumerated().map { ($1.id, $0) })

        var state = molecule.cxState ?? CDKCxSmilesState()

        let generatedCoordinateBonds = Dictionary(
            grouping: orderedBonds.compactMap { bond -> (Int, Int)? in
                guard let referenceAtomID = bond.coordinateBondReferenceAtomID,
                    let atomIndex = atomIndexByID[referenceAtomID],
                    let bondIndex = bondIndexByID[bond.id]
                else { return nil }
                return (atomIndex, bondIndex)
            }, by: \.0
        )
        .mapValues { entries in entries.map(\.1).sorted() }
        if !generatedCoordinateBonds.isEmpty {
            state.coordinateBonds = generatedCoordinateBonds
        }

        if shouldEmitAtomLabels {
            state.atomLabels = [:]
            for atomID in atomIDs {
                guard let atom = atomByID[atomID],
                    let outputIndex = atomIndexByID[atomID],
                    let label = cxAtomLabel(for: atom)
                else { continue }
                state.atomLabels[outputIndex] = label
            }
        }

        if shouldEmitAtomValues {
            state.atomValues = [:]
            for atomID in atomIDs {
                guard let atom = atomByID[atomID],
                    let outputIndex = atomIndexByID[atomID],
                    let value = atom.atomValue
                else { continue }
                state.atomValues[outputIndex] = value
            }
        }

        if shouldEmitCoordinates(for: molecule) {
            state.atomCoordinates = atomIDs.compactMap { atomID in
                guard let atom = atomByID[atomID] else { return nil }
                return CxCoordinate(
                    x: atom.position.x,
                    y: atom.position.y,
                    z: atom.zPosition ?? 0)
            }
            state.has3DCoordinates = state.atomCoordinates.contains { abs($0.z) > 0.000_000_1 }
        }

        if shouldEmitRadicals {
            state.atomRadicals = [:]
            for atomID in atomIDs {
                guard let atom = atomByID[atomID],
                    let outputIndex = atomIndexByID[atomID],
                    let radicalType = atom.radicalType ?? radicalType(for: atom.radical)
                else { continue }
                state.atomRadicals[outputIndex] = radicalType
            }
        }

        if shouldEmitLigandOrdering {
            state.ligandOrdering = [:]
            for atomID in atomIDs {
                guard let atom = atomByID[atomID],
                    let outputIndex = atomIndexByID[atomID],
                    let orderingIDs = atom.ligandOrderingAtomIDs
                else { continue }
                let orderingIndices = orderingIDs.compactMap { atomIndexByID[$0] }
                guard !orderingIndices.isEmpty else { continue }
                state.ligandOrdering[outputIndex] = orderingIndices
            }
        }

        if shouldEmitEnhancedStereo {
            state.stereoGroups = [:]
            for atomID in atomIDs {
                guard let atom = atomByID[atomID],
                    let outputIndex = atomIndexByID[atomID],
                    let group = atom.cxStereoGroup
                else { continue }
                state.stereoGroups[outputIndex] = group
            }

            if !state.stereoGroups.isEmpty && !state.racemic {
                let groupedIndices = Set(state.stereoGroups.keys)
                for atomID in atomIDs {
                    guard let atom = atomByID[atomID],
                        atom.chirality != .none,
                        let outputIndex = atomIndexByID[atomID],
                        !groupedIndices.contains(outputIndex)
                    else { continue }
                    state.stereoGroups[outputIndex] = CDKCxSmilesParser.encodeStereoGroup(
                        kind: "abs", number: 0)
                }
            }
        }

        if shouldEmitMulticenter {
            state.positionalVariations = [:]
            for sgroup in molecule.sgroups where sgroup.kind == .extMulticenter {
                guard let beginID = sgroup.atomIDs.first,
                    let beginIndex = atomIndexByID[beginID]
                else { continue }
                let endpointIndices = sgroup.atomIDs.dropFirst().compactMap { atomIndexByID[$0] }
                guard !endpointIndices.isEmpty else { continue }
                state.positionalVariations[beginIndex] = endpointIndices
            }
        }

        if shouldEmitPolymer || shouldEmitDataSgroups {
            var generatedSgroups: [CDKCxSmilesState.SgroupDefinition] = []
            for sgroup in molecule.sgroups
            where sgroup.kind == .structureRepeatUnit || sgroup.kind == .polymer || sgroup.kind == .data {
                switch sgroup.kind {
                case .data:
                    if shouldEmitDataSgroups {
                        generatedSgroups.append(
                            .init(
                                kind: .data,
                                atomIndices: sgroup.atomIDs.compactMap { atomIndexByID[$0] },
                                dataFieldName: sgroup.dataFieldName,
                                dataValue: sgroup.dataValue,
                                dataOperator: sgroup.dataOperator,
                                dataUnit: sgroup.dataUnit,
                                dataTag: sgroup.dataTag)
                        )
                    }
                case .structureRepeatUnit, .polymer:
                    generatedSgroups.append(
                        .init(
                            kind: .polymer,
                            keyword: sgroup.keyword ?? defaultPolymerKeyword(for: sgroup),
                            atomIndices: sgroup.atomIDs.compactMap { atomIndexByID[$0] },
                            subscriptText: sgroup.subscriptText,
                            superscriptText: sgroup.connectivity ?? sgroup.superscriptText,
                            childIndices: sgroup.childGroupIndices)
                    )
                case .extMulticenter, .generic:
                    break
                }
            }
            state.sgroups = generatedSgroups
        }

        return createCxLayer(from: state, atomOutputOrder: atomIDs)
    }

    private func createCxLayer(from state: CDKCxSmilesState, atomOutputOrder: [Int]) -> String? {
        var layers: [String] = []

        if shouldEmitFragmentGroups, !state.fragmentGroups.isEmpty {
            let groups = state.fragmentGroups.map { $0.map(String.init).joined(separator: ".") }.joined(
                separator: ",")
            layers.append("f:\(groups)")
        }

        if shouldEmitAtomLabels,
            let layer = serializeAtomMetadata(state.atomLabels, prefix: "$", suffix: "$")
        {
            layers.append(layer)
        }

        if shouldEmitAtomValues,
            let layer = serializeAtomMetadata(state.atomValues, prefix: "$_AV:", suffix: "$")
        {
            layers.append(layer)
        }

        if shouldEmitEnhancedStereo {
            if state.racemic {
                layers.append("r")
            } else {
                if !state.racemicFragments.isEmpty {
                    layers.append("r:" + state.racemicFragments.map(String.init).joined(separator: ","))
                }
                layers.append(contentsOf: serializeStereoGroups(state.stereoGroups))
            }
        } else if state.racemic {
            layers.append("r")
        } else if !state.racemicFragments.isEmpty {
            layers.append("r:" + state.racemicFragments.map(String.init).joined(separator: ","))
        }

        if shouldEmitCoordinates, let layer = serializeCoordinates(state.atomCoordinates) {
            layers.append(layer)
        }

        if shouldEmitCoordinateBonds,
            let coordinateBonds = state.coordinateBonds,
            !coordinateBonds.isEmpty
        {
            let entries = coordinateBonds.keys.sorted().flatMap { atomIndex in
                coordinateBonds[atomIndex, default: []].sorted().map { "\(atomIndex).\($0)" }
            }
            if !entries.isEmpty {
                layers.append("C:" + entries.joined(separator: ","))
            }
        }

        if shouldEmitMulticenter, !state.positionalVariations.isEmpty {
            let entries = state.positionalVariations.keys.sorted().compactMap { key -> String? in
                guard let values = state.positionalVariations[key], !values.isEmpty else { return nil }
                return "\(key):" + values.map(String.init).joined(separator: ".")
            }
            if !entries.isEmpty {
                layers.append("m:" + entries.joined(separator: ","))
            }
        }

        if shouldEmitLigandOrdering, !state.ligandOrdering.isEmpty {
            let entries = state.ligandOrdering.keys.sorted().compactMap { key -> String? in
                guard let values = state.ligandOrdering[key], !values.isEmpty else { return nil }
                return "\(key):" + values.map(String.init).joined(separator: ".")
            }
            if !entries.isEmpty {
                layers.append("LO:" + entries.joined(separator: ","))
            }
        }

        if shouldEmitPolymer || shouldEmitDataSgroups, !state.sgroups.isEmpty {
            for sgroup in state.sgroups {
                switch sgroup.kind {
                case .polymer:
                    guard shouldEmitPolymer else { continue }
                    let keyword = sgroup.keyword ?? "n"
                    let atomList = sgroup.atomIndices.map(String.init).joined(separator: ",")
                    let subscriptText = escapeCxText(sgroup.subscriptText ?? keyword)
                    let superscriptText = escapeCxText(sgroup.superscriptText ?? "")
                    layers.append("Sg:\(keyword):\(atomList):\(subscriptText):\(superscriptText)")
                case .data:
                    guard shouldEmitDataSgroups else { continue }
                    let atomList = sgroup.atomIndices.map(String.init).joined(separator: ",")
                    let field = escapeCxText(sgroup.dataFieldName ?? "")
                    let value = escapeCxText(sgroup.dataValue ?? "")
                    let op = escapeCxText(sgroup.dataOperator ?? "")
                    let unit = escapeCxText(sgroup.dataUnit ?? "")
                    var layer = "SgD:\(atomList):\(field):\(value)"
                    if !op.isEmpty || !unit.isEmpty {
                        layer += ":\(op)"
                        if !unit.isEmpty {
                            layer += ":\(unit)"
                        }
                    } else {
                        layer += "::"
                    }
                    layers.append(layer)
                }
            }

            let hierarchyEntries = state.sgroups.enumerated().compactMap { index, sgroup -> String? in
                guard !sgroup.childIndices.isEmpty else { return nil }
                return "\(index):" + sgroup.childIndices.sorted().map(String.init).joined(separator: ".")
            }
            if !hierarchyEntries.isEmpty {
                layers.append("SgH:" + hierarchyEntries.joined(separator: ","))
            }
        }

        if shouldEmitRadicals, !state.atomRadicals.isEmpty {
            let grouped = Dictionary(grouping: state.atomRadicals) { $0.value }
            for radicalType in grouped.keys.sorted(by: { $0.rawValue < $1.rawValue }) {
                let atomIndices = grouped[radicalType, default: []]
                    .map(\.key)
                    .sorted()
                    .map(String.init)
                    .joined(separator: ",")
                layers.append("^\(radicalType.rawValue):\(atomIndices)")
            }
        }

        if shouldEmitRGroups, !state.rGroupDefinitions.isEmpty {
            let orderedLabels =
                !state.rGroupOrder.isEmpty ? state.rGroupOrder : state.rGroupDefinitions.keys.sorted()
            let definitions = orderedLabels.compactMap { label -> String? in
                guard let values = state.rGroupDefinitions[label], !values.isEmpty else { return nil }
                return "_\(label)=" + values.map { "{\($0)}" }.joined(separator: ",")
            }
            if !definitions.isEmpty {
                layers.append("RG:" + definitions.joined(separator: ","))
            }
        }

        guard !layers.isEmpty else { return nil }
        return " |" + layers.joined(separator: ",") + "|"
    }

    private var shouldEmitAtomLabels: Bool {
        flavor.contains(.cxsmiles) || flavor.contains(.cxAtomLabel) || flavor.contains(.cxRGroups)
    }

    private var shouldEmitAtomValues: Bool {
        flavor.contains(.cxsmiles) || flavor.contains(.cxAtomValue)
    }

    private var shouldEmitCoordinates: Bool {
        flavor.contains(.cxsmiles) || flavor.contains(.cxCoordinates)
    }

    private func shouldEmitCoordinates(for molecule: Molecule) -> Bool {
        if molecule.cxState?.atomCoordinates.isEmpty == false {
            return true
        }
        if flavor.contains(.cxCoordinates) {
            return molecule.atoms.contains { $0.zPosition != nil }
        }
        return false
    }

    private var shouldEmitMulticenter: Bool {
        flavor.contains(.cxsmiles) || flavor.contains(.cxMulticenter)
    }

    private var shouldEmitPolymer: Bool {
        flavor.contains(.cxsmiles) || flavor.contains(.cxPolymer)
    }

    private var shouldEmitRadicals: Bool {
        flavor.contains(.cxsmiles) || flavor.contains(.cxRadical)
    }

    private var shouldEmitFragmentGroups: Bool {
        flavor.contains(.cxsmiles) || flavor.contains(.cxFragmentGroup)
    }

    private var shouldEmitLigandOrdering: Bool {
        flavor.contains(.cxsmiles) || flavor.contains(.cxLigandOrder)
    }

    private var shouldEmitEnhancedStereo: Bool {
        flavor.contains(.cxsmiles) || flavor.contains(.cxEnhancedStereo)
    }

    private var shouldEmitDataSgroups: Bool {
        flavor.contains(.cxsmiles) || flavor.contains(.cxDataSgroups)
    }

    private var shouldEmitRGroups: Bool {
        flavor.contains(.cxsmiles) || flavor.contains(.cxRGroups)
    }

    private var shouldEmitCoordinateBonds: Bool {
        flavor.contains(.cxsmiles) || flavor.contains(.cxCoordinateBonds)
    }

    private struct EdgeKey: Hashable, Comparable {
        let a: Int
        let b: Int

        init(_ u: Int, _ v: Int) {
            a = min(u, v)
            b = max(u, v)
        }

        static func < (lhs: EdgeKey, rhs: EdgeKey) -> Bool {
            if lhs.a != rhs.a { return lhs.a < rhs.a }
            return lhs.b < rhs.b
        }
    }

    private struct RingClosure {
        let number: Int
        let bondSymbol: String
        let partner: Int
    }

    private struct TreeBuildFrame {
        let atomID: Int
        let neighbors: [Int]
        var nextIndex: Int
    }

    private func buildAdjacency(_ molecule: Molecule) -> [Int: [Int]] {
        var map: [Int: Set<Int>] = [:]
        for atom in molecule.atoms {
            map[atom.id] = []
        }
        for bond in molecule.bonds {
            map[bond.a1, default: []].insert(bond.a2)
            map[bond.a2, default: []].insert(bond.a1)
        }
        return map.mapValues { Array($0).sorted() }
    }

    private func buildBondLookup(_ molecule: Molecule) -> [EdgeKey: Bond] {
        var lookup: [EdgeKey: Bond] = [:]
        for bond in molecule.bonds {
            lookup[EdgeKey(bond.a1, bond.a2)] = bond
        }
        return lookup
    }

    private func connectedComponents(
        _ molecule: Molecule,
        adjacency: [Int: [Int]]
    ) -> [[Int]] {
        var visited: Set<Int> = []
        var components: [[Int]] = []

        for atomID in molecule.atoms.map(\.id).sorted() {
            if visited.contains(atomID) { continue }

            var stack = [atomID]
            var component: [Int] = []
            visited.insert(atomID)

            while let node = stack.popLast() {
                component.append(node)
                for neighbor in adjacency[node, default: []] where !visited.contains(neighbor) {
                    visited.insert(neighbor)
                    stack.append(neighbor)
                }
            }

            component.sort()
            components.append(component)
        }

        components.sort { lhs, rhs in
            (lhs.first ?? Int.max) < (rhs.first ?? Int.max)
        }

        return components
    }

    private func buildTree(
        component: Set<Int>,
        root: Int,
        adjacency: [Int: [Int]]
    ) -> (adjacency: [Int: [Int]], edges: Set<EdgeKey>) {
        var treeAdjacency: [Int: [Int]] = [:]
        var treeEdges: Set<EdgeKey> = []
        var visited: Set<Int> = []

        func componentNeighbors(of atomID: Int) -> [Int] {
            adjacency[atomID, default: []]
                .filter { component.contains($0) }
                .sorted()
        }

        visited.insert(root)
        var stack: [TreeBuildFrame] = [
            TreeBuildFrame(atomID: root, neighbors: componentNeighbors(of: root), nextIndex: 0)
        ]

        while var frame = stack.popLast() {
            guard frame.nextIndex < frame.neighbors.count else {
                continue
            }

            let neighbor = frame.neighbors[frame.nextIndex]
            frame.nextIndex += 1
            if frame.nextIndex < frame.neighbors.count {
                stack.append(frame)
            }

            guard !visited.contains(neighbor) else {
                continue
            }

            visited.insert(neighbor)
            treeAdjacency[frame.atomID, default: []].append(neighbor)
            treeAdjacency[neighbor, default: []].append(frame.atomID)
            treeEdges.insert(EdgeKey(frame.atomID, neighbor))
            stack.append(
                TreeBuildFrame(
                    atomID: neighbor,
                    neighbors: componentNeighbors(of: neighbor),
                    nextIndex: 0))
        }

        for atomID in component {
            treeAdjacency[atomID, default: []].sort()
        }

        return (treeAdjacency, treeEdges)
    }

    private func assignRingNumbers(to ringEdges: Set<EdgeKey>) -> [EdgeKey: Int] {
        var assignments: [EdgeKey: Int] = [:]
        var next = 1

        for edge in ringEdges.sorted() {
            assignments[edge] = next
            next += 1
            if next > 99 { next = 1 }
        }

        return assignments
    }

    private func renderTreeAtom(
        _ atomID: Int,
        parent: Int?,
        treeAdjacency: [Int: [Int]],
        ringClosuresByAtom: [Int: [RingClosure]],
        bondByEdge: [EdgeKey: Bond],
        atomByID: [Int: Atom],
        aromaticBondIDs: Set<Int>,
        atomOutputOrder: inout [Int],
        visited: inout Set<Int>
    ) -> String {
        var out = ""
        var stack: [RenderEvent] = [.atom(id: atomID, parent: parent)]

        while let event = stack.popLast() {
            switch event {
            case .text(let text):
                out += text

            case .atom(let atomID, let parent):
                visited.insert(atomID)
                atomOutputOrder.append(atomID)

                guard let atom = atomByID[atomID] else { continue }

                out += atomToken(atom)

                for closure in ringClosuresByAtom[atomID] ?? [] {
                    out += closure.bondSymbol + ringToken(closure.number)
                }

                let children = treeAdjacency[atomID, default: []]
                    .filter { $0 != parent }
                    .sorted()

                guard let mainChild = children.last else {
                    continue
                }

                var pending: [RenderEvent] = []
                for childID in children.dropLast() {
                    let edge = EdgeKey(atomID, childID)
                    guard let bond = bondByEdge[edge] else { continue }

                    pending.append(.text("("))
                    pending.append(
                        .text(bondToken(bond, atomByID: atomByID, aromaticBondIDs: aromaticBondIDs)))
                    pending.append(.atom(id: childID, parent: atomID))
                    pending.append(.text(")"))
                }

                let edge = EdgeKey(atomID, mainChild)
                if let bond = bondByEdge[edge] {
                    pending.append(
                        .text(bondToken(bond, atomByID: atomByID, aromaticBondIDs: aromaticBondIDs)))
                    pending.append(.atom(id: mainChild, parent: atomID))
                }

                stack.append(contentsOf: pending.reversed())
            }
        }

        return out
    }

    private func ringToken(_ number: Int) -> String {
        if number < 10 {
            return "\(number)"
        }
        return "%\(number)"
    }

    private func atomToken(_ atom: Atom) -> String {
        if cxAtomLabel(for: atom) != nil {
            return "*"
        }

        if atom.queryType != nil || atom.atomList != nil {
            return "*"
        }

        let element = atom.element.trimmingCharacters(in: .whitespacesAndNewlines)
        if element == "*" {
            return "*"
        }

        let plain = plainSymbol(for: atom)
        if !requiresBracket(atom), let plain {
            return plain
        }

        var token = "["

        if flavor.contains(.isomeric), let isotope = atom.isotopeMassNumber {
            token += "\(isotope)"
        }

        token += bracketElementSymbol(for: atom)

        if flavor.contains(.isomeric) {
            switch atom.chirality {
            case .clockwise:
                token += "@"
            case .anticlockwise:
                token += "@@"
            case .none:
                break
            }
        }

        if let h = atom.explicitHydrogenCount {
            token += hydrogenToken(h)
        }

        if atom.charge != 0 {
            token += chargeToken(atom.charge)
        }

        token += "]"
        return token
    }

    private func cxAtomLabel(for atom: Atom) -> String? {
        if let attachmentPoint = atom.attachmentPoint {
            return "_AP\(attachmentPoint)"
        }

        if let rGroupLabel = atom.rGroupLabel {
            if rGroupLabel == 0,
                atom.element.trimmingCharacters(in: .whitespacesAndNewlines).uppercased() == "R"
            {
                return "R"
            }
            return "R\(rGroupLabel)"
        }

        let element = atom.element.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !element.isEmpty, element != "*" else { return nil }

        let uppercase = element.uppercased()
        if plainSymbol(for: atom) != nil {
            return nil
        }
        if uppercase == "R" || uppercase == "R#" {
            return element
        }
        return element
    }

    private func radicalType(for radicalCount: Int?) -> CxRadicalType? {
        switch radicalCount {
        case 1:
            return .monovalent
        case 2:
            return .divalent
        case 3:
            return .trivalent
        default:
            return nil
        }
    }

    private func defaultPolymerKeyword(for sgroup: MoleculeSgroup) -> String {
        switch sgroup.kind {
        case .structureRepeatUnit:
            return "n"
        case .polymer:
            return "gen"
        case .data, .extMulticenter, .generic:
            return "gen"
        }
    }

    private func serializeAtomMetadata(
        _ values: [Int: String],
        prefix: String,
        suffix: String
    ) -> String? {
        guard !values.isEmpty else { return nil }
        let lastIndex = values.keys.max() ?? -1
        guard lastIndex >= 0 else { return nil }
        let entries = (0...lastIndex).map { index in
            values[index].map(escapeCxText) ?? ""
        }
        return prefix + entries.joined(separator: ";") + suffix
    }

    private func serializeCoordinates(_ coordinates: [CxCoordinate]) -> String? {
        guard !coordinates.isEmpty else { return nil }
        let tuples = coordinates.map { coord in
            let x = formatCxNumber(coord.x)
            let y = formatCxNumber(coord.y)
            let z = formatCxNumber(coord.z)
            return "\(x),\(y),\(z)"
        }
        return "(" + tuples.joined(separator: ";") + ")"
    }

    private func serializeStereoGroups(_ stereoGroups: [Int: Int]) -> [String] {
        guard !stereoGroups.isEmpty else { return [] }
        let grouped = Dictionary(grouping: stereoGroups) { $0.value }
        let sortedEntries = grouped.compactMap { encoded, pairs -> (kind: String, atoms: [Int])? in
            guard let decoded = CDKCxSmilesParser.decodeStereoGroup(encoded) else { return nil }
            return (decoded.kind, pairs.map(\.key).sorted())
        }

        guard !sortedEntries.isEmpty else { return [] }

        if sortedEntries.allSatisfy({ $0.kind == "abs" }) {
            return []
        }

        if sortedEntries.count == 1,
            sortedEntries.first?.kind == "and"
        {
            return ["r"]
        }

        let ordered = sortedEntries.sorted { lhs, rhs in
            (lhs.atoms.first ?? Int.max) < (rhs.atoms.first ?? Int.max)
        }

        var nextAndNumber = 1
        var nextOrNumber = 1
        return ordered.compactMap { entry in
            let atoms = entry.atoms.map(String.init).joined(separator: ",")
            switch entry.kind {
            case "abs":
                return atoms.isEmpty ? nil : "a:\(atoms)"
            case "or":
                defer { nextOrNumber += 1 }
                return atoms.isEmpty ? nil : "o\(nextOrNumber):\(atoms)"
            case "and":
                defer { nextAndNumber += 1 }
                return atoms.isEmpty ? nil : "&\(nextAndNumber):\(atoms)"
            default:
                return nil
            }
        }
    }

    private func escapeCxText(_ text: String) -> String {
        var out = ""
        for scalar in text.unicodeScalars {
            switch scalar {
            case ":":
                out += "&#58;"
            case "\n":
                out += "&#10;"
            case "$":
                out += "&#36;"
            case "'":
                out += "&#39;"
            default:
                out += String(scalar)
            }
        }
        return out
    }

    private func formatCxNumber(_ value: Double) -> String {
        if abs(value) < 0.000_000_1 {
            return ""
        }
        let rounded = (value * 1000).rounded() / 1000
        let string = String(rounded)
        if string.hasSuffix(".0") {
            return String(string.dropLast(2))
        }
        return string
    }

    private func requiresBracket(_ atom: Atom) -> Bool {
        if atom.queryType != nil || atom.atomList != nil { return true }
        if atom.charge != 0 { return true }
        if atom.explicitHydrogenCount != nil { return true }
        if atom.radical != nil || atom.radicalType != nil { return true }
        if flavor.contains(.isomeric) {
            if atom.isotopeMassNumber != nil { return true }
            if atom.chirality != .none { return true }
        }
        return plainSymbol(for: atom) == nil
    }

    private func plainSymbol(for atom: Atom) -> String? {
        let upper = atom.element.uppercased()
        if upper == "*" {
            return "*"
        }

        let aromaticEnabled = flavor.contains(.useAromaticSymbols)
        if aromaticEnabled && atom.aromatic {
            switch upper {
            case "B", "C", "N", "O", "P", "S":
                return upper.lowercased()
            case "SE":
                return "se"
            case "AS":
                return "as"
            default:
                return nil
            }
        }

        switch upper {
        case "B", "C", "N", "O", "P", "S", "F", "CL", "BR", "I":
            return normalizedElementSymbol(upper)
        default:
            return nil
        }
    }

    private func bracketElementSymbol(for atom: Atom) -> String {
        let upper = atom.element.uppercased()
        if upper == "*" { return "*" }

        if flavor.contains(.useAromaticSymbols), atom.aromatic {
            switch upper {
            case "B", "C", "N", "O", "P", "S":
                return upper.lowercased()
            case "SE":
                return "se"
            case "AS":
                return "as"
            default:
                break
            }
        }

        return normalizedElementSymbol(upper)
    }

    private func normalizedElementSymbol(_ upper: String) -> String {
        guard let first = upper.first else { return upper }
        let head = String(first)
        let tail = String(upper.dropFirst()).lowercased()
        return head + tail
    }

    private func hydrogenToken(_ count: Int) -> String {
        if count <= 0 {
            return "H0"
        }
        if count == 1 {
            return "H"
        }
        return "H\(count)"
    }

    private func chargeToken(_ charge: Int) -> String {
        if charge == 1 { return "+" }
        if charge == -1 { return "-" }

        let sign = charge > 0 ? "+" : "-"
        let magnitude = abs(charge)
        if magnitude <= 3 {
            return String(repeating: sign, count: magnitude)
        }
        return "\(sign)\(magnitude)"
    }

    private func bondToken(
        _ bond: Bond,
        atomByID: [Int: Atom],
        aromaticBondIDs: Set<Int>
    ) -> String {
        if bond.queryType != nil {
            return "-"
        }

        switch bond.order {
        case .single:
            if flavor.contains(.isomeric) {
                switch bond.stereo {
                case .up, .upReversed:
                    return "/"
                case .down, .downReversed:
                    return "\\"
                case .either, .none:
                    break
                }
            }

            if flavor.contains(.useAromaticSymbols), aromaticBondIDs.contains(bond.id) {
                let leftAromatic = atomByID[bond.a1]?.aromatic == true
                let rightAromatic = atomByID[bond.a2]?.aromatic == true
                if leftAromatic && rightAromatic {
                    return ""
                }
            }

            return ""

        case .double:
            return "="

        case .triple:
            return "#"

        case .aromatic:
            if flavor.contains(.useAromaticSymbols) {
                let leftAromatic = atomByID[bond.a1]?.aromatic == true
                let rightAromatic = atomByID[bond.a2]?.aromatic == true
                if leftAromatic && rightAromatic {
                    return ""
                }
            }
            return ":"
        }
    }
}

/// Factory facade mirroring CDK-style generator construction patterns.
public final class CDKSmilesGeneratorFactory: @unchecked Sendable {
    public static let shared = CDKSmilesGeneratorFactory()

    private init() {}

    public func newSmilesGenerator(flavor: CDKSmiFlavor = .cdkDefault) -> CDKSmilesGenerator {
        CDKSmilesGenerator(flavor: flavor)
    }
}
