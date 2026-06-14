import Foundation
#if canImport(CoreGraphics)
import CoreGraphics
#endif

public struct CDKModelBuilder3DOptions: Hashable, Sendable {
    public var targetBondLength: Double
    public var ringPuckerAmplitude: Double
    public var stereoDepth: Double
    public var tetrahedralBranchDepth: Double
    public var coordinateTolerance: Double

    public init(
        targetBondLength: Double = 1.50,
        ringPuckerAmplitude: Double = 0.42,
        stereoDepth: Double = 0.72,
        tetrahedralBranchDepth: Double = 0.36,
        coordinateTolerance: Double = 0.000001
    ) {
        self.targetBondLength = targetBondLength
        self.ringPuckerAmplitude = ringPuckerAmplitude
        self.stereoDepth = stereoDepth
        self.tetrahedralBranchDepth = tetrahedralBranchDepth
        self.coordinateTolerance = coordinateTolerance
    }
}

/// CDK-inspired 3D coordinate generation for viewer scene conversion.
///
/// The upstream Java CDK `ModelBuilder3D` creates a 3D model with MM2-style
/// parameters, ring templates, Z-matrix chain placement, and tetrahedral ligand
/// placement. This Swift port keeps the same default rendering contract:
/// non-flat coordinates are preserved, while flat 2D structures are converted
/// into a generated 3D coordinate set before rendering.
public enum CDKModelBuilder3D {
    public static func hasMeaningful3DCoordinates(
        _ molecule: Molecule,
        tolerance: Double = CDKModelBuilder3DOptions().coordinateTolerance
    ) -> Bool {
        molecule.atoms.contains { atom in
            guard let z = atom.zPosition else { return false }
            return abs(z) > tolerance
        }
    }

    public static func generate3DCoordinates(
        molecule: Molecule,
        options: CDKModelBuilder3DOptions = CDKModelBuilder3DOptions()
    ) -> Molecule {
        guard !molecule.atoms.isEmpty else { return molecule }

        var result = molecule
        let graph = MolecularGraph(molecule: molecule)
        let atomByID = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0) })
        let scale = averageBondScale(molecule: molecule, options: options)
        let componentCenters = Dictionary(uniqueKeysWithValues: graph.connectedComponents().map { component in
            (component, centroid(for: component, atoms: atomByID))
        })
        var zByAtomID = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, atomZ($0)) })

        applyRingPuckering(
            molecule: molecule,
            graph: graph,
            options: options,
            zByAtomID: &zByAtomID)
        applyTetrahedralBranchDepth(
            molecule: molecule,
            graph: graph,
            options: options,
            zByAtomID: &zByAtomID)
        applyStereochemicalDepth(
            molecule: molecule,
            graph: graph,
            options: options,
            zByAtomID: &zByAtomID)

        centerZByConnectedComponent(graph: graph, zByAtomID: &zByAtomID)

        for index in result.atoms.indices {
            let atom = result.atoms[index]
            let component = graph.component(containing: atom.id)
            let center = component.flatMap { componentCenters[$0] } ?? .zero
            let position = generatedXYPosition(
                for: atom,
                graph: graph,
                componentCenter: center,
                scale: scale)
            result.atoms[index].position = position
            result.atoms[index].zPosition = zByAtomID[atom.id] ?? 0
        }

        return result
    }
}

private struct ModelBuilderPoint: Hashable {
    var x: Double
    var y: Double

    static let zero = ModelBuilderPoint(x: 0, y: 0)
}

private struct ModelBuilderEdgeKey: Hashable {
    let a: Int
    let b: Int

    init(_ lhs: Int, _ rhs: Int) {
        a = min(lhs, rhs)
        b = max(lhs, rhs)
    }
}

private struct MolecularGraph {
    let atomIDs: [Int]
    let neighborsByAtomID: [Int: [Int]]
    let bondByEdge: [ModelBuilderEdgeKey: Bond]
    private let componentByAtomID: [Int: Set<Int>]

    init(molecule: Molecule) {
        atomIDs = molecule.atoms.map(\.id).sorted()
        var neighborSets: [Int: Set<Int>] = [:]
        var bondByEdge: [ModelBuilderEdgeKey: Bond] = [:]
        for atomID in atomIDs {
            neighborSets[atomID] = []
        }
        for bond in molecule.bonds {
            neighborSets[bond.a1, default: []].insert(bond.a2)
            neighborSets[bond.a2, default: []].insert(bond.a1)
            bondByEdge[ModelBuilderEdgeKey(bond.a1, bond.a2)] = bond
        }
        neighborsByAtomID = neighborSets.mapValues { $0.sorted() }
        self.bondByEdge = bondByEdge

        let components = MolecularGraph.makeConnectedComponents(atomIDs: atomIDs, neighborsByAtomID: neighborsByAtomID)
        var componentByAtomID: [Int: Set<Int>] = [:]
        for component in components {
            for atomID in component {
                componentByAtomID[atomID] = component
            }
        }
        self.componentByAtomID = componentByAtomID
    }

    func neighbors(of atomID: Int) -> [Int] {
        neighborsByAtomID[atomID] ?? []
    }

    func bond(between lhs: Int, and rhs: Int) -> Bond? {
        bondByEdge[ModelBuilderEdgeKey(lhs, rhs)]
    }

    func connectedComponents() -> [Set<Int>] {
        MolecularGraph.makeConnectedComponents(atomIDs: atomIDs, neighborsByAtomID: neighborsByAtomID)
    }

    func component(containing atomID: Int) -> Set<Int>? {
        componentByAtomID[atomID]
    }

    func isConnected(_ lhs: Int, to rhs: Int, removing removedEdge: ModelBuilderEdgeKey) -> Bool {
        guard lhs != rhs else { return true }
        var seen: Set<Int> = [lhs]
        var queue: [Int] = [lhs]
        var index = 0
        while index < queue.count {
            let current = queue[index]
            index += 1
            for next in neighbors(of: current) {
                if ModelBuilderEdgeKey(current, next) == removedEdge { continue }
                if next == rhs { return true }
                guard seen.insert(next).inserted else { continue }
                queue.append(next)
            }
        }
        return false
    }

    func fragment(from seed: Int, excluding stop: Int, removedEdge: ModelBuilderEdgeKey) -> Set<Int> {
        var seen: Set<Int> = [seed]
        var queue: [Int] = [seed]
        var index = 0
        while index < queue.count {
            let current = queue[index]
            index += 1
            for next in neighbors(of: current) {
                if next == stop { continue }
                if ModelBuilderEdgeKey(current, next) == removedEdge { continue }
                guard seen.insert(next).inserted else { continue }
                queue.append(next)
            }
        }
        return seen
    }

    private static func makeConnectedComponents(
        atomIDs: [Int],
        neighborsByAtomID: [Int: [Int]]
    ) -> [Set<Int>] {
        var components: [Set<Int>] = []
        var seen: Set<Int> = []
        for seed in atomIDs where !seen.contains(seed) {
            var component: Set<Int> = [seed]
            var queue: [Int] = [seed]
            var index = 0
            seen.insert(seed)
            while index < queue.count {
                let current = queue[index]
                index += 1
                for next in neighborsByAtomID[current, default: []] where !seen.contains(next) {
                    seen.insert(next)
                    component.insert(next)
                    queue.append(next)
                }
            }
            components.append(component)
        }
        return components
    }
}

private extension CDKModelBuilder3D {
    static func atomZ(_ atom: Atom) -> Double {
        guard let z = atom.zPosition, z.isFinite else { return 0 }
        return z
    }

    static func averageBondScale(
        molecule: Molecule,
        options: CDKModelBuilder3DOptions
    ) -> Double {
        let atomByID = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0) })
        let distances = molecule.bonds.compactMap { bond -> Double? in
            guard let lhs = atomByID[bond.a1], let rhs = atomByID[bond.a2] else { return nil }
            let dx = Double(rhs.position.x - lhs.position.x)
            let dy = Double(rhs.position.y - lhs.position.y)
            let distance = sqrt(dx * dx + dy * dy)
            return distance > options.coordinateTolerance ? distance : nil
        }
        guard !distances.isEmpty else { return 1 }
        let average = distances.reduce(0, +) / Double(distances.count)
        guard average > options.coordinateTolerance else { return 1 }
        return options.targetBondLength / average
    }

    static func centroid(for component: Set<Int>, atoms: [Int: Atom]) -> ModelBuilderPoint {
        let points = component.compactMap { atoms[$0]?.position }
        guard !points.isEmpty else { return .zero }
        let sum = points.reduce(CGPoint.zero) { partial, point in
            CGPoint(x: partial.x + point.x, y: partial.y + point.y)
        }
        return ModelBuilderPoint(
            x: Double(sum.x) / Double(points.count),
            y: Double(sum.y) / Double(points.count))
    }

    static func generatedXYPosition(
        for atom: Atom,
        graph: MolecularGraph,
        componentCenter: ModelBuilderPoint,
        scale: Double
    ) -> CGPoint {
        let source = ModelBuilderPoint(x: Double(atom.position.x), y: Double(atom.position.y))
        let dx = (source.x - componentCenter.x) * scale
        let dy = (source.y - componentCenter.y) * scale
        if abs(dx) > 0.000001 || abs(dy) > 0.000001 {
            return CGPoint(x: dx + componentCenter.x, y: dy + componentCenter.y)
        }

        let component = graph.component(containing: atom.id) ?? [atom.id]
        guard component.count > 1 else {
            return atom.position
        }

        let ordered = component.sorted()
        guard let offset = ordered.firstIndex(of: atom.id) else { return atom.position }
        let angle = (Double(offset) / Double(max(1, ordered.count))) * 2.0 * Double.pi
        let radius = max(1.0, optionsFallbackBondLength(graph: graph, component: component))
        return CGPoint(
            x: componentCenter.x + cos(angle) * radius,
            y: componentCenter.y + sin(angle) * radius)
    }

    static func optionsFallbackBondLength(graph: MolecularGraph, component: Set<Int>) -> Double {
        max(1.0, Double(component.count) * 0.18 + Double(graph.bondByEdge.count) * 0.01)
    }

    static func applyRingPuckering(
        molecule: Molecule,
        graph: MolecularGraph,
        options: CDKModelBuilder3DOptions,
        zByAtomID: inout [Int: Double]
    ) {
        let rings = molecule.simpleCycles(maxSize: 12)
        guard !rings.isEmpty else { return }
        var contributionByAtomID: [Int: [Double]] = [:]
        for ring in rings where ring.count >= 5 {
            guard !isPlanarPiRing(ring, molecule: molecule, graph: graph) else { continue }
            let amplitude = options.ringPuckerAmplitude * ringPuckerScale(ringSize: ring.count)
            for index in ring.indices {
                let phase = (2.0 * Double.pi * Double(index)) / Double(ring.count)
                let z = sin(phase) * amplitude
                contributionByAtomID[ring[index], default: []].append(z)
            }
        }

        for (atomID, contributions) in contributionByAtomID {
            guard !contributions.isEmpty else { continue }
            let average = contributions.reduce(0, +) / Double(contributions.count)
            addDepth(average, to: atomID, zByAtomID: &zByAtomID, tolerance: options.coordinateTolerance)
        }
    }

    static func isPlanarPiRing(_ ring: [Int], molecule: Molecule, graph: MolecularGraph) -> Bool {
        let atomByID = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0) })
        let aromaticAtoms = ring.allSatisfy { atomByID[$0]?.aromatic ?? false }
        if aromaticAtoms { return true }

        var orders: [BondOrder] = []
        for index in ring.indices {
            let lhs = ring[index]
            let rhs = ring[(index + 1) % ring.count]
            guard let bond = graph.bond(between: lhs, and: rhs) else { return false }
            if bond.order == .aromatic { return true }
            orders.append(bond.order)
        }
        guard orders.count == ring.count, ring.count % 2 == 0 else { return false }
        let startsSingle = orders.enumerated().allSatisfy { index, order in
            order == (index % 2 == 0 ? .single : .double)
        }
        let startsDouble = orders.enumerated().allSatisfy { index, order in
            order == (index % 2 == 0 ? .double : .single)
        }
        let doubleCount = orders.filter { $0 == .double }.count
        return doubleCount >= ring.count / 2 && (startsSingle || startsDouble)
    }

    static func ringPuckerScale(ringSize: Int) -> Double {
        switch ringSize {
        case 5: 0.70
        case 6: 1.00
        case 7: 0.85
        default: 0.55
        }
    }

    static func applyTetrahedralBranchDepth(
        molecule: Molecule,
        graph: MolecularGraph,
        options: CDKModelBuilder3DOptions,
        zByAtomID: inout [Int: Double]
    ) {
        let atomByID = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0) })
        for atom in molecule.atoms {
            let neighbors = graph.neighbors(of: atom.id)
            guard neighbors.count >= 3, likelyTetrahedral(atom: atom, molecule: molecule, graph: graph) else { continue }
            let sortedNeighbors = neighbors.sorted { lhs, rhs in
                angle(from: atom.id, to: lhs, atomByID: atomByID) < angle(from: atom.id, to: rhs, atomByID: atomByID)
            }
            let offsets = tetrahedralOffsets(count: sortedNeighbors.count, amplitude: options.tetrahedralBranchDepth)
            for (neighborID, offset) in zip(sortedNeighbors, offsets) {
                guard shouldApplyBranchDepth(atomID: neighborID, centerID: atom.id, graph: graph) else { continue }
                addDepth(offset, to: neighborID, zByAtomID: &zByAtomID, tolerance: options.coordinateTolerance)
            }
        }
    }

    static func likelyTetrahedral(atom: Atom, molecule: Molecule, graph: MolecularGraph) -> Bool {
        if atom.aromatic { return false }
        let bonds = graph.neighbors(of: atom.id).compactMap { graph.bond(between: atom.id, and: $0) }
        guard !bonds.isEmpty else { return false }
        if bonds.contains(where: { $0.order == .double || $0.order == .triple || $0.order == .aromatic }) {
            return false
        }
        let element = atom.element.trimmingCharacters(in: .whitespacesAndNewlines).uppercased()
        _ = molecule
        return ["C", "N", "P", "S", "SI"].contains(element)
    }

    static func angle(from centerID: Int, to neighborID: Int, atomByID: [Int: Atom]) -> Double {
        guard let center = atomByID[centerID]?.position, let neighbor = atomByID[neighborID]?.position else {
            return 0
        }
        return atan2(Double(neighbor.y - center.y), Double(neighbor.x - center.x))
    }

    static func tetrahedralOffsets(count: Int, amplitude: Double) -> [Double] {
        switch count {
        case 0: []
        case 1: [amplitude]
        case 2: [amplitude, -amplitude]
        case 3: [amplitude, -amplitude, 0]
        default: [amplitude, -amplitude, amplitude * 0.5, -amplitude * 0.5] + Array(repeating: 0, count: count - 4)
        }
    }

    static func shouldApplyBranchDepth(atomID: Int, centerID: Int, graph: MolecularGraph) -> Bool {
        let removedEdge = ModelBuilderEdgeKey(atomID, centerID)
        return !graph.isConnected(atomID, to: centerID, removing: removedEdge)
    }

    static func applyStereochemicalDepth(
        molecule: Molecule,
        graph: MolecularGraph,
        options: CDKModelBuilder3DOptions,
        zByAtomID: inout [Int: Double]
    ) {
        for bond in molecule.bonds {
            guard let depth = stereochemicalDepth(for: bond, options: options) else { continue }
            let anchorID = stereochemicalAnchor(for: bond)
            let ligandID = stereochemicalLigand(for: bond)
            let removedEdge = ModelBuilderEdgeKey(anchorID, ligandID)
            let targetIDs: Set<Int>
            if graph.isConnected(ligandID, to: anchorID, removing: removedEdge) {
                targetIDs = [ligandID]
            } else {
                targetIDs = graph.fragment(from: ligandID, excluding: anchorID, removedEdge: removedEdge)
            }
            for atomID in targetIDs {
                addDepth(depth, to: atomID, zByAtomID: &zByAtomID, tolerance: options.coordinateTolerance)
            }
        }

        for atom in molecule.atoms where atom.chirality != .none {
            let hasDirectionalBond = graph.neighbors(of: atom.id)
                .compactMap { graph.bond(between: atom.id, and: $0) }
                .contains { directionalStereo($0.stereo) }
            guard !hasDirectionalBond else { continue }
            let neighbors = (atom.ligandOrderingAtomIDs ?? graph.neighbors(of: atom.id)).filter {
                graph.neighbors(of: atom.id).contains($0)
            }
            guard neighbors.count >= 3 else { continue }
            let firstSign: Double = atom.chirality == .clockwise ? 1 : -1
            addDepth(firstSign * options.stereoDepth, to: neighbors[0], zByAtomID: &zByAtomID, tolerance: options.coordinateTolerance)
            addDepth(-firstSign * options.stereoDepth, to: neighbors[1], zByAtomID: &zByAtomID, tolerance: options.coordinateTolerance)
        }
    }

    static func stereochemicalDepth(for bond: Bond, options: CDKModelBuilder3DOptions) -> Double? {
        switch bond.stereo {
        case .up, .upReversed:
            return options.stereoDepth
        case .down, .downReversed:
            return -options.stereoDepth
        case .none, .either:
            return nil
        }
    }

    static func stereochemicalAnchor(for bond: Bond) -> Int {
        switch bond.stereo {
        case .up, .down, .either, .none:
            return bond.a1
        case .upReversed, .downReversed:
            return bond.a2
        }
    }

    static func stereochemicalLigand(for bond: Bond) -> Int {
        switch bond.stereo {
        case .up, .down, .either, .none:
            return bond.a2
        case .upReversed, .downReversed:
            return bond.a1
        }
    }

    static func directionalStereo(_ stereo: BondStereo) -> Bool {
        switch stereo {
        case .up, .down, .upReversed, .downReversed:
            return true
        case .none, .either:
            return false
        }
    }

    static func addDepth(
        _ depth: Double,
        to atomID: Int,
        zByAtomID: inout [Int: Double],
        tolerance: Double
    ) {
        guard abs(depth) > tolerance else { return }
        let current = zByAtomID[atomID] ?? 0
        if abs(current) <= tolerance {
            zByAtomID[atomID] = depth
        } else if current.sign == depth.sign {
            zByAtomID[atomID] = (current + depth) * 0.5
        } else if abs(depth) > abs(current) {
            zByAtomID[atomID] = depth
        }
    }

    static func centerZByConnectedComponent(graph: MolecularGraph, zByAtomID: inout [Int: Double]) {
        for component in graph.connectedComponents() {
            let values = component.map { zByAtomID[$0] ?? 0 }
            guard !values.isEmpty else { continue }
            let average = values.reduce(0, +) / Double(values.count)
            guard abs(average) > 0.000001 else { continue }
            for atomID in component {
                zByAtomID[atomID] = (zByAtomID[atomID] ?? 0) - average
            }
        }
    }
}
