import Foundation

/// ABS centres retain their drawn parity. An OR group permits joint inversion of
/// all its centres; independent OR groups vary independently. AND/racemic mixtures
/// are not individual compounds and are deliberately not flattened into a family.
struct MarkushStereoPlan {
    let source: Molecule
    let groups: [[Int]]
    let groupNumbers: [Int]
    let hasMetadata: Bool

    init(_ molecule: Molecule, concrete: Bool = false) throws {
        guard molecule.cxState?.racemic != true, molecule.cxState?.racemicFragments.isEmpty != false else {
            throw CDKRGroupEnumerationError.unsupportedFeature("Racemic mixtures require mixture semantics.")
        }
        var source = molecule
        // Reconcile CX array indices with atom metadata before clearing either representation.
        for (index, encoded) in molecule.cxState?.stereoGroups ?? [:] {
            guard molecule.atoms.indices.contains(index),
                molecule.atoms[index].cxStereoGroup == nil || molecule.atoms[index].cxStereoGroup == encoded
            else { throw CDKRGroupEnumerationError.invalidQuery("Inconsistent enhanced stereo metadata.") }
            source.atoms[index].cxStereoGroup = encoded
        }
        hasMetadata = source.atoms.contains { $0.cxStereoGroup != nil }
        var members: [Int: [Int]] = [:]
        for index in source.atoms.indices {
            let atom = source.atoms[index]
            guard let encoded = atom.cxStereoGroup else { continue }
            guard let group = CDKCxSmilesParser.decodeStereoGroup(encoded),
                atom.chirality == .clockwise || atom.chirality == .anticlockwise
            else {
                throw CDKRGroupEnumerationError.unsupportedFeature(
                    "An enhanced stereo centre is not a supported tetrahedron.")
            }
            guard group.kind == "abs" || (group.kind == "or" && !concrete) else {
                throw CDKRGroupEnumerationError.unsupportedFeature(
                    "The compound must specify one stereoisomer; mixture groups are unsupported.")
            }
            if group.kind == "or" { members[group.number, default: []].append(atom.id) }
            source.atoms[index].cxStereoGroup = nil
        }
        guard members.count <= 16 else {
            throw CDKRGroupEnumerationError.unsupportedFeature("More than 16 independent stereo groups.")
        }
        source.cxState?.stereoGroups = [:]
        self.source = source
        groupNumbers = members.keys.sorted()
        groups = groupNumbers.map { members[$0]! }
    }

    func enumerate(checkpoint: () throws -> Void, _ visit: (Molecule, [Int]) throws -> Bool) throws {
        func expand(_ index: Int, _ graph: Molecule, _ inversions: [Int]) throws -> Bool {
            try checkpoint()
            if index == groups.count { return try visit(graph, inversions) }
            if try !expand(index + 1, graph, inversions) { return false }
            var inverted = graph
            let ids = Set(groups[index])
            for i in inverted.atoms.indices where ids.contains(inverted.atoms[i].id) {
                inverted.atoms[i].chirality = inverted.atoms[i].chirality == .clockwise ? .anticlockwise : .clockwise
            }
            return try expand(index + 1, inverted, inversions + [groupNumbers[index]])
        }
        _ = try expand(0, source, [])
    }
}
