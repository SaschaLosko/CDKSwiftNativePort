import Foundation

// CDXML/CDX double-bond stereo is geometric. Translate the model's directional
// SMILES markers and explicit cis/trans references into 2D geometry before writing.
enum ChemDrawStereo {
    static func prepare(_ source: Molecule) throws -> Molecule {
        var result = source
        var directional = Set<Int>()
        for double in source.bonds where double.order == .double {
            var references = double.stereoReferenceAtomIDs
            var configuration = double.doubleBondStereo
            if configuration == nil, double.stereo == .either {
                let left = source.bonds.first {
                    $0.order == .single && $0.stereo != .none && $0.stereo != .either
                        && ($0.a1 == double.a1 || $0.a2 == double.a1)
                }
                let right = source.bonds.first {
                    $0.order == .single && $0.stereo != .none && $0.stereo != .either
                        && ($0.a1 == double.a2 || $0.a2 == double.a2)
                }
                if let left, let right {
                    references = [
                        left.a1 == double.a1 ? left.a2 : left.a1, double.a1, double.a2,
                        right.a1 == double.a2 ? right.a2 : right.a1,
                    ]
                    configuration =
                        direction(left, from: double.a1) == direction(right, from: double.a2) ? .cis : .trans
                    directional.formUnion([left.id, right.id])
                }
            }
            guard let configuration else {
                if source.coordinatesAreGenerated == true,
                    let index = result.bonds.firstIndex(where: { $0.id == double.id })
                {
                    result.bonds[index].stereo = .either
                }
                continue
            }
            guard let refs = references, refs.count == 4,
                let first = result.atom(id: refs[0]), let a = result.atom(id: refs[1]),
                let b = result.atom(id: refs[2]), let last = result.atom(id: refs[3]),
                Set([refs[1], refs[2]]) == Set([double.a1, double.a2])
            else {
                throw ChemError.unsupported("ChemDraw export requires valid double-bond stereo references.")
            }
            let dx = b.position.x - a.position.x, dy = b.position.y - a.position.y
            let crossA = dx * (first.position.y - a.position.y) - dy * (first.position.x - a.position.x)
            let crossB = dx * (last.position.y - a.position.y) - dy * (last.position.x - a.position.x)
            guard abs(crossA) > 0.00001, abs(crossB) > 0.00001 else {
                throw ChemError.unsupported(
                    "ChemDraw export cannot preserve double-bond stereo with collinear coordinates. Generate a new 2D layout first."
                )
            }
            if (crossA * crossB > 0) != (configuration == .cis) {
                var branch: Set<Int> = [b.id]
                var queue = [b.id]
                while let current = queue.popLast() {
                    for edge in result.bonds where edge.id != double.id && (edge.a1 == current || edge.a2 == current) {
                        let other = edge.a1 == current ? edge.a2 : edge.a1
                        if branch.insert(other).inserted { queue.append(other) }
                    }
                }
                guard !branch.contains(a.id) else {
                    throw ChemError.unsupported(
                        "ChemDraw export cannot reconcile this ring's double-bond stereo and drawing coordinates.")
                }
                let denominator = dx * dx + dy * dy
                for index in result.atoms.indices where branch.contains(result.atoms[index].id) {
                    let px = result.atoms[index].position.x - a.position.x
                    let py = result.atoms[index].position.y - a.position.y
                    let projection = (px * dx + py * dy) / denominator
                    result.atoms[index].position.x = a.position.x + 2 * projection * dx - px
                    result.atoms[index].position.y = a.position.y + 2 * projection * dy - py
                }
                for index in result.bonds.indices
                where branch.contains(result.bonds[index].a1) && branch.contains(result.bonds[index].a2) {
                    result.bonds[index].stereo = inverted(result.bonds[index].stereo)
                }
            }
            if let index = result.bonds.firstIndex(where: { $0.id == double.id }) { result.bonds[index].stereo = .none }
        }
        for index in result.bonds.indices where directional.contains(result.bonds[index].id) {
            result.bonds[index].stereo = .none
        }
        try assignTetrahedralStereo(to: &result)
        return result
    }

    private static func assignTetrahedralStereo(to molecule: inout Molecule) throws {
        let centers = molecule.atoms.filter { $0.chirality != .none && $0.ligandOrderingAtomIDs?.count == 4 }
        let centerIDs = Set(centers.map(\.id))
        // Preserve drawn stereo when no explicit parity ordering exists. Otherwise
        // choose an outgoing bond and calculate its wedge from ordered ligand vectors.
        for center in centers {
            for index in molecule.bonds.indices {
                let bond = molecule.bonds[index]
                let owner: Int? =
                    switch bond.stereo {
                    case .up, .down: bond.a1
                    case .upReversed, .downReversed: bond.a2
                    case .none, .either: nil
                    }
                if owner == center.id { molecule.bonds[index].stereo = .none }
            }
        }
        for center in centers {
            let candidates = molecule.bonds.indices.filter {
                let bond = molecule.bonds[$0]
                return bond.order == .single && bond.stereo == .none
                    && (bond.a1 == center.id || bond.a2 == center.id)
            }.sorted { lhs, rhs in
                let l = molecule.bonds[lhs], r = molecule.bonds[rhs]
                let ln = l.a1 == center.id ? l.a2 : l.a1
                let rn = r.a1 == center.id ? r.a2 : r.a1
                if centerIDs.contains(ln) != centerIDs.contains(rn) { return !centerIDs.contains(ln) }
                return l.id < r.id
            }
            var assigned = false
            for index in candidates {
                let bond = molecule.bonds[index]
                let neighbor = bond.a1 == center.id ? bond.a2 : bond.a1
                let vectors = center.ligandOrderingAtomIDs!.compactMap { id -> SIMD3<Double>? in
                    guard let atom = molecule.atom(id: id) else { return nil }
                    return SIMD3(
                        Double(atom.position.x - center.position.x),
                        Double(atom.position.y - center.position.y), id == neighbor ? 1 : 0)
                }
                guard vectors.count == 4 else { continue }
                let a = vectors[0] - vectors[3], b = vectors[1] - vectors[3], c = vectors[2] - vectors[3]
                let determinant =
                    a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x) + a.z * (b.x * c.y - b.y * c.x)
                guard abs(determinant) > 0.00001 else { continue }
                let up = (determinant > 0) == (center.chirality == .clockwise)
                molecule.bonds[index].stereo =
                    bond.a1 == center.id ? (up ? .up : .down) : (up ? .upReversed : .downReversed)
                assigned = true
                break
            }
            guard assigned else {
                throw ChemError.unsupported(
                    "ChemDraw export cannot preserve this tetrahedral stereocenter in the current 2D layout.")
            }
        }
    }

    private static func direction(_ bond: Bond, from atom: Int) -> Int {
        let sign =
            switch bond.stereo {
            case .up, .downReversed: 1;
            default: -1
            }
        return bond.a1 == atom ? sign : -sign
    }

    private static func inverted(_ stereo: BondStereo) -> BondStereo {
        switch stereo {
        case .up: .down
        case .down: .up
        case .upReversed: .downReversed
        case .downReversed: .upReversed
        case .none, .either: stereo
        }
    }
}
