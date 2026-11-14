import Foundation

/// CTfile parity uses output atom order (hydrogen last), independently of bond
/// order. Double-bond geometry is explicit on export; synthetic coordinates
/// must not turn an unspecified alkene into a specified isomer.
enum CDKMDLStereo {
    static func attachParityCarriers(to molecule: inout Molecule) {
        for index in molecule.atoms.indices where molecule.atoms[index].chirality != .none {
            let atom = molecule.atoms[index]
            let neighbors = molecule.neighbors(of: atom.id).sorted {
                let leftH = molecule.atom(id: $0)?.element == "H"
                let rightH = molecule.atom(id: $1)?.element == "H"
                return leftH == rightH ? $0 < $1 : !leftH
            }
            if neighbors.count == 3 {
                molecule.atoms[index].ligandOrderingAtomIDs = neighbors + [atom.id]
            } else if neighbors.count == 4 {
                molecule.atoms[index].ligandOrderingAtomIDs = neighbors
            }
        }
    }

    /// Read defined 2D alkene geometry before generating any layout. Crossed
    /// double bonds and adjacent wavy bonds explicitly suppress assignment.
    static func attachDoubleBondGeometry(to molecule: inout Molecule) {
        guard molecule.atoms.allSatisfy({ abs($0.zPosition ?? 0) < 1e-8 }) else { return }
        for index in molecule.bonds.indices {
            var bond = molecule.bonds[index]
            guard bond.order == .double, bond.stereo == .none,
                bond.doubleBondStereo == nil,
                !molecule.bonds(forAtom: bond.a1).contains(where: { $0.stereo == .either }),
                !molecule.bonds(forAtom: bond.a2).contains(where: { $0.stereo == .either }),
                let left = molecule.neighbors(of: bond.a1).filter({ $0 != bond.a2 }).sorted().first,
                let right = molecule.neighbors(of: bond.a2).filter({ $0 != bond.a1 }).sorted().first
            else { continue }
            bond.stereoReferenceAtomIDs = [left, bond.a1, bond.a2, right]
            guard let value = geometry(bond, in: molecule) else { continue }
            bond.doubleBondStereo = value
            molecule.bonds[index] = bond
        }
    }

    static func parity(_ atom: Atom, in molecule: Molecule, outputIDs: [Int]) -> Int {
        guard atom.chirality != .none else { return 0 }
        var parity = atom.chirality == .clockwise ? 1 : 2
        guard let references = atom.ligandOrderingAtomIDs, references.count == 4 else { return parity }
        let order = Dictionary(uniqueKeysWithValues: outputIDs.enumerated().map { ($1, $0) })
        let ranks = references.map { id in
            id == atom.id || molecule.atom(id: id)?.element == "H" ? outputIDs.count + 1 : order[id] ?? -1
        }
        for i in 0..<4 { for j in (i + 1)..<4 where ranks[i] > ranks[j] { parity ^= 3 } }
        return parity
    }

    static func preparedForWriting(_ source: Molecule) throws -> Molecule {
        var molecule = source
        let specified = molecule.bonds.filter { $0.doubleBondStereo != nil }
        if !specified.isEmpty, specified.contains(where: { geometry($0, in: molecule) == nil }) {
            molecule = Depiction2DGenerator.generate(for: molecule)
            molecule.coordinatesAreGenerated = true
        }
        for bond in specified {
            guard let actual = geometry(bond, in: molecule), let desired = bond.doubleBondStereo else {
                throw ChemError.unsupported(
                    "Specified alkene stereochemistry cannot be represented by this CTfile layout.")
            }
            if actual == desired { continue }
            // Reflect the whole component beyond a bridge, preserving its internal
            // shape. Cyclic or degenerate layouts require a separate constraint solver.
            var component: Set<Int> = [bond.a2]
            var pending = [bond.a2]
            while let current = pending.popLast() {
                for edge in molecule.bonds(forAtom: current) where edge.id != bond.id {
                    let next = edge.a1 == current ? edge.a2 : edge.a1
                    if component.insert(next).inserted { pending.append(next) }
                }
            }
            guard !component.contains(bond.a1),
                let a = molecule.atom(id: bond.a1)?.position,
                let b = molecule.atom(id: bond.a2)?.position
            else {
                throw ChemError.unsupported("Inconsistent cyclic alkene geometry cannot be exported losslessly.")
            }
            let dx = b.x - a.x, dy = b.y - a.y
            let length2 = dx * dx + dy * dy
            guard length2 > 1e-12 else { throw ChemError.parseFailed("Degenerate alkene coordinates.") }
            for index in molecule.atoms.indices where component.contains(molecule.atoms[index].id) {
                let p = molecule.atoms[index].position
                let projection = ((p.x - a.x) * dx + (p.y - a.y) * dy) / length2
                molecule.atoms[index].position = CGPoint(
                    x: 2 * (a.x + projection * dx) - p.x,
                    y: 2 * (a.y + projection * dy) - p.y)
            }
            for index in molecule.bonds.indices
            where component.contains(molecule.bonds[index].a1) && component.contains(molecule.bonds[index].a2) {
                switch molecule.bonds[index].stereo {
                case .up: molecule.bonds[index].stereo = .down
                case .down: molecule.bonds[index].stereo = .up
                case .upReversed: molecule.bonds[index].stereo = .downReversed
                case .downReversed: molecule.bonds[index].stereo = .upReversed
                default: break
                }
            }
        }
        for bond in specified where geometry(bond, in: molecule) != bond.doubleBondStereo {
            throw ChemError.unsupported("CTfile layout does not preserve all specified alkene configurations.")
        }
        for index in molecule.bonds.indices where molecule.bonds[index].doubleBondStereo != nil {
            // SMILES uses `either` as a depiction hint for ordered E/Z. CTfile
            // uses it for genuinely unknown geometry, so it must be cleared.
            molecule.bonds[index].stereo = .none
        }
        if molecule.coordinatesAreGenerated == true {
            for index in molecule.bonds.indices {
                let bond = molecule.bonds[index]
                if bond.order == .double, bond.doubleBondStereo == nil,
                    molecule.neighbors(of: bond.a1).count > 1, molecule.neighbors(of: bond.a2).count > 1
                {
                    molecule.bonds[index].stereo = .either
                }
            }
        }
        return molecule
    }

    private static func geometry(_ bond: Bond, in molecule: Molecule) -> DoubleBondStereo? {
        guard let ids = bond.stereoReferenceAtomIDs, ids.count == 4, Set(ids).count == 4,
            Set([ids[1], ids[2]]) == Set([bond.a1, bond.a2]),
            molecule.bond(between: ids[0], and: ids[1]) != nil,
            molecule.bond(between: ids[2], and: ids[3]) != nil
        else { return nil }
        let points = ids.compactMap { molecule.atom(id: $0)?.position }
        guard points.count == 4 else { return nil }
        let a = points[1], b = points[2]
        let left = (b.x - a.x) * (points[0].y - a.y) - (b.y - a.y) * (points[0].x - a.x)
        let right = (b.x - a.x) * (points[3].y - a.y) - (b.y - a.y) * (points[3].x - a.x)
        guard abs(left) > 1e-8, abs(right) > 1e-8 else { return nil }
        return left * right > 0 ? .cis : .trans
    }
}
