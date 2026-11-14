import Foundation

struct MCSStereo {
    let queryTetrahedra: [CDKQueryAtomStereo]
    let targetTetrahedra: [CDKQueryAtomStereo]
    let queryAlkenes: [QueryStereoSupport.Alkene]
    let targetAlkenes: [QueryStereoSupport.Alkene]

    init(query: Molecule, target: Molecule) throws {
        queryTetrahedra = try Self.tetrahedra(query)
        targetTetrahedra = try Self.tetrahedra(target)
        queryAlkenes = try QueryStereoSupport.alkenes(query)
        targetAlkenes = try QueryStereoSupport.alkenes(target)
    }

    private static func tetrahedra(_ source: Molecule) throws -> [CDKQueryAtomStereo] {
        var molecule = source
        // MCS compares relative configurations, not family/mixture membership.
        // Global racemic metadata permits one joint inversion; explicit groups override it.
        if let state = molecule.cxState {
            var next = -1
            var component: [Int: Int] = [:]
            var index = 0
            for atom in molecule.atoms where component[atom.id] == nil {
                var queue = [atom.id]; component[atom.id] = index; var cursor = 0
                while cursor < queue.count {
                    for neighbor in molecule.neighbors(of: queue[cursor]) where component[neighbor] == nil {
                        component[neighbor] = index; queue.append(neighbor)
                    }
                    cursor += 1
                }
                index += 1
            }
            var groupByFragment: [Int: Int] = [:]
            for fragment in state.racemicFragments.sorted() { groupByFragment[fragment] = next; next -= 1 }
            let global = next
            // Use arbitrary distinct group IDs after extracting validated ordinary elements.
            molecule.cxState?.racemic = false
            molecule.cxState?.racemicFragments = []
            var elements = try QueryStereoSupport.elements(molecule, explicit: [])
            for i in elements.indices where elements[i].relativeGroup == nil {
                if let encoded = molecule.atom(id: elements[i].focus)?.cxStereoGroup,
                    CDKCxSmilesParser.decodeStereoGroup(encoded)?.kind == "abs"
                {
                    continue
                }
                if let fragment = component[elements[i].focus], let group = groupByFragment[fragment] {
                    elements[i].relativeGroup = group
                } else if state.racemic {
                    elements[i].relativeGroup = global
                }
            }
            return elements
        }
        return try QueryStereoSupport.elements(molecule, explicit: [])
    }

    func accepts(mapping: [Int: Int], commonBonds: [Int: Int], query q: MCSGraph, target t: MCSGraph) throws -> Bool {
        var left: [CDKQueryAtomStereo] = []
        var right: [CDKQueryAtomStereo] = []
        for original in queryTetrahedra {
            guard let focus = mapping[original.focus],
                let other = targetTetrahedra.first(where: { $0.focus == focus })
            else { continue }
            let neighbors = q.edges[original.focus, default: [:]].filter { neighbor, bond in
                commonBonds[bond.id] != nil && mapping[neighbor].map { t.edges[focus]?[$0] != nil } == true
            }.keys
            guard neighbors.count >= 3 else { continue }
            // A fourth unmapped ligand has the one remaining position, including an implicit H.
            var map = mapping
            let fixed = Set(neighbors.compactMap { map[$0] })
            let remainingTarget = other.carriers.filter { !fixed.contains($0) }
            let remainingQuery = original.carriers.filter { !neighbors.contains($0) }
            if remainingQuery.count == 1 && remainingTarget.count == 1 {
                map[remainingQuery[0]] = remainingTarget[0]
            }
            guard original.carriers.allSatisfy({ map[$0] != nil }) else { continue }
            var translated = original
            translated.focus = focus
            translated.carriers = original.carriers.map { map[$0]! }
            left.append(translated); right.append(other)
        }
        let ids = Set(left.flatMap { $0.carriers + [$0.focus] })
        let identity = Dictionary(uniqueKeysWithValues: ids.map { ($0, $0) })
        if try QueryStereoSupport.compare(query: left, target: right, mapping: identity, targetMolecule: t.molecule)
            == .no
        {
            return false
        }
        for expected in queryAlkenes {
            guard let a = mapping[expected.a], let b = mapping[expected.b],
                let doubleBond = q.edges[expected.a]?[expected.b], commonBonds[doubleBond.id] != nil,
                let actual = targetAlkenes.first(where: { Set([$0.a, $0.b]) == Set([a, b]) })
            else { continue }
            func mappedSide(_ focus: Int, excluding other: Int) -> Int? {
                q.edges[focus, default: [:]].keys.sorted().first {
                    $0 != other && mapping[$0] != nil && commonBonds[q.edges[focus]![$0]!.id] != nil
                }
            }
            guard let l = mappedSide(expected.a, excluding: expected.b),
                let r = mappedSide(expected.b, excluding: expected.a)
            else { continue }
            let targetLeft = actual.a == a ? actual.left : actual.right
            let targetRight = actual.a == a ? actual.right : actual.left
            let queryCis = expected.cis != ((l != expected.left) != (r != expected.right))
            let targetCis = actual.cis != ((mapping[l] != targetLeft) != (mapping[r] != targetRight))
            if queryCis != targetCis { return false }
        }
        return true
    }
}
