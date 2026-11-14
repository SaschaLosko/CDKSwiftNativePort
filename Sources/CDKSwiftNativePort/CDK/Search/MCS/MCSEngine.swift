import Foundation

struct MCSEngine {
    let query: Molecule
    let target: Molecule
    let pattern: CDKQueryPattern?
    let options: CDKMCSOptions
    let budget: MarkushMatchBudget
    var q: MCSGraph
    var t: MCSGraph
    var predicateMatcher: QuerySearch?
    var stereo: MCSStereo?
    var ids: [Int] = []
    var candidates: [Int: [Int]] = [:]
    var mapping: [Int: Int] = [:]
    var used: Set<Int> = []
    var best: [CDKSubstructureMatch] = []
    var bestAtoms = 0
    var bestBonds = 0
    var bestKeys: Set<[Int?]> = []
    var truncated = false

    init(
        query: Molecule, target: Molecule, pattern: CDKQueryPattern?, options: CDKMCSOptions, budget: MarkushMatchBudget
    ) throws {
        self.query = query; self.target = target; self.pattern = pattern
        self.options = options; self.budget = budget
        q = try MCSGraph(query); t = try MCSGraph(target)
    }

    mutating func prepare() throws {
        try budget.check()
        func requireGraph(_ molecule: Molecule, allowPredicates: Bool) throws {
            guard molecule.rGroupLogicDefinitions.isEmpty, molecule.cxState?.linkNodes.isEmpty != false,
                !molecule.sgroups.contains(where: {
                    [.structureRepeatUnit, .polymer, .extMulticenter].contains($0.kind)
                }),
                molecule.atoms.allSatisfy({ $0.rGroupMembership == nil && $0.attachmentPoint == nil }),
                allowPredicates
                    || (molecule.atoms.allSatisfy {
                        $0.queryType == nil && $0.atomList == nil && $0.substitutionCount == nil
                            && $0.unsaturated == nil && $0.ringBondCount == nil
                    }
                        && molecule.bonds.allSatisfy { $0.queryType == nil && $0.topology == nil })
            else {
                throw CDKMCSError.unsupported("Expand family/variable definitions and use explicit query predicates.")
            }
        }
        try requireGraph(query, allowPredicates: pattern != nil)
        try requireGraph(target, allowPredicates: false)
        if query.atoms.isEmpty || target.atoms.isEmpty { return }
        if options.completeRings || options.matchRingTopology {
            try q.prepareRings(complete: options.completeRings, budget: budget)
            try t.prepareRings(complete: options.completeRings, budget: budget)
        }
        if let pattern {
            let facts = try QueryGraphFacts(target, budget: budget)
            let matcher = try QuerySearch(
                query: pattern, target: facts, stereo: QueryStereoSupport.elements(target, explicit: []),
                mode: .substructure, budget: budget)
            guard matcher.positions.placeholders.isEmpty else {
                throw CDKMCSError.unsupported("MCS positional-variation predicates require expansion.")
            }
            predicateMatcher = matcher
        } else if options.matchStereochemistry {
            stereo = try MCSStereo(query: query, target: target)
        }
        for atom in query.atoms {
            for other in target.atoms {
                try budget.check()
                if try atomsMatch(atom, other) { candidates[atom.id, default: []].append(other.id) }
            }
            candidates[atom.id]?.sort()
        }
        ids = query.atoms.map(\.id).sorted { a, b in
            let ac = candidates[a, default: []].count, bc = candidates[b, default: []].count
            if ac != bc { return ac < bc }
            let ad = q.edges[a, default: [:]].count, bd = q.edges[b, default: [:]].count
            return ad == bd ? a < b : ad > bd
        }
    }

    func atomsMatch(_ a: Atom, _ b: Atom) throws -> Bool {
        if options.completeRings || options.matchRingTopology,
            q.ringAtoms.contains(a.id) != t.ringAtoms.contains(b.id)
        {
            return false
        }
        if let predicateMatcher, let pattern {
            let truth = try predicateMatcher.expression(pattern.atoms[a.id]!, at: b.id, level: 0)
            guard truth != .unknown else {
                throw CDKMCSError.unsupported("Unresolved recursive query stereochemistry.")
            }
            return truth == .yes
        }
        switch options.atomMatching {
        case .symbol: if a.element != b.element { return false }
        case .element:
            if let number = QueryElements.number(a.element), number > 0,
                number != QueryElements.number(b.element)
            {
                return false
            }
        case .any: break
        }
        if options.matchFormalCharge && a.charge != b.charge { return false }
        if options.matchIsotopes, let isotope = a.isotopeMassNumber, isotope != b.isotopeMassNumber { return false }
        return true
    }

    func bondsMatch(_ a: Bond, _ b: Bond) throws -> Bool {
        if options.completeRings || options.matchRingTopology,
            q.ringBonds.contains(a.id) != t.ringBonds.contains(b.id)
        {
            return false
        }
        if let predicateMatcher, let pattern {
            return try predicateMatcher.target.bond(pattern.bonds[a.id]!, b, budget: budget)
        }
        return options.ignoreBondOrder || a.order == b.order
    }

    mutating func visit(_ index: Int) throws {
        try budget.check()
        try consider()
        guard index < ids.count else { return }
        let possibleAtoms = mapping.count + min(ids.count - index, target.atoms.count - used.count)
        if options.minimumBondsPerFragment == nil && possibleAtoms < bestAtoms { return }
        let id = ids[index]
        // Try extensions sharing the most existing common edges first, then stable IDs.
        let ordered = candidates[id, default: []].filter { !used.contains($0) }.sorted { a, b in
            func links(_ candidate: Int) -> Int {
                q.edges[id, default: [:]].keys.filter { neighbor in
                    mapping[neighbor].map { t.edges[candidate]?[$0] != nil } ?? false
                }.count
            }
            let ac = links(a), bc = links(b)
            return ac == bc ? a < b : ac > bc
        }
        for candidate in ordered {
            mapping[id] = candidate; used.insert(candidate)
            try visit(index + 1)
            mapping[id] = nil; used.remove(candidate)
        }
        try visit(index + 1)  // Exclude this query atom; no connectivity assumption prunes valid solutions.
    }

    mutating func consider() throws {
        guard !mapping.isEmpty else { return }
        var bonds: [Int: Int] = [:]
        var adjacency: [Int: Set<Int>] = [:]
        for bond in query.bonds {
            guard let a = mapping[bond.a1], let b = mapping[bond.a2], let other = t.edges[a]?[b],
                try bondsMatch(bond, other)
            else { continue }
            bonds[bond.id] = other.id
            adjacency[bond.a1, default: []].insert(bond.a2)
            adjacency[bond.a2, default: []].insert(bond.a1)
        }
        let atoms = mapping.count, count = bonds.count
        let improves: Bool
        let ties = atoms == bestAtoms && count == bestBonds
        if options.minimumBondsPerFragment != nil {
            improves = count > bestBonds || (count == bestBonds && atoms > bestAtoms)
        } else {
            improves = atoms > bestAtoms || (atoms == bestAtoms && count > bestBonds)
        }
        guard improves || ties else { return }
        var remaining = Set(mapping.keys)
        var components = 0
        while let seed = remaining.min() {
            components += 1
            var queue = [seed]; remaining.remove(seed); var cursor = 0; var degreeSum = 0
            while cursor < queue.count {
                let neighbors = adjacency[queue[cursor], default: []]
                degreeSum += neighbors.count
                for neighbor in neighbors where remaining.remove(neighbor) != nil { queue.append(neighbor) }
                cursor += 1
            }
            if let minimum = options.minimumBondsPerFragment {
                if degreeSum / 2 < minimum { return }
            } else if components > 1 {
                return
            }
        }
        if options.completeRings {
            guard q.ringsComplete(atomIDs: Set(mapping.keys), bondIDs: Set(bonds.keys)),
                t.ringsComplete(atomIDs: Set(mapping.values), bondIDs: Set(bonds.values))
            else { return }
        }
        if let stereo, try !stereo.accepts(mapping: mapping, commonBonds: bonds, query: q, target: t) { return }
        let key = query.atoms.map { mapping[$0.id] }
        if improves {
            bestAtoms = atoms; bestBonds = count; best = []; bestKeys = []; truncated = false
        }
        guard !bestKeys.contains(key) else { return }
        if best.count < options.maximumMappings {
            best.append(.init(atomIDByQueryAtomID: mapping, bondIDByQueryBondID: bonds)); bestKeys.insert(key)
        } else {
            truncated = true
        }
    }
}
