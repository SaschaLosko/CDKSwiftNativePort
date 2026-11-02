import Foundation

/// Finite single-bond position variation with an explicit unlabelled dummy endpoint.
/// Carbon placeholders and coordination/multicentre bonds are ambiguous and remain unsupported.
struct MarkushPositionPlan {
    struct Site { let dummy: Int; let bond: Bond; let targets: [Int] }
    let source: Molecule
    let sites: [Site]

    init(_ molecule: Molecule) throws {
        var sites: [Site] = []
        let groups = molecule.sgroups.filter { $0.kind == .extMulticenter }
        for (index, endpoints) in molecule.cxState?.positionalVariations ?? [:] {
            guard molecule.atoms.indices.contains(index), endpoints.allSatisfy(molecule.atoms.indices.contains),
                groups.contains(where: {
                    $0.atomIDs == [molecule.atoms[index].id] + endpoints.map { molecule.atoms[$0].id }
                })
            else {
                throw CDKRGroupEnumerationError.invalidQuery("Inconsistent positional attachment metadata.")
            }
        }
        for group in groups {
            guard group.keyword == "m", let dummy = group.atomIDs.first,
                let atom = molecule.atom(id: dummy), atom.element == "*",
                atom.aliasLabel == nil, atom.rGroupLabel == nil, atom.rGroupMembership == nil,
                atom.charge == 0, atom.chirality == .none,
                molecule.bonds(forAtom: dummy).count == 1,
                let bond = molecule.bonds(forAtom: dummy).first, bond.order == .single,
                bond.coordinateBondReferenceAtomID == nil, group.atomIDs.count > 1
            else {
                throw CDKRGroupEnumerationError.unsupportedFeature(
                    "Position variation requires an explicit single-bond dummy endpoint.")
            }
            let neighbor = bond.a1 == dummy ? bond.a2 : bond.a1
            let targets = Array(group.atomIDs.dropFirst())
            guard Set(targets).count == targets.count,
                targets.allSatisfy({ id in
                    guard id != dummy, id != neighbor, let target = molecule.atom(id: id) else { return false }
                    return target.rGroupMembership == nil && target.rGroupLabel == nil
                        && target.chirality == .none && !target.isUnresolvedPseudoAtom
                        && CDKDescriptorSupport.implicitHydrogenCount(on: id, in: molecule) > 0
                        && !molecule.neighbors(of: neighbor).contains(id)
                })
            else {
                throw CDKRGroupEnumerationError.unsupportedFeature(
                    "A variable attachment has an invalid, stereogenic or fully substituted target.")
            }
            sites.append(.init(dummy: dummy, bond: bond, targets: targets))
        }
        guard Set(sites.map(\.dummy)).count == sites.count, sites.count <= 16,
            sites.allSatisfy({ site in sites.allSatisfy { other in !site.targets.contains(other.dummy) } })
        else {
            throw CDKRGroupEnumerationError.unsupportedFeature("Variable attachment sites overlap or exceed the limit.")
        }
        if molecule.cxState?.positionalVariations.isEmpty == false && sites.isEmpty {
            throw CDKRGroupEnumerationError.invalidQuery("Position metadata has no corresponding attachment groups.")
        }
        var source = molecule
        source.sgroups.removeAll { $0.kind == .extMulticenter }
        source.cxState?.positionalVariations = [:]
        self.source = source
        self.sites = sites
    }

    func enumerate(checkpoint: () throws -> Void, _ visit: (Molecule, [Int: Int]) throws -> Bool) throws {
        func expand(_ index: Int, _ graph: Molecule, _ targets: [Int: Int]) throws -> Bool {
            try checkpoint()
            guard index < sites.count else { return try visit(graph, targets) }
            let site = sites[index]
            for target in site.targets {
                // Multiple independently varying substituents cannot occupy an exhausted site.
                guard CDKDescriptorSupport.implicitHydrogenCount(on: target, in: graph) > 0 else { continue }
                var member = graph
                member.atoms.removeAll { $0.id == site.dummy }
                for i in member.atoms.indices {
                    if member.atoms[i].id == target, let h = member.atoms[i].explicitHydrogenCount {
                        member.atoms[i].explicitHydrogenCount = max(0, h - 1)
                    }
                    member.atoms[i].ligandOrderingAtomIDs = member.atoms[i].ligandOrderingAtomIDs?.map {
                        $0 == site.dummy ? target : $0
                    }
                }
                guard let b = member.bonds.firstIndex(where: { $0.id == site.bond.id }) else {
                    throw CDKRGroupEnumerationError.invalidQuery("A variable attachment bond is missing.")
                }
                member.bonds[b] = RGroupEnumerationPlan.copyBond(
                    site.bond, id: site.bond.id,
                    a1: site.bond.a1 == site.dummy ? target : site.bond.a1,
                    a2: site.bond.a2 == site.dummy ? target : site.bond.a2,
                    mapReference: { $0 == site.dummy ? target : $0 })
                var selection = targets; selection[site.dummy] = target
                if try !expand(index + 1, member, selection) { return false }
            }
            return true
        }
        _ = try expand(0, source, [:])
    }
}
