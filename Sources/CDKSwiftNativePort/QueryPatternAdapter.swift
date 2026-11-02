import Foundation

extension CDKQueryPattern {
    /// Adapts imported MDL/CX atom/bond predicates and positional-variation metadata.
    /// R-group definitions must be expanded first; this does not interpret text labels.
    public init(molecule: Molecule) throws {
        try QueryGraphFacts.validate(molecule)
        var atoms: [Int: CDKAtomExpression] = [:]
        var bonds: [Int: CDKBondExpression] = [:]
        let aromatic = CDKDescriptorRingSupport.aromaticBondIDs(in: molecule, perceiveAromaticity: true)
        let aromaticAtoms = Set(molecule.bonds.filter { aromatic.contains($0.id) }.flatMap { [$0.a1, $0.a2] })
        for atom in molecule.atoms {
            guard atom.rGroupLabel == nil, atom.rGroupMembership == nil, !atom.isUnresolvedPseudoAtom,
                atom.attachmentPoint == nil
            else { throw CDKQueryError.unsupported("Unexpanded variable or attachment atom.") }
            var expression = CDKAtomExpression.any
            if let list = atom.atomList {
                let numbers = try list.map { symbol -> Int in
                    guard let number = QueryElements.number(symbol) else {
                        throw CDKQueryError.invalid("Unknown atom-list element.")
                    }
                    return number
                }
                expression = numbers.reduce(.none) { .or($0, .element($1)) }
                if atom.atomListIsNegated { expression = .not(expression) }
            } else if let type = atom.queryType {
                switch type {
                case .anyAtom: expression = .any
                case .anyNonHydrogen: expression = .not(.element(1))
                case .anyHetero: expression = .and(.not(.element(1)), .not(.element(6)))
                }
            } else if atom.element != "*" {
                guard let number = QueryElements.number(atom.element) else {
                    throw CDKQueryError.unsupported("Undefined query atom label.")
                }
                expression = .and(.element(number), .aromatic(atom.aromatic || aromaticAtoms.contains(atom.id)))
            }
            expression = .and(expression, .charge(atom.charge))
            if let isotope = atom.isotopeMassNumber { expression = .and(expression, .isotope(isotope)) }
            if let hydrogen = atom.explicitHydrogenCount { expression = .and(expression, .hydrogens(hydrogen)) }
            if let count = atom.substitutionCount {
                guard count >= 0 else { throw CDKQueryError.unsupported("Special substitution-count query.") }
                expression = .and(expression, .heavyDegree(count))
            }
            if let count = atom.ringBondCount {
                guard count >= 0 else { throw CDKQueryError.unsupported("Special ring-bond-count query.") }
                expression = .and(expression, .ringBonds(count))
            }
            if atom.unsaturated != nil {
                throw CDKQueryError.unsupported("Unsaturation query needs an explicit expression.")
            }
            if atom.rGroupRestH == true {
                expression = .and(
                    expression,
                    .heavyDegree(
                        molecule.neighbors(of: atom.id).filter {
                            !["H", "D", "T"].contains(molecule.atom(id: $0)!.element)
                        }.count))
            }
            // A variable-attachment dummy is a mapping placeholder, not a charge constraint.
            let isPlaceholder = molecule.sgroups.contains { $0.kind == .extMulticenter && $0.atomIDs.first == atom.id }
            if isPlaceholder {
                guard atom.element == "*", atom.charge == 0, atom.queryType == nil,
                    atom.aliasLabel == nil, atom.isotopeMassNumber == nil, atom.explicitHydrogenCount == nil,
                    atom.atomList == nil, atom.substitutionCount == nil, atom.ringBondCount == nil,
                    atom.chirality == .none
                else {
                    throw CDKQueryError.unsupported("A positional placeholder cannot carry chemical predicates.")
                }
                atoms[atom.id] = .any
            } else {
                atoms[atom.id] = expression
            }
        }
        for bond in molecule.bonds {
            var expression: CDKBondExpression
            switch bond.queryType {
            case .any: expression = .any
            case .singleOrDouble: expression = .or(.order(.single), .order(.double))
            case .singleOrAromatic: expression = .or(.order(.single), .order(.aromatic))
            case .doubleOrAromatic: expression = .or(.order(.double), .order(.aromatic))
            case nil: expression = .order(aromatic.contains(bond.id) ? .aromatic : bond.order)
            }
            if let topology = bond.topology {
                expression = .and(expression, topology == .ring ? .inRing : .not(.inRing))
            }
            bonds[bond.id] = expression
        }
        self.init(molecule: molecule, atoms: atoms, bonds: bonds)
    }
}
