import Foundation

/// CML's ordered stereo references carry chemistry independently of drawing coordinates.
/// See https://www.xml-cml.org/convention/molecular (atomParity and bondStereo).
enum CMLStereo {
    static func validate(_ molecule: Molecule) throws {
        let ids = Set(molecule.atoms.map(\.id))
        for atom in molecule.atoms where atom.chirality != .none {
            // Older coordinate-based models may have only wedge/hatch information.
            // Ordered metadata, when supplied, must never be partially serialized.
            if let refs = atom.ligandOrderingAtomIDs {
                guard refs.count == 4, Set(refs).count == 4, Set(refs).isSubset(of: ids) else {
                    throw ChemError.parseFailed("CML atom parity requires four distinct, existing atom references.")
                }
            }
        }
        for bond in molecule.bonds where bond.doubleBondStereo != nil {
            guard bond.order == .double, let refs = bond.stereoReferenceAtomIDs,
                refs.count == 4, Set(refs).count == 4, Set(refs).isSubset(of: ids),
                Set(refs[1...2]) == Set([bond.a1, bond.a2])
            else {
                throw ChemError.parseFailed("CML double-bond stereo requires four ordered atom references.")
            }
        }
    }

    static func atomParity(_ atom: Atom, reference: (Int) -> String) -> String? {
        guard atom.chirality != .none, let refs = atom.ligandOrderingAtomIDs, refs.count == 4 else { return nil }
        let parity = atom.chirality == .anticlockwise ? 1 : -1
        return "<atomParity atomRefs4=\"\(refs.map(reference).joined(separator: " "))\">\(parity)</atomParity>"
    }

    static func bondStereo(_ bond: Bond, reference: (Int) -> String) -> String? {
        if let stereo = bond.doubleBondStereo, let refs = bond.stereoReferenceAtomIDs {
            return
                "<bondStereo atomRefs4=\"\(refs.map(reference).joined(separator: " "))\">\(stereo == .cis ? "C" : "T")</bondStereo>"
        }
        let symbol: String
        let refs: [Int]
        switch bond.stereo {
        case .up: symbol = "W"; refs = [bond.a1, bond.a2]
        case .down: symbol = "H"; refs = [bond.a1, bond.a2]
        case .upReversed: symbol = "W"; refs = [bond.a2, bond.a1]
        case .downReversed: symbol = "H"; refs = [bond.a2, bond.a1]
        default: return nil
        }
        return
            "<bondStereo atomRefs2=\"\(refs.map(reference).joined(separator: " "))\" dictRef=\"cml:\(symbol)\">\(symbol)</bondStereo>"
    }

    static func parity(_ text: String, attributes: [String: String]) throws -> (AtomChirality, [String]?) {
        guard let value = Double(text.trimmingCharacters(in: .whitespacesAndNewlines)), value.isFinite else {
            throw ChemError.parseFailed("Invalid CML atom parity value.")
        }
        if value == 0 { return (.none, nil) }
        let refs = attributes["atomRefs4"]?.split(whereSeparator: \.isWhitespace).map(String.init) ?? []
        guard refs.count == 4, Set(refs).count == 4 else {
            throw ChemError.parseFailed("CML atom parity requires four distinct atom references.")
        }
        return (value > 0 ? .anticlockwise : .clockwise, refs)
    }

    static func wedge(_ value: BondStereo, attributes: [String: String], first: String, second: String) throws
        -> BondStereo
    {
        guard let raw = attributes["atomRefs2"] else { return value }  // Older CML documents omit this attribute.
        let refs = raw.split(whereSeparator: \.isWhitespace).map(String.init)
        if refs == [first, second] { return value }
        if refs == [second, first] { return value == .up ? .upReversed : .downReversed }
        throw ChemError.parseFailed("CML wedge references do not identify the parent bond.")
    }
}
