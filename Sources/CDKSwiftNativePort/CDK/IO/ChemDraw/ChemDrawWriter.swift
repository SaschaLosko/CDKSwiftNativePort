import Foundation

/// Native Swift ChemDraw interchange writers. These export editable chemical graphs,
/// not rendered images. Unsupported query/polymer semantics fail explicitly.
public enum CDKChemDrawWriter {
    public static func cdxml(molecules: [Molecule]) throws -> String {
        var builder = ChemDrawBuilder()
        return try xml(builder.document(molecules: molecules))
    }

    public static func cdx(molecules: [Molecule]) throws -> Data {
        var builder = ChemDrawBuilder()
        return try binary(builder.document(molecules: molecules))
    }

    public static func cdxml(reactions: [CDKReaction]) throws -> String {
        var builder = ChemDrawBuilder()
        return try xml(builder.document(reactions: reactions))
    }

    public static func cdx(reactions: [CDKReaction]) throws -> Data {
        var builder = ChemDrawBuilder()
        return try binary(builder.document(reactions: reactions))
    }

    private static func xml(_ document: ChemDrawObject) -> String {
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" + document.xml()
    }

    private static func binary(_ document: ChemDrawObject) throws -> Data {
        var data = Data("VjCD0100".utf8)
        data.append(contentsOf: [4, 3, 2, 1])
        data.append(Data(repeating: 0, count: 10))
        try document.appendBinary(to: &data)
        data.append(littleEndian(UInt16(0)))
        return data
    }
}

struct ChemDrawBuilder {
    private var nextID: UInt32 = 1

    mutating func object(_ tag: UInt16, _ name: String) throws -> ChemDrawObject {
        guard nextID < UInt32.max else { throw ChemError.unsupported("Too many ChemDraw objects.") }
        defer { nextID += 1 }
        return ChemDrawObject(tag: tag, name: name, id: nextID)
    }

    private mutating func root() throws -> ChemDrawObject {
        var root = try object(0x8000, "CDXML")
        // Platform 0, one UTF-8 font, ID 3, five-byte name Arial.
        var fonts = Data()
        for value: UInt16 in [0, 1, 3, 65001, 5] { fonts.append(littleEndian(value)) }
        fonts.append(contentsOf: "Arial".utf8)
        root.properties = [
            ChemDrawProperty(tag: 0x0100, name: "", text: "", bytes: fonts),
            .integer(0x0805, "BondLength", Int32(30 * 65536), text: "30"),
            .integer(0x081A, "LabelFont", UInt16(3)),
        ]
        return root
    }

    mutating func document(molecules: [Molecule]) throws -> ChemDrawObject {
        guard !molecules.isEmpty else { throw ChemError.emptyInput }
        var root = try root()
        for molecule in molecules {
            var page = try object(0x8001, "page")
            page.children.append(try fragment(molecule, x: 40, y: 60).object)
            page.children.append(try text(molecule.name, x: 40, y: 30))
            root.children.append(page)
        }
        return root
    }

    mutating func document(reactions: [CDKReaction]) throws -> ChemDrawObject {
        guard !reactions.isEmpty else { throw ChemError.emptyInput }
        var root = try root()
        for reaction in reactions { root.children.append(try reactionPage(reaction)) }
        return root
    }

    mutating func text(_ value: String, x: Double, y: Double) throws -> ChemDrawObject {
        try validateText(value)
        var text = try object(0x8006, "t")
        text.properties = [
            try .coordinates(0x0200, "p", [x, y], binaryOrder: [1, 0]),
            .string(0x0700, "", value),
        ]
        return text
    }

    mutating func fragment(_ source: Molecule, x: Double, y: Double)
        throws -> (object: ChemDrawObject, width: Double, height: Double, maps: [Int: UInt32])
    {
        try validate(source)
        var molecule = source
        let positions = Dictionary(uniqueKeysWithValues: source.atoms.map { ($0.id, $0.position) })
        let lengths = source.bonds.map { bond in
            let a = positions[bond.a1]!, b = positions[bond.a2]!
            return hypot(Double(a.x - b.x), Double(a.y - b.y))
        }.filter { $0 > 0.0001 }
        if lengths.count != source.bonds.count || (source.atoms.count > 1 && source.boundingBox()?.size == .zero) {
            molecule = Depiction2DGenerator.generate(for: source)
        }
        molecule = try ChemDrawStereo.prepare(molecule)
        let bounds = molecule.boundingBox()!
        let byID = Dictionary(uniqueKeysWithValues: molecule.atoms.map { ($0.id, $0.position) })
        let normalizedLengths = molecule.bonds.map { bond in
            let a = byID[bond.a1]!, b = byID[bond.a2]!
            return hypot(Double(a.x - b.x), Double(a.y - b.y))
        }.filter { $0 > 0.0001 }.sorted()
        let scale = normalizedLengths.isEmpty ? 1 : 30 / normalizedLengths[normalizedLengths.count / 2]
        var fragment = try object(0x8003, "fragment")
        fragment.properties = [.string(0x0008, "Name", molecule.name)]
        var ids: [Int: UInt32] = [:]
        var maps: [Int: UInt32] = [:]
        for atom in molecule.atoms {
            var node = try object(0x8004, "n")
            ids[atom.id] = node.id
            if let map = atom.atomMapNumber, map > 0 {
                guard maps.updateValue(node.id, forKey: map) == nil else {
                    throw ChemError.unsupported("ChemDraw export requires unique atom map numbers per participant.")
                }
            }
            let atomicNumber = Self.elements.firstIndex(of: atom.element)!
            node.properties = [
                try .coordinates(
                    0x0200, "p",
                    [
                        x + Double(atom.position.x - bounds.minX) * scale,
                        y + Double(bounds.maxY - atom.position.y) * scale,
                    ], binaryOrder: [1, 0]),
                .integer(0x0402, "Element", UInt16(atomicNumber)),
                .integer(0x0421, "Charge", Int8(atom.charge)),
            ]
            if let isotope = atom.isotopeMassNumber {
                node.properties.append(.integer(0x0420, "Isotope", UInt16(isotope)))
            }
            if let hydrogens = atom.explicitHydrogenCount {
                node.properties.append(.integer(0x042B, "NumHydrogens", UInt16(hydrogens)))
            }
            if let map = atom.atomMapNumber { node.properties.append(.string(0x0439, "AtomNumber", String(map))) }
            fragment.children.append(node)
        }
        for bond in molecule.bonds {
            var edge = try object(0x8005, "b")
            let order: UInt16 =
                switch bond.order {
                case .single: 1;
                case .double: 2;
                case .triple: 4;
                case .aromatic: 128
                }
            edge.properties = [
                .integer(0x0604, "B", ids[bond.a1]!), .integer(0x0605, "E", ids[bond.a2]!),
                .integer(0x0600, "Order", order, text: bond.order == .aromatic ? "1.5" : String(bond.order.rawValue)),
            ]
            let display: (UInt16, String) =
                switch bond.stereo {
                case .none: (0, "Solid")
                case .up: (6, "WedgeBegin")
                case .down: (3, "WedgedHashBegin")
                case .upReversed: (7, "WedgeEnd")
                case .downReversed: (4, "WedgedHashEnd")
                case .either: (8, "Wavy")
                }
            edge.properties.append(.integer(0x0601, "Display", display.0, text: display.1))
            fragment.children.append(edge)
        }
        return (fragment, max(30, Double(bounds.width) * scale), max(30, Double(bounds.height) * scale), maps)
    }

    private func validateText(_ value: String) throws {
        guard
            value.unicodeScalars.allSatisfy({
                (32...0xD7FF).contains($0.value) || (0xE000...0xFFFD).contains($0.value)
                    || (0x10000...0x10FFFF).contains($0.value) || [9, 10, 13].contains($0.value)
            })
        else {
            throw ChemError.unsupported("ChemDraw text contains unsupported control characters.")
        }
    }

    private func validate(_ molecule: Molecule) throws {
        guard !molecule.atoms.isEmpty else { throw ChemError.emptyInput }
        try validateText(molecule.name)
        let ids = Set(molecule.atoms.map(\.id))
        guard ids.count == molecule.atoms.count,
            Set(molecule.bonds.map(\.id)).count == molecule.bonds.count,
            molecule.bonds.allSatisfy({ ids.contains($0.a1) && ids.contains($0.a2) && $0.a1 != $0.a2 })
        else {
            throw ChemError.unsupported("ChemDraw export requires a valid graph with unique atom and bond IDs.")
        }
        guard molecule.sgroups.isEmpty, molecule.rGroupLogicDefinitions.isEmpty,
            molecule.cxState?.racemic != true,
            molecule.cxState?.racemicFragments.isEmpty != false,
            molecule.cxState?.stereoGroups.isEmpty != false
        else {
            throw ChemError.unsupported(
                "ChemDraw export does not yet support polymer groups or R-group definitions. Use MOL V3000 or RGfile.")
        }
        for atom in molecule.atoms {
            guard atom.zPosition == nil || atom.zPosition == 0 else {
                throw ChemError.unsupported(
                    "ChemDraw export currently supports 2D structures. Use MOL V3000 to preserve 3D coordinates.")
            }
            guard Self.elements.contains(atom.element), atom.queryType == nil, atom.atomList == nil,
                atom.rGroupLabel == nil, atom.rGroupMembership == nil, atom.attachmentPoint == nil,
                atom.substitutionCount == nil, atom.unsaturated == nil, atom.ringBondCount == nil,
                atom.cxStereoGroup == nil, atom.radical == nil, atom.radicalType == nil,
                atom.aliasLabel == nil, atom.valenceOverride == nil
            else {
                throw ChemError.unsupported(
                    "ChemDraw export does not yet support this atom's query, radical, alias, or enhanced stereo features. Use MOL V3000."
                )
            }
            guard atom.position.x.isFinite, atom.position.y.isFinite, (-128...127).contains(atom.charge),
                (atom.isotopeMassNumber.map { (1...32767).contains($0) } ?? true),
                (atom.explicitHydrogenCount.map { (0...65535).contains($0) } ?? true)
            else {
                throw ChemError.unsupported(
                    "Invalid atom coordinates, charge, isotope, or hydrogen count for ChemDraw export.")
            }
        }
        guard
            molecule.bonds.allSatisfy({
                $0.queryType == nil && $0.topology == nil && $0.coordinateBondReferenceAtomID == nil
                    && $0.reactingCenterStatus == nil
            })
        else {
            throw ChemError.unsupported(
                "ChemDraw export does not yet support query, coordinate, or reaction-center bond annotations. Use MOL V3000."
            )
        }
    }

    private static let elements =
        ("_ H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
        + "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg "
        + "Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og")
        .split(separator: " ").map(String.init)
}
