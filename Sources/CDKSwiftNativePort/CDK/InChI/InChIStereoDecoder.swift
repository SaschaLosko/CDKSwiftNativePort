import Foundation
import IUPACInChI

/// Decode ordered stereo with the bundled reference library, independently of
/// the source-identifier cache. Unknown parity never becomes a defined isomer.
enum CDKInChIStereoDecoder {
    static func apply(inchi: String, to molecule: inout Molecule) throws -> [String] {
        CDKInChIOfficialLibraryGenerator.lock.lock()
        defer { CDKInChIOfficialLibraryGenerator.lock.unlock() }
        var result = CDKInChIStereoResult()
        let status = inchi.withCString { cdk_inchi_read_stereo($0, &result) }
        defer { cdk_inchi_free_stereo(&result) }
        guard status == 0, Int(result.atom_count) >= molecule.atoms.count else {
            throw ChemError.parseFailed("Could not decode InChI stereochemistry.")
        }
        let atoms = molecule.atoms.sorted { $0.id < $1.id }
        for index in atoms.indices {
            guard let symbols = result.atom_symbols,
                String(cString: symbols.advanced(by: index * 8)) == atoms[index].element
            else {
                throw ChemError.unsupported("InChI stereo atom ordering does not match the decoded graph.")
            }
        }
        var hydrogenCounts = atoms.indices.map { Int(result.atom_hydrogens?[$0] ?? 0) }
        for index in atoms.count..<Int(result.atom_count) {
            let parent = Int(result.hydrogen_parents?[index] ?? -1)
            guard hydrogenCounts.indices.contains(parent) else {
                throw ChemError.unsupported("Additional isotope or disconnected atoms in InChI stereo decoding.")
            }
            hydrogenCounts[parent] += 1
        }
        for index in atoms.indices {
            let target = molecule.atoms.firstIndex { $0.id == atoms[index].id }!
            molecule.atoms[target].explicitHydrogenCount = hydrogenCounts[index]
            molecule.atoms[target].charge = Int(result.atom_charges?[index] ?? 0)
        }
        for index in 0..<Int(result.bond_count) {
            guard let bonds = result.bonds else { break }
            let a = Int(bonds[index * 3]), b = Int(bonds[index * 3 + 1])
            guard a >= 0, b >= 0 else { throw ChemError.parseFailed("Invalid InChI bond endpoints.") }
            if a >= atoms.count || b >= atoms.count { continue }
            guard
                let target = molecule.bonds.firstIndex(where: {
                    Set([$0.a1, $0.a2]) == Set([atoms[a].id, atoms[b].id])
                }), let order = BondOrder(rawValue: Int(bonds[index * 3 + 2]))
            else {
                throw ChemError.unsupported("InChI stereo graph connectivity differs.")
            }
            molecule.bonds[target].order = order
        }
        var ignored: [String] = []
        for index in 0..<Int(result.count) {
            guard let records = result.records else { break }
            let record = records.advanced(by: index * 7)
            var parity = record[6]
            guard parity == 1 || parity == 2 else { continue }
            let indices = (2..<6).map { Int(record[$0]) }
            guard indices.allSatisfy({ $0 >= 0 && $0 < Int(result.atom_count) }) else {
                throw ChemError.parseFailed("Invalid InChI stereo atom references.")
            }
            var ids = indices.map { $0 < atoms.count ? atoms[$0].id : -1 }
            switch record[0] {
            case 1:
                guard indices[1] < atoms.count, indices[2] < atoms.count else {
                    throw ChemError.unsupported("Extended InChI stereo endpoints.")
                }
                for (carrier, center, other) in [(0, 1, 2), (3, 2, 1)] where ids[carrier] == -1 {
                    guard result.hydrogen_parents?[indices[carrier]] == Int32(indices[center]),
                        molecule.implicitHydrogenCount(for: ids[center]) > 0,
                        let replacement = molecule.neighbors(of: ids[center]).sorted().first(where: { $0 != ids[other] }
                        )
                    else { throw ChemError.unsupported("Unmapped InChI stereo hydrogen.") }
                    ids[carrier] = replacement
                    parity = parity == 1 ? 2 : 1
                }
                guard
                    let bondIndex = molecule.bonds.firstIndex(where: {
                        Set([$0.a1, $0.a2]) == Set([ids[1], ids[2]]) && $0.order == .double
                    }), molecule.neighbors(of: ids[1]).contains(ids[0]),
                    molecule.neighbors(of: ids[2]).contains(ids[3])
                else {
                    ignored.append("extended double-bond stereo")
                    continue
                }
                molecule.bonds[bondIndex].doubleBondStereo = parity == 1 ? .cis : .trans
                molecule.bonds[bondIndex].stereoReferenceAtomIDs = ids
                molecule.bonds[bondIndex].stereo = .none
            case 2:
                let center = Int(record[1])
                guard center >= 0, center < atoms.count,
                    let atomIndex = molecule.atoms.firstIndex(where: { $0.id == atoms[center].id })
                else {
                    throw ChemError.parseFailed("Invalid InChI tetrahedral center.")
                }
                for index in ids.indices where ids[index] == -1 {
                    guard result.hydrogen_parents?[indices[index]] == Int32(center),
                        molecule.implicitHydrogenCount(for: atoms[center].id) > 0
                    else {
                        throw ChemError.unsupported("Unmapped InChI tetrahedral hydrogen.")
                    }
                    ids[index] = atoms[center].id
                }
                molecule.atoms[atomIndex].chirality = parity == 1 ? .clockwise : .anticlockwise
                molecule.atoms[atomIndex].ligandOrderingAtomIDs = ids
            default:
                ignored.append("unsupported InChI stereo geometry")
            }
        }
        return ignored
    }
}
