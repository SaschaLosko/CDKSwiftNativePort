import Foundation
import IUPACInChI

enum CDKInChIOfficialLibraryGenerator {
    enum Mode {
        case standard
        case fixedH
    }

    private static let lock = NSLock()

    static func generate(
        for molecule: Molecule,
        mode: Mode = .standard
    ) throws -> CDKInChINativeGenerationResult {
        let normalized = normalizeInput(molecule)
        if normalized.atoms.contains(where: { $0.ligandOrderingAtomIDs?.count == 4 })
            || normalized.bonds.contains(where: { $0.doubleBondStereo != nil })
        {
            return try generateFromAtoms(normalized, mode: mode)
        }

        let molfile = try CDKMDLV2000Writer.write(
            normalized,
            options: .init(programName: "CDKSwiftNativePort")
        )

        lock.lock()
        defer { lock.unlock() }

        var result = CDKInChIResult()
        let bridgeStatus = molfile.withCString { cMolfile in
            switch mode {
            case .standard:
                cdk_inchi_generate_standard_molfile(cMolfile, &result)
            case .fixedH:
                cdk_inchi_generate_fixed_h_molfile(cMolfile, &result)
            }
        }
        defer {
            cdk_inchi_free_result(&result)
        }

        let message = string(from: result.message)
        let log = string(from: result.log)
        guard bridgeStatus == 0,
            let inchiPointer = result.inchi,
            let keyPointer = result.inchi_key
        else {
            let detail =
                [message, log]
                .map { $0.trimmingCharacters(in: .whitespacesAndNewlines) }
                .first { !$0.isEmpty }
                ?? "Official InChI generation failed."
            throw ChemError.parseFailed(detail)
        }

        return CDKInChINativeGenerationResult(
            inchi: String(cString: inchiPointer),
            inchiKey: String(cString: keyPointer),
            auxInfo: string(from: result.aux_info),
            status: .success,
            message: message
        )
    }

    private static func generateFromAtoms(
        _ molecule: Molecule,
        mode: Mode
    ) throws -> CDKInChINativeGenerationResult {
        let atoms = molecule.atoms.sorted { $0.id < $1.id }
        let atomIndexByID = Dictionary(uniqueKeysWithValues: atoms.enumerated().map { ($1.id, Int32($0)) })
        let bonds = molecule.bonds.sorted { $0.id < $1.id }
        let symbolStride = 8

        var symbols = Array(repeating: CChar(0), count: atoms.count * symbolStride)
        var charges: [Int32] = []
        var isotopes: [Int32] = []
        var implicitHydrogens: [Int32] = []

        let hasExplicitHydrogenAtoms = atoms.contains {
            CDKDescriptorSupport.canonicalElementSymbol($0.element).uppercased() == "H"
        }
        for (atomIndex, atom) in atoms.enumerated() {
            let symbol = CDKDescriptorSupport.canonicalElementSymbol(atom.element)
            let encoded = Array(symbol.utf8.prefix(symbolStride - 1))
            for (byteIndex, byte) in encoded.enumerated() {
                symbols[atomIndex * symbolStride + byteIndex] = CChar(bitPattern: byte)
            }
            charges.append(Int32(atom.charge))
            isotopes.append(Int32(atom.isotopeMassNumber ?? 0))

            if CDKDescriptorSupport.canonicalElementSymbol(atom.element).uppercased() == "H" || hasExplicitHydrogenAtoms
            {
                implicitHydrogens.append(0)
            } else if let explicitHydrogenCount = atom.explicitHydrogenCount {
                implicitHydrogens.append(Int32(max(0, explicitHydrogenCount)))
            } else {
                implicitHydrogens.append(-1)
            }
        }

        var bondFrom: [Int32] = []
        var bondTo: [Int32] = []
        var bondOrder: [Int32] = []
        bondFrom.reserveCapacity(bonds.count)
        bondTo.reserveCapacity(bonds.count)
        bondOrder.reserveCapacity(bonds.count)

        for bond in bonds {
            guard let from = atomIndexByID[bond.a1],
                let to = atomIndexByID[bond.a2]
            else {
                throw ChemError.parseFailed("Bond references unknown atom while generating InChI.")
            }
            bondFrom.append(from)
            bondTo.append(to)
            switch bond.order {
            case .single:
                bondOrder.append(1)
            case .double:
                bondOrder.append(2)
            case .triple:
                bondOrder.append(3)
            case .aromatic:
                bondOrder.append(4)
            }
        }

        var stereoCenters: [Int32] = []
        var stereoTypes: [Int32] = []
        var stereoNeighbors: [Int32] = []
        var stereoParities: [Int32] = []

        for atom in atoms where atom.chirality != .none {
            guard let ligandOrdering = atom.ligandOrderingAtomIDs,
                ligandOrdering.count == 4,
                let center = atomIndexByID[atom.id]
            else {
                continue
            }

            let mappedNeighbors = ligandOrdering.compactMap { atomIndexByID[$0] }
            guard mappedNeighbors.count == 4 else { continue }

            stereoCenters.append(center)
            stereoTypes.append(2)
            stereoNeighbors.append(contentsOf: mappedNeighbors)
            stereoParities.append(atom.chirality == .clockwise ? 1 : 2)
        }

        for bond in bonds where bond.doubleBondStereo != nil {
            guard bond.order == .double,
                let doubleBondStereo = bond.doubleBondStereo,
                let referenceAtomIDs = bond.stereoReferenceAtomIDs,
                referenceAtomIDs.count == 4,
                Set(referenceAtomIDs[1...2]) == Set([bond.a1, bond.a2])
            else {
                throw ChemError.parseFailed("Double-bond stereochemistry did not define four valid atom references.")
            }
            let mappedNeighbors = referenceAtomIDs.compactMap { atomIndexByID[$0] }
            guard mappedNeighbors.count == 4 else {
                throw ChemError.parseFailed("Double-bond stereochemistry references an unknown atom.")
            }

            stereoCenters.append(-1)
            stereoTypes.append(1)
            stereoNeighbors.append(contentsOf: mappedNeighbors)
            stereoParities.append(doubleBondStereo == .cis ? 1 : 2)
        }

        lock.lock()
        defer { lock.unlock() }

        var result = CDKInChIResult()
        let bridgeStatus = symbols.withUnsafeBufferPointer { symbolPointer in
            charges.withUnsafeBufferPointer { chargePointer in
                isotopes.withUnsafeBufferPointer { isotopePointer in
                    implicitHydrogens.withUnsafeBufferPointer { hydrogenPointer in
                        bondFrom.withUnsafeBufferPointer { bondFromPointer in
                            bondTo.withUnsafeBufferPointer { bondToPointer in
                                bondOrder.withUnsafeBufferPointer { bondOrderPointer in
                                    stereoTypes.withUnsafeBufferPointer { stereoTypePointer in
                                        stereoCenters.withUnsafeBufferPointer { stereoCenterPointer in
                                            stereoNeighbors.withUnsafeBufferPointer { stereoNeighborPointer in
                                                stereoParities.withUnsafeBufferPointer { stereoParityPointer in
                                                    switch mode {
                                                    case .standard:
                                                        cdk_inchi_generate_standard_atoms_with_stereo_types(
                                                            Int32(atoms.count),
                                                            symbolPointer.baseAddress,
                                                            Int32(symbolStride),
                                                            chargePointer.baseAddress,
                                                            isotopePointer.baseAddress,
                                                            hydrogenPointer.baseAddress,
                                                            Int32(bonds.count),
                                                            bondFromPointer.baseAddress,
                                                            bondToPointer.baseAddress,
                                                            bondOrderPointer.baseAddress,
                                                            Int32(stereoCenters.count),
                                                            stereoTypePointer.baseAddress,
                                                            stereoCenterPointer.baseAddress,
                                                            stereoNeighborPointer.baseAddress,
                                                            stereoParityPointer.baseAddress,
                                                            &result
                                                        )
                                                    case .fixedH:
                                                        cdk_inchi_generate_fixed_h_atoms_with_stereo_types(
                                                            Int32(atoms.count),
                                                            symbolPointer.baseAddress,
                                                            Int32(symbolStride),
                                                            chargePointer.baseAddress,
                                                            isotopePointer.baseAddress,
                                                            hydrogenPointer.baseAddress,
                                                            Int32(bonds.count),
                                                            bondFromPointer.baseAddress,
                                                            bondToPointer.baseAddress,
                                                            bondOrderPointer.baseAddress,
                                                            Int32(stereoCenters.count),
                                                            stereoTypePointer.baseAddress,
                                                            stereoCenterPointer.baseAddress,
                                                            stereoNeighborPointer.baseAddress,
                                                            stereoParityPointer.baseAddress,
                                                            &result
                                                        )
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        defer {
            cdk_inchi_free_result(&result)
        }

        let message = string(from: result.message)
        let log = string(from: result.log)
        guard bridgeStatus == 0,
            let inchiPointer = result.inchi,
            let keyPointer = result.inchi_key
        else {
            let detail =
                [message, log]
                .map { $0.trimmingCharacters(in: .whitespacesAndNewlines) }
                .first { !$0.isEmpty }
                ?? "Official InChI generation failed."
            throw ChemError.parseFailed(detail)
        }

        return CDKInChINativeGenerationResult(
            inchi: String(cString: inchiPointer),
            inchiKey: String(cString: keyPointer),
            auxInfo: string(from: result.aux_info),
            status: .success,
            message: message
        )
    }

    static func inchiKey(for inchi: String) throws -> String {
        lock.lock()
        defer { lock.unlock() }

        var keyBuffer = Array(repeating: CChar(0), count: 28)
        let status = inchi.withCString { cInChI in
            keyBuffer.withUnsafeMutableBufferPointer { keyPointer in
                cdk_inchi_key_from_standard_inchi(cInChI, keyPointer.baseAddress, Int32(keyPointer.count))
            }
        }
        guard status == 0 else {
            throw ChemError.parseFailed("Official InChIKey generation failed with status \(status).")
        }
        let endIndex = keyBuffer.firstIndex(of: 0) ?? keyBuffer.endIndex
        let bytes = keyBuffer[..<endIndex].map { UInt8(bitPattern: $0) }
        return String(decoding: bytes, as: UTF8.self)
    }

    static var version: String {
        guard let pointer = cdk_inchi_version() else {
            return "unknown"
        }
        return String(cString: pointer)
    }

    private static func string(from pointer: UnsafeMutablePointer<CChar>?) -> String {
        guard let pointer else { return "" }
        return String(cString: pointer)
    }

    private static func normalizeInput(_ molecule: Molecule) -> Molecule {
        var copy = molecule
        for index in copy.atoms.indices {
            var atom = copy.atoms[index]
            let rawElement = atom.element
            var normalized = normalizedElementSymbol(rawElement)

            if atom.charge == 0 {
                atom.charge += inferredCharge(fromRawElement: rawElement)
            }

            switch normalized.uppercased() {
            case "D":
                normalized = "H"
                if atom.isotopeMassNumber == nil {
                    atom.isotopeMassNumber = 2
                }
            case "T":
                normalized = "H"
                if atom.isotopeMassNumber == nil {
                    atom.isotopeMassNumber = 3
                }
            default:
                break
            }

            if isSupportedElementSymbol(normalized) {
                atom.element = normalized
            }
            copy.atoms[index] = atom
        }
        return copy
    }

    private static func normalizedElementSymbol(_ raw: String) -> String {
        let trimmed = raw.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { return "" }

        let letters = trimmed.prefix { $0.isLetter }
        guard !letters.isEmpty else { return trimmed }

        let source = String(letters)
        let first = String(source.prefix(1)).uppercased()
        if source.count == 1 {
            return first
        }

        let second = String(source.dropFirst().prefix(1)).lowercased()
        let candidate = first + second
        return candidate
    }

    private static func inferredCharge(fromRawElement raw: String) -> Int {
        let trimmed = raw.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { return 0 }

        let letterPrefixLength = trimmed.prefix { $0.isLetter }.count
        guard letterPrefixLength < trimmed.count else { return 0 }
        let suffix = String(trimmed.dropFirst(letterPrefixLength))
        guard !suffix.isEmpty else { return 0 }

        var index = suffix.startIndex
        var total = 0
        while index < suffix.endIndex {
            let signCharacter = suffix[index]
            guard signCharacter == "+" || signCharacter == "-" else {
                index = suffix.index(after: index)
                continue
            }

            let sign = signCharacter == "+" ? 1 : -1
            index = suffix.index(after: index)

            let digitStart = index
            while index < suffix.endIndex, suffix[index].isNumber {
                index = suffix.index(after: index)
            }
            if digitStart != index {
                total += sign * (Int(suffix[digitStart..<index]) ?? 1)
                continue
            }

            var repeatedSigns = 1
            while index < suffix.endIndex, suffix[index] == signCharacter {
                repeatedSigns += 1
                index = suffix.index(after: index)
            }
            total += sign * repeatedSigns
        }
        return total
    }

    private static func isSupportedElementSymbol(_ symbol: String) -> Bool {
        let canonical = CDKDescriptorSupport.canonicalElementSymbol(symbol)
        let uppercased = canonical.uppercased()
        if uppercased == "D" || uppercased == "T" {
            return true
        }
        return CDKDescriptorSupport.averageAtomicMass(forElementSymbol: canonical) > 0
    }
}
