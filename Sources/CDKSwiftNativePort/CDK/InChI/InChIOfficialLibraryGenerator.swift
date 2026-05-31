import Foundation
import IUPACInChI

enum CDKInChIOfficialLibraryGenerator {
    private static let lock = NSLock()

    static func generate(for molecule: Molecule) throws -> CDKInChINativeGenerationResult {
        let normalized = normalizeInput(molecule)
        let molfile = try CDKMDLV2000Writer.write(
            normalized,
            options: .init(programName: "CDKSwiftNativePort")
        )

        lock.lock()
        defer { lock.unlock() }

        var result = CDKInChIResult()
        let bridgeStatus = molfile.withCString { cMolfile in
            cdk_inchi_generate_standard_molfile(cMolfile, &result)
        }
        defer {
            cdk_inchi_free_result(&result)
        }

        let message = string(from: result.message)
        let log = string(from: result.log)
        guard bridgeStatus == 0,
              let inchiPointer = result.inchi,
              let keyPointer = result.inchi_key else {
            let detail = [message, log]
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
