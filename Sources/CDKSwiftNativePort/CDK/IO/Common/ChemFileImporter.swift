import Foundation

public struct CDKFileImporterFormat: Hashable {
    public let displayName: String
    public let fileExtensions: [String]
    public let utiIdentifiers: [String]

    public init(displayName: String,
                fileExtensions: [String],
                utiIdentifiers: [String]) {
        self.displayName = displayName
        self.fileExtensions = fileExtensions
        self.utiIdentifiers = utiIdentifiers
    }
}

/// Unified file importer dispatch for CDKSwiftNativePort-backed formats.
public enum CDKFileImporter {
    public static let formats: [CDKFileImporterFormat] = [
        CDKFileImporterFormat(displayName: "MDL Molfile",
                              fileExtensions: ["mol"],
                              utiIdentifiers: ["chemical/x-mdl-molfile", "net.sourceforge.openbabel.mdl"]),
        CDKFileImporterFormat(displayName: "MDL SDFile",
                              fileExtensions: ["sdf", "sd"],
                              utiIdentifiers: ["chemical/x-mdl-sdfile", "net.sourceforge.openbabel.mdl"]),
        CDKFileImporterFormat(displayName: "SMILES",
                              fileExtensions: ["smi", "smiles", "ism", "can"],
                              utiIdentifiers: ["chemical/x-daylight-smiles", "chemical/x-smiles"]),
        CDKFileImporterFormat(displayName: "InChI",
                              fileExtensions: ["inchi", "ich"],
                              utiIdentifiers: ["chemical/x-inchi"]),
        CDKFileImporterFormat(displayName: "Tripos MOL2",
                              fileExtensions: ["mol2"],
                              utiIdentifiers: ["chemical/x-mol2"]),
        CDKFileImporterFormat(displayName: "Protein Data Bank",
                              fileExtensions: ["pdb", "ent"],
                              utiIdentifiers: ["chemical/x-pdb"]),
        CDKFileImporterFormat(displayName: "XYZ Coordinates",
                              fileExtensions: ["xyz"],
                              utiIdentifiers: ["chemical/x-xyz"]),
        CDKFileImporterFormat(displayName: "Chemical Markup Language",
                              fileExtensions: ["cml"],
                              utiIdentifiers: ["chemical/x-cml"]),
        CDKFileImporterFormat(displayName: "MDL RXN",
                              fileExtensions: ["rxn"],
                              utiIdentifiers: ["chemical/x-mdl-rxnfile"]),
        CDKFileImporterFormat(displayName: "MDL RDF",
                              fileExtensions: ["rdf"],
                              utiIdentifiers: ["chemical/x-mdl-rdfile"])
    ]

    public static var supportedFileExtensions: [String] {
        Array(Set(formats.flatMap(\.fileExtensions))).sorted()
    }

    public static var supportedUTIIdentifiers: [String] {
        Array(Set(formats.flatMap(\.utiIdentifiers))).sorted()
    }

    public static func supports(fileExtension: String) -> Bool {
        supportedFileExtensions.contains(fileExtension.lowercased()) || fileExtension.lowercased() == "txt"
    }

    public static func preferredInputFormat(forFileExtension ext: String,
                                            text: String? = nil) -> ChemFormat {
        let lower = ext.lowercased()
        if ["smi", "smiles", "ism", "can"].contains(lower) {
            return .smiles
        }
        if ["inchi", "ich"].contains(lower) {
            return .inchi
        }
        if lower == "txt", let text, looksLikeInChI(text) {
            return .inchi
        }
        return .sdf
    }

    public static func readMolecules(from url: URL) throws -> [Molecule] {
        let ext = url.pathExtension.lowercased()
        let shouldStopAccess: Bool
        if url.isFileURL {
            shouldStopAccess = url.startAccessingSecurityScopedResource()
        } else {
            shouldStopAccess = false
        }

        defer {
            if shouldStopAccess {
                url.stopAccessingSecurityScopedResource()
            }
        }

        let text = try decodeText(from: url)
        return try readMolecules(text: text, fileExtension: ext)
    }

    public static func readMolecules(text: String,
                                     fileExtension: String?) throws -> [Molecule] {
        let ext = (fileExtension ?? "").lowercased()

        switch ext {
        case "sdf", "sd":
            return try CDKIteratingSDFReader.read(text: text)
        case "mol":
            return [try CDKMDLReader.read(text: text)]
        case "smi", "smiles", "ism", "can":
            return try CDKSMILESReader.read(text: text)
        case "inchi", "ich":
            return try CDKInChIReader.read(text: text)
        case "xyz":
            return try CDKXYZReader.read(text: text)
        case "pdb", "ent":
            return try CDKPDBReader.read(text: text)
        case "mol2":
            return try CDKMol2Reader.read(text: text)
        case "cml":
            return try CDKCMLReader.read(text: text)
        case "rxn":
            return try CDKRXNReader.read(text: text)
        case "rdf":
            return try CDKRDFReader.read(text: text)
        case "txt":
            return try readTextWithAutoDetection(text)
        default:
            return try readTextWithAutoDetection(text)
        }
    }

    private static func readTextWithAutoDetection(_ text: String) throws -> [Molecule] {
        if looksLikeInChI(text) {
            return try CDKInChIReader.read(text: text)
        }
        if looksLikeRXN(text) {
            return try CDKRXNReader.read(text: text)
        }
        if looksLikeRDF(text) {
            return try CDKRDFReader.read(text: text)
        }
        if looksLikeMol2(text) {
            return try CDKMol2Reader.read(text: text)
        }
        if looksLikeCML(text) {
            return try CDKCMLReader.read(text: text)
        }
        if looksLikePDB(text) {
            return try CDKPDBReader.read(text: text)
        }

        let parsers: [(String, () throws -> [Molecule])] = [
            ("SDF", { try CDKIteratingSDFReader.read(text: text) }),
            ("MDL", { [try CDKMDLReader.read(text: text)] }),
            ("SMILES", { try CDKSMILESReader.read(text: text) }),
            ("InChI", { try CDKInChIReader.read(text: text) }),
            ("XYZ", { try CDKXYZReader.read(text: text) })
        ]

        for (_, parser) in parsers {
            if let molecules = try? parser(), !molecules.isEmpty {
                return molecules
            }
        }

        throw ChemError.unsupported("Unable to detect a supported molecule/reaction format.")
    }

    private static func decodeText(from url: URL) throws -> String {
        let data = try Data(contentsOf: url)
        if let utf8 = String(data: data, encoding: .utf8) {
            return utf8
        }
        return String(decoding: data, as: UTF8.self)
    }

    private static func firstMeaningfulLine(in text: String) -> String? {
        for raw in text.components(separatedBy: .newlines) {
            let line = raw.trimmingCharacters(in: .whitespacesAndNewlines)
            if line.isEmpty || line.hasPrefix("#") || line.hasPrefix("//") {
                continue
            }
            return line
        }
        return nil
    }

    private static func looksLikeInChI(_ text: String) -> Bool {
        firstMeaningfulLine(in: text)?.hasPrefix("InChI=") == true
    }

    private static func looksLikeMol2(_ text: String) -> Bool {
        text.localizedCaseInsensitiveContains("@<TRIPOS>MOLECULE")
    }

    private static func looksLikeCML(_ text: String) -> Bool {
        text.localizedCaseInsensitiveContains("<cml") || text.localizedCaseInsensitiveContains("<molecule")
    }

    private static func looksLikeRXN(_ text: String) -> Bool {
        text.components(separatedBy: .newlines).contains {
            $0.trimmingCharacters(in: .whitespacesAndNewlines) == "$RXN"
        }
    }

    private static func looksLikeRDF(_ text: String) -> Bool {
        text.localizedCaseInsensitiveContains("$RDFILE") && looksLikeRXN(text)
    }

    private static func looksLikePDB(_ text: String) -> Bool {
        var checked = 0
        for raw in text.components(separatedBy: .newlines) {
            let line = raw.trimmingCharacters(in: .whitespacesAndNewlines)
            if line.isEmpty { continue }
            checked += 1
            if line.hasPrefix("ATOM") || line.hasPrefix("HETATM") || line.hasPrefix("HEADER") {
                return true
            }
            if checked >= 8 { break }
        }
        return false
    }
}
