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
        CDKFileImporterFormat(displayName: "MDL RGfile",
                              fileExtensions: ["rgf"],
                              utiIdentifiers: ["chemical/x-mdl-rgfile"]),
        CDKFileImporterFormat(displayName: "MDL SDFile",
                              fileExtensions: ["sdf", "sd"],
                              utiIdentifiers: ["chemical/x-mdl-sdfile", "net.sourceforge.openbabel.mdl"]),
        CDKFileImporterFormat(displayName: "CXSMILES",
                              fileExtensions: ["cxsmiles"],
                              utiIdentifiers: []),
        CDKFileImporterFormat(displayName: "SMILES",
                              fileExtensions: ["smi", "smiles", "ism", "can"],
                              utiIdentifiers: ["chemical/x-daylight-smiles", "chemical/x-smiles"]),
        CDKFileImporterFormat(displayName: "Reaction SMILES",
                              fileExtensions: ["rsmi"],
                              utiIdentifiers: ["chemical/x-reaction-smiles"]),
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
        if ["smi", "smiles", "ism", "can", "cxsmiles"].contains(lower) {
            return .smiles
        }
        if lower == "rsmi" {
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

    public static func readMolecules(from url: URL,
                                     coordinateAccess: Bool = true) throws -> [Molecule] {
        let ext = url.pathExtension.lowercased()
        let text = try decodeText(from: url, coordinateAccess: coordinateAccess)
        return try readMolecules(text: text, fileExtension: ext)
    }

    public static func readReaction(from url: URL,
                                    coordinateAccess: Bool = true) throws -> CDKReaction {
        let ext = url.pathExtension.lowercased()
        let text = try decodeText(from: url, coordinateAccess: coordinateAccess)
        return try readReaction(text: text, fileExtension: ext)
    }

    public static func readReactionHierarchy(from url: URL,
                                             coordinateAccess: Bool = true) throws -> CDKReactionHierarchy {
        let ext = url.pathExtension.lowercased()
        let text = try decodeText(from: url, coordinateAccess: coordinateAccess)
        return try readReactionHierarchy(text: text, fileExtension: ext)
    }

    public static func readReactionSet(from url: URL,
                                       coordinateAccess: Bool = true) throws -> CDKReactionSet {
        try readReactionHierarchy(from: url, coordinateAccess: coordinateAccess).asSet
    }

    public static func readMolecules(text: String,
                                     fileExtension: String?) throws -> [Molecule] {
        let ext = (fileExtension ?? "").lowercased()

        switch ext {
        case "sdf", "sd":
            return try CDKIteratingSDFReader.read(text: text)
        case "mol":
            if looksLikeRGFile(text) {
                return [try CDKRGFileReader.readFlattenedMolecule(text: text)]
            }
            return [try CDKMDLReader.read(text: text)]
        case "rgf":
            return [try CDKRGFileReader.readFlattenedMolecule(text: text)]
        case "smi", "smiles", "ism", "can":
            return try CDKSMILESReader.read(text: text)
        case "rsmi":
            return try molecules(from: parseReactionSmilesRecords(text))
        case "inchi", "ich":
            return try CDKInChIReader.read(text: text)
        case "xyz":
            return try CDKXYZReader.read(text: text)
        case "pdb", "ent":
            return try CDKPDBReader.read(text: text)
        case "mol2":
            return try CDKMol2Reader.read(text: text)
        case "cml":
            if CDKCMLReactionReader.containsReactionMarkup(text) {
                return try molecules(from: CDKCMLReactionReader.readReactions(text: text))
            }
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

    public static func readReaction(text: String,
                                    fileExtension: String?) throws -> CDKReaction {
        let ext = (fileExtension ?? "").lowercased()

        switch ext {
        case "cml":
            return try CDKCMLReactionReader.readReaction(text: text)
        case "rxn":
            return try CDKRXNReader.readReaction(text: text)
        case "rdf":
            return try CDKRDFReader.readReaction(text: text)
        case "smi", "smiles", "ism", "can":
            return try firstReaction(from: parseReactionSmilesRecords(text))
        case "rsmi":
            return try firstReaction(from: parseReactionSmilesRecords(text))
        case "txt":
            return try readTextReactionWithAutoDetection(text)
        default:
            return try readTextReactionWithAutoDetection(text)
        }
    }

    public static func readReactionHierarchy(text: String,
                                             fileExtension: String?) throws -> CDKReactionHierarchy {
        let ext = (fileExtension ?? "").lowercased()

        switch ext {
        case "cml":
            return try CDKCMLReactionReader.readHierarchy(text: text)
        case "rxn":
            return hierarchy(from: try CDKRXNReader.readReactions(text: text))
        case "rdf":
            return hierarchy(from: try CDKRDFReader.readReactions(text: text))
        case "smi", "smiles", "ism", "can", "rsmi":
            return hierarchy(from: try parseReactionSmilesRecords(text))
        case "txt":
            return try readTextReactionHierarchyWithAutoDetection(text)
        default:
            return try readTextReactionHierarchyWithAutoDetection(text)
        }
    }

    public static func readReactionSet(text: String,
                                       fileExtension: String?) throws -> CDKReactionSet {
        try readReactionHierarchy(text: text, fileExtension: fileExtension).asSet
    }

    private static func readTextWithAutoDetection(_ text: String) throws -> [Molecule] {
        if looksLikeInChI(text) {
            return try CDKInChIReader.read(text: text)
        }
        if looksLikeRGFile(text) {
            return [try CDKRGFileReader.readFlattenedMolecule(text: text)]
        }
        if looksLikeRXN(text) {
            return try CDKRXNReader.read(text: text)
        }
        if looksLikeRDF(text) {
            return try CDKRDFReader.read(text: text)
        }
        if looksLikeReactionSmiles(text) {
            return try molecules(from: parseReactionSmilesRecords(text))
        }
        if looksLikeMol2(text) {
            return try CDKMol2Reader.read(text: text)
        }
        if looksLikeCML(text) {
            if CDKCMLReactionReader.containsReactionMarkup(text) {
                return try molecules(from: CDKCMLReactionReader.readReactions(text: text))
            }
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

    private static func readTextReactionWithAutoDetection(_ text: String) throws -> CDKReaction {
        let hierarchy = try readTextReactionHierarchyWithAutoDetection(text)
        guard let first = hierarchy.flattenedReactions.first else {
            throw ChemError.unsupported("Unable to detect a supported reaction format.")
        }
        return first
    }

    private static func readTextReactionHierarchyWithAutoDetection(_ text: String) throws -> CDKReactionHierarchy {
        if looksLikeRDF(text) {
            return hierarchy(from: try CDKRDFReader.readReactions(text: text))
        }
        if looksLikeRXN(text) {
            return hierarchy(from: try CDKRXNReader.readReactions(text: text))
        }
        if looksLikeCML(text), CDKCMLReactionReader.containsReactionMarkup(text) {
            return try CDKCMLReactionReader.readHierarchy(text: text)
        }
        if looksLikeReactionSmiles(text) {
            return hierarchy(from: try parseReactionSmilesRecords(text))
        }
        throw ChemError.unsupported("Unable to detect a supported reaction format.")
    }

    public static func looksLikeReaction(text: String, fileExtension: String?) -> Bool {
        let ext = (fileExtension ?? "").lowercased()
        if ext == "rxn" || ext == "rdf" || ext == "rsmi" {
            return true
        }
        if ext == "cml" && CDKCMLReactionReader.containsReactionMarkup(text) {
            return true
        }
        if looksLikeRDF(text) || looksLikeRXN(text) {
            return true
        }
        if looksLikeCML(text), CDKCMLReactionReader.containsReactionMarkup(text) {
            return true
        }
        return looksLikeReactionSmiles(text)
    }

    private static func parseReactionSmilesRecords(_ text: String) throws -> [CDKReaction] {
        let lines = meaningfulLines(in: text)
        guard !lines.isEmpty else {
            throw ChemError.emptyInput
        }
        let parser = CDKSmilesParser(flavor: .cdkDefault)
        return try lines.map(parser.parseReactionSmiles)
    }

    private static func firstReaction(from reactions: [CDKReaction]) throws -> CDKReaction {
        guard let first = reactions.first else {
            throw ChemError.emptyInput
        }
        return first
    }

    private static func molecules(from reaction: CDKReaction) -> [Molecule] {
        reaction.reactants + reaction.agents + reaction.products
    }

    private static func molecules(from reactions: [CDKReaction]) -> [Molecule] {
        reactions.flatMap(molecules(from:))
    }

    private static func hierarchy(from reactions: [CDKReaction]) -> CDKReactionHierarchy {
        if reactions.count == 1, let reaction = reactions.first {
            return .reaction(reaction)
        }
        return .set(CDKReactionSet(reactions: reactions))
    }

    private static func decodeText(from url: URL,
                                   coordinateAccess: Bool = true) throws -> String {
        try CDKFileAccess.decodeText(from: url,
                                     preferredEncodings: [.utf8],
                                     coordinateAccess: coordinateAccess)
    }

    private static func firstMeaningfulLine(in text: String) -> String? {
        meaningfulLines(in: text).first
    }

    private static func meaningfulLines(in text: String) -> [String] {
        text.components(separatedBy: .newlines).compactMap { raw in
            let line = raw.trimmingCharacters(in: .whitespacesAndNewlines)
            if line.isEmpty || line.hasPrefix("#") || line.hasPrefix("//") {
                return nil
            }
            return line
        }
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
            $0.trimmingCharacters(in: .whitespacesAndNewlines).uppercased().hasPrefix("$RXN")
        }
    }

    private static func looksLikeRGFile(_ text: String) -> Bool {
        let meaningfulLines = text.components(separatedBy: .newlines)
            .map { $0.trimmingCharacters(in: .whitespacesAndNewlines) }
            .filter { !$0.isEmpty }
        guard meaningfulLines.count >= 2 else { return false }
        return meaningfulLines[0].hasPrefix("$MDL") && meaningfulLines[1] == "$MOL"
    }

    private static func looksLikeRDF(_ text: String) -> Bool {
        text.localizedCaseInsensitiveContains("$RDFILE") && looksLikeRXN(text)
    }

    private static func looksLikeReactionSmiles(_ text: String) -> Bool {
        guard let line = firstMeaningfulLine(in: text) else { return false }
        if line.hasPrefix("InChI=") { return false }
        let separatorCount = line.filter { $0 == ">" }.count
        return separatorCount == 2
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
