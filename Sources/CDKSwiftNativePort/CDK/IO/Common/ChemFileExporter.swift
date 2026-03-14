import Foundation
import CoreGraphics

public enum CDKFileExportFormat: String, CaseIterable, Identifiable {
    case mol
    case molV3000
    case sdf
    case smiles
    case isomericSmiles
    case inchi
    case mol2
    case pdb
    case xyz
    case cml
    case rxn
    case rdf
    case svg

    public var id: String { rawValue }
}

public struct CDKFileExporterFormat: Hashable, Identifiable {
    public let format: CDKFileExportFormat
    public let displayName: String
    public let fileExtensions: [String]
    public let utiIdentifiers: [String]
    public let supportsMultipleMolecules: Bool

    public var id: CDKFileExportFormat { format }
    public var primaryFileExtension: String { fileExtensions.first ?? "txt" }

    public init(format: CDKFileExportFormat,
                displayName: String,
                fileExtensions: [String],
                utiIdentifiers: [String],
                supportsMultipleMolecules: Bool = true) {
        self.format = format
        self.displayName = displayName
        self.fileExtensions = fileExtensions
        self.utiIdentifiers = utiIdentifiers
        self.supportsMultipleMolecules = supportsMultipleMolecules
    }
}

public struct CDKFileExportOptions {
    public var smilesFlavor: CDKSmiFlavor
    public var isomericSmilesFlavor: CDKSmiFlavor
    public var sdfOptions: CDKSDFWriterOptions
    public var renderStyle: RenderStyle
    public var svgCanvasSize: CGSize
    public var svgIncludeBackground: Bool

    public init(smilesFlavor: CDKSmiFlavor = [.useAromaticSymbols, .strict],
                isomericSmilesFlavor: CDKSmiFlavor = [.useAromaticSymbols, .isomeric, .strict],
                sdfOptions: CDKSDFWriterOptions = CDKSDFWriterOptions(),
                renderStyle: RenderStyle = RenderStyle(),
                svgCanvasSize: CGSize = CGSize(width: 1400, height: 920),
                svgIncludeBackground: Bool = true) {
        self.smilesFlavor = smilesFlavor
        self.isomericSmilesFlavor = isomericSmilesFlavor
        self.sdfOptions = sdfOptions
        self.renderStyle = renderStyle
        self.svgCanvasSize = svgCanvasSize
        self.svgIncludeBackground = svgIncludeBackground
    }
}

/// Unified file exporter dispatch for CDKSwiftNativePort-backed formats.
public enum CDKFileExporter {
    public static let formats: [CDKFileExporterFormat] = [
        CDKFileExporterFormat(format: .mol,
                              displayName: "MDL Molfile (V2000)",
                              fileExtensions: ["mol"],
                              utiIdentifiers: ["chemical/x-mdl-molfile", "net.sourceforge.openbabel.mdl"],
                              supportsMultipleMolecules: false),
        CDKFileExporterFormat(format: .molV3000,
                              displayName: "MDL Molfile (V3000)",
                              fileExtensions: ["mol"],
                              utiIdentifiers: ["chemical/x-mdl-molfile", "net.sourceforge.openbabel.mdl"],
                              supportsMultipleMolecules: false),
        CDKFileExporterFormat(format: .sdf,
                              displayName: "MDL SDFile",
                              fileExtensions: ["sdf", "sd"],
                              utiIdentifiers: ["chemical/x-mdl-sdfile", "net.sourceforge.openbabel.mdl"]),
        CDKFileExporterFormat(format: .smiles,
                              displayName: "SMILES",
                              fileExtensions: ["smi", "smiles", "can"],
                              utiIdentifiers: ["chemical/x-daylight-smiles", "chemical/x-smiles"]),
        CDKFileExporterFormat(format: .isomericSmiles,
                              displayName: "SMILES (Isomeric)",
                              fileExtensions: ["ism"],
                              utiIdentifiers: ["chemical/x-daylight-smiles", "chemical/x-smiles"]),
        CDKFileExporterFormat(format: .inchi,
                              displayName: "InChI",
                              fileExtensions: ["inchi", "ich"],
                              utiIdentifiers: ["chemical/x-inchi"]),
        CDKFileExporterFormat(format: .mol2,
                              displayName: "Tripos MOL2",
                              fileExtensions: ["mol2"],
                              utiIdentifiers: ["chemical/x-mol2"]),
        CDKFileExporterFormat(format: .pdb,
                              displayName: "Protein Data Bank",
                              fileExtensions: ["pdb", "ent"],
                              utiIdentifiers: ["chemical/x-pdb"]),
        CDKFileExporterFormat(format: .xyz,
                              displayName: "XYZ Coordinates",
                              fileExtensions: ["xyz"],
                              utiIdentifiers: ["chemical/x-xyz"]),
        CDKFileExporterFormat(format: .cml,
                              displayName: "Chemical Markup Language",
                              fileExtensions: ["cml"],
                              utiIdentifiers: ["chemical/x-cml"]),
        CDKFileExporterFormat(format: .rxn,
                              displayName: "MDL RXN",
                              fileExtensions: ["rxn"],
                              utiIdentifiers: ["chemical/x-mdl-rxnfile"]),
        CDKFileExporterFormat(format: .rdf,
                              displayName: "MDL RDF",
                              fileExtensions: ["rdf"],
                              utiIdentifiers: ["chemical/x-mdl-rdfile"]),
        CDKFileExporterFormat(format: .svg,
                              displayName: "SVG (Depiction)",
                              fileExtensions: ["svg"],
                              utiIdentifiers: ["public.svg-image"],
                              supportsMultipleMolecules: false)
    ]

    public static var supportedFileExtensions: [String] {
        Array(Set(formats.flatMap(\.fileExtensions))).sorted()
    }

    public static var supportedUTIIdentifiers: [String] {
        Array(Set(formats.flatMap(\.utiIdentifiers))).sorted()
    }

    public static func format(forFileExtension ext: String) -> CDKFileExportFormat? {
        let lower = ext.lowercased()
        for format in formats where format.fileExtensions.contains(lower) {
            return format.format
        }
        return nil
    }

    public static func format(forUTIIdentifier uti: String) -> CDKFileExportFormat? {
        for format in formats where format.utiIdentifiers.contains(uti) {
            return format.format
        }
        return nil
    }

    public static func write(molecule: Molecule,
                             as format: CDKFileExportFormat,
                             options: CDKFileExportOptions = CDKFileExportOptions()) throws -> String {
        try write(molecules: [molecule], as: format, options: options)
    }

    public static func write(molecules: [Molecule],
                             as format: CDKFileExportFormat,
                             options: CDKFileExportOptions = CDKFileExportOptions()) throws -> String {
        guard !molecules.isEmpty else { throw ChemError.emptyInput }

        switch format {
        case .mol:
            guard let first = molecules.first, molecules.count == 1 else {
                throw ChemError.unsupported("Molfile export supports a single molecule only. Use SDF for multiple molecules.")
            }
            return try CDKMDLV2000Writer.write(first)
        case .molV3000:
            guard let first = molecules.first, molecules.count == 1 else {
                throw ChemError.unsupported("V3000 Molfile export supports a single molecule only. Use SDF for multiple molecules.")
            }
            return try CDKMDLV3000Writer.write(first, options: CDKMDLV3000Writer.Options(includeDataFields: false))
        case .sdf:
            return try CDKSDFWriter.write(molecules, options: options.sdfOptions)
        case .smiles:
            return try CDKSMILESWriter.write(molecules, flavor: options.smilesFlavor)
        case .isomericSmiles:
            return try CDKSMILESWriter.write(molecules, flavor: options.isomericSmilesFlavor)
        case .inchi:
            return try CDKInChIWriter.write(molecules)
        case .mol2:
            return try CDKMol2Writer.write(molecules)
        case .pdb:
            return try CDKPDBWriter.write(molecules)
        case .xyz:
            return try CDKXYZWriter.write(molecules)
        case .cml:
            return try CDKCMLWriter.write(molecules)
        case .rxn:
            return try CDKRXNWriter.write(reactants: molecules)
        case .rdf:
            return try CDKRDFWriter.write(reactants: molecules)
        case .svg:
            guard let first = molecules.first, molecules.count == 1 else {
                throw ChemError.unsupported("SVG depiction export supports a single molecule only.")
            }
            return CDKDepictionGenerator.toSVG(molecule: first,
                                               style: options.renderStyle,
                                               canvasSize: options.svgCanvasSize,
                                               includeBackground: options.svgIncludeBackground)
        }
    }

    public static func write(reaction: CDKReaction,
                             as format: CDKFileExportFormat,
                             options: CDKFileExportOptions = CDKFileExportOptions()) throws -> String {
        try write(reactions: [reaction], as: format, options: options)
    }

    public static func write(reactions: [CDKReaction],
                             as format: CDKFileExportFormat,
                             options: CDKFileExportOptions = CDKFileExportOptions()) throws -> String {
        guard !reactions.isEmpty else { throw ChemError.emptyInput }

        switch format {
        case .cml:
            return try reactions.count == 1
                ? CDKCMLReactionWriter.write(reactions[0])
                : CDKCMLReactionWriter.write(reactions)
        case .rxn:
            guard reactions.count == 1, let reaction = reactions.first else {
                throw ChemError.unsupported("RXN export supports a single reaction only.")
            }
            return try CDKRXNWriter.write(reactants: reaction.reactants,
                                          products: reaction.products,
                                          agents: reaction.agents,
                                          reactionName: reaction.name ?? reaction.id ?? "CDKSwiftNativePort Reaction")
        case .rdf:
            if reactions.count == 1, let reaction = reactions.first {
                return try CDKRDFWriter.write(reactants: reaction.reactants,
                                              products: reaction.products,
                                              agents: reaction.agents,
                                              reactionName: reaction.name ?? reaction.id ?? "CDKSwiftNativePort RDF")
            }

            return try reactions.enumerated().map { index, reaction in
                try CDKRDFWriter.write(reactants: reaction.reactants,
                                       products: reaction.products,
                                       agents: reaction.agents,
                                       reactionName: reaction.name ?? reaction.id ?? "CDKSwiftNativePort RDF \(index + 1)")
                    .trimmingCharacters(in: .whitespacesAndNewlines)
            }.joined(separator: "\n")
                + "\n"
        case .smiles, .isomericSmiles:
            let generator = CDKSmilesGenerator(flavor: format == .smiles ? options.smilesFlavor : options.isomericSmilesFlavor)
            return reactions.map(generator.create).joined(separator: "\n") + "\n"
        default:
            throw ChemError.unsupported("Format \(format.rawValue) does not support reaction export.")
        }
    }

    public static func write(molecule: Molecule,
                             to url: URL,
                             as format: CDKFileExportFormat? = nil,
                             options: CDKFileExportOptions = CDKFileExportOptions()) throws {
        try write(molecules: [molecule], to: url, as: format, options: options)
    }

    public static func write(molecules: [Molecule],
                             to url: URL,
                             as format: CDKFileExportFormat? = nil,
                             options: CDKFileExportOptions = CDKFileExportOptions()) throws {
        let chosenFormat: CDKFileExportFormat
        if let format {
            chosenFormat = format
        } else if let inferred = self.format(forFileExtension: url.pathExtension) {
            chosenFormat = inferred
        } else {
            throw ChemError.unsupported("Unable to infer export format from extension '\(url.pathExtension)'.")
        }

        let text = try write(molecules: molecules, as: chosenFormat, options: options)
        try text.write(to: url, atomically: true, encoding: .utf8)
    }
}
