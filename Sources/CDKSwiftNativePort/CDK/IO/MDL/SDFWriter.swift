import Foundation

public struct CDKSDFWriterOptions {
    public var alwaysV3000: Bool
    public var writeDataFields: Bool
    public var acceptedDataFieldNames: Set<String>?
    public var truncateLongDataFields: Bool
    public var programName: String?

    public init(alwaysV3000: Bool = false,
                writeDataFields: Bool = true,
                acceptedDataFieldNames: Set<String>? = nil,
                truncateLongDataFields: Bool = false,
                programName: String? = nil) {
        self.alwaysV3000 = alwaysV3000
        self.writeDataFields = writeDataFields
        self.acceptedDataFieldNames = acceptedDataFieldNames
        self.truncateLongDataFields = truncateLongDataFields
        self.programName = programName
    }
}

/// CDK-style SDFile writer built on the V2000 mol writer.
public enum CDKSDFWriter {
    public static func write(_ molecule: Molecule,
                             options: CDKSDFWriterOptions = CDKSDFWriterOptions()) throws -> String {
        try write([molecule], options: options)
    }

    public static func write(_ molecules: [Molecule],
                             options: CDKSDFWriterOptions = CDKSDFWriterOptions()) throws -> String {
        guard !molecules.isEmpty else { throw ChemError.emptyInput }

        var chunks: [String] = []
        chunks.reserveCapacity(molecules.count)
        for molecule in molecules {
            let mol = try molfileRecord(for: molecule, options: options)
                .trimmingCharacters(in: .whitespacesAndNewlines)
            let dataFields = options.writeDataFields
                ? CDKSDFDataFieldSupport.serializeDataFields(for: molecule,
                                                             acceptedFieldNames: options.acceptedDataFieldNames,
                                                             truncateLongValues: options.truncateLongDataFields)
                : ""
            let chunk = dataFields.isEmpty
                ? "\(mol)\n$$$$"
                : "\(mol)\n\(dataFields)\n$$$$"
            chunks.append(chunk)
        }

        return chunks.joined(separator: "\n") + "\n"
    }

    private static func molfileRecord(for molecule: Molecule,
                                      options: CDKSDFWriterOptions) throws -> String {
        if shouldWriteV3000(for: molecule, options: options) {
            return try CDKMDLV3000Writer.write(molecule,
                                               options: CDKMDLV3000Writer.Options(programName: options.programName,
                                                                                  includeDataFields: false))
        }

        return try CDKMDLV2000Writer.write(molecule,
                                           options: CDKMDLV2000Writer.Options(programName: options.programName))
    }

    private static func shouldWriteV3000(for molecule: Molecule,
                                         options: CDKSDFWriterOptions) -> Bool {
        if options.alwaysV3000 { return true }
        if molecule.atomCount > 999 || molecule.bondCount > 999 { return true }

        let stereoGroups = molecule.atoms
            .filter { $0.chirality != .none }
            .map { $0.cxStereoGroup ?? CDKCxSmilesParser.encodeStereoGroup(kind: "abs", number: 0) }
        if let firstGroup = stereoGroups.first {
            if stereoGroups.contains(where: { $0 != firstGroup }) {
                return true
            }
            if let decoded = CDKCxSmilesParser.decodeStereoGroup(firstGroup), decoded.kind == "and" {
                return true
            }
        }

        if molecule.sgroups.contains(where: { $0.kind == .extMulticenter }) {
            return true
        }

        return false
    }
}
