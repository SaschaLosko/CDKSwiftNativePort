import Foundation

/// CDK-style MDL RXN writer with V2000 and V3000 support.
public enum CDKRXNWriter {
    public struct Options {
        public var alwaysV3000: Bool
        public var programName: String?

        public init(alwaysV3000: Bool = false,
                    programName: String? = nil) {
            self.alwaysV3000 = alwaysV3000
            self.programName = programName
        }
    }

    public static func write(reactants: [Molecule],
                             products: [Molecule] = [],
                             agents: [Molecule] = [],
                             reactionName: String = "CDKSwiftNativePort Reaction",
                             options: Options = Options()) throws -> String {
        let reaction = CDKReaction(reactants: reactants,
                                   agents: agents,
                                   products: products,
                                   name: reactionName)
        return try write(reaction: reaction, options: options)
    }

    public static func write(reaction: CDKReaction,
                             options: Options = Options()) throws -> String {
        try validateParticipants(in: reaction)
        if shouldWriteV3000(for: reaction, options: options) {
            return try writeV3000(reaction: reaction, options: options)
        }
        return try writeV2000(reaction: reaction, options: options)
    }

    public static func write(reactions: [CDKReaction],
                             options: Options = Options()) throws -> String {
        guard !reactions.isEmpty else { throw ChemError.emptyInput }
        let blocks = try reactions.map { reaction in
            try write(reaction: reaction, options: options).trimmingCharacters(in: .whitespacesAndNewlines)
        }
        return blocks.joined(separator: "\n") + "\n"
    }

    private static func writeV2000(reaction: CDKReaction,
                                   options: Options) throws -> String {
        var lines: [String] = []
        lines.append("$RXN")
        lines.append(normalizedName(reaction.name, fallback: "CDKSwiftNativePort Reaction"))
        lines.append(headerProgramLine(programName: options.programName))
        lines.append("")
        lines.append(String(format: "%3d%3d%3d",
                            min(999, reaction.reactantParticipants.count),
                            min(999, reaction.productParticipants.count),
                            min(999, reaction.agentParticipants.count)))

        for participant in reaction.reactantParticipants + reaction.productParticipants + reaction.agentParticipants {
            lines.append(molHeader(for: participant))
            let mol = try CDKMDLV2000Writer.write(participant.molecule,
                                                  options: CDKMDLV2000Writer.Options(programName: options.programName))
                .trimmingCharacters(in: .whitespacesAndNewlines)
            lines.append(mol)
        }

        return lines.joined(separator: "\n") + "\n"
    }

    private static func writeV3000(reaction: CDKReaction,
                                   options: Options) throws -> String {
        var lines: [String] = []
        lines.append("$RXN V3000")
        lines.append(normalizedName(reaction.name, fallback: "CDKSwiftNativePort Reaction"))
        lines.append(headerProgramLine(programName: options.programName))
        lines.append("")
        lines.append("M  V30 COUNTS \(reaction.reactantParticipants.count) \(reaction.productParticipants.count) \(reaction.agentParticipants.count)")
        lines.append(contentsOf: try sectionLines(name: "REACTANT",
                                                  participants: reaction.reactantParticipants,
                                                  options: options))
        lines.append(contentsOf: try sectionLines(name: "PRODUCT",
                                                  participants: reaction.productParticipants,
                                                  options: options))
        if !reaction.agentParticipants.isEmpty {
            lines.append(contentsOf: try sectionLines(name: "AGENT",
                                                      participants: reaction.agentParticipants,
                                                      options: options))
        }
        lines.append("M  END")
        return lines.joined(separator: "\n") + "\n"
    }

    private static func sectionLines(name: String,
                                     participants: [CDKReactionParticipant],
                                     options: Options) throws -> [String] {
        var lines = ["M  V30 BEGIN \(name)"]
        for participant in participants {
            lines.append(contentsOf: try ctabLines(for: participant.molecule, options: options))
        }
        lines.append("M  V30 END \(name)")
        return lines
    }

    private static func ctabLines(for molecule: Molecule,
                                  options: Options) throws -> [String] {
        let molfile = try CDKMDLV3000Writer.write(molecule,
                                                  options: CDKMDLV3000Writer.Options(programName: options.programName,
                                                                                     includeDataFields: false))
        let lines = molfile.components(separatedBy: .newlines)
        guard let begin = lines.firstIndex(where: { $0.trimmingCharacters(in: .whitespacesAndNewlines).uppercased() == "M  V30 BEGIN CTAB" }),
              let end = lines.firstIndex(where: { $0.trimmingCharacters(in: .whitespacesAndNewlines).uppercased() == "M  V30 END CTAB" }),
              end >= begin else {
            throw ChemError.parseFailed("Failed to extract V3000 CTAB while writing RXN.")
        }
        return Array(lines[begin...end])
    }

    private static func shouldWriteV3000(for reaction: CDKReaction,
                                         options: Options) -> Bool {
        if options.alwaysV3000 { return true }
        return (reaction.reactantParticipants + reaction.agentParticipants + reaction.productParticipants).contains {
            requiresV3000($0.molecule)
        }
    }

    private static func requiresV3000(_ molecule: Molecule) -> Bool {
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

    private static func validateParticipants(in reaction: CDKReaction) throws {
        let total = reaction.reactantParticipants.count + reaction.productParticipants.count + reaction.agentParticipants.count
        guard total > 0 else {
            throw ChemError.emptyInput
        }
        for participant in reaction.reactantParticipants + reaction.productParticipants + reaction.agentParticipants
            where participant.molecule.atoms.isEmpty {
            throw ChemError.parseFailed("Cannot write RXN with an empty molecule block.")
        }
    }

    private static func molHeader(for participant: CDKReactionParticipant) -> String {
        if let stoichiometry = participant.stoichiometry {
            let normalized = stoichiometry.rounded(.towardZero) == stoichiometry
                ? String(Int(stoichiometry))
                : String(stoichiometry)
            return "$MOL STOICH \(normalized)"
        }
        return "$MOL"
    }

    private static func normalizedName(_ raw: String?, fallback: String) -> String {
        let cleaned = (raw ?? "").replacingOccurrences(of: "\r", with: " ")
            .replacingOccurrences(of: "\n", with: " ")
            .trimmingCharacters(in: .whitespacesAndNewlines)
        return cleaned.isEmpty ? fallback : cleaned
    }

    private static func headerProgramLine(programName: String?) -> String {
        let trimmed = (programName ?? "CDKSwiftNativePort")
            .trimmingCharacters(in: .whitespacesAndNewlines)
        if trimmed.isEmpty {
            return "  CDKSwiftNativePort"
        }
        return "  \(String(trimmed.prefix(16)))"
    }
}
