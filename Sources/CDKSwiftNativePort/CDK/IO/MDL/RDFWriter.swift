import Foundation

/// CDK-style MDL RDF writer wrapping embedded RXN payloads.
public enum CDKRDFWriter {
    public struct Options {
        public var rxnOptions: CDKRXNWriter.Options
        public var fixedTimestamp: Date?

        public init(rxnOptions: CDKRXNWriter.Options = CDKRXNWriter.Options(),
                    fixedTimestamp: Date? = nil) {
            self.rxnOptions = rxnOptions
            self.fixedTimestamp = fixedTimestamp
        }
    }

    public static func write(reactants: [Molecule],
                             products: [Molecule] = [],
                             agents: [Molecule] = [],
                             reactionName: String = "CDKSwiftNativePort RDF",
                             options: Options = Options()) throws -> String {
        let reaction = CDKReaction(reactants: reactants,
                                   agents: agents,
                                   products: products,
                                   name: reactionName)
        return try write(reaction: reaction, options: options)
    }

    public static func write(reaction: CDKReaction,
                             options: Options = Options()) throws -> String {
        let rxnPayload = try CDKRXNWriter.write(reaction: reaction, options: options.rxnOptions)
            .trimmingCharacters(in: .whitespacesAndNewlines)
        return wrapRXNPayload(rxnPayload, date: options.fixedTimestamp ?? Date())
    }

    public static func write(reactions: [CDKReaction],
                             options: Options = Options()) throws -> String {
        guard !reactions.isEmpty else { throw ChemError.emptyInput }
        let blocks = try reactions.map { reaction in
            try write(reaction: reaction, options: options).trimmingCharacters(in: .whitespacesAndNewlines)
        }
        return blocks.joined(separator: "\n") + "\n"
    }

    private static func wrapRXNPayload(_ rxnPayload: String, date: Date) -> String {
        let timestamp = rdfTimestamp(date)
        var lines: [String] = []
        lines.append("$RDFILE 1")
        lines.append("$DATM    \(timestamp)")
        lines.append("$RFMT")
        lines.append(rxnPayload)
        return lines.joined(separator: "\n") + "\n"
    }

    private static func rdfTimestamp(_ date: Date) -> String {
        let formatter = DateFormatter()
        formatter.locale = Locale(identifier: "en_US_POSIX")
        formatter.dateFormat = "MM/dd/yyHH:mm"
        return formatter.string(from: date)
    }
}
