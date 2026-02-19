import Foundation

/// CDK-style MDL RDF writer wrapping an embedded RXN payload.
public enum CDKRDFWriter {
    public static func write(reactants: [Molecule],
                             products: [Molecule] = [],
                             agents: [Molecule] = [],
                             reactionName: String = "CDKSwiftNativePort RDF") throws -> String {
        let rxnPayload = try CDKRXNWriter.write(reactants: reactants,
                                                products: products,
                                                agents: agents,
                                                reactionName: reactionName)
            .trimmingCharacters(in: .whitespacesAndNewlines)

        let timestamp = rdfTimestamp(Date())
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
