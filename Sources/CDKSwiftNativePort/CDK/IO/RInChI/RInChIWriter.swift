import Foundation

/// Line-oriented RInChI writer mirroring the text-based reaction writer behavior.
public enum CDKRInChIWriter {
    public static func write(_ reactions: [CDKReaction]) throws -> String {
        guard !reactions.isEmpty else { throw ChemError.emptyInput }

        var lines: [String] = []
        lines.reserveCapacity(reactions.count)

        for reaction in reactions {
            let generator = CDKRInChIGenerator().generate(reaction)
            let rinchi = generator.getRInChI()?.trimmingCharacters(in: .whitespacesAndNewlines) ?? ""
            if rinchi.isEmpty {
                let message = generator.getMessages().joined(separator: " ").trimmingCharacters(in: .whitespacesAndNewlines)
                if message.isEmpty {
                    throw ChemError.parseFailed("RInChI generator produced an empty line.")
                }
                throw ChemError.unsupported("RInChI export failed: \(message)")
            }

            let name = normalizedName(reaction.name ?? "")
            if name.isEmpty {
                lines.append(rinchi)
            } else {
                lines.append("\(rinchi) \(name)")
            }
        }

        return lines.joined(separator: "\n") + "\n"
    }

    private static func normalizedName(_ raw: String) -> String {
        raw.replacingOccurrences(of: "\r", with: " ")
            .replacingOccurrences(of: "\n", with: " ")
            .trimmingCharacters(in: .whitespacesAndNewlines)
    }
}
