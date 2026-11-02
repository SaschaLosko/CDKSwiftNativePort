import Foundation

extension CDKFiniteMarkushDefinition {
    private enum CodingKeys: String, CodingKey {
        case schemaVersion, sourceIdentifier, backbone, reviewed, reviewedDefinition, unresolvedConditions, variables
    }

    public init(from decoder: Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        let version = try container.decodeIfPresent(Int.self, forKey: .schemaVersion) ?? 1
        guard version == 1 || version == 2 else {
            throw DecodingError.dataCorruptedError(
                forKey: .schemaVersion, in: container, debugDescription: "Unsupported definition schema.")
        }
        sourceIdentifier = try container.decode(String.self, forKey: .sourceIdentifier)
        backbone = try container.decode(String.self, forKey: .backbone)
        reviewed = try container.decode(Bool.self, forKey: version == 2 ? .reviewedDefinition : .reviewed)
        unresolvedConditions = try container.decode([String].self, forKey: .unresolvedConditions)
        variables = try container.decode([Variable].self, forKey: .variables)
        if version == 1 && hasExtendedDefinitions {
            throw DecodingError.dataCorruptedError(
                forKey: .schemaVersion, in: container, debugDescription: "Extended definitions require schema 2.")
        }
    }

    public func encode(to encoder: Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)
        try container.encode(sourceIdentifier, forKey: .sourceIdentifier)
        try container.encode(backbone, forKey: .backbone)
        try container.encode(unresolvedConditions, forKey: .unresolvedConditions)
        try container.encode(variables, forKey: .variables)
        if hasExtendedDefinitions {
            try container.encode(2, forKey: .schemaVersion)
            // Legacy readers ignore unknown keys. Keep their known review gate closed,
            // rather than allowing them to execute a silently narrowed/widened definition.
            try container.encode(false, forKey: .reviewed)
            try container.encode(reviewed, forKey: .reviewedDefinition)
        } else {
            try container.encode(reviewed, forKey: .reviewed)
        }
    }

    private var hasExtendedDefinitions: Bool {
        variables.contains { !($0.chemicalClasses ?? []).isEmpty || !($0.localDefinitions ?? []).isEmpty }
    }
}
