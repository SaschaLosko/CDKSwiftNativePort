import Foundation

public struct CDKReactionSet: Equatable {
    public var id: String?
    public var name: String?
    public var members: [CDKReactionSetMember]
    public var properties: [String: String]

    public var flattenedReactions: [CDKReaction] {
        members.flatMap(\.flattenedReactions)
    }

    public init(id: String? = nil,
                name: String? = nil,
                members: [CDKReactionSetMember],
                properties: [String: String] = [:]) {
        self.id = id
        self.name = name
        self.members = members
        self.properties = properties
    }

    public init(reactions: [CDKReaction],
                id: String? = nil,
                name: String? = nil,
                properties: [String: String] = [:]) {
        self.init(id: id,
                  name: name,
                  members: reactions.map(CDKReactionSetMember.reaction),
                  properties: properties)
    }
}

public indirect enum CDKReactionSetMember: Equatable {
    case reaction(CDKReaction)
    case list(CDKReactionList)
    case scheme(CDKReactionScheme)

    public var flattenedReactions: [CDKReaction] {
        switch self {
        case .reaction(let reaction):
            return [reaction]
        case .list(let list):
            return list.flattenedReactions
        case .scheme(let scheme):
            return scheme.flattenedReactions
        }
    }
}

public struct CDKReactionList: Equatable {
    public var id: String?
    public var name: String?
    public var entries: [CDKReactionListEntry]
    public var properties: [String: String]
    public var isStepList: Bool

    public var reactions: [CDKReaction] {
        get { entries.flatMap(\.flattenedReactions) }
        set { entries = newValue.map(CDKReactionListEntry.reaction) }
    }

    public var flattenedReactions: [CDKReaction] {
        entries.flatMap(\.flattenedReactions)
    }

    public init(id: String? = nil,
                name: String? = nil,
                entries: [CDKReactionListEntry],
                properties: [String: String] = [:],
                isStepList: Bool = false) {
        self.id = id
        self.name = name
        self.entries = entries
        self.properties = properties
        self.isStepList = isStepList
    }

    public init(id: String? = nil,
                name: String? = nil,
                reactions: [CDKReaction],
                properties: [String: String] = [:],
                isStepList: Bool = false) {
        self.init(id: id,
                  name: name,
                  entries: reactions.map(CDKReactionListEntry.reaction),
                  properties: properties,
                  isStepList: isStepList)
    }
}

public indirect enum CDKReactionListEntry: Equatable {
    case reaction(CDKReaction)
    case list(CDKReactionList)
    case scheme(CDKReactionScheme)

    public var flattenedReactions: [CDKReaction] {
        switch self {
        case .reaction(let reaction):
            return [reaction]
        case .list(let list):
            return list.flattenedReactions
        case .scheme(let scheme):
            return scheme.flattenedReactions
        }
    }
}

public struct CDKReactionScheme: Equatable {
    public var id: String?
    public var name: String?
    public var entries: [CDKReactionSchemeEntry]
    public var properties: [String: String]

    public var flattenedReactions: [CDKReaction] {
        entries.flatMap(\.flattenedReactions)
    }

    public init(id: String? = nil,
                name: String? = nil,
                entries: [CDKReactionSchemeEntry],
                properties: [String: String] = [:]) {
        self.id = id
        self.name = name
        self.entries = entries
        self.properties = properties
    }
}

public indirect enum CDKReactionSchemeEntry: Equatable {
    case reaction(CDKReaction)
    case list(CDKReactionList)
    case scheme(CDKReactionScheme)

    public var flattenedReactions: [CDKReaction] {
        switch self {
        case .reaction(let reaction):
            return [reaction]
        case .list(let list):
            return list.flattenedReactions
        case .scheme(let scheme):
            return scheme.flattenedReactions
        }
    }
}

public enum CDKReactionHierarchy: Equatable {
    case reaction(CDKReaction)
    case list(CDKReactionList)
    case scheme(CDKReactionScheme)
    case set(CDKReactionSet)

    public var flattenedReactions: [CDKReaction] {
        switch self {
        case .reaction(let reaction):
            return [reaction]
        case .list(let list):
            return list.flattenedReactions
        case .scheme(let scheme):
            return scheme.flattenedReactions
        case .set(let set):
            return set.flattenedReactions
        }
    }

    public var asSet: CDKReactionSet {
        switch self {
        case .reaction(let reaction):
            return CDKReactionSet(reactions: [reaction])
        case .list(let list):
            return CDKReactionSet(members: [.list(list)])
        case .scheme(let scheme):
            return CDKReactionSet(members: [.scheme(scheme)])
        case .set(let set):
            return set
        }
    }
}
