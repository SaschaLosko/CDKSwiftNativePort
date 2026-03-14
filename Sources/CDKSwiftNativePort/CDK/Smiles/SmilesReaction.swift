import Foundation

/// Lightweight Swift counterpart to CDK's reaction container for SMILES parsing.
public enum CDKReactionRole: String, CaseIterable, Codable, Hashable, Sendable {
    case reactant
    case agent
    case product
}

public struct CDKReactionParticipant: Equatable, Hashable, Codable, Sendable {
    public var molecule: Molecule
    public var role: CDKReactionRole
    public var stoichiometry: Double?

    public init(molecule: Molecule, role: CDKReactionRole, stoichiometry: Double? = nil) {
        self.molecule = molecule
        self.role = role
        self.stoichiometry = stoichiometry
    }
}

public struct CDKReaction: Equatable {
    public var id: String? = nil
    public var reactantParticipants: [CDKReactionParticipant]
    public var agentParticipants: [CDKReactionParticipant]
    public var productParticipants: [CDKReactionParticipant]
    public var name: String? = nil
    public var properties: [String: String] = [:]
    public var cxState: CDKCxSmilesState? = nil

    public var reactants: [Molecule] {
        get { reactantParticipants.map(\.molecule) }
        set {
            reactantParticipants = newValue.map { molecule in
                CDKReactionParticipant(molecule: molecule, role: .reactant)
            }
        }
    }

    public var agents: [Molecule] {
        get { agentParticipants.map(\.molecule) }
        set {
            agentParticipants = newValue.map { molecule in
                CDKReactionParticipant(molecule: molecule, role: .agent)
            }
        }
    }

    public var products: [Molecule] {
        get { productParticipants.map(\.molecule) }
        set {
            productParticipants = newValue.map { molecule in
                CDKReactionParticipant(molecule: molecule, role: .product)
            }
        }
    }

    public var reactantCount: Int { reactants.count }
    public var agentCount: Int { agents.count }
    public var productCount: Int { products.count }
    public var participants: [CDKReactionParticipant] {
        reactantParticipants + agentParticipants + productParticipants
    }

    public init(reactants: [Molecule],
                agents: [Molecule],
                products: [Molecule],
                id: String? = nil,
                name: String? = nil,
                properties: [String: String] = [:],
                cxState: CDKCxSmilesState? = nil) {
        self.id = id
        self.reactantParticipants = reactants.map { molecule in
            CDKReactionParticipant(molecule: molecule, role: .reactant)
        }
        self.agentParticipants = agents.map { molecule in
            CDKReactionParticipant(molecule: molecule, role: .agent)
        }
        self.productParticipants = products.map { molecule in
            CDKReactionParticipant(molecule: molecule, role: .product)
        }
        self.name = name
        self.properties = properties
        self.cxState = cxState
    }

    public init(reactantParticipants: [CDKReactionParticipant],
                agentParticipants: [CDKReactionParticipant],
                productParticipants: [CDKReactionParticipant],
                id: String? = nil,
                name: String? = nil,
                properties: [String: String] = [:],
                cxState: CDKCxSmilesState? = nil) {
        self.id = id
        self.reactantParticipants = Self.normalizedParticipants(reactantParticipants, role: .reactant)
        self.agentParticipants = Self.normalizedParticipants(agentParticipants, role: .agent)
        self.productParticipants = Self.normalizedParticipants(productParticipants, role: .product)
        self.name = name
        self.properties = properties
        self.cxState = cxState
    }

    public init(participants: [CDKReactionParticipant],
                id: String? = nil,
                name: String? = nil,
                properties: [String: String] = [:],
                cxState: CDKCxSmilesState? = nil) {
        self.id = id
        self.reactantParticipants = participants
            .filter { $0.role == .reactant }
        self.agentParticipants = participants
            .filter { $0.role == .agent }
        self.productParticipants = participants
            .filter { $0.role == .product }
        self.name = name
        self.properties = properties
        self.cxState = cxState
    }

    private static func normalizedParticipants(_ participants: [CDKReactionParticipant],
                                               role: CDKReactionRole) -> [CDKReactionParticipant] {
        participants.map { participant in
            CDKReactionParticipant(molecule: participant.molecule,
                                   role: role,
                                   stoichiometry: participant.stoichiometry)
        }
    }
}
