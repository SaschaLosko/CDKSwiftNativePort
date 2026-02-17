import Foundation

/// Lightweight Swift counterpart to CDK's reaction container for SMILES parsing.
public struct CDKReaction: Equatable {
    public var reactants: [Molecule]
    public var agents: [Molecule]
    public var products: [Molecule]
    public var cxState: CDKCxSmilesState? = nil

    public var reactantCount: Int { reactants.count }
    public var agentCount: Int { agents.count }
    public var productCount: Int { products.count }

    public init(reactants: [Molecule], agents: [Molecule], products: [Molecule], cxState: CDKCxSmilesState? = nil) {
        self.reactants = reactants
        self.agents = agents
        self.products = products
        self.cxState = cxState
    }
}
