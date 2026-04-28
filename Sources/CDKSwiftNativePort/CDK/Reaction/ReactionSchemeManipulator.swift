import Foundation

public enum CDKReactionSchemeManipulator {
    public static func getAllAtomContainers(_ scheme: CDKReactionScheme) -> [Molecule] {
        getAllAtomContainers(scheme, accumulatingInto: [])
    }

    public static func getAllAtomContainers(_ scheme: CDKReactionScheme,
                                            accumulatingInto molecules: [Molecule]) -> [Molecule] {
        var collected = molecules
        appendMolecules(from: scheme, into: &collected)
        return collected
    }

    public static func getAllIDs(_ scheme: CDKReactionScheme) -> [String] {
        var ids: [String] = []
        appendIDs(from: scheme, into: &ids)
        return ids
    }

    public static func getAllReactions(_ scheme: CDKReactionScheme) -> CDKReactionSet {
        CDKReactionSet(id: scheme.id,
                       name: scheme.name,
                       members: scheme.flattenedReactions.map(CDKReactionSetMember.reaction),
                       properties: scheme.properties)
    }

    public static func createReactionScheme(_ reactionSet: CDKReactionSet) -> CDKReactionScheme {
        let reactions = reactionSet.flattenedReactions
        let topReactions = reactions.filter { directPrecursors(of: $0, in: reactions).isEmpty }

        var entries: [CDKReactionSchemeEntry] = []
        for reaction in topReactions {
            entries.append(.reaction(reaction))
            let nested = buildScheme(after: reaction, in: reactions)
            if !nested.entries.isEmpty {
                entries.append(.scheme(nested))
            }
        }

        return CDKReactionScheme(id: reactionSet.id,
                                 name: reactionSet.name,
                                 entries: entries,
                                 properties: reactionSet.properties)
    }

    public static func extractTopReactions(_ scheme: CDKReactionScheme) -> CDKReactionSet {
        let allReactions = scheme.flattenedReactions
        let top = allReactions.filter { directPrecursors(of: $0, in: allReactions).isEmpty }
        return CDKReactionSet(reactions: top)
    }

    public static func getAtomContainerSet(originMol: Molecule,
                                           finalMol: Molecule,
                                           reactionScheme: CDKReactionScheme) -> [[Molecule]] {
        let reactions = reactionScheme.flattenedReactions
        guard reactions.contains(where: { $0.reactants.contains(originMol) }) else { return [] }

        var paths: [[Molecule]] = []

        func dfs(current: Molecule, path: [Molecule], visitedReactionIndices: Set<Int>) {
            if current == finalMol {
                paths.append(path)
                return
            }

            for (index, reaction) in reactions.enumerated() where !visitedReactionIndices.contains(index) {
                guard reaction.reactants.contains(current) else { continue }

                var nextVisited = visitedReactionIndices
                nextVisited.insert(index)

                for product in reaction.products {
                    if path.contains(product) { continue }
                    dfs(current: product, path: path + [product], visitedReactionIndices: nextVisited)
                }
            }
        }

        dfs(current: originMol, path: [originMol], visitedReactionIndices: [])
        return paths
    }

    private static func appendMolecules(from scheme: CDKReactionScheme, into molecules: inout [Molecule]) {
        for entry in scheme.entries {
            appendMolecules(from: entry, into: &molecules)
        }
    }

    private static func appendMolecules(from list: CDKReactionList, into molecules: inout [Molecule]) {
        for entry in list.entries {
            appendMolecules(from: entry, into: &molecules)
        }
    }

    private static func appendMolecules(from entry: CDKReactionSchemeEntry, into molecules: inout [Molecule]) {
        switch entry {
        case .reaction(let reaction):
            for molecule in CDKReactionManipulator.getAllAtomContainers(reaction) where !molecules.contains(molecule) {
                molecules.append(molecule)
            }
        case .list(let list):
            appendMolecules(from: list, into: &molecules)
        case .scheme(let scheme):
            appendMolecules(from: scheme, into: &molecules)
        }
    }

    private static func appendMolecules(from entry: CDKReactionListEntry, into molecules: inout [Molecule]) {
        switch entry {
        case .reaction(let reaction):
            for molecule in CDKReactionManipulator.getAllAtomContainers(reaction) where !molecules.contains(molecule) {
                molecules.append(molecule)
            }
        case .list(let list):
            appendMolecules(from: list, into: &molecules)
        case .scheme(let scheme):
            appendMolecules(from: scheme, into: &molecules)
        }
    }

    private static func appendIDs(from scheme: CDKReactionScheme, into ids: inout [String]) {
        if let id = scheme.id {
            ids.append(id)
        }
        for entry in scheme.entries {
            appendIDs(from: entry, into: &ids)
        }
    }

    private static func appendIDs(from list: CDKReactionList, into ids: inout [String]) {
        if let id = list.id {
            ids.append(id)
        }
        for entry in list.entries {
            appendIDs(from: entry, into: &ids)
        }
    }

    private static func appendIDs(from entry: CDKReactionSchemeEntry, into ids: inout [String]) {
        switch entry {
        case .reaction(let reaction):
            ids.append(contentsOf: CDKReactionManipulator.getAllIDs(reaction))
        case .list(let list):
            appendIDs(from: list, into: &ids)
        case .scheme(let scheme):
            appendIDs(from: scheme, into: &ids)
        }
    }

    private static func appendIDs(from entry: CDKReactionListEntry, into ids: inout [String]) {
        switch entry {
        case .reaction(let reaction):
            ids.append(contentsOf: CDKReactionManipulator.getAllIDs(reaction))
        case .list(let list):
            appendIDs(from: list, into: &ids)
        case .scheme(let scheme):
            appendIDs(from: scheme, into: &ids)
        }
    }

    private static func buildScheme(after reaction: CDKReaction,
                                    in reactions: [CDKReaction]) -> CDKReactionScheme {
        let direct = directSubsequents(of: reaction, in: reactions)
        guard !direct.isEmpty else {
            return CDKReactionScheme(entries: [])
        }

        var entries: [CDKReactionSchemeEntry] = []
        for subsequent in direct {
            entries.append(.reaction(subsequent))
            let nested = buildScheme(after: subsequent, in: reactions)
            if !nested.entries.isEmpty {
                entries.append(.scheme(nested))
            }
        }
        return CDKReactionScheme(entries: entries)
    }

    private static func directPrecursors(of reaction: CDKReaction,
                                         in reactions: [CDKReaction]) -> [CDKReaction] {
        reactions.filter { candidate in
            guard candidate != reaction else { return false }
            return candidate.products.contains { product in
                reaction.reactants.contains(product)
            }
        }
    }

    private static func directSubsequents(of reaction: CDKReaction,
                                          in reactions: [CDKReaction]) -> [CDKReaction] {
        reactions.filter { candidate in
            guard candidate != reaction else { return false }
            return candidate.reactants.contains { reactant in
                reaction.products.contains(reactant)
            }
        }
    }
}
