import Foundation
#if canImport(FoundationXML)
import FoundationXML
#endif
#if canImport(CoreGraphics)
import CoreGraphics
#endif

public enum CDKCMLReactionReader {
    public static func readReaction(text: String) throws -> CDKReaction {
        let reactions = try readReactions(text: text)
        guard let first = reactions.first else {
            throw ChemError.parseFailed("CML did not contain any reactions.")
        }
        return first
    }

    public static func readReactions(text: String) throws -> [CDKReaction] {
        try readHierarchy(text: text).flattenedReactions
    }

    public static func readHierarchy(text: String) throws -> CDKReactionHierarchy {
        let data = Data(text.utf8)
        let delegate = CMLReactionParserDelegate()

        let parser = XMLParser(data: data)
        parser.delegate = delegate
        guard parser.parse() else {
            let message = parser.parserError?.localizedDescription ?? "Malformed CML document."
            throw ChemError.parseFailed(message)
        }

        guard let hierarchy = delegate.buildHierarchy() else {
            throw ChemError.parseFailed("CML did not contain any reactions.")
        }
        return hierarchy
    }

    public static func readReactionList(text: String) throws -> CDKReactionList {
        let hierarchy = try readHierarchy(text: text)
        guard case .list(let list) = hierarchy else {
            throw ChemError.parseFailed("CML did not contain a top-level reaction list.")
        }
        return list
    }

    public static func readReactionScheme(text: String) throws -> CDKReactionScheme {
        let hierarchy = try readHierarchy(text: text)
        guard case .scheme(let scheme) = hierarchy else {
            throw ChemError.parseFailed("CML did not contain a top-level reaction scheme.")
        }
        return scheme
    }

    public static func readReactionSet(text: String) throws -> CDKReactionSet {
        try readHierarchy(text: text).asSet
    }

    public static func containsReactionMarkup(_ text: String) -> Bool {
        let normalized = text.lowercased()
        return normalized.contains("<reaction")
            || normalized.contains("<reactionlist")
            || normalized.contains("<reactionscheme")
            || normalized.contains("<reactionstep")
    }
}

private final class CMLReactionParserDelegate: NSObject, XMLParserDelegate {
    private final class PendingMolecule {
        var id: String?
        var title: String?
        var formulaValues: [String] = []
        var atoms: [PendingAtom] = []
        var bonds: [PendingBond] = []

        init(id: String? = nil, title: String? = nil) {
            self.id = id
            self.title = title
        }

        func appendFormula(_ raw: String?) {
            guard let raw else { return }
            let trimmed = raw.trimmingCharacters(in: .whitespacesAndNewlines)
            guard !trimmed.isEmpty else { return }
            if !formulaValues.contains(trimmed) {
                formulaValues.append(trimmed)
            }
        }
    }

    private struct PendingAtom {
        let xmlID: String
        var element: String
        var position: CGPoint
        var zPosition: Double?
        var charge: Int
        var isotope: Int?
        var hydrogenCount: Int?
        var aromatic: Bool
        var aliasLabel: String?
    }

    private struct PendingBond {
        var xmlID: String?
        let refA: String
        let refB: String
        var order: BondOrder
        var stereo: BondStereo
        var aromatic: Bool
        var hadExplicitOrder: Bool
    }

    private struct PendingReaction {
        var id: String?
        var name: String?
        var properties: [String: String] = [:]
        var reactants: [PendingMolecule] = []
        var products: [PendingMolecule] = []
        var agents: [PendingMolecule] = []
    }

    private final class PendingReactionList {
        let isStepList: Bool
        var id: String?
        var name: String?
        var properties: [String: String] = [:]
        var entries: [PendingHierarchyNode] = []

        init(id: String?, name: String?, isStepList: Bool) {
            self.id = id
            self.name = name
            self.isStepList = isStepList
        }
    }

    private final class PendingReactionScheme {
        var id: String?
        var name: String?
        var properties: [String: String] = [:]
        var entries: [PendingHierarchyNode] = []

        init(id: String?, name: String?) {
            self.id = id
            self.name = name
        }
    }

    private enum PendingHierarchyNode {
        case reaction(PendingReaction)
        case list(PendingReactionList)
        case scheme(PendingReactionScheme)
    }

    private enum PendingContainerRef {
        case list(PendingReactionList)
        case scheme(PendingReactionScheme)
    }

    private struct PendingParticipant {
        let role: CDKReactionRole
        var molecule: PendingMolecule
    }

    private enum CaptureContext {
        case atom
        case bond
        case reaction
    }

    private struct TextCapture {
        let elementName: String
        let attributes: [String: String]
        let context: CaptureContext
        var text: String
    }

    private struct ParsedCoordinates {
        let position: CGPoint
        let zPosition: Double?
    }

    private var moleculesByID: [String: PendingMolecule] = [:]
    private var orderedMolecules: [PendingMolecule] = []
    private var moleculeStack: [PendingMolecule] = []
    private var participantStack: [PendingParticipant] = []
    private var reactionStack: [PendingReaction] = []
    private var containerStack: [PendingContainerRef] = []
    private var rootNodes: [PendingHierarchyNode] = []
    private var currentAtom: PendingAtom?
    private var currentBond: PendingBond?
    private var textCaptures: [TextCapture] = []
    private var generatedAtomCounter = 1
    private var formulaDepth = 0

    func parser(_ parser: XMLParser,
                didStartElement elementName: String,
                namespaceURI: String?,
                qualifiedName qName: String?,
                attributes attributeDict: [String: String] = [:]) {
        let lower = localName(elementName).lowercased()

        switch lower {
        case "reaction":
            reactionStack.append(PendingReaction(id: normalizedLabel(attributeDict["id"]),
                                                 name: normalizedLabel(attributeDict["title"] ?? attributeDict["name"])))

        case "reactionlist":
            containerStack.append(.list(PendingReactionList(id: normalizedLabel(attributeDict["id"]),
                                                            name: normalizedLabel(attributeDict["title"] ?? attributeDict["name"]),
                                                            isStepList: false)))

        case "reactionsteplist":
            containerStack.append(.list(PendingReactionList(id: normalizedLabel(attributeDict["id"]),
                                                            name: normalizedLabel(attributeDict["title"] ?? attributeDict["name"]),
                                                            isStepList: true)))

        case "reactionscheme":
            containerStack.append(.scheme(PendingReactionScheme(id: normalizedLabel(attributeDict["id"]),
                                                                name: normalizedLabel(attributeDict["title"] ?? attributeDict["name"]))))

        case "reactant":
            ensureCurrentReaction()
            participantStack.append(PendingParticipant(role: .reactant,
                                                      molecule: participantSeedMolecule(attributes: attributeDict)))

        case "product":
            ensureCurrentReaction()
            participantStack.append(PendingParticipant(role: .product,
                                                      molecule: participantSeedMolecule(attributes: attributeDict)))

        case "substance":
            ensureCurrentReaction()
            participantStack.append(PendingParticipant(role: .agent,
                                                      molecule: participantSeedMolecule(attributes: attributeDict)))

        case "molecule":
            let molecule = resolveMolecule(attributes: attributeDict)
            moleculeStack.append(molecule)

        case "formula":
            formulaDepth += 1
            if let molecule = moleculeStack.last {
                molecule.appendFormula(attributeDict["concise"] ?? attributeDict["formula"])
            }

        case "atom":
            guard formulaDepth == 0 else { return }
            ensureCurrentMolecule()
            currentAtom = makeAtom(attributes: attributeDict)

        case "bond":
            guard formulaDepth == 0 else { return }
            ensureCurrentMolecule()
            currentBond = makeBond(attributes: attributeDict)

        case "atomarray":
            guard formulaDepth == 0 else { return }
            ensureCurrentMolecule()
            guard let molecule = moleculeStack.last else { return }
            applyVectorAtomArray(attributes: attributeDict, to: molecule)

        case "bondarray":
            guard formulaDepth == 0 else { return }
            ensureCurrentMolecule()
            guard let molecule = moleculeStack.last else { return }
            applyVectorBondArray(attributes: attributeDict, to: molecule)

        case "bondtype":
            guard currentBond != nil else { return }
            if isAromaticDictRef(attributeDict["dictRef"]) {
                currentBond?.aromatic = true
            }

        case "bondstereo":
            if currentBond != nil {
                textCaptures.append(TextCapture(elementName: lower,
                                                attributes: attributeDict,
                                                context: .bond,
                                                text: ""))
            }

        case "label":
            if currentAtom != nil {
                textCaptures.append(TextCapture(elementName: lower,
                                                attributes: attributeDict,
                                                context: .atom,
                                                text: ""))
            }

        case "scalar", "string", "float", "integer":
            if currentAtom != nil {
                textCaptures.append(TextCapture(elementName: lower,
                                                attributes: attributeDict,
                                                context: .atom,
                                                text: ""))
            } else if currentBond != nil {
                textCaptures.append(TextCapture(elementName: lower,
                                                attributes: attributeDict,
                                                context: .bond,
                                                text: ""))
            } else if !reactionStack.isEmpty || !containerStack.isEmpty {
                textCaptures.append(TextCapture(elementName: lower,
                                                attributes: attributeDict,
                                                context: .reaction,
                                                text: ""))
            }

        default:
            break
        }
    }

    func parser(_ parser: XMLParser, foundCharacters string: String) {
        guard !textCaptures.isEmpty else { return }
        textCaptures[textCaptures.count - 1].text.append(string)
    }

    func parser(_ parser: XMLParser,
                didEndElement elementName: String,
                namespaceURI: String?,
                qualifiedName qName: String?) {
        let lower = localName(elementName).lowercased()

        switch lower {
        case "atom":
            guard formulaDepth == 0, let molecule = moleculeStack.last, let atom = currentAtom else { return }
            upsert(atom: atom, into: molecule)
            currentAtom = nil

        case "bond":
            guard formulaDepth == 0, let molecule = moleculeStack.last, let bond = currentBond else { return }
            upsert(bond: bond, into: molecule)
            currentBond = nil

        case "molecule":
            _ = moleculeStack.popLast()

        case "reactant", "product", "substance":
            guard let participant = participantStack.popLast(), !reactionStack.isEmpty else { return }
            switch participant.role {
            case .reactant:
                reactionStack[reactionStack.count - 1].reactants.append(participant.molecule)
            case .product:
                reactionStack[reactionStack.count - 1].products.append(participant.molecule)
            case .agent:
                reactionStack[reactionStack.count - 1].agents.append(participant.molecule)
            }

        case "reaction":
            guard let reaction = reactionStack.popLast() else { return }
            append(reaction: reaction)

        case "reactionlist", "reactionsteplist":
            guard let container = containerStack.popLast() else { return }
            guard case .list(let list) = container else { return }
            append(list: list)

        case "reactionscheme":
            guard let container = containerStack.popLast() else { return }
            guard case .scheme(let scheme) = container else { return }
            append(scheme: scheme)

        case "formula":
            formulaDepth = max(0, formulaDepth - 1)

        case "bondstereo", "label", "scalar", "string", "float", "integer":
            guard let capture = popCapture(for: lower) else { return }
            applyCapture(capture)

        default:
            break
        }
    }

    func buildHierarchy() -> CDKReactionHierarchy? {
        let hierarchyNodes = rootNodes.map(buildHierarchyNode)
        guard !hierarchyNodes.isEmpty else { return nil }
        if hierarchyNodes.count == 1 {
            return hierarchyNodes[0]
        }
        return .set(CDKReactionSet(members: hierarchyNodes.map(asSetMember)))
    }

    private func buildHierarchyNode(_ node: PendingHierarchyNode) -> CDKReactionHierarchy {
        switch node {
        case .reaction(let reaction):
            return .reaction(buildReaction(from: reaction, index: 1))
        case .list(let list):
            return .list(buildReactionList(from: list))
        case .scheme(let scheme):
            return .scheme(buildReactionScheme(from: scheme))
        }
    }

    private func buildReaction(from pending: PendingReaction, index: Int) -> CDKReaction {
        let reactants = pending.reactants.enumerated().map {
            buildMolecule(from: $0.element,
                          allowEmpty: true,
                          fallbackName: "Reactant \($0.offset + 1)")
        }
        let products = pending.products.enumerated().map {
            buildMolecule(from: $0.element,
                          allowEmpty: true,
                          fallbackName: "Product \($0.offset + 1)")
        }
        let agents = pending.agents.enumerated().map {
            buildMolecule(from: $0.element,
                          allowEmpty: true,
                          fallbackName: "Agent \($0.offset + 1)")
        }
        return CDKReaction(reactants: reactants,
                           agents: agents,
                           products: products,
                           id: pending.id,
                           name: pending.name,
                           properties: pending.properties)
    }

    private func buildReactionList(from pending: PendingReactionList) -> CDKReactionList {
        let entries = pending.entries.enumerated().map { offset, entry in
            buildReactionListEntry(entry, defaultIndex: offset + 1)
        }
        return CDKReactionList(id: pending.id,
                               name: pending.name,
                               entries: entries,
                               properties: pending.properties,
                               isStepList: pending.isStepList)
    }

    private func buildReactionScheme(from pending: PendingReactionScheme) -> CDKReactionScheme {
        let entries = pending.entries.enumerated().map { offset, entry in
            buildReactionSchemeEntry(entry, defaultIndex: offset + 1)
        }
        return CDKReactionScheme(id: pending.id,
                                 name: pending.name,
                                 entries: entries,
                                 properties: pending.properties)
    }

    private func buildReactionListEntry(_ entry: PendingHierarchyNode,
                                        defaultIndex: Int) -> CDKReactionListEntry {
        switch entry {
        case .reaction(let reaction):
            return .reaction(buildReaction(from: reaction, index: defaultIndex))
        case .list(let list):
            return .list(buildReactionList(from: list))
        case .scheme(let scheme):
            return .scheme(buildReactionScheme(from: scheme))
        }
    }

    private func buildReactionSchemeEntry(_ entry: PendingHierarchyNode,
                                          defaultIndex: Int) -> CDKReactionSchemeEntry {
        switch entry {
        case .reaction(let reaction):
            return .reaction(buildReaction(from: reaction, index: defaultIndex))
        case .list(let list):
            return .list(buildReactionList(from: list))
        case .scheme(let scheme):
            return .scheme(buildReactionScheme(from: scheme))
        }
    }

    private func buildMolecule(from source: PendingMolecule,
                               allowEmpty: Bool,
                               fallbackName: String) -> Molecule {
        let preferredName = source.title?.trimmingCharacters(in: .whitespacesAndNewlines)
        let moleculeName = (preferredName?.isEmpty == false)
            ? preferredName!
            : (source.id?.isEmpty == false ? source.id! : fallbackName)

        var atoms: [Atom] = []
        var atomIDByXMLID: [String: Int] = [:]
        atoms.reserveCapacity(source.atoms.count)

        for (atomIndex, atomSource) in source.atoms.enumerated() {
            let atomID = atomIndex + 1
            atomIDByXMLID[atomSource.xmlID] = atomID
            atoms.append(Atom(id: atomID,
                              externalID: atomSource.xmlID,
                              element: atomSource.element,
                              position: atomSource.position,
                              zPosition: atomSource.zPosition,
                              charge: atomSource.charge,
                              isotopeMassNumber: atomSource.isotope,
                              aromatic: atomSource.aromatic,
                              explicitHydrogenCount: atomSource.hydrogenCount,
                              aliasLabel: atomSource.aliasLabel))
        }

        var bonds: [Bond] = []
        var aromaticAtomIDs = Set(atoms.filter(\.aromatic).map(\.id))
        var nextBondID = 1
        for bondSource in source.bonds {
            guard let a1 = atomIDByXMLID[bondSource.refA],
                  let a2 = atomIDByXMLID[bondSource.refB],
                  a1 != a2 else {
                continue
            }
            bonds.append(Bond(id: nextBondID,
                              externalID: bondSource.xmlID,
                              a1: a1,
                              a2: a2,
                              order: bondSource.order,
                              stereo: bondSource.stereo))
            if bondSource.aromatic || bondSource.order == .aromatic {
                aromaticAtomIDs.insert(a1)
                aromaticAtomIDs.insert(a2)
            }
            nextBondID += 1
        }

        if bonds.isEmpty, !atoms.isEmpty {
            bonds = CDKBondPerception.inferSingleBonds(for: atoms)
        }

        for index in atoms.indices where aromaticAtomIDs.contains(atoms[index].id) {
            atoms[index].aromatic = true
        }

        var molecule = Molecule(name: moleculeName,
                                externalID: source.id,
                                atoms: atoms,
                                bonds: bonds)
        for formula in source.formulaValues {
            molecule.appendDataFieldValue(formula, named: "Formula")
        }

        if !atoms.isEmpty,
           let box = molecule.boundingBox(),
           box.width <= 0.0001 && box.height <= 0.0001 {
            molecule = Depiction2DGenerator.generate(for: molecule)
        }

        if atoms.isEmpty && !allowEmpty && molecule.dataFields.isEmpty {
            return Molecule(name: moleculeName, externalID: source.id)
        }

        return molecule
    }

    private func ensureCurrentReaction() {
        if reactionStack.isEmpty {
            reactionStack.append(PendingReaction())
        }
    }

    private func append(reaction: PendingReaction) {
        if let lastContainer = containerStack.last {
            switch lastContainer {
            case .list(let list):
                list.entries.append(.reaction(reaction))
            case .scheme(let scheme):
                scheme.entries.append(.reaction(reaction))
            }
            return
        }
        rootNodes.append(.reaction(reaction))
    }

    private func append(list: PendingReactionList) {
        if let lastContainer = containerStack.last {
            switch lastContainer {
            case .list(let parentList):
                parentList.entries.append(.list(list))
            case .scheme(let scheme):
                scheme.entries.append(.list(list))
            }
            return
        }
        rootNodes.append(.list(list))
    }

    private func append(scheme: PendingReactionScheme) {
        if let lastContainer = containerStack.last {
            switch lastContainer {
            case .list(let list):
                list.entries.append(.scheme(scheme))
            case .scheme(let parentScheme):
                parentScheme.entries.append(.scheme(scheme))
            }
            return
        }
        rootNodes.append(.scheme(scheme))
    }

    private func participantSeedMolecule(attributes: [String: String]) -> PendingMolecule {
        if let identifier = normalizedLabel(attributes["ref"] ?? attributes["id"]) {
            return molecule(forID: identifier)
        }
        return PendingMolecule(id: nil, title: nil)
    }

    private func resolveMolecule(attributes: [String: String]) -> PendingMolecule {
        let explicitID = normalizedLabel(attributes["id"])
        let referenceID = normalizedLabel(attributes["ref"])

        let molecule: PendingMolecule
        if let referenceID {
            molecule = self.molecule(forID: referenceID)
        } else if let explicitID {
            molecule = self.molecule(forID: explicitID)
        } else if let participant = participantStack.last {
            molecule = participant.molecule
        } else if let current = moleculeStack.last {
            molecule = current
        } else {
            molecule = PendingMolecule()
        }

        if molecule.id == nil {
            molecule.id = explicitID ?? referenceID
        }
        if let title = normalizedLabel(attributes["title"] ?? attributes["name"]) {
            molecule.title = title
        }
        molecule.appendFormula(attributes["formula"])
        registerMoleculeIfNeeded(molecule)

        if !participantStack.isEmpty {
            participantStack[participantStack.count - 1].molecule = molecule
        }
        return molecule
    }

    private func ensureCurrentMolecule() {
        if moleculeStack.isEmpty {
            let molecule = participantStack.last?.molecule ?? PendingMolecule()
            registerMoleculeIfNeeded(molecule)
            moleculeStack.append(molecule)
        }
    }

    private func molecule(forID id: String) -> PendingMolecule {
        if let existing = moleculesByID[id] {
            return existing
        }
        let molecule = PendingMolecule(id: id, title: nil)
        moleculesByID[id] = molecule
        orderedMolecules.append(molecule)
        return molecule
    }

    private func registerMoleculeIfNeeded(_ molecule: PendingMolecule) {
        if let id = molecule.id, moleculesByID[id] == nil {
            moleculesByID[id] = molecule
        }
        let identity = ObjectIdentifier(molecule)
        if !orderedMolecules.contains(where: { ObjectIdentifier($0) == identity }) {
            orderedMolecules.append(molecule)
        }
    }

    private func makeAtom(attributes: [String: String]) -> PendingAtom {
        let rawID = attributes["id"]?.trimmingCharacters(in: .whitespacesAndNewlines)
        let xmlID: String
        if let rawID, !rawID.isEmpty {
            xmlID = rawID
        } else {
            xmlID = "a\(generatedAtomCounter)"
            generatedAtomCounter += 1
        }

        let title = normalizedLabel(attributes["title"])
        let rawElement = attributes["elementType"] ?? attributes["element"] ?? "C"
        let element = normalizedElement(rawElement, aliasLabel: title)
        let coords = coordinates(from: attributes)
        let charge = Int(attributes["formalCharge"] ?? "") ?? 0
        let isotope = Int(attributes["isotopeNumber"] ?? "")
        let hydrogenCount = Int(attributes["hydrogenCount"] ?? "")
        let aromatic = isTruthy(attributes["aromatic"])
            || (attributes["atomType"] ?? "").lowercased().contains(".ar")

        let aliasLabel: String?
        if isPseudoElement(rawElement) {
            aliasLabel = title ?? "R"
        } else {
            aliasLabel = title
        }

        return PendingAtom(xmlID: xmlID,
                           element: element,
                           position: coords.position,
                           zPosition: coords.zPosition,
                           charge: charge,
                           isotope: isotope,
                           hydrogenCount: hydrogenCount,
                           aromatic: aromatic,
                           aliasLabel: aliasLabel)
    }

    private func makeBond(attributes: [String: String]) -> PendingBond? {
        let refs: [String]
        if let atomRefs2 = attributes["atomRefs2"] {
            refs = atomRefs2.split(whereSeparator: \.isWhitespace).map(String.init)
        } else if let atomRefs = attributes["atomRefs"] {
            refs = atomRefs.split(whereSeparator: \.isWhitespace).map(String.init)
        } else if let refA = attributes["atomRef1"], let refB = attributes["atomRef2"] {
            refs = [refA, refB]
        } else {
            refs = []
        }

        guard refs.count >= 2 else { return nil }
        let rawOrder = attributes["order"]
        let order = mapBondOrder(rawOrder)
        return PendingBond(xmlID: normalizedLabel(attributes["id"]),
                           refA: refs[0],
                           refB: refs[1],
                           order: order,
                           stereo: .none,
                           aromatic: order == .aromatic || isTruthy(attributes["aromatic"]),
                           hadExplicitOrder: rawOrder?.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty == false)
    }

    private func upsert(atom: PendingAtom, into molecule: PendingMolecule) {
        if let existingIndex = molecule.atoms.firstIndex(where: { $0.xmlID == atom.xmlID }) {
            molecule.atoms[existingIndex] = atom
        } else {
            molecule.atoms.append(atom)
        }
    }

    private func upsert(bond: PendingBond, into molecule: PendingMolecule) {
        if let xmlID = bond.xmlID,
           let existingIndex = molecule.bonds.firstIndex(where: { $0.xmlID == xmlID }) {
            molecule.bonds[existingIndex] = bond
            return
        }
        if let existingIndex = molecule.bonds.firstIndex(where: {
            ($0.refA == bond.refA && $0.refB == bond.refB) || ($0.refA == bond.refB && $0.refB == bond.refA)
        }) {
            molecule.bonds[existingIndex] = bond
        } else {
            molecule.bonds.append(bond)
        }
    }

    private func applyVectorAtomArray(attributes: [String: String], to molecule: PendingMolecule) {
        guard let atomIDsRaw = attributes["atomID"] else { return }

        let atomIDs = atomIDsRaw.split(whereSeparator: \.isWhitespace).map(String.init)
        guard !atomIDs.isEmpty else { return }

        let elementTypes = splitTokens(attributes["elementType"])
        let x2 = parseDoubles(attributes["x2"])
        let y2 = parseDoubles(attributes["y2"])
        let x3 = parseDoubles(attributes["x3"])
        let y3 = parseDoubles(attributes["y3"])
        let z3 = parseDoubles(attributes["z3"])
        let xy2 = parseCoordinatePairs(attributes["xy2"])
        let xyz3 = parseCoordinateTriplets(attributes["xyz3"])
        let charges = parseInts(attributes["formalCharge"])
        let isotopes = parseInts(attributes["isotopeNumber"])
        let hydrogens = parseInts(attributes["hydrogenCount"])
        let aromaticFlags = parseBooleans(attributes["aromatic"])

        for idx in atomIDs.indices {
            let xmlID = atomIDs[idx]
            if molecule.atoms.contains(where: { $0.xmlID == xmlID }) {
                continue
            }

            let rawElement = idx < elementTypes.count ? elementTypes[idx] : "C"
            let coords = coordinates(at: idx,
                                     xy2: xy2,
                                     xyz3: xyz3,
                                     x2: x2,
                                     y2: y2,
                                     x3: x3,
                                     y3: y3,
                                     z3: z3)

            molecule.atoms.append(PendingAtom(xmlID: xmlID,
                                              element: normalizedElement(rawElement, aliasLabel: nil),
                                              position: coords.position,
                                              zPosition: coords.zPosition,
                                              charge: idx < charges.count ? charges[idx] : 0,
                                              isotope: idx < isotopes.count ? isotopes[idx] : nil,
                                              hydrogenCount: idx < hydrogens.count ? hydrogens[idx] : nil,
                                              aromatic: idx < aromaticFlags.count ? aromaticFlags[idx] : false,
                                              aliasLabel: nil))
        }
    }

    private func applyVectorBondArray(attributes: [String: String], to molecule: PendingMolecule) {
        let orders = splitTokens(attributes["order"])
        let ids = splitTokens(attributes["bondID"])

        if let atomRef1 = attributes["atomRef1"], let atomRef2 = attributes["atomRef2"] {
            let refs1 = atomRef1.split(whereSeparator: \.isWhitespace).map(String.init)
            let refs2 = atomRef2.split(whereSeparator: \.isWhitespace).map(String.init)
            let count = min(refs1.count, refs2.count)
            guard count > 0 else { return }

            for idx in 0..<count {
                let rawOrder = idx < orders.count ? orders[idx] : nil
                let order = mapBondOrder(rawOrder)
                molecule.bonds.append(PendingBond(xmlID: idx < ids.count ? ids[idx] : nil,
                                                  refA: refs1[idx],
                                                  refB: refs2[idx],
                                                  order: order,
                                                  stereo: .none,
                                                  aromatic: order == .aromatic,
                                                  hadExplicitOrder: rawOrder?.isEmpty == false))
            }
            return
        }

        if let atomRefs2 = attributes["atomRefs2"] ?? attributes["atomRefs"] {
            let tokens = atomRefs2.split(whereSeparator: \.isWhitespace).map(String.init)
            guard tokens.count >= 2 else { return }

            var pairIndex = 0
            var idx = 0
            while idx + 1 < tokens.count {
                let rawOrder = pairIndex < orders.count ? orders[pairIndex] : nil
                let order = mapBondOrder(rawOrder)
                molecule.bonds.append(PendingBond(xmlID: pairIndex < ids.count ? ids[pairIndex] : nil,
                                                  refA: tokens[idx],
                                                  refB: tokens[idx + 1],
                                                  order: order,
                                                  stereo: .none,
                                                  aromatic: order == .aromatic,
                                                  hadExplicitOrder: rawOrder?.isEmpty == false))
                pairIndex += 1
                idx += 2
            }
        }
    }

    private func applyCapture(_ capture: TextCapture) {
        switch capture.context {
        case .atom:
            applyAtomCapture(capture)
        case .bond:
            applyBondCapture(capture)
        case .reaction:
            applyReactionCapture(capture)
        }
    }

    private func applyAtomCapture(_ capture: TextCapture) {
        switch capture.elementName {
        case "bondstereo":
            return

        case "label":
            applyAtomLabel(capture.text)

        case "scalar", "string", "float", "integer":
            applyAtomScalarCapture(capture)

        default:
            break
        }
    }

    private func applyBondCapture(_ capture: TextCapture) {
        switch capture.elementName {
        case "bondstereo":
            if let stereo = parseBondStereo(attributes: capture.attributes, text: capture.text) {
                currentBond?.stereo = stereo
            }

        case "scalar", "string", "float", "integer":
            applyBondScalarCapture(capture)

        default:
            break
        }
    }

    private func applyReactionCapture(_ capture: TextCapture) {
        let dictRef = normalizedDictRef(capture.attributes["dictRef"])
        guard dictRef?.hasSuffix("reactionproperty") == true else { return }
        guard let title = normalizedLabel(capture.attributes["title"]) else { return }
        let value = capture.text.trimmingCharacters(in: .whitespacesAndNewlines)
        if !reactionStack.isEmpty {
            reactionStack[reactionStack.count - 1].properties[title] = value
        } else if !containerStack.isEmpty {
            switch containerStack[containerStack.count - 1] {
            case .list(let list):
                list.properties[title] = value
            case .scheme(let scheme):
                scheme.properties[title] = value
            }
        }
    }

    private func asSetMember(_ hierarchy: CDKReactionHierarchy) -> CDKReactionSetMember {
        switch hierarchy {
        case .reaction(let reaction):
            return .reaction(reaction)
        case .list(let list):
            return .list(list)
        case .scheme(let scheme):
            return .scheme(scheme)
        case .set(let set):
            return .list(CDKReactionList(entries: set.members.map(asListEntry)))
        }
    }

    private func asListEntry(_ member: CDKReactionSetMember) -> CDKReactionListEntry {
        switch member {
        case .reaction(let reaction):
            return .reaction(reaction)
        case .list(let list):
            return .list(list)
        case .scheme(let scheme):
            return .scheme(scheme)
        }
    }

    private func applyAtomLabel(_ rawText: String) {
        guard let label = normalizedLabel(rawText), !label.isEmpty else { return }
        currentAtom?.aliasLabel = label
        if let atom = currentAtom, isPseudoElement(atom.element) || atom.element == "*" || atom.element.isEmpty {
            currentAtom?.element = "R"
        }
    }

    private func applyAtomScalarCapture(_ capture: TextCapture) {
        let builtin = scalarBuiltin(capture.attributes)
        let dictRef = normalizedDictRef(capture.attributes["dictRef"])
        let text = capture.text.trimmingCharacters(in: .whitespacesAndNewlines)

        switch builtin ?? dictRef {
        case "label":
            applyAtomLabel(text)

        case "elementtype":
            let aliasLabel = currentAtom?.aliasLabel
            currentAtom?.element = normalizedElement(text, aliasLabel: aliasLabel)

        case "x2":
            if let value = Double(text) { currentAtom?.position.x = value }
        case "y2":
            if let value = Double(text) { currentAtom?.position.y = value }
        case "x3":
            if let value = Double(text) { currentAtom?.position.x = value }
        case "y3":
            if let value = Double(text) { currentAtom?.position.y = value }
        case "z3":
            if let value = Double(text) { currentAtom?.zPosition = value }
        case "xy2":
            if let pair = parseCoordinatePair(text) {
                currentAtom?.position = CGPoint(x: pair.0, y: pair.1)
            }
        case "xyz3":
            if let triplet = parseCoordinateTriplet(text) {
                currentAtom?.position = CGPoint(x: triplet.0, y: triplet.1)
                currentAtom?.zPosition = triplet.2
            }
        case "formalcharge":
            if let charge = Int(text) { currentAtom?.charge = charge }
        case "hydrogencount":
            currentAtom?.hydrogenCount = Int(text)
        case "isotopenumber":
            currentAtom?.isotope = Int(text)
        case "cdk:aromaticatom", "aromaticatom":
            currentAtom?.aromatic = true
        default:
            break
        }
    }

    private func applyBondScalarCapture(_ capture: TextCapture) {
        let builtin = scalarBuiltin(capture.attributes)
        let dictRef = normalizedDictRef(capture.attributes["dictRef"])
        let text = capture.text.trimmingCharacters(in: .whitespacesAndNewlines)

        switch builtin ?? dictRef {
        case "order", "cdk:bondorder", "bondorder":
            if !text.isEmpty {
                currentBond?.order = mapBondOrder(text)
                currentBond?.hadExplicitOrder = true
            }
        case "cdk:aromaticbond", "aromaticbond":
            currentBond?.aromatic = true
            if currentBond?.hadExplicitOrder == false {
                currentBond?.order = .aromatic
            }
        default:
            if isAromaticDictRef(capture.attributes["dictRef"]) {
                currentBond?.aromatic = true
                if currentBond?.hadExplicitOrder == false {
                    currentBond?.order = .aromatic
                }
            }
        }
    }

    private func localName(_ elementName: String) -> String {
        String(elementName.split(separator: ":").last ?? Substring(elementName))
    }

    private func popCapture(for elementName: String) -> TextCapture? {
        guard let last = textCaptures.last, last.elementName == elementName else { return nil }
        return textCaptures.popLast()
    }

    private func coordinates(from attributes: [String: String]) -> ParsedCoordinates {
        if let xyz3 = parseCoordinateTriplet(attributes["xyz3"]) {
            return ParsedCoordinates(position: CGPoint(x: xyz3.0, y: xyz3.1), zPosition: xyz3.2)
        }
        if let xy2 = parseCoordinatePair(attributes["xy2"]) {
            return ParsedCoordinates(position: CGPoint(x: xy2.0, y: xy2.1), zPosition: nil)
        }

        let x2 = Double(attributes["x2"] ?? "")
        let y2 = Double(attributes["y2"] ?? "")
        let x3 = Double(attributes["x3"] ?? "")
        let y3 = Double(attributes["y3"] ?? "")
        let z3 = Double(attributes["z3"] ?? "")

        let x = x2 ?? x3 ?? 0
        let y = y2 ?? y3 ?? 0
        return ParsedCoordinates(position: CGPoint(x: x, y: y), zPosition: z3)
    }

    private func coordinates(at index: Int,
                             xy2: [(Double, Double)],
                             xyz3: [(Double, Double, Double)],
                             x2: [Double],
                             y2: [Double],
                             x3: [Double],
                             y3: [Double],
                             z3: [Double]) -> ParsedCoordinates {
        if index < xyz3.count {
            let triplet = xyz3[index]
            return ParsedCoordinates(position: CGPoint(x: triplet.0, y: triplet.1), zPosition: triplet.2)
        }
        if index < xy2.count {
            let pair = xy2[index]
            return ParsedCoordinates(position: CGPoint(x: pair.0, y: pair.1), zPosition: nil)
        }

        let x = index < x2.count ? x2[index] : (index < x3.count ? x3[index] : 0)
        let y = index < y2.count ? y2[index] : (index < y3.count ? y3[index] : 0)
        let z = index < z3.count ? z3[index] : nil
        return ParsedCoordinates(position: CGPoint(x: x, y: y), zPosition: z)
    }

    private func normalizedElement(_ raw: String, aliasLabel: String?) -> String {
        let trimmed = raw.trimmingCharacters(in: .whitespacesAndNewlines)
        let normalized = trimmed.uppercased()
        if normalized == "DU" || normalized == "DUMMY" || normalized == "R" || normalized == "R#" {
            return "R"
        }
        let canonical = CDKDescriptorSupport.canonicalElementSymbol(trimmed)
        if canonical.isEmpty {
            return aliasLabel == nil ? "C" : "R"
        }
        return canonical
    }

    private func isPseudoElement(_ raw: String) -> Bool {
        let normalized = raw.trimmingCharacters(in: .whitespacesAndNewlines).uppercased()
        return normalized == "DU" || normalized == "DUMMY" || normalized == "R" || normalized == "R#"
    }

    private func normalizedLabel(_ raw: String?) -> String? {
        guard let raw else { return nil }
        let trimmed = raw.trimmingCharacters(in: .whitespacesAndNewlines)
        return trimmed.isEmpty ? nil : trimmed
    }

    private func scalarBuiltin(_ attributes: [String: String]) -> String? {
        if let builtin = attributes["builtin"]?.trimmingCharacters(in: .whitespacesAndNewlines).lowercased(),
           !builtin.isEmpty {
            return builtin
        }
        if let dictRef = normalizedDictRef(attributes["dictRef"]),
           dictRef.hasPrefix("cml:") {
            return String(dictRef.dropFirst(4))
        }
        return normalizedDictRef(attributes["dictRef"])
    }

    private func normalizedDictRef(_ raw: String?) -> String? {
        guard let raw else { return nil }
        let trimmed = raw.trimmingCharacters(in: .whitespacesAndNewlines)
        return trimmed.isEmpty ? nil : trimmed.lowercased()
    }

    private func isAromaticDictRef(_ raw: String?) -> Bool {
        guard let normalized = normalizedDictRef(raw) else { return false }
        return normalized.hasSuffix("aromaticbond")
    }

    private func parseBondStereo(attributes: [String: String], text: String) -> BondStereo? {
        if let dictRef = attributes["dictRef"],
           let stereo = bondStereoTokenToValue(String(dictRef.split(separator: ":").last ?? "")) {
            return stereo
        }
        return bondStereoTokenToValue(text)
    }

    private func bondStereoTokenToValue(_ raw: String) -> BondStereo? {
        let token = raw.trimmingCharacters(in: .whitespacesAndNewlines).uppercased()
        switch token {
        case "W", "UP", "WEDGEBEGIN":
            return .up
        case "H", "DOWN", "HATCH":
            return .down
        default:
            return nil
        }
    }

    private func mapBondOrder(_ raw: String?) -> BondOrder {
        let token = (raw ?? "").trimmingCharacters(in: .whitespacesAndNewlines).uppercased()
        switch token {
        case "2", "D", "DOUBLE":
            return .double
        case "3", "T", "TRIPLE":
            return .triple
        case "A", "AR", "AROMATIC", "4":
            return .aromatic
        default:
            return .single
        }
    }

    private func splitTokens(_ raw: String?) -> [String] {
        guard let raw else { return [] }
        return raw.split(whereSeparator: \.isWhitespace).map(String.init)
    }

    private func parseInts(_ raw: String?) -> [Int] {
        splitTokens(raw).compactMap(Int.init)
    }

    private func parseBooleans(_ raw: String?) -> [Bool] {
        splitTokens(raw).map(isTruthy)
    }

    private func parseDoubles(_ raw: String?) -> [Double] {
        splitTokens(raw).compactMap(Double.init)
    }

    private func parseCoordinatePair(_ raw: String?) -> (Double, Double)? {
        guard let raw else { return nil }
        let values = splitTokens(raw).compactMap(Double.init)
        guard values.count >= 2 else { return nil }
        return (values[0], values[1])
    }

    private func parseCoordinateTriplet(_ raw: String?) -> (Double, Double, Double)? {
        guard let raw else { return nil }
        let values = splitTokens(raw).compactMap(Double.init)
        guard values.count >= 3 else { return nil }
        return (values[0], values[1], values[2])
    }

    private func parseCoordinatePairs(_ raw: String?) -> [(Double, Double)] {
        let values = parseDoubles(raw)
        guard values.count >= 2 else { return [] }
        var result: [(Double, Double)] = []
        var index = 0
        while index + 1 < values.count {
            result.append((values[index], values[index + 1]))
            index += 2
        }
        return result
    }

    private func parseCoordinateTriplets(_ raw: String?) -> [(Double, Double, Double)] {
        let values = parseDoubles(raw)
        guard values.count >= 3 else { return [] }
        var result: [(Double, Double, Double)] = []
        var index = 0
        while index + 2 < values.count {
            result.append((values[index], values[index + 1], values[index + 2]))
            index += 3
        }
        return result
    }

    private func isTruthy(_ raw: String?) -> Bool {
        guard let raw else { return false }
        switch raw.trimmingCharacters(in: .whitespacesAndNewlines).lowercased() {
        case "true", "1", "yes":
            return true
        default:
            return false
        }
    }
}
