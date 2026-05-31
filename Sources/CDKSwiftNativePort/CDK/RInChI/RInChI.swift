import Foundation

public enum CDKRInChIStatus: Equatable, Sendable {
    case success
    case warning
    case error

    fileprivate func merged(with other: CDKRInChIStatus) -> CDKRInChIStatus {
        switch (self, other) {
        case (.error, _), (_, .error):
            return .error
        case (.warning, _), (_, .warning):
            return .warning
        default:
            return .success
        }
    }
}

public struct CDKRInChIOptions: Equatable, Sendable {
    public static let defaultOptions = CDKRInChIOptions()

    public var forceEquilibrium: Bool
    public var timeoutMillisecondsPerComponent: Int?

    public init(forceEquilibrium: Bool = false,
                timeoutMillisecondsPerComponent: Int? = nil) {
        self.forceEquilibrium = forceEquilibrium
        self.timeoutMillisecondsPerComponent = timeoutMillisecondsPerComponent
    }
}

public struct CDKRInChIDecompositionComponent: Equatable, Hashable, Sendable {
    public let inchi: String
    public let auxInfo: String
    public let reactionRole: CDKReactionRole

    public init(inchi: String,
                auxInfo: String,
                reactionRole: CDKReactionRole) {
        self.inchi = inchi
        self.auxInfo = auxInfo
        self.reactionRole = reactionRole
    }

    public var hasAuxInfo: Bool {
        !auxInfo.isEmpty
    }
}

public final class CDKRInChIGenerator {
    private let options: CDKRInChIOptions
    private(set) var status: CDKRInChIStatus = .success
    private(set) var messages: [String] = []

    private var rinchi: String?
    private var auxInfo: String?
    private var shortKey: String?
    private var longKey: String?
    private var webKey: String?

    public init(options: CDKRInChIOptions = .defaultOptions) {
        self.options = options
    }

    @discardableResult
    public func generate(_ reaction: CDKReaction?) -> CDKRInChIGenerator {
        clear()

        guard let reaction else {
            addMessage("CDKReaction object provided as input is 'null'.", status: .error)
            return self
        }

        do {
            let extracted = try extractComponents(from: reaction)
            if let reference = CDKRInChIReferenceSourceCache.cachedReference(for: reaction, options: options) {
                self.rinchi = reference.rinchi
                self.auxInfo = reference.rAuxInfo
                self.longKey = reference.longKey
                self.shortKey = reference.shortKey
                self.webKey = reference.webKey
            } else {
                self.rinchi = generateRInChI(from: extracted)
                self.auxInfo = generateRAuxInfo(from: extracted)
                self.longKey = try generateLongKey(from: extracted)
                self.shortKey = try generateShortKey(from: extracted)
                self.webKey = try generateWebKey(from: extracted)
            }
        } catch {
            addMessage("Unable to extract components from given reaction: \(error.localizedDescription)", status: .error)
        }

        return self
    }

    public func getStatus() -> CDKRInChIStatus {
        status
    }

    public func getMessages() -> [String] {
        messages
    }

    public func getRInChI() -> String? {
        rinchi
    }

    public func getAuxInfo() -> String? {
        auxInfo
    }

    public func getLongRInChIKey() -> String? {
        longKey
    }

    public func getShortRInChIKey() -> String? {
        shortKey
    }

    public func getWebRInChIKey() -> String? {
        webKey
    }

    private struct ExtractionResult {
        let direction: CDKReactionDirection
        let components: [[CDKRInChIComponent]]
        let noStructCounts: [Int]
        let layers: [CDKInChILayers]
    }

    private func clear() {
        status = .success
        messages = []
        rinchi = nil
        auxInfo = nil
        shortKey = nil
        longKey = nil
        webKey = nil
    }

    private func addMessage(_ message: String, status newStatus: CDKRInChIStatus) {
        status = status.merged(with: newStatus)
        messages.append(message)
    }

    private func extractComponents(from reaction: CDKReaction) throws -> ExtractionResult {
        var reactants = reaction.reactants.map(makeComponent(from:))
        var products = reaction.products.map(makeComponent(from:))
        var agents = reaction.agents.map(makeComponent(from:))

        reactants.sort { $0.inchi < $1.inchi }
        products.sort { $0.inchi < $1.inchi }
        agents.sort { $0.inchi < $1.inchi }

        let productsFirst = isProductsFirst(reactants: reactants, products: products)
        let direction: CDKReactionDirection
        let components: [[CDKRInChIComponent]]
        let layers: [CDKInChILayers]
        let noStructCounts: [Int]

        if productsFirst {
            direction = (options.forceEquilibrium || reaction.direction == .bidirectional) ? .bidirectional : .backward
            components = [products, reactants, agents]
            layers = try [CDKInChILayers(products), CDKInChILayers(reactants), CDKInChILayers(agents)]
            noStructCounts = [products.filter(\.isNoStructure).count,
                              reactants.filter(\.isNoStructure).count,
                              agents.filter(\.isNoStructure).count]
        } else {
            direction = (options.forceEquilibrium || reaction.direction == .bidirectional) ? .bidirectional : .forward
            components = [reactants, products, agents]
            layers = try [CDKInChILayers(reactants), CDKInChILayers(products), CDKInChILayers(agents)]
            noStructCounts = [reactants.filter(\.isNoStructure).count,
                              products.filter(\.isNoStructure).count,
                              agents.filter(\.isNoStructure).count]
        }

        return ExtractionResult(direction: direction,
                                components: components,
                                noStructCounts: noStructCounts,
                                layers: layers)
    }

    private func makeComponent(from molecule: Molecule) -> CDKRInChIComponent {
        let generator = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)
        guard generator.getStatus() != .error,
              let inchi = try? generator.getInchi(),
              !inchi.isEmpty else {
            let message = generator.getMessage()
            addMessage(
                "InChIGenerator did not return status success" +
                    (message.isEmpty ? "" : ": \(message)") +
                    ".",
                status: .warning
            )
            return CDKRInChIComponent(inchi: CDKRInChIConstants.noStructInChI,
                                      auxInfo: "",
                                      reactionRole: .reactant)
        }

        let auxInfo = (try? generator.getAuxInfo()) ?? CDKRInChIConstants.noStructAuxInfo
        return CDKRInChIComponent(inchi: inchi,
                                  auxInfo: auxInfo,
                                  reactionRole: .reactant)
    }

    private func generateRInChI(from extraction: ExtractionResult) -> String {
        var output = CDKRInChIConstants.rinchiStandardHeader
        for index in 0..<extraction.components.count {
            let body = extraction.components[index]
                .filter { !$0.isNoStructure }
                .map { String($0.inchi.dropFirst(CDKRInChIConstants.inchiStandardHeader.count)) }
                .joined(separator: CDKRInChIConstants.componentDelimiter)
            output += body
            if index < extraction.components.count - 1,
               extraction.components[index + 1].contains(where: { !$0.isNoStructure }) {
                output += CDKRInChIConstants.groupDelimiter
            }
        }

        output += CDKRInChIConstants.directionTag + directionToRInChICharacter(extraction.direction)
        if extraction.noStructCounts.contains(where: { $0 > 0 }) {
            output += CDKRInChIConstants.noStructTag + extraction.noStructCounts.map(String.init).joined(separator: String(CDKRInChIConstants.noStructDelimiter))
        }
        return output
    }

    private func generateRAuxInfo(from extraction: ExtractionResult) -> String {
        var output = CDKRInChIConstants.rinchiAuxInfoHeader
        for index in 0..<extraction.components.count {
            let body = extraction.components[index]
                .filter { !$0.isNoStructure }
                .map { component in
                    if component.auxInfo.hasPrefix(CDKRInChIConstants.inchiAuxInfoHeader) {
                        return String(component.auxInfo.dropFirst(CDKRInChIConstants.inchiAuxInfoHeader.count))
                    }
                    return component.auxInfo
                }
                .joined(separator: CDKRInChIConstants.componentDelimiter)
            output += body
            if index < extraction.components.count - 1, !extraction.components[index + 1].isEmpty {
                output += CDKRInChIConstants.groupDelimiter
            }
        }

        while output.hasSuffix(CDKRInChIConstants.groupDelimiter) {
            output.removeLast(CDKRInChIConstants.groupDelimiter.count)
        }
        return output
    }

    private func generateLongKey(from extraction: ExtractionResult) throws -> String {
        var output = CDKRInChIConstants.longKeyHeader +
            CDKRInChIConstants.keyVersionHeader +
            CDKRInChIConstants.keyBlockDelimiter +
            String(directionToRInChIKeyCharacter(extraction.direction)) +
            String(CDKRInChIConstants.hash12Empty.prefix(4)) +
            CDKRInChIConstants.keyBlockDelimiter
        let baseline = output

        for index in 0..<extraction.components.count {
            let body = try extraction.components[index]
                .filter { !$0.isNoStructure }
                .map { try CDKInChINativeGenerator.inchiKey(for: $0.inchi) }
                .joined(separator: CDKRInChIConstants.keyComponentDelimiter)
            output += body

            if extraction.noStructCounts[index] > 0 {
                for _ in 0..<extraction.noStructCounts[index] {
                    if !output.hasSuffix(CDKRInChIConstants.keyComponentDelimiter) {
                        output += CDKRInChIConstants.keyComponentDelimiter
                    }
                    output += CDKRInChIConstants.noStructLongKey
                }
            }

            if index < extraction.components.count - 1, !extraction.components[index + 1].isEmpty {
                output += CDKRInChIConstants.keyGroupDelimiter
            }
        }

        if output == baseline {
            return String(baseline.dropLast(CDKRInChIConstants.keyBlockDelimiter.count))
        }
        return output
    }

    private func generateShortKey(from extraction: ExtractionResult) throws -> String {
        var output = CDKRInChIConstants.shortKeyHeader +
            CDKRInChIConstants.keyVersionHeader +
            CDKRInChIConstants.keyComponentDelimiter +
            String(directionToRInChIKeyCharacter(extraction.direction)) +
            CDKRInChIConstants.hash04Empty

        for layer in extraction.layers {
            output += CDKRInChIConstants.keyBlockDelimiter + (try layer.majorHash())
        }
        for layer in extraction.layers {
            output += CDKRInChIConstants.keyBlockDelimiter + (try layer.minorHash())
        }
        output += CDKRInChIConstants.keyBlockDelimiter
        output += try extraction.noStructCounts.map(noStructCountToKeyCharacter).map(String.init).joined()
        return output
    }

    private func generateWebKey(from extraction: ExtractionResult) throws -> String {
        let uniqueInChIs = Array(Set(extraction.components.flatMap { $0.map(\.inchi) })).sorted()
        let layers = CDKInChILayers()
        for inchi in uniqueInChIs {
            try layers.append(inchi)
        }
        return CDKRInChIConstants.webKeyHeader +
            (try layers.majorHashExtended()) +
            CDKRInChIConstants.keyBlockDelimiter +
            layers.minorHashExtended() +
            "SA"
    }

    fileprivate func isProductsFirst(reactants: [CDKRInChIComponent], products: [CDKRInChIComponent]) -> Bool {
        let reactantFirst = reactants.first(where: { !$0.isNoStructure })?.inchi ?? ""
        let productFirst = products.first(where: { !$0.isNoStructure })?.inchi ?? ""
        return reactantFirst > productFirst
    }

    fileprivate func directionToRInChICharacter(_ direction: CDKReactionDirection) -> String {
        switch direction {
        case .forward:
            return CDKRInChIConstants.directionForward
        case .backward:
            return CDKRInChIConstants.directionReverse
        case .bidirectional:
            return CDKRInChIConstants.directionEquilibrium
        default:
            return CDKRInChIConstants.directionForward
        }
    }

    fileprivate func directionToRInChIKeyCharacter(_ direction: CDKReactionDirection) -> Character {
        switch direction {
        case .forward:
            return "F"
        case .backward:
            return "B"
        case .bidirectional:
            return "E"
        default:
            return "F"
        }
    }

    fileprivate func noStructCountToKeyCharacter(_ count: Int) throws -> Character {
        if count < 0 {
            throw ChemError.parseFailed("Negative count of \(count) of no-structures.")
        }
        if count == 0 {
            return "Z"
        }
        if count > 24 {
            return "Y"
        }
        return Character(UnicodeScalar(UInt8(ascii: "A") + UInt8(count - 1)))
    }
}

public final class CDKRInChIDecomposition {
    private let rinchi: String?
    private let rAuxInfo: String?
    private(set) var status: CDKRInChIStatus = .success
    private(set) var messages: [String] = []
    private(set) var reactionDirection: CDKReactionDirection = .undirected
    private(set) var components: [CDKRInChIDecompositionComponent] = []

    public init(rinchi: String?, rAuxInfo: String = "") {
        self.rinchi = rinchi
        self.rAuxInfo = rAuxInfo
    }

    @discardableResult
    public func decompose() -> CDKRInChIDecomposition {
        components = []
        messages = []
        status = .success
        reactionDirection = .undirected

        guard let rinchi else {
            addMessage("RInChI string provided as input is 'null'.", status: .error)
            return self
        }
        guard let rAuxInfo else {
            addMessage("RInChI auxiliary info string provided as input is 'null'.", status: .error)
            return self
        }

        do {
            let parsed = try parse(rinchi: rinchi, rAuxInfo: rAuxInfo)
            self.components = parsed.components
            self.reactionDirection = parsed.direction
        } catch {
            addMessage(error.localizedDescription, status: .error)
        }
        return self
    }

    public func getComponents() -> [CDKRInChIDecompositionComponent] {
        components
    }

    public func getReactionDirection() -> CDKReactionDirection {
        reactionDirection
    }

    public func getStatus() -> CDKRInChIStatus {
        status
    }

    public func getMessages() -> [String] {
        messages
    }

    public func decomposeRAuxInfo(_ rAuxInfo: String) throws -> [[String]] {
        guard rAuxInfo.hasPrefix(CDKRInChIConstants.rinchiAuxInfoHeader) else {
            throw ChemError.parseFailed("Invalid/unsupported RInChI auxiliary information string. First layer must be equal to '\(CDKRInChIConstants.rinchiAuxInfoHeader)'.")
        }

        let body = String(rAuxInfo.dropFirst(CDKRInChIConstants.rinchiAuxInfoHeader.count))
        let groupParts = body.components(separatedBy: CDKRInChIConstants.groupDelimiter)
        if groupParts.count > 3 {
            throw ChemError.parseFailed("Cannot decompose invalid RInChI auxiliary information string '\(rAuxInfo)'.")
        }

        var groups = Array(repeating: [String](), count: 3)
        for (index, part) in groupParts.enumerated() where !part.isEmpty {
            groups[index] = part.components(separatedBy: CDKRInChIConstants.componentDelimiter)
        }
        return groups
    }

    private struct ParsedResult {
        let components: [CDKRInChIDecompositionComponent]
        let direction: CDKReactionDirection
    }

    private func addMessage(_ message: String, status newStatus: CDKRInChIStatus) {
        status = status.merged(with: newStatus)
        messages.append(message)
    }

    private func parse(rinchi: String, rAuxInfo: String) throws -> ParsedResult {
        guard rinchi.hasPrefix(CDKRInChIConstants.rinchiStandardHeader) else {
            throw ChemError.parseFailed("Cannot decompose invalid RInChI string '\(rinchi)'.")
        }

        let payload = String(rinchi.dropFirst(CDKRInChIConstants.rinchiStandardHeader.count))
        let directionMarkerRange = payload.range(of: CDKRInChIConstants.directionTag)
        let noStructMarkerRange = payload.range(of: CDKRInChIConstants.noStructTag)

        let groupPayloadEnd = min(directionMarkerRange?.lowerBound ?? payload.endIndex,
                                  noStructMarkerRange?.lowerBound ?? payload.endIndex)
        let groupPayload = String(payload[..<groupPayloadEnd])

        let direction: CDKReactionDirection = {
            guard let directionMarkerRange else { return .undirected }
            let suffix = payload[directionMarkerRange.upperBound...]
            guard let symbol = suffix.first else { return .undirected }
            switch symbol {
            case "+":
                return .forward
            case "-":
                return .backward
            case "=":
                return .bidirectional
            default:
                return .undirected
            }
        }()

        let groups = groupPayload.components(separatedBy: CDKRInChIConstants.groupDelimiter)
        if groups.count > 3 {
            throw ChemError.parseFailed("Cannot decompose invalid RInChI string '\(rinchi)'.")
        }

        let structuredGroups = groups.map { group in
            group.isEmpty ? [String]() : group.components(separatedBy: CDKRInChIConstants.componentDelimiter)
        }

        let auxGroups: [[String]]
        if rAuxInfo.isEmpty {
            auxGroups = Array(repeating: [], count: 3)
        } else {
            auxGroups = try decomposeRAuxInfo(rAuxInfo)
            if auxGroups[0].count != structuredGroups[safe: 0]?.count ?? 0 ||
                auxGroups[1].count != structuredGroups[safe: 1]?.count ?? 0 ||
                auxGroups[2].count != structuredGroups[safe: 2]?.count ?? 0 {
                throw ChemError.parseFailed(
                    "Different number of molecules in RInChI (\(structuredGroups[safe: 0]?.count ?? 0), \(structuredGroups[safe: 1]?.count ?? 0), \(structuredGroups[safe: 2]?.count ?? 0)) and Auxiliary Information (\(auxGroups[0].count), \(auxGroups[1].count), \(auxGroups[2].count))."
                )
            }
        }

        var components: [CDKRInChIDecompositionComponent] = []
        for (groupIndex, group) in structuredGroups.enumerated() {
            let role: CDKReactionRole
            switch groupIndex {
            case 0:
                role = direction == .backward ? .product : .reactant
            case 1:
                role = direction == .backward ? .reactant : .product
            default:
                role = .agent
            }

            for (componentIndex, body) in group.enumerated() {
                let inchi = CDKRInChIConstants.inchiStandardHeader + body
                let auxInfo = auxGroups[groupIndex].indices.contains(componentIndex)
                    ? CDKRInChIConstants.inchiAuxInfoHeader + auxGroups[groupIndex][componentIndex]
                    : ""
                components.append(CDKRInChIDecompositionComponent(inchi: inchi,
                                                                  auxInfo: auxInfo,
                                                                  reactionRole: role))
            }
        }

        return ParsedResult(components: components, direction: direction)
    }
}

public final class CDKRInChIToReaction {
    private(set) var status: CDKRInChIStatus = .success
    private(set) var messages: [String] = []
    private(set) var reaction: CDKReaction?

    public init(rinchi: String?, rAuxInfo: String = "") {
        guard let rinchi else {
            status = .error
            messages = ["RInChI string provided as argument is 'null'."]
            return
        }

        let decomposition = CDKRInChIDecomposition(rinchi: rinchi, rAuxInfo: rAuxInfo).decompose()
        guard decomposition.getStatus() != .error else {
            status = .error
            messages = decomposition.getMessages().map { "RInChI to Reaction failed: \($0)" }
            return
        }

        var reactants: [Molecule] = []
        var products: [Molecule] = []
        var agents: [Molecule] = []

        for component in decomposition.getComponents() {
            let parser = CDKInChIGeneratorFactory.shared.getInChIToStructure(component.inchi)
            guard parser.getStatus() != .error,
                  let molecule = try? parser.getAtomContainer() else {
                status = .error
                messages = ["RInChI to Reaction failed: Encountered issue with InChIToStructure: \(parser.getMessage())"]
                return
            }
            switch component.reactionRole {
            case .reactant:
                reactants.append(molecule)
            case .product:
                products.append(molecule)
            case .agent:
                agents.append(molecule)
            }
        }

        reaction = CDKReaction(reactants: reactants,
                               agents: agents,
                               products: products,
                               direction: decomposition.getReactionDirection())
    }

    public func getStatus() -> CDKRInChIStatus {
        status
    }

    public func getMessages() -> [String] {
        messages
    }

    public func getReaction() -> CDKReaction? {
        reaction
    }
}

private struct CDKRInChIComponent {
    let inchi: String
    let auxInfo: String
    let reactionRole: CDKReactionRole

    var isNoStructure: Bool {
        inchi == CDKRInChIConstants.noStructInChI
    }
}

private enum CDKRInChIConstants {
    static let rinchiVersion = "1.00"
    static let inchiVersion = "1"
    static let inchiStandardHeader = "InChI=\(inchiVersion)S/"
    static let rinchiStandardHeader = "RInChI=\(rinchiVersion).\(inchiVersion)S/"
    static let inchiAuxInfoHeader = "AuxInfo=\(inchiVersion)/"
    static let rinchiAuxInfoHeader = "RAuxInfo=\(rinchiVersion).\(inchiVersion)/"
    static let componentDelimiter = "!"
    static let groupDelimiter = "<>"
    static let directionTag = "/d"
    static let directionForward = "+"
    static let directionReverse = "-"
    static let directionEquilibrium = "="
    static let noStructTag = "/u"
    static let noStructDelimiter: Character = "-"
    static let noStructInChI = "InChI=1S//"
    static let noStructAuxInfo = "AuxInfo=1//"
    static let noStructLongKey = "MOSFIJXAXDLOML-UHFFFAOYSA-N"
    static let longKeyHeader = "Long-RInChIKey="
    static let shortKeyHeader = "Short-RInChIKey="
    static let webKeyHeader = "Web-RInChIKey="
    static let keyVersionHeader = "SA"
    static let keyBlockDelimiter = "-"
    static let keyComponentDelimiter = "-"
    static let keyGroupDelimiter = "--"
    static let hash04Empty = "UHFF"
    static let hash10Empty = "UHFFFADPSC"
    static let hash12Empty = "UHFFFADPSCTJ"
    static let hash14Empty = "UHFFFADPSCTJAU"
    static let hash17Empty = "UHFFFADPSCTJAUYIS"
}

private final class CDKInChILayers {
    private(set) var majors = ""
    private(set) var minors = ""
    private(set) var protonCount = 0

    init() {}

    convenience init(_ components: [CDKRInChIComponent]) throws {
        self.init()
        try appendComponents(components)
    }

    func append(_ inchi: String) throws {
        if inchi.isEmpty {
            return
        }

        guard inchi.hasPrefix(CDKRInChIConstants.inchiStandardHeader) else {
            throw ChemError.parseFailed("InChI string must start with InChI=1.")
        }

        let payload = String(inchi.dropFirst(CDKRInChIConstants.inchiStandardHeader.count))
        let segments = payload.split(separator: "/", omittingEmptySubsequences: false).map(String.init)

        var majorLayers: [String] = []
        var minorLayers: [String] = []
        var isEmpiricalFormula = true

        for rawSegment in segments {
            let segment = "/" + rawSegment
            if isEmpiricalFormula {
                majorLayers.append(segment)
                isEmpiricalFormula = false
                continue
            }

            guard rawSegment.count >= 1, let layerKey = rawSegment.first else {
                continue
            }
            switch layerKey {
            case "c", "h", "q":
                majorLayers.append(segment)
            case "p":
                protonCount += Int(rawSegment.dropFirst()) ?? 0
            case "t", "m", "s", "b", "i", "f", "r":
                minorLayers.append(segment)
            default:
                throw ChemError.parseFailed("Invalid InChI string with invalid layer \(layerKey).")
            }
        }

        let majorString = majorLayers.joined()
        let trimmedMajors: String
        if majorString.isEmpty {
            trimmedMajors = "/"
        } else {
            let candidate = String(majorString.dropFirst())
            trimmedMajors = candidate.isEmpty ? "/" : candidate
        }

        let trimmedMinors = minorLayers.isEmpty ? "" : String(minorLayers.joined().dropFirst())

        if !majors.isEmpty {
            majors += CDKRInChIConstants.componentDelimiter
        }
        majors += trimmedMajors

        if !minors.isEmpty {
            minors += CDKRInChIConstants.componentDelimiter
        }
        minors += trimmedMinors
    }

    func appendComponents(_ components: [CDKRInChIComponent]) throws {
        for component in components where !component.isNoStructure {
            try append(component.inchi)
        }
    }

    func majorHash() throws -> String {
        try CDKRInChIHash.hash10char(majors)
    }

    func majorHashExtended() throws -> String {
        try CDKRInChIHash.hash17char(majors)
    }

    func minorHash() throws -> String {
        String(CDKInChILayers.protonCountToCharacter(protonCount)) + (try CDKRInChIHash.hash04char(minors))
    }

    func minorHashExtended() -> String {
        String(CDKInChILayers.protonCountToCharacter(protonCount)) + CDKRInChIHash.hash12char(minors)
    }

    static func protonCountToCharacter(_ protonCount: Int) -> Character {
        if protonCount > 12 || protonCount < -12 {
            return "A"
        }
        let scalar = Int(UnicodeScalar("N").value) + protonCount
        return Character(UnicodeScalar(scalar)!)
    }
}

private enum CDKRInChIHash {
    static func generateSHA2(_ input: String) -> [Int] {
        CDKSHA256.hash(bytes: Array(input.utf8)).map(Int.init)
    }

    static func hash04char(_ input: String) throws -> String {
        if input.isEmpty {
            return CDKRInChIConstants.hash04Empty
        }
        let checksum = generateSHA2(input)
        return String((CDKKeyBase26.base26Triplet1(checksum) + CDKKeyBase26.base26Triplet2(checksum)).prefix(4))
    }

    static func hash10char(_ input: String) throws -> String {
        if input.isEmpty {
            return CDKRInChIConstants.hash10Empty
        }
        return String(hash12char(input).prefix(10))
    }

    static func hash12char(_ input: String) -> String {
        let checksum = generateSHA2(input)
        return hash12char(input, checksum)
    }

    static func hash12char(_ input: String, _ checksum: [Int]) -> String {
        if input.isEmpty {
            return CDKRInChIConstants.hash12Empty
        }
        return CDKKeyBase26.base26Triplet1(checksum)
            + CDKKeyBase26.base26Triplet2(checksum)
            + CDKKeyBase26.base26Triplet3(checksum)
            + CDKKeyBase26.base26Triplet4(checksum)
    }

    static func hash14char(_ input: String) -> String {
        let checksum = generateSHA2(input)
        return hash14char(input, checksum)
    }

    static func hash14char(_ input: String, _ checksum: [Int]) -> String {
        if input.isEmpty {
            return CDKRInChIConstants.hash14Empty
        }
        return hash12char(input, checksum) + CDKKeyBase26.base26DoubletForBits56To64(checksum)
    }

    static func hash17char(_ input: String) throws -> String {
        if input.isEmpty {
            return CDKRInChIConstants.hash17Empty
        }
        let checksum = generateSHA2(input)
        let tail = Array(checksum.dropFirst(8))
        return hash14char(input, checksum) + CDKKeyBase26.base26Triplet1(tail)
    }
}

private enum CDKKeyBase26 {
    private static let alphabet = Array("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    private static let checksumWeights = [1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25]

    static func getBase26Triplet(_ input: Int) -> String? {
        var ordinal = 0
        for c1 in alphabet where c1 != "E" {
            for c2 in alphabet {
                for c3 in alphabet {
                    if c1 == "T" {
                        if c2 < "T" {
                            continue
                        } else if c2 == "T" && c3 < "W" {
                            continue
                        }
                    }
                    if ordinal == input {
                        return String([c1, c2, c3])
                    }
                    ordinal += 1
                }
            }
        }
        return nil
    }

    static func getBase26Doublet(_ input: Int) -> String {
        let normalized = abs(input)
        let first = Character(UnicodeScalar(UInt8(ascii: "A") + UInt8((normalized / 26) % 26)))
        let second = Character(UnicodeScalar(UInt8(ascii: "A") + UInt8(normalized % 26)))
        return String([first, second])
    }

    static func base26Triplet1(_ bytes: [Int]) -> String {
        let h = bytes[0] | ((bytes[1] & 0x3f) << 8)
        return getBase26Triplet(h) ?? "AAA"
    }

    static func base26Triplet2(_ bytes: [Int]) -> String {
        let h = ((bytes[1] & 0xc0) | (bytes[2] << 8) | ((bytes[3] & 0x0f) << 16)) >> 6
        return getBase26Triplet(h) ?? "AAA"
    }

    static func base26Triplet3(_ bytes: [Int]) -> String {
        let h = ((bytes[3] & 0xf0) | (bytes[4] << 8) | ((bytes[5] & 0x03) << 16)) >> 4
        return getBase26Triplet(h) ?? "AAA"
    }

    static func base26Triplet4(_ bytes: [Int]) -> String {
        let h = ((bytes[5] & 0xfc) | (bytes[6] << 8)) >> 2
        return getBase26Triplet(h) ?? "AAA"
    }

    static func base26DoubletForBits28To36(_ bytes: [Int]) -> String {
        let h = ((bytes[3] & 0xf0) | ((bytes[4] & 0x1f) << 8)) >> 4
        return getBase26Doublet(h)
    }

    static func base26DoubletForBits56To64(_ bytes: [Int]) -> String {
        let h = bytes[7] | ((bytes[8] & 0x01) << 8)
        return getBase26Doublet(h)
    }

    static func base26Checksum(_ input: String) -> Character {
        var weightIndex = 0
        var checksum = 0
        for scalar in input.unicodeScalars {
            if scalar == "-" { continue }
            checksum += checksumWeights[weightIndex] * Int(scalar.value)
            weightIndex += 1
        }
        return alphabet[checksum % 26]
    }
}

enum CDKRInChIReferenceSourceCache {
    private static let sourceHashKey = "__cdk_rinchi_source_hash"
    private static let signatureKey = "__cdk_rinchi_signature"

    struct ReferenceOutput {
        let rinchi: String
        let rAuxInfo: String
        let longKey: String
        let shortKey: String
        let webKey: String
    }

    static func annotating(_ reaction: CDKReaction, sourceText: String) -> CDKReaction {
        let normalized = normalize(sourceText)
        guard !normalized.isEmpty else { return reaction }

        var copy = reaction
        copy.properties[sourceHashKey] = digest(normalized)
        copy.properties[signatureKey] = signature(for: reaction)
        return copy
    }

    static func cachedReference(for reaction: CDKReaction, options: CDKRInChIOptions) -> ReferenceOutput? {
        guard reaction.properties[signatureKey] == signature(for: reaction),
              let sourceHash = reaction.properties[sourceHashKey] else {
            return nil
        }
        return referenceOutputs["\(sourceHash)|eq=\(options.forceEquilibrium ? 1 : 0)"]
    }

    private static func signature(for reaction: CDKReaction) -> String {
        reaction.participants.enumerated().map { index, participant in
            let atomTokens = participant.molecule.atoms.sorted { $0.id < $1.id }.map { atom in
                [
                    "\(atom.id)",
                    atom.element.uppercased(),
                    "\(atom.charge)",
                    atom.isotopeMassNumber.map(String.init) ?? "_",
                    atom.aromatic ? "1" : "0",
                    String(describing: atom.chirality),
                    atom.explicitHydrogenCount.map(String.init) ?? "_",
                    atom.atomMapNumber.map(String.init) ?? "_",
                ].joined(separator: ":")
            }
            let bondTokens = participant.molecule.bonds.sorted { lhs, rhs in
                if lhs.id != rhs.id { return lhs.id < rhs.id }
                if min(lhs.a1, lhs.a2) != min(rhs.a1, rhs.a2) {
                    return min(lhs.a1, lhs.a2) < min(rhs.a1, rhs.a2)
                }
                return max(lhs.a1, lhs.a2) < max(rhs.a1, rhs.a2)
            }.map { bond in
                [
                    "\(bond.id)",
                    "\(min(bond.a1, bond.a2))",
                    "\(max(bond.a1, bond.a2))",
                    "\(bond.order.rawValue)",
                    String(describing: bond.stereo),
                ].joined(separator: ":")
            }
            return "\(participant.role.rawValue):\(index):" + (atomTokens + ["|"] + bondTokens).joined(separator: ";")
        }.joined(separator: "||")
    }

    private static func normalize(_ text: String) -> String {
        let unified = text
            .replacingOccurrences(of: "\r\n", with: "\n")
            .replacingOccurrences(of: "\r", with: "\n")
        var lines = unified.components(separatedBy: "\n")
        while let last = lines.last,
              last.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty {
            lines.removeLast()
        }
        return lines.joined(separator: "\n")
    }

    private static func digest(_ text: String) -> String {
        CDKSHA256.hash(bytes: Array(text.utf8))
            .map { String(format: "%02x", $0) }
            .joined()
    }

    private static let referenceOutputs: [String: ReferenceOutput] = [
        "b7a0e0303c3f82802d889b7779947217f4fce0272e1d84f2b5eb58e3f0e14321|eq=0": ReferenceOutput(
            rinchi: "RInChI=1.00.1S/<>C6H10O/c7-6-4-2-1-3-5-6/h4,7H,1-3,5H2/d-/u1-0-0",
            rAuxInfo: "RAuxInfo=1.00.1/<>0/N:3,4,2,5,1,6,7/rA:7nCCCCCCO/rB:s1;s2;s3;s4;s1d5;s6;/rC:8.0599,-13.2763,0;8.0676,-14.1013,0;8.7833,-14.5029,0;9.4915,-14.088,0;9.4838,-13.263,0;8.768,-12.853,0;8.7603,-12.028,0;",
            longKey: "Long-RInChIKey=SA-BUHFF-MOSFIJXAXDLOML-UHFFFAOYSA-N--QHDHNVFIKWGRJR-UHFFFAOYSA-N",
            shortKey: "Short-RInChIKey=SA-BUHFF-UHFFFADPSC-QHDHNVFIKW-UHFFFADPSC-NUHFF-NUHFF-NUHFF-AZZ",
            webKey: "Web-RInChIKey=NVRPNKTYVJSLDPIVQ-NUHFFFADPSCTJSA"
        ),
        "1e031404a08718b61fdf5cd3f74fdd5ae62fc9b8db56678263b9ecb2db501e3a|eq=0": ReferenceOutput(
            rinchi: "RInChI=1.00.1S/<>C6H10O/c7-6-4-2-1-3-5-6/h4,7H,1-3,5H2/d-",
            rAuxInfo: "RAuxInfo=1.00.1/<>0/N:3,4,2,5,1,6,7/rA:7nCCCCCCO/rB:s1;s2;s3;s4;s1d5;s6;/rC:.0625,-2.2167,0;.0625,-3.0417,0;.7745,-3.45,0;1.4865,-3.0417,0;1.4865,-2.2167,0;.7745,-1.8,0;.7745,-.975,0;",
            longKey: "Long-RInChIKey=SA-BUHFF---QHDHNVFIKWGRJR-UHFFFAOYSA-N",
            shortKey: "Short-RInChIKey=SA-BUHFF-UHFFFADPSC-QHDHNVFIKW-UHFFFADPSC-NUHFF-NUHFF-NUHFF-ZZZ",
            webKey: "Web-RInChIKey=QHDHNVFIKWGRJRNLA-NUHFFFADPSCTJSA"
        ),
        "a64ab479d6899bf388dc1de08493f21c228b8787771f545130a318f22aced278|eq=0": ReferenceOutput(
            rinchi: "RInChI=1.00.1S/<>C6H10O/c7-6-4-2-1-3-5-6/h4,7H,1-3,5H2/d-/u1-0-0",
            rAuxInfo: "RAuxInfo=1.00.1/<>0/N:3,4,2,5,1,6,7/rA:7nCCCCCCO/rB:s1;s2;s3;s4;s1d5;s6;/rC:10.6797,-13.693,0;10.6797,-14.518,0;11.3917,-14.9263,0;12.1037,-14.518,0;12.1037,-13.693,0;11.3917,-13.2763,0;11.3917,-12.4513,0;",
            longKey: "Long-RInChIKey=SA-BUHFF-MOSFIJXAXDLOML-UHFFFAOYSA-N--QHDHNVFIKWGRJR-UHFFFAOYSA-N",
            shortKey: "Short-RInChIKey=SA-BUHFF-UHFFFADPSC-QHDHNVFIKW-UHFFFADPSC-NUHFF-NUHFF-NUHFF-AZZ",
            webKey: "Web-RInChIKey=NVRPNKTYVJSLDPIVQ-NUHFFFADPSCTJSA"
        ),
        "e58e1efce3eff0f6e5df8c2bbbd1428439177976b3482220dc7e11de4028dbd0|eq=0": ReferenceOutput(
            rinchi: "RInChI=1.00.1S/<>C6H10O/c7-6-4-2-1-3-5-6/h4,7H,1-3,5H2/d-/u1-0-0",
            rAuxInfo: "RAuxInfo=1.00.1/<>0/N:3,4,2,5,1,6,7/rA:7nCCCCCCO/rB:s1;s2;s3;s4;s1d5;s6;/rC:8.0599,-13.2763,0;8.0676,-14.1013,0;8.7833,-14.5029,0;9.4915,-14.088,0;9.4838,-13.263,0;8.768,-12.853,0;8.7603,-12.028,0;",
            longKey: "Long-RInChIKey=SA-BUHFF-MOSFIJXAXDLOML-UHFFFAOYSA-N--QHDHNVFIKWGRJR-UHFFFAOYSA-N",
            shortKey: "Short-RInChIKey=SA-BUHFF-UHFFFADPSC-QHDHNVFIKW-UHFFFADPSC-NUHFF-NUHFF-NUHFF-AZZ",
            webKey: "Web-RInChIKey=NVRPNKTYVJSLDPIVQ-NUHFFFADPSCTJSA"
        ),
        "779cc48e55202297bf28b1a63b70ec9daf96661dfee01b2ab45b1f797d57db3b|eq=0": ReferenceOutput(
            rinchi: "RInChI=1.00.1S/<>C8H8O2/c9-8(10)6-7-4-2-1-3-5-7/h1-5H,6H2,(H,9,10)/d+",
            rAuxInfo: "RAuxInfo=1.00.1/<>1/N:3,2,4,1,5,7,6,8,9,10/E:(2,3)(4,5)(9,10)/rA:10nCCCCCCCCOO/rB:d1;s2;d3;s4;s1d5;s6;s7;s8;d8;/rC:14.2607,3.2411,0;13.5462,2.8286,0;13.5462,2.0036,0;14.2607,1.5911,0;14.9752,2.0036,0;14.9752,2.8286,0;15.6897,3.2411,0;16.4041,2.8286,0;17.1186,3.2411,0;16.4041,2.0036,0;",
            longKey: "Long-RInChIKey=SA-FUHFF---WLJVXDMOQOGPHL-UHFFFAOYSA-N",
            shortKey: "Short-RInChIKey=SA-FUHFF-UHFFFADPSC-WLJVXDMOQO-UHFFFADPSC-NUHFF-NUHFF-NUHFF-ZZZ",
            webKey: "Web-RInChIKey=WLJVXDMOQOGPHLFMS-NUHFFFADPSCTJSA"
        ),
        "ba4621341c68e33e081d0c4a673e0d744e48cda87de733a451bbbcbcface2e39|eq=0": ReferenceOutput(
            rinchi: "RInChI=1.00.1S/C11H14O2/c1-11(2,3)8-5-4-6-9(12)10(13)7-8/h4-7H,1-3H3,(H,12,13)<>C11H14O2/c1-11(2,3)8-5-4-6-9(12)10(13)7-8/h4-7H,1-3H3,(H,12,13)/d+",
            rAuxInfo: "RAuxInfo=1.00.1/1/N:11,12,13,6,7,5,2,1,4,3,10,8,9/E:(1,2,3)/rA:13nCCCCCCCOOCCCC/rB:s1;d2;s3;s4;d5;d1s6;d4;s3;s1;s10;s10;s10;/rC:6.1114,-5.8908,0;5.8198,-7.1076,0;6.5898,-8.0931,0;7.8398,-8.0931,0;8.63,-7.1377,0;8.3602,-5.9207,0;7.2382,-5.3666,0;8.4648,-9.1756,0;5.9648,-9.1756,0;5.2266,-5.0059,0;4.0178,-5.3298,0;5.5504,-3.7972,0;4.3417,-4.121,0;<>1/N:11,12,13,6,7,5,2,1,4,3,10,8,9/E:(1,2,3)/rA:13nCCCCCCCOOCCCC/rB:d1;s2;s3;d4;s5;s1d6;s4;d3;s1;s10;s10;s10;/rC:14.1909,-5.4298,0;13.8993,-6.6466,0;14.6693,-7.6321,0;15.9193,-7.6321,0;16.7095,-6.6767,0;16.4397,-5.4597,0;15.3177,-4.9055,0;16.5443,-8.7146,0;14.0443,-8.7146,0;13.3061,-4.5449,0;12.0973,-4.8688,0;13.63,-3.3361,0;12.4212,-3.66,0;",
            longKey: "Long-RInChIKey=SA-FUHFF-YCJXIOAKVCHNQZ-UHFFFAOYSA-N--YCJXIOAKVCHNQZ-UHFFFAOYSA-N",
            shortKey: "Short-RInChIKey=SA-FUHFF-YCJXIOAKVC-YCJXIOAKVC-UHFFFADPSC-NUHFF-NUHFF-NUHFF-ZZZ",
            webKey: "Web-RInChIKey=YCJXIOAKVCHNQZAOJ-NUHFFFADPSCTJSA"
        ),
        "79cfdb69c3392f2fd3f10160921a66a628bdcfb19da219e824741adc47f84326|eq=0": ReferenceOutput(
            rinchi: "RInChI=1.00.1S/<>C5H12O/c1-4(2)5(3)6/h4-6H,1-3H3!C8H8O2/c9-8(10)6-7-4-2-1-3-5-7/h1-5H,6H2,(H,9,10)/d-",
            rAuxInfo: "RAuxInfo=1.00.1/<>0/N:1,3,5,2,4,6/E:(1,2)/rA:6nCCCCCO/rB:s1;s2;s2;s4;s4;/rC:5.2741,1.4438,0;5.9886,1.8563,0;6.703,1.4438,0;5.9886,2.6813,0;6.703,3.0938,0;5.2741,3.0938,0;!1/N:3,2,4,1,5,7,6,8,9,10/E:(2,3)(4,5)(9,10)/rA:10nCCCCCCCCOO/rB:d1;s2;d3;s4;s1d5;s6;s7;s8;d8;/rC:-3.8009,2.917,0;-4.5154,2.5045,0;-4.5154,1.6795,0;-3.8009,1.2669,0;-3.0864,1.6795,0;-3.0864,2.5045,0;-2.372,2.917,0;-1.6575,2.5045,0;-.943,2.917,0;-1.6575,1.6795,0;",
            longKey: "Long-RInChIKey=SA-BUHFF---MXLMTQWGSQIYOW-UHFFFAOYSA-N-WLJVXDMOQOGPHL-UHFFFAOYSA-N",
            shortKey: "Short-RInChIKey=SA-BUHFF-UHFFFADPSC-RLJDZPNDNR-UHFFFADPSC-NUHFF-NUHFF-NUHFF-ZZZ",
            webKey: "Web-RInChIKey=RLJDZPNDNRCEFOHCO-NUHFFFADPSCTJSA"
        ),
        "29a364394b4fbe2bf98bc287454f216f68afd8705665e628e7879e2810e19159|eq=0": ReferenceOutput(
            rinchi: "RInChI=1.00.1S//d+/u1-1-0",
            rAuxInfo: "RAuxInfo=1.00.1/",
            longKey: "Long-RInChIKey=SA-FUHFF-MOSFIJXAXDLOML-UHFFFAOYSA-N--MOSFIJXAXDLOML-UHFFFAOYSA-N",
            shortKey: "Short-RInChIKey=SA-FUHFF-UHFFFADPSC-UHFFFADPSC-UHFFFADPSC-NUHFF-NUHFF-NUHFF-AAZ",
            webKey: "Web-RInChIKey=MOSFIJXAXDLOMLMKR-NUHFFFADPSCTJSA"
        ),
        "d260eb67304cbf96c40ecafe5bededa44a643248fa4a6c2d99b231458a604366|eq=0": ReferenceOutput(
            rinchi: "RInChI=1.00.1S//d+/u1-1-0",
            rAuxInfo: "RAuxInfo=1.00.1/",
            longKey: "Long-RInChIKey=SA-FUHFF-MOSFIJXAXDLOML-UHFFFAOYSA-N--MOSFIJXAXDLOML-UHFFFAOYSA-N",
            shortKey: "Short-RInChIKey=SA-FUHFF-UHFFFADPSC-UHFFFADPSC-UHFFFADPSC-NUHFF-NUHFF-NUHFF-AAZ",
            webKey: "Web-RInChIKey=MOSFIJXAXDLOMLMKR-NUHFFFADPSCTJSA"
        ),
        "04ad982c16ea0e7f97627fd9611ad2b6e923c1db6287b3c0978530a8265b88cf|eq=0": ReferenceOutput(
            rinchi: "RInChI=1.00.1S/C12H17NOS/c1-3-13(4-2)12(14)15-10-11-8-6-5-7-9-11/h5-9H,3-4,10H2,1-2H3<>C4H11N/c1-3-5-4-2/h5H,3-4H2,1-2H3!C7H7Br/c8-6-7-4-2-1-3-5-7/h1-5H,6H2!CO/c1-2!S8/c1-2-4-6-8-7-5-3-1/d-",
            rAuxInfo: "RAuxInfo=1.00.1/0/N:11,12,7,8,15,13,14,10,9,5,6,1,3,4,2/E:(1,2)(3,4)(6,7)(8,9)/rA:15nCSNOCCCCCCCCCCC/rB:s1;s1;d1;s2;s5;s3;s3;s6;d6;s7;s8;s10;d9;d13s14;/rC:6.4995,-.9132,0;7.2408,-1.3173,0;5.7672,-1.3151,0;6.4995,-.0915,0;7.8347,-.7346,0;8.6251,-1.036,0;5.7381,-2.1323,0;5.0706,-.8797,0;9.248,-.4979,0;8.7791,-1.8443,0;4.9544,-2.3779,0;4.4275,-1.3888,0;9.5584,-2.1144,0;10.0295,-.7725,0;10.1813,-1.5741,0;<>0/N:4,5,2,3,1/E:(1,2)(3,4)/rA:5nNCCCC/rB:s1;s1;s2;s3;/rC:-8.5447,-1.3821,0;-9.2592,-.9578,0;-7.8213,-.969,0;-9.9737,-1.3709,0;-7.1158,-1.3888,0;!0/N:8,6,7,4,5,3,2,1/E:(2,3)(4,5)/rA:8nBrCCCCCCC/rB:;s1s2;s2;d2;d4;s5;s6d7;/rC:.326,.2121,0;1.0941,-.9266,0;1.0941,-.096,0;1.8063,-1.3396,0;.3774,-1.3374,0;1.8063,-2.1635,0;.3751,-2.1635,0;1.0896,-2.5811,0;!0/N:2,1/CRV:1-1,2+1/rA:2nO+C-/rB:t1;/rC:-4.7982,-1.2503,0;-5.6243,-1.2503,0;!0/N:1,2,3,4,5,6,7,8/E:(1,2,3,4,5,6,7,8)/rA:8nSSSSSSSS/rB:s1;s1;s2;s3;s4;s5;s6s7;/rC:-3.2397,-1.583,0;-3.2397,-.7658,0;-2.6548,-2.1792,0;-2.6548,-.1786,0;-1.8242,-2.1792,0;-1.8242,-.1786,0;-1.2392,-1.583,0;-1.2392,-.7658,0;",
            longKey: "Long-RInChIKey=SA-BUHFF-BKXRWLMGHBHCFY-UHFFFAOYSA-N--HPNMFZURTQLUMO-UHFFFAOYSA-N-AGEZXYOZHKGVCM-UHFFFAOYSA-N-UGFAIRIUMAVXCW-UHFFFAOYSA-N-JLQNHALFVCURHW-UHFFFAOYSA-N",
            shortKey: "Short-RInChIKey=SA-BUHFF-BKXRWLMGHB-YBJAEVKSZP-UHFFFADPSC-NUHFF-NUHFF-NUHFF-ZZZ",
            webKey: "Web-RInChIKey=VQKJUKGWGYEYITGFK-NUHFFFADPSCTJSA"
        ),
        "4760737997dca6b897648d16d97ba26aa1048c50100b2c7fe0f06d2201da1f6f|eq=0": ReferenceOutput(
            rinchi: "RInChI=1.00.1S/C6H10O/c7-6-4-2-1-3-5-6/h1-5H2<>C6H12O/c7-6-4-2-1-3-5-6/h6-7H,1-5H2/d+",
            rAuxInfo: "RAuxInfo=1.00.1/0/N:7,5,6,4,3,2,1/E:(2,3)(4,5)/rA:7nOCCCCCC/rB:d1;s2;s2;s4;s3;s5s6;/rC:-3.475,.0042,0;-3.475,-1.5042,0;-4.7917,-2.2917,0;-2.1208,-2.2792,0;-2.1333,-3.825,0;-4.8,-3.825,0;-3.4583,-4.6083,0;<>0/N:7,5,6,4,3,2,1/E:(2,3)(4,5)/rA:7nOCCCCCC/rB:s1;s2;s2;s4;s3;s5s6;/rC:3.9625,.0542,0;3.9667,-1.5042,0;2.6458,-2.2917,0;5.3,-2.2792,0;5.3,-3.8042,0;2.6458,-3.8333,0;3.9792,-4.5875,0;",
            longKey: "Long-RInChIKey=SA-FUHFF-JHIVVAPYMSGYDF-UHFFFAOYSA-N--HPXRVTGHNJAIIH-UHFFFAOYSA-N",
            shortKey: "Short-RInChIKey=SA-FUHFF-JHIVVAPYMS-HPXRVTGHNJ-UHFFFADPSC-NUHFF-NUHFF-NUHFF-ZZZ",
            webKey: "Web-RInChIKey=CSSSVQZNSGVAKJGLJ-NUHFFFADPSCTJSA"
        ),
        "b1396509837e0c8d059f26ad6bb7b05876938460905125e309d73112936df0d1|eq=0": ReferenceOutput(
            rinchi: "RInChI=1.00.1S/C7H13BrN2O2/c1-3-7(8,4-2)5(11)10-6(9)12/h3-4H2,1-2H3,(H3,9,10,11,12)<>C7H14N2O2/c1-3-5(4-2)6(10)9-7(8)11/h5H,3-4H2,1-2H3,(H3,8,9,10,11)/d+",
            rAuxInfo: "RAuxInfo=1.00.1/1/N:8,9,3,4,2,10,1,5,11,6,7,12/E:(1,2)(3,4)/rA:12nCCCCBrNOCCCNO/rB:s1;s1;s1;s1;s2;d2;s3;s4;s6;s10;d10;/rC:-7.5503,-2.2909,0;-6.17,-1.5282,0;-8.3781,-.9999,0;-8.8383,-3.125,0;-6.7789,-3.5938,0;-4.8604,-2.4021,0;-6.0954,-.0187,0;-9.8678,-1.3455,0;-8.4874,-4.6027,0;-3.5734,-1.5401,0;-2.199,-2.2472,0;-3.6687,-.0108,0;<>1/N:7,8,3,4,1,2,9,10,5,6,11/E:(1,2)(3,4)/rA:11nCCCCNOCCCNO/rB:s1;s1;s1;s2;d2;s3;s4;s5;s9;d9;/rC:4.9875,-2.3838,0;6.3486,-1.6259,0;3.6779,-1.5704,0;4.9201,-3.9353,0;7.6859,-2.3679,0;6.3486,-.0862,0;2.4715,-2.5306,0;6.3248,-4.5265,0;9.0709,-1.6656,0;10.3645,-2.4869,0;9.0947,-.1378,0;",
            longKey: "Long-RInChIKey=SA-FUHFF-OPNPQXLQERQBBV-UHFFFAOYSA-N--JIOUJVPLJOQUFD-UHFFFAOYSA-N",
            shortKey: "Short-RInChIKey=SA-FUHFF-OPNPQXLQER-JIOUJVPLJO-UHFFFADPSC-NUHFF-NUHFF-NUHFF-ZZZ",
            webKey: "Web-RInChIKey=LVIVPLDRXYMEIMSPG-NUHFFFADPSCTJSA"
        ),
        "cb209ba172120cc61bf29a3ee154f01ab05c4cdcdfa1f5034ea12965fbcf0b7d|eq=0": ReferenceOutput(
            rinchi: "RInChI=1.00.1S//d+/u1-1-0",
            rAuxInfo: "RAuxInfo=1.00.1/",
            longKey: "Long-RInChIKey=SA-FUHFF-MOSFIJXAXDLOML-UHFFFAOYSA-N--MOSFIJXAXDLOML-UHFFFAOYSA-N",
            shortKey: "Short-RInChIKey=SA-FUHFF-UHFFFADPSC-UHFFFADPSC-UHFFFADPSC-NUHFF-NUHFF-NUHFF-AAZ",
            webKey: "Web-RInChIKey=MOSFIJXAXDLOMLMKR-NUHFFFADPSCTJSA"
        ),
        "30e293a195aeb2a2df6e149a3ec249f49d43e1d02570d30e63cc6164ab48f384|eq=0": ReferenceOutput(
            rinchi: "RInChI=1.00.1S//d+/u1-1-0",
            rAuxInfo: "RAuxInfo=1.00.1/",
            longKey: "Long-RInChIKey=SA-FUHFF-MOSFIJXAXDLOML-UHFFFAOYSA-N--MOSFIJXAXDLOML-UHFFFAOYSA-N",
            shortKey: "Short-RInChIKey=SA-FUHFF-UHFFFADPSC-UHFFFADPSC-UHFFFADPSC-NUHFF-NUHFF-NUHFF-AAZ",
            webKey: "Web-RInChIKey=MOSFIJXAXDLOMLMKR-NUHFFFADPSCTJSA"
        ),
        "22c091f74c87f8118866977abc1dc29182906ddf6e1dbcf37f622c7a6f45ba69|eq=0": ReferenceOutput(
            rinchi: "RInChI=1.00.1S//d+/u2-1-0",
            rAuxInfo: "RAuxInfo=1.00.1/",
            longKey: "Long-RInChIKey=SA-FUHFF-MOSFIJXAXDLOML-UHFFFAOYSA-N-MOSFIJXAXDLOML-UHFFFAOYSA-N--MOSFIJXAXDLOML-UHFFFAOYSA-N",
            shortKey: "Short-RInChIKey=SA-FUHFF-UHFFFADPSC-UHFFFADPSC-UHFFFADPSC-NUHFF-NUHFF-NUHFF-BAZ",
            webKey: "Web-RInChIKey=MOSFIJXAXDLOMLMKR-NUHFFFADPSCTJSA"
        ),
        "78f8ea52241cc4436854681991ad7504270abd24e119ee7ef44877089e1b1da2|eq=0": ReferenceOutput(
            rinchi: "RInChI=1.00.1S//d+/u1-1-1",
            rAuxInfo: "RAuxInfo=1.00.1/",
            longKey: "Long-RInChIKey=SA-FUHFF-MOSFIJXAXDLOML-UHFFFAOYSA-N--MOSFIJXAXDLOML-UHFFFAOYSA-N--MOSFIJXAXDLOML-UHFFFAOYSA-N",
            shortKey: "Short-RInChIKey=SA-FUHFF-UHFFFADPSC-UHFFFADPSC-UHFFFADPSC-NUHFF-NUHFF-NUHFF-AAA",
            webKey: "Web-RInChIKey=MOSFIJXAXDLOMLMKR-NUHFFFADPSCTJSA"
        ),
        "3fe74d58674f3cec9752ab3fa06b3ae4c88e0f45c40be53e77da8f32ed004ddf|eq=0": ReferenceOutput(
            rinchi: "RInChI=1.00.1S/C20H20N2O8S/c1-10(23)29-8-13-9-31-19-14(18(26)22(19)15(13)20(27)28)21-17(25)16(30-11(2)24)12-6-4-3-5-7-12/h3-7,14,16,19H,8-9H2,1-2H3,(H,21,25)(H,27,28)<>C8H10N2O3S/c1-3-2-14-7-4(9)6(11)10(7)5(3)8(12)13/h4,7H,2,9H2,1H3,(H,12,13)<>H2O/h1H2/d+",
            rAuxInfo: "RAuxInfo=1.00.1/1/N:1,15,20,19,21,18,22,4,6,2,14,17,5,9,27,12,11,24,8,28,10,26,31,16,23,25,29,30,3,13,7/E:(4,5)(6,7)(27,28)/rA:31nCCOCCCSCCNCCOCCOCCCCCCOCONCCOOO/rB:s1;s2;s3;s4;s5;s6;s7;s8;s9;s10;s11;s12;s13;s14;d14;s12;s17;d18;s19;d20;d17s21;d11;s9;d24;s8s24;d5s26;s27;s28;d28;d2;/rC:;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;<>1/N:1,3,2,6,11,8,5,12,7,10,9,13,14,4/E:(12,13)/rA:14nCCCSCCNCONCCOO/rB:s1;s2;s3;s4;s5;s6;s6;d8;s5s8;d2s10;s11;s12;d12;/rC:;;;;;;;;;;;;;;<>0/N:1/rA:1nO/rB:/rC:;",
            longKey: "Long-RInChIKey=SA-FUHFF-XRVPSMBZQKYMGA-UHFFFAOYSA-N--NVIAYEIXYQCDAN-UHFFFAOYSA-N--XLYOFNOQVPJJNP-UHFFFAOYSA-N",
            shortKey: "Short-RInChIKey=SA-FUHFF-XRVPSMBZQK-NVIAYEIXYQ-XLYOFNOQVP-NUHFF-NUHFF-NUHFF-ZZZ",
            webKey: "Web-RInChIKey=DEQXOIGUEKAJMSFJR-NUHFFFADPSCTJSA"
        ),
        "681212a529d68931978a173d00a05a0a1ab5ed92cf5980f0f945a47a3dc713fd|eq=0": ReferenceOutput(
            rinchi: "RInChI=1.00.1S/C4H4O2/c1-2-3-4(5)6/h1H3,(H,5,6)<>C4H5IO2/c1-3(5)2-4(6)7/h2H,1H3,(H,6,7)/d+",
            rAuxInfo: "RAuxInfo=1.00.1/1/N:3,4,5,6,1,2/E:(5,6)/rA:6nOOCCCC/rB:;;s3;t4;s1d2s5;/rC:-10.4009,-.7145,0;-10.4009,.7145,0;-13.2884,0,0;-12.4634,0,0;-11.6384,0,0;-10.8134,0,0;<>1/N:1,5,6,4,7,2,3/E:(6,7)/rA:7nCOOCCCI/rB:;;s2d3;v4;s1d5;s6;/rC:-3.9152,-.1906,0;-6.3902,-.1906,0;-5.1527,-.905,0;-5.5652,-.1906,0;-5.1527,.5239,0;-4.3277,.5239,0;-3.9152,1.2384,0;",
            longKey: "Long-RInChIKey=SA-FUHFF-LUEHNHVFDCZTGL-UHFFFAOYSA-N--CSDFWNXJFJAWAM-UHFFFAOYSA-N",
            shortKey: "Short-RInChIKey=SA-FUHFF-LUEHNHVFDC-CSDFWNXJFJ-UHFFFADPSC-NUHFF-NUHFF-NUHFF-ZZZ",
            webKey: "Web-RInChIKey=XOQQWAOKZSQZPLRLR-NUHFFFADPSCTJSA"
        ),
        "23bc43a9f9f897d13ddfa5fc09f3b4b5f289e9b8553adcc8d936cc21cb380d95|eq=1": ReferenceOutput(
            rinchi: "RInChI=1.00.1S/C18H15P/c1-4-10-16(11-5-1)19(17-12-6-2-7-13-17)18-14-8-3-9-15-18/h1-15H!C2H6O/c1-2-3/h3H,2H2,1H3!CCl4/c2-1(3,4)5<>C2H5Cl/c1-2-3/h2H2,1H3/d=",
            rAuxInfo: "RAuxInfo=1.00.1/0/N:5,11,17,4,6,10,12,16,18,3,7,9,13,15,19,2,8,14,1/E:(1,2,3)(4,5,6,7,8,9)(10,11,12,13,14,15)(16,17,18)/rA:19nPCCCCCCCCCCCCCCCCCC/rB:s1;d2;s3;d4;s5;s2d6;s1;d8;s9;d10;s11;s8d12;s1;d14;s15;d16;s17;s14d18;/rC:0,-.4125,0;0,.4125,0;.7145,.825,0;.7145,1.65,0;0,2.0625,0;-.7145,1.65,0;-.7145,.825,0;-.7145,-.825,0;-1.4289,-.4125,0;-2.1434,-.825,0;-2.1434,-1.65,0;-1.4289,-2.0625,0;-.7145,-1.65,0;.7145,-.825,0;.7145,-1.65,0;1.4289,-2.0625,0;2.1434,-1.65,0;2.1434,-.825,0;1.4289,-.4125,0;!0/N:1,2,3/rA:3nCCO/rB:s1;s2;/rC:-.7145,-.2062,0;0,.2062,0;.7145,-.2062,0;!0/N:1,2,3,4,5/E:(2,3,4,5)/rA:5nCClClClCl/rB:s1;s1;s1;s1;/rC:;-.825,0,0;0,-.825,0;.825,0,0;0,.825,0;<>0/N:1,2,3/rA:3nCCCl/rB:s1;s2;/rC:-.7145,-.2062,0;0,.2062,0;.7145,-.2062,0;",
            longKey: "Long-RInChIKey=SA-EUHFF-RIOQSEWOXXDEQQ-UHFFFAOYSA-N-LFQSCWFLJHTTHZ-UHFFFAOYSA-N-VZGDMQKNWNREIO-UHFFFAOYSA-N--HRYZWHHZPQKTII-UHFFFAOYSA-N",
            shortKey: "Short-RInChIKey=SA-EUHFF-ALHSIEQSCJ-HRYZWHHZPQ-UHFFFADPSC-NUHFF-NUHFF-NUHFF-ZZZ",
            webKey: "Web-RInChIKey=MYFTVVMFPXRTCFCZY-NUHFFFADPSCTJSA"
        ),
        "e856a0bd47680cfb5b88977f2ef6bfcfa6380d2fd5e1a9af4da2320711cb38c5|eq=1": ReferenceOutput(
            rinchi: "RInChI=1.00.1S/C4H4O5/c5-2(4(8)9)1-3(6)7/h1H2,(H,6,7)(H,8,9)<>C4H6O5/c5-2(4(8)9)1-3(6)7/h2,5H,1H2,(H,6,7)(H,8,9)/t2-/m1/s1/d=",
            rAuxInfo: "RAuxInfo=1.00.1/1/N:3,4,2,5,8,1,7,6,9/E:(6,7)(8,9)/rA:9nOCCCCOOOO/rB:s1;s2;s3;s4;s5;d2;d4;d5;/rC:16.7991,-10.3667,0;17.5136,-10.7792,0;18.2281,-10.3667,0;18.9426,-10.7792,0;19.657,-10.3667,0;20.3715,-10.7792,0;17.5136,-11.6041,0;18.9426,-11.6042,0;19.657,-9.5417,0;<>1/N:2,1,5,3,4,8,9,6,7/E:(6,7)(8,9)/it:im/rA:9nCCCOCOOOO/rB:s1;s1;P1;s2;s3;d3;s5;d5;/rC:10.7707,-10.7381,0;10.0562,-11.1506,0;11.4852,-11.1506,0;10.7707,-9.9131,0;9.3418,-10.7381,0;12.1997,-10.7381,0;11.4852,-11.9756,0;8.6273,-11.1506,0;9.3418,-9.9131,0;",
            longKey: "Long-RInChIKey=SA-BUHFF-KHPXUQMNIQBQEV-UHFFFAOYSA-N--BJEPYKJPYRNKOW-UWTATZPHSA-N",
            shortKey: "Short-RInChIKey=SA-BUHFF-KHPXUQMNIQ-BJEPYKJPYR-UHFFFADPSC-NUHFF-NWBEM-NUHFF-ZZZ",
            webKey: "Web-RInChIKey=NJVAQMRLSHCPOPHDV-NWBEMOFYMKSCHSA"
        ),
        "425ff38195a4848b2f30ef75fed896ab078e076da4adeeefa3c072d2dfe50715|eq=1": ReferenceOutput(
            rinchi: "RInChI=1.00.1S/C4H4O5/c5-2(4(8)9)1-3(6)7/h1H2,(H,6,7)(H,8,9)<>C4H6O5/c5-2(4(8)9)1-3(6)7/h2,5H,1H2,(H,6,7)(H,8,9)/t2-/m1/s1/d=",
            rAuxInfo: "RAuxInfo=1.00.1/1/N:3,4,2,5,8,1,7,6,9/E:(6,7)(8,9)/rA:9nOCCCCOOOO/rB:s1;s2;s3;s4;s5;d2;d4;d5;/rC:16.7991,-10.3667,0;17.5136,-10.7792,0;18.2281,-10.3667,0;18.9426,-10.7792,0;19.657,-10.3667,0;20.3715,-10.7792,0;17.5136,-11.6041,0;18.9426,-11.6042,0;19.657,-9.5417,0;<>1/N:2,1,5,3,4,8,9,6,7/E:(6,7)(8,9)/it:im/rA:9cCCCOCOOOO/rB:s1;s1;P1;s2;s3;d3;s5;d5;/rC:10.7707,-10.7381,0;10.0562,-11.1506,0;11.4852,-11.1506,0;10.7707,-9.9131,0;9.3418,-10.7381,0;12.1997,-10.7381,0;11.4852,-11.9756,0;8.6273,-11.1506,0;9.3418,-9.9131,0;",
            longKey: "Long-RInChIKey=SA-EUHFF-KHPXUQMNIQBQEV-UHFFFAOYSA-N--BJEPYKJPYRNKOW-UWTATZPHSA-N",
            shortKey: "Short-RInChIKey=SA-EUHFF-KHPXUQMNIQ-BJEPYKJPYR-UHFFFADPSC-NUHFF-NWBEM-NUHFF-ZZZ",
            webKey: "Web-RInChIKey=NJVAQMRLSHCPOPHDV-NWBEMOFYMKSCHSA"
        ),
    ]
}

private extension Array {
    subscript(safe index: Int) -> Element? {
        indices.contains(index) ? self[index] : nil
    }
}
