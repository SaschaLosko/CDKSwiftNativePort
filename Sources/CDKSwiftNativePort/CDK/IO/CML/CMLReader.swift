import Foundation
#if canImport(FoundationXML)
import FoundationXML
#endif
#if canImport(CoreGraphics)
import CoreGraphics
#endif

/// CDK-style CML reader supporting common atomArray/bondArray variants emitted by CDK and other toolkits.
public enum CDKCMLReader {
    public static func read(text: String) throws -> [Molecule] {
        let data = Data(text.utf8)
        let delegate = CMLParserDelegate()

        let parser = XMLParser(data: data)
        parser.delegate = delegate
        guard parser.parse() else {
            let message = parser.parserError?.localizedDescription ?? "Malformed CML document."
            throw ChemError.parseFailed(message)
        }

        let molecules = delegate.buildMolecules()
        guard !molecules.isEmpty else {
            throw ChemError.parseFailed("CML did not contain any atoms.")
        }
        return molecules
    }
}

private final class CMLParserDelegate: NSObject, XMLParserDelegate {
    private final class PendingMolecule {
        var id: String?
        var title: String?
        var atoms: [PendingAtom] = []
        var bonds: [PendingBond] = []

        init(id: String?, title: String?) {
            self.id = id
            self.title = title
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
        var chirality: AtomChirality
        var ligandOrderingRefs: [String]?
    }

    private struct PendingBond {
        let refA: String
        let refB: String
        var order: BondOrder
        var stereo: BondStereo
        var doubleBondStereo: DoubleBondStereo?
        var stereoReferenceAtomRefs: [String]?
        var aromatic: Bool
        var hadExplicitOrder: Bool
    }

    private struct TextCapture {
        let elementName: String
        let attributes: [String: String]
        var text: String
    }

    private struct ParsedCoordinates {
        let position: CGPoint
        let zPosition: Double?
    }

    private var moleculeStack: [PendingMolecule] = []
    private var parsedMolecules: [PendingMolecule] = []
    private var currentAtom: PendingAtom?
    private var currentBond: PendingBond?
    private var textCaptures: [TextCapture] = []
    private var generatedAtomCounter = 1

    func parser(_ parser: XMLParser,
                didStartElement elementName: String,
                namespaceURI: String?,
                qualifiedName qName: String?,
                attributes attributeDict: [String: String] = [:]) {
        let lower = localName(elementName).lowercased()

        switch lower {
        case "molecule":
            let molecule = PendingMolecule(id: attributeDict["id"],
                                           title: attributeDict["title"] ?? attributeDict["name"])
            moleculeStack.append(molecule)

        case "atom":
            ensureCurrentMolecule()
            currentAtom = makeAtom(attributes: attributeDict)

        case "bond":
            ensureCurrentMolecule()
            currentBond = makeBond(attributes: attributeDict)

        case "atomarray":
            ensureCurrentMolecule()
            guard let molecule = moleculeStack.last else { return }
            applyVectorAtomArray(attributes: attributeDict, to: molecule)

        case "bondarray":
            ensureCurrentMolecule()
            guard let molecule = moleculeStack.last else { return }
            applyVectorBondArray(attributes: attributeDict, to: molecule)

        case "bondtype":
            guard currentBond != nil else { return }
            if isAromaticDictRef(attributeDict["dictRef"]) {
                currentBond?.aromatic = true
            }

        case "bondstereo", "atomparity", "label", "scalar", "string", "float", "integer":
            guard currentAtom != nil || currentBond != nil else { return }
            textCaptures.append(TextCapture(elementName: lower, attributes: attributeDict, text: ""))

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
            guard let molecule = moleculeStack.last, let atom = currentAtom else { return }
            molecule.atoms.append(atom)
            currentAtom = nil

        case "bond":
            guard let molecule = moleculeStack.last, let bond = currentBond else { return }
            molecule.bonds.append(bond)
            currentBond = nil

        case "molecule":
            guard let molecule = moleculeStack.popLast() else { return }
            if !molecule.atoms.isEmpty {
                parsedMolecules.append(molecule)
            }

        case "bondstereo", "atomparity", "label", "scalar", "string", "float", "integer":
            guard let capture = popCapture(for: lower) else { return }
            applyCapture(capture)

        default:
            break
        }
    }

    func buildMolecules() -> [Molecule] {
        var pending = parsedMolecules
        for trailing in moleculeStack where !trailing.atoms.isEmpty {
            pending.append(trailing)
        }

        var output: [Molecule] = []
        output.reserveCapacity(pending.count)

        for (idx, source) in pending.enumerated() {
            let fallbackName = "CML Molecule \(idx + 1)"
            let name = source.title?.trimmingCharacters(in: .whitespacesAndNewlines)
            let moleculeName = (name?.isEmpty == false) ? name! : (source.id?.isEmpty == false ? source.id! : fallbackName)

            var atoms: [Atom] = []
            var atomIDByXMLID: [String: Int] = [:]
            atoms.reserveCapacity(source.atoms.count)

            for (atomIndex, atomSource) in source.atoms.enumerated() {
                atomIDByXMLID[atomSource.xmlID] = atomIndex + 1
            }

            for (atomIndex, atomSource) in source.atoms.enumerated() {
                let atomID = atomIndex + 1
                atoms.append(Atom(id: atomID,
                                  element: atomSource.element,
                                  position: atomSource.position,
                                  zPosition: atomSource.zPosition,
                                  charge: atomSource.charge,
                                  isotopeMassNumber: atomSource.isotope,
                                  aromatic: atomSource.aromatic,
                                  chirality: atomSource.chirality,
                                  explicitHydrogenCount: atomSource.hydrogenCount,
                                  ligandOrderingAtomIDs: atomSource.ligandOrderingRefs?.compactMap { atomIDByXMLID[$0] },
                                  aliasLabel: atomSource.aliasLabel))
            }

            var bonds: [Bond] = []
            var aromaticIDs = Set(atoms.filter(\.aromatic).map(\.id))
            var nextBondID = 1
            for bondSource in source.bonds {
                guard let a1 = atomIDByXMLID[bondSource.refA],
                      let a2 = atomIDByXMLID[bondSource.refB],
                      a1 != a2 else {
                    continue
                }
                bonds.append(Bond(id: nextBondID,
                                  a1: a1,
                                  a2: a2,
                                  order: bondSource.order,
                                  stereo: bondSource.stereo,
                                  doubleBondStereo: bondSource.doubleBondStereo,
                                  stereoReferenceAtomIDs: bondSource.stereoReferenceAtomRefs?.compactMap {
                                      atomIDByXMLID[$0]
                                  }))
                if bondSource.aromatic || bondSource.order == .aromatic {
                    aromaticIDs.insert(a1)
                    aromaticIDs.insert(a2)
                }
                nextBondID += 1
            }

            if bonds.isEmpty {
                bonds = CDKBondPerception.inferSingleBonds(for: atoms)
            }

            for index in atoms.indices where aromaticIDs.contains(atoms[index].id) {
                atoms[index].aromatic = true
            }

            var molecule = Molecule(name: moleculeName, atoms: atoms, bonds: bonds)
            if let box = molecule.boundingBox(), box.width <= 0.0001 && box.height <= 0.0001 {
                molecule = Depiction2DGenerator.generate(for: molecule)
            }
            molecule.assignWedgeHashFromChiralCenters()
            output.append(molecule)
        }

        return output
    }

    private func ensureCurrentMolecule() {
        if moleculeStack.isEmpty {
            moleculeStack.append(PendingMolecule(id: nil, title: nil))
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
                           aliasLabel: aliasLabel,
                           chirality: .none,
                           ligandOrderingRefs: nil)
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
        return PendingBond(refA: refs[0],
                           refB: refs[1],
                           order: order,
                           stereo: .none,
                           doubleBondStereo: nil,
                           stereoReferenceAtomRefs: nil,
                           aromatic: order == .aromatic || isTruthy(attributes["aromatic"]),
                           hadExplicitOrder: rawOrder?.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty == false)
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
                                              aliasLabel: nil,
                                              chirality: .none,
                                              ligandOrderingRefs: nil))
        }
    }

    private func applyVectorBondArray(attributes: [String: String], to molecule: PendingMolecule) {
        let orders = splitTokens(attributes["order"])

        if let atomRef1 = attributes["atomRef1"], let atomRef2 = attributes["atomRef2"] {
            let refs1 = atomRef1.split(whereSeparator: \.isWhitespace).map(String.init)
            let refs2 = atomRef2.split(whereSeparator: \.isWhitespace).map(String.init)
            let count = min(refs1.count, refs2.count)
            guard count > 0 else { return }

            for idx in 0..<count {
                let rawOrder = idx < orders.count ? orders[idx] : nil
                let order = mapBondOrder(rawOrder)
                molecule.bonds.append(PendingBond(refA: refs1[idx],
                                                  refB: refs2[idx],
                                                  order: order,
                                                  stereo: .none,
                                                  doubleBondStereo: nil,
                                                  stereoReferenceAtomRefs: nil,
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
                molecule.bonds.append(PendingBond(refA: tokens[idx],
                                                  refB: tokens[idx + 1],
                                                  order: order,
                                                  stereo: .none,
                                                  doubleBondStereo: nil,
                                                  stereoReferenceAtomRefs: nil,
                                                  aromatic: order == .aromatic,
                                                  hadExplicitOrder: rawOrder?.isEmpty == false))
                pairIndex += 1
                idx += 2
            }
        }
    }

    private func applyCapture(_ capture: TextCapture) {
        switch capture.elementName {
        case "bondstereo":
            guard currentBond != nil else { return }
            if let doubleBondStereo = parseDoubleBondStereo(attributes: capture.attributes, text: capture.text) {
                currentBond?.doubleBondStereo = doubleBondStereo
                let refs = capture.attributes["atomRefs4"]?
                    .split(whereSeparator: \.isWhitespace)
                    .map(String.init)
                currentBond?.stereoReferenceAtomRefs = refs?.count == 4 ? refs : nil
            }
            if let stereo = parseBondStereo(attributes: capture.attributes, text: capture.text) {
                currentBond?.stereo = stereo
            }

        case "atomparity":
            guard currentAtom != nil else { return }
            applyAtomParityCapture(capture)

        case "label":
            guard currentAtom != nil else { return }
            applyAtomLabel(capture.text)

        case "scalar", "string", "float", "integer":
            if currentAtom != nil {
                applyAtomScalarCapture(capture)
            } else if currentBond != nil {
                applyBondScalarCapture(capture)
            }

        default:
            break
        }
    }

    private func applyAtomLabel(_ rawText: String) {
        guard let label = normalizedLabel(rawText), !label.isEmpty else { return }
        currentAtom?.aliasLabel = label
        if let atom = currentAtom, isPseudoElement(atom.element) || atom.element == "*" || atom.element.isEmpty {
            currentAtom?.element = "R"
        }
    }

    private func applyAtomParityCapture(_ capture: TextCapture) {
        let parityText = capture.text.trimmingCharacters(in: .whitespacesAndNewlines)
        let parity = Int(parityText) ?? 0
        switch parity.signum() {
        case 1:
            currentAtom?.chirality = .anticlockwise
        case -1:
            currentAtom?.chirality = .clockwise
        default:
            currentAtom?.chirality = .none
        }

        if let refs = capture.attributes["atomRefs4"]?
            .split(whereSeparator: \.isWhitespace)
            .map(String.init),
           refs.count == 4 {
            currentAtom?.ligandOrderingRefs = refs
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
            if let value = Double(text) {
                currentAtom?.position.x = value
            }
        case "y2":
            if let value = Double(text) {
                currentAtom?.position.y = value
            }
        case "x3":
            if let value = Double(text) {
                currentAtom?.position.x = value
            }
        case "y3":
            if let value = Double(text) {
                currentAtom?.position.y = value
            }
        case "z3":
            if let value = Double(text) {
                currentAtom?.zPosition = value
            }
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
            if let charge = Int(text) {
                currentAtom?.charge = charge
            }
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
        case "label":
            break
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
        if normalized == "DU" || normalized == "DUMMY" {
            return "R"
        }
        if normalized == "R" || normalized == "R#" {
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

    private func parseDoubleBondStereo(attributes: [String: String], text: String) -> DoubleBondStereo? {
        let raw = attributes["dictRef"].map { String($0.split(separator: ":").last ?? "") } ?? text
        switch raw.trimmingCharacters(in: .whitespacesAndNewlines).uppercased() {
        case "C", "CIS":
            return .cis
        case "T", "TRANS":
            return .trans
        default:
            return nil
        }
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
