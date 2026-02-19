import Foundation
import CoreGraphics

/// CDK-style CML reader supporting atomArray/bondArray forms used by common CML exports.
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
        let element: String
        let position: CGPoint
        let charge: Int
        let isotope: Int?
        let hydrogenCount: Int?
        let aromatic: Bool
    }

    private struct PendingBond {
        let refA: String
        let refB: String
        let order: BondOrder
    }

    private var moleculeStack: [PendingMolecule] = []
    private var parsedMolecules: [PendingMolecule] = []
    private var generatedAtomCounter = 1

    func parser(_ parser: XMLParser,
                didStartElement elementName: String,
                namespaceURI: String?,
                qualifiedName qName: String?,
                attributes attributeDict: [String: String] = [:]) {
        let lower = elementName.lowercased()

        switch lower {
        case "molecule":
            let molecule = PendingMolecule(id: attributeDict["id"],
                                           title: attributeDict["title"] ?? attributeDict["name"])
            moleculeStack.append(molecule)

        case "atom":
            ensureCurrentMolecule()
            guard let molecule = moleculeStack.last,
                  let atom = makeAtom(attributes: attributeDict) else { return }
            molecule.atoms.append(atom)

        case "bond":
            ensureCurrentMolecule()
            guard let molecule = moleculeStack.last,
                  let bond = makeBond(attributes: attributeDict) else { return }
            molecule.bonds.append(bond)

        case "atomarray":
            ensureCurrentMolecule()
            guard let molecule = moleculeStack.last else { return }
            applyVectorAtomArray(attributes: attributeDict, to: molecule)

        case "bondarray":
            ensureCurrentMolecule()
            guard let molecule = moleculeStack.last else { return }
            applyVectorBondArray(attributes: attributeDict, to: molecule)

        default:
            break
        }
    }

    func parser(_ parser: XMLParser,
                didEndElement elementName: String,
                namespaceURI: String?,
                qualifiedName qName: String?) {
        let lower = elementName.lowercased()
        guard lower == "molecule" else { return }

        guard let molecule = moleculeStack.popLast() else { return }
        if !molecule.atoms.isEmpty {
            parsedMolecules.append(molecule)
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
                let atomID = atomIndex + 1
                atomIDByXMLID[atomSource.xmlID] = atomID
                atoms.append(Atom(id: atomID,
                                  element: atomSource.element,
                                  position: atomSource.position,
                                  charge: atomSource.charge,
                                  isotopeMassNumber: atomSource.isotope,
                                  aromatic: atomSource.aromatic,
                                  explicitHydrogenCount: atomSource.hydrogenCount))
            }

            var bonds: [Bond] = []
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
                                  stereo: .none))
                nextBondID += 1
            }

            if bonds.isEmpty {
                bonds = CDKBondPerception.inferSingleBonds(for: atoms)
            }

            if bonds.contains(where: { $0.order == .aromatic }) {
                var aromaticIDs = Set<Int>()
                for bond in bonds where bond.order == .aromatic {
                    aromaticIDs.insert(bond.a1)
                    aromaticIDs.insert(bond.a2)
                }
                for idx in atoms.indices where aromaticIDs.contains(atoms[idx].id) {
                    atoms[idx].aromatic = true
                }
            }

            var molecule = Molecule(name: moleculeName, atoms: atoms, bonds: bonds)
            if let box = molecule.boundingBox(), box.width <= 0.0001 && box.height <= 0.0001 {
                molecule = Depiction2DGenerator.generate(for: molecule)
            }
            output.append(molecule)
        }

        return output
    }

    private func ensureCurrentMolecule() {
        if moleculeStack.isEmpty {
            moleculeStack.append(PendingMolecule(id: nil, title: nil))
        }
    }

    private func makeAtom(attributes: [String: String]) -> PendingAtom? {
        let rawID = attributes["id"]?.trimmingCharacters(in: .whitespacesAndNewlines)
        let xmlID: String
        if let rawID, !rawID.isEmpty {
            xmlID = rawID
        } else {
            xmlID = "a\(generatedAtomCounter)"
            generatedAtomCounter += 1
        }

        let element = CDKDescriptorSupport.canonicalElementSymbol(attributes["elementType"] ?? attributes["element"] ?? "C")

        let x2 = Double(attributes["x2"] ?? "")
        let y2 = Double(attributes["y2"] ?? "")
        let x3 = Double(attributes["x3"] ?? "")
        let y3 = Double(attributes["y3"] ?? "")

        let x = x2 ?? x3 ?? 0
        let y = y2 ?? y3 ?? 0

        let charge = Int(attributes["formalCharge"] ?? "") ?? 0
        let isotope = Int(attributes["isotopeNumber"] ?? "")
        let hydrogenCount = Int(attributes["hydrogenCount"] ?? "")
        let aromatic = (attributes["aromatic"] ?? "").lowercased() == "true"
            || (attributes["atomType"] ?? "").lowercased().contains(".ar")

        return PendingAtom(xmlID: xmlID,
                           element: element,
                           position: CGPoint(x: x, y: y),
                           charge: charge,
                           isotope: isotope,
                           hydrogenCount: hydrogenCount,
                           aromatic: aromatic)
    }

    private func makeBond(attributes: [String: String]) -> PendingBond? {
        if let refs = attributes["atomRefs2"] {
            let tokens = refs.split(whereSeparator: \.isWhitespace).map(String.init)
            guard tokens.count >= 2 else { return nil }
            return PendingBond(refA: tokens[0],
                               refB: tokens[1],
                               order: mapBondOrder(attributes["order"]))
        }

        if let refA = attributes["atomRef1"], let refB = attributes["atomRef2"] {
            return PendingBond(refA: refA,
                               refB: refB,
                               order: mapBondOrder(attributes["order"]))
        }

        return nil
    }

    private func applyVectorAtomArray(attributes: [String: String], to molecule: PendingMolecule) {
        guard let atomIDsRaw = attributes["atomID"] else { return }

        let atomIDs = atomIDsRaw.split(whereSeparator: \.isWhitespace).map(String.init)
        guard !atomIDs.isEmpty else { return }

        let elementTypes = (attributes["elementType"] ?? "")
            .split(whereSeparator: \.isWhitespace)
            .map(String.init)
        let x2 = parseDoubles(attributes["x2"])
        let y2 = parseDoubles(attributes["y2"])
        let x3 = parseDoubles(attributes["x3"])
        let y3 = parseDoubles(attributes["y3"])

        for idx in atomIDs.indices {
            let xmlID = atomIDs[idx]
            if molecule.atoms.contains(where: { $0.xmlID == xmlID }) {
                continue
            }

            let element = idx < elementTypes.count
                ? CDKDescriptorSupport.canonicalElementSymbol(elementTypes[idx])
                : "C"
            let x = idx < x2.count ? x2[idx] : (idx < x3.count ? x3[idx] : 0)
            let y = idx < y2.count ? y2[idx] : (idx < y3.count ? y3[idx] : 0)

            molecule.atoms.append(PendingAtom(xmlID: xmlID,
                                              element: element,
                                              position: CGPoint(x: x, y: y),
                                              charge: 0,
                                              isotope: nil,
                                              hydrogenCount: nil,
                                              aromatic: false))
        }
    }

    private func applyVectorBondArray(attributes: [String: String], to molecule: PendingMolecule) {
        let orders = (attributes["order"] ?? "")
            .split(whereSeparator: \.isWhitespace)
            .map(String.init)

        if let atomRef1 = attributes["atomRef1"], let atomRef2 = attributes["atomRef2"] {
            let refs1 = atomRef1.split(whereSeparator: \.isWhitespace).map(String.init)
            let refs2 = atomRef2.split(whereSeparator: \.isWhitespace).map(String.init)
            let count = min(refs1.count, refs2.count)
            guard count > 0 else { return }

            for idx in 0..<count {
                molecule.bonds.append(PendingBond(refA: refs1[idx],
                                                  refB: refs2[idx],
                                                  order: mapBondOrder(idx < orders.count ? orders[idx] : nil)))
            }
            return
        }

        if let refs = attributes["atomRefs2"] {
            let tokens = refs.split(whereSeparator: \.isWhitespace).map(String.init)
            guard tokens.count >= 2 else { return }

            var pairIndex = 0
            var idx = 0
            while idx + 1 < tokens.count {
                molecule.bonds.append(PendingBond(refA: tokens[idx],
                                                  refB: tokens[idx + 1],
                                                  order: mapBondOrder(pairIndex < orders.count ? orders[pairIndex] : nil)))
                pairIndex += 1
                idx += 2
            }
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

    private func parseDoubles(_ raw: String?) -> [Double] {
        guard let raw else { return [] }
        return raw
            .split(whereSeparator: \.isWhitespace)
            .compactMap { Double($0) }
    }
}
