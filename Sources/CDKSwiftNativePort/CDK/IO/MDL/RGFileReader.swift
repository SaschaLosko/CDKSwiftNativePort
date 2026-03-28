import Foundation

/// Reader for MDL RGfiles carrying advanced R-group logic (`M  LOG`) and definition blocks.
public enum CDKRGFileReader {
    public static func read(text: String) throws -> CDKRGroupQuery {
        let normalized = text
            .replacingOccurrences(of: "\r\n", with: "\n")
            .replacingOccurrences(of: "\r", with: "\n")
        let lines = normalized.components(separatedBy: "\n")
        return try read(lines: lines)
    }

    public static func readFlattenedMolecule(text: String) throws -> Molecule {
        CDKRGroupQueryManipulator.toFlatMolecule(try read(text: text))
    }

    static func read(lines: [String]) throws -> CDKRGroupQuery {
        var index = 0

        func advance(_ expectedPrefix: String) throws {
            guard index < lines.count else {
                throw ChemError.parseFailed("Unexpected end of RGfile while reading \(expectedPrefix).")
            }
            guard lines[index].trimmingCharacters(in: .whitespaces).hasPrefix(expectedPrefix) else {
                throw ChemError.parseFailed("Expected \(expectedPrefix) while reading RGfile.")
            }
            index += 1
        }

        try advance("$MDL")
        try advance("$MOL")
        try advance("$HDR")

        while index < lines.count,
              lines[index].trimmingCharacters(in: .whitespaces) != "$END HDR" {
            index += 1
        }
        try advance("$END HDR")
        try advance("$CTAB")

        let rootCTAB = try collectBlock(lines: lines, index: &index, endMarker: "$END CTAB")
        var root = try parseMolfileCTAB(rootCTAB, title: "Root structure")
        root.rGroupLogicDefinitions = parseLogicDefinitions(from: rootCTAB)
        applyAttachmentOrderLines(from: rootCTAB, to: &root)

        let rootLabels = Set(root.atoms.compactMap(\.rGroupLabel)).sorted()
        var definitions: [Int: CDKRGroupList] = [:]
        for label in rootLabels {
            let logic = root.rGroupLogicDefinitions[label] ?? MoleculeRGroupLogic()
            definitions[label] = try CDKRGroupList(rGroupNumber: label,
                                                   restH: logic.restH,
                                                   occurrence: logic.occurrence,
                                                   requiredRGroupNumber: logic.requiredRGroupNumber)
        }

        while index < lines.count {
            let line = lines[index].trimmingCharacters(in: .whitespaces)
            if line.isEmpty {
                index += 1
                continue
            }
            if line == "$END MOL" {
                break
            }
            guard line == "$RGP" else {
                throw ChemError.parseFailed("Unexpected RGfile block marker '\(line)'.")
            }
            index += 1
            guard index < lines.count, let rGroupNumber = Int(lines[index].trimmingCharacters(in: .whitespaces)) else {
                throw ChemError.parseFailed("Missing RGfile R-group number after $RGP.")
            }
            index += 1

            guard definitions[rGroupNumber] != nil else {
                throw ChemError.parseFailed("R\(rGroupNumber) is defined in $RGP but missing from root structure.")
            }

            var groups = definitions[rGroupNumber]?.rGroups ?? []
            while index < lines.count {
                let marker = lines[index].trimmingCharacters(in: .whitespaces)
                if marker == "$END RGP" {
                    index += 1
                    break
                }
                guard marker == "$CTAB" else {
                    throw ChemError.parseFailed("Expected $CTAB inside $RGP block.")
                }
                index += 1
                let groupCTAB = try collectBlock(lines: lines, index: &index, endMarker: "$END CTAB")
                let group = try parseMolfileCTAB(groupCTAB, title: "R\(rGroupNumber)")
                groups.append(CDKRGroup(group: group,
                                        firstAttachmentPointAtomID: attachmentPointAtomID(in: group, values: [1, 3]),
                                        secondAttachmentPointAtomID: attachmentPointAtomID(in: group, values: [2, 3])))
            }

            guard var definition = definitions[rGroupNumber] else {
                throw ChemError.parseFailed("Missing RGfile definition for R\(rGroupNumber).")
            }
            definition.rGroups = groups
            definitions[rGroupNumber] = definition
        }

        return CDKRGroupQuery(rootStructure: root, rGroupDefinitions: definitions)
    }

    private static func collectBlock(lines: [String], index: inout Int, endMarker: String) throws -> [String] {
        let start = index
        while index < lines.count,
              lines[index].trimmingCharacters(in: .whitespaces) != endMarker {
            index += 1
        }
        guard index < lines.count else {
            throw ChemError.parseFailed("Unterminated RGfile block, expected \(endMarker).")
        }
        let block = Array(lines[start..<index])
        index += 1
        return block
    }

    private static func parseMolfileCTAB(_ ctabLines: [String], title: String) throws -> Molecule {
        let molfileLines = [title, "", ""] + ctabLines
        return try CDKMDLV2000Reader.read(lines: molfileLines)
    }

    private static func parseLogicDefinitions(from ctabLines: [String]) -> [Int: MoleculeRGroupLogic] {
        var definitions: [Int: MoleculeRGroupLogic] = [:]

        for line in ctabLines {
            let trimmed = line.trimmingCharacters(in: .whitespaces)
            guard trimmed.hasPrefix("M  LOG") else { continue }
            let fields = trimmed.split(whereSeparator: \.isWhitespace).map(String.init)
            guard fields.count >= 6,
                  let rGroupNumber = Int(fields[3]),
                  let required = Int(fields[4]),
                  let restHValue = Int(fields[5]) else {
                continue
            }
            let occurrence = fields.count > 6 ? fields.dropFirst(6).joined() : CDKRGroupList.defaultOccurrence
            definitions[rGroupNumber] = MoleculeRGroupLogic(occurrence: occurrence,
                                                            requiredRGroupNumber: required,
                                                            restH: restHValue == 1)
        }

        return definitions
    }

    private static func applyAttachmentOrderLines(from ctabLines: [String], to molecule: inout Molecule) {
        for line in ctabLines {
            let trimmed = line.trimmingCharacters(in: .whitespaces)
            guard trimmed.hasPrefix("M  AAL") else { continue }
            let fields = trimmed.split(whereSeparator: \.isWhitespace).map(String.init)
            guard fields.count >= 4,
                  let rootAtomID = Int(fields[2]),
                  let pairCount = Int(fields[3]),
                  let atomIndex = molecule.indexOfAtom(id: rootAtomID) else {
                continue
            }

            var orderedPartners: [(Int, Int)] = []
            for pairIndex in 0..<pairCount {
                let base = 4 + pairIndex * 2
                guard base + 1 < fields.count,
                      let partnerID = Int(fields[base]),
                      let order = Int(fields[base + 1]) else {
                    continue
                }
                orderedPartners.append((order, partnerID))
            }

            molecule.atoms[atomIndex].ligandOrderingAtomIDs = orderedPartners
                .sorted(by: { $0.0 < $1.0 })
                .map(\.1)
        }
    }

    private static func attachmentPointAtomID(in molecule: Molecule, values: Set<Int>) -> Int? {
        molecule.atoms.first(where: { atom in
            guard let attachmentPoint = atom.attachmentPoint else { return false }
            return values.contains(attachmentPoint)
        })?.id
    }
}
