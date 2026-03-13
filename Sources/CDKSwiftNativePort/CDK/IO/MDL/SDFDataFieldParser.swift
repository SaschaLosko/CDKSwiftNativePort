import Foundation

enum CDKSDFDataFieldParser {

    static func applyParsedFields(from lines: [String], to molecule: inout Molecule) {
        let (dataFields, fieldOrder) = parse(lines: lines)
        guard !fieldOrder.isEmpty else { return }

        for fieldName in fieldOrder {
            molecule.ensureDataField(named: fieldName)
            for value in dataFields[fieldName] ?? [] {
                molecule.appendDataFieldValue(value, named: fieldName)
            }
        }
    }

    static func parse(lines: [String]) -> (dataFields: [String: [String]], fieldOrder: [String]) {
        var dataFields: [String: [String]] = [:]
        var fieldOrder: [String] = []
        var didReachStructureEnd = false
        var index = 0

        while index < lines.count {
            let trimmed = lines[index].trimmingCharacters(in: .whitespacesAndNewlines)

            if trimmed == "$$$$" {
                break
            }

            if !didReachStructureEnd {
                if trimmed == "M  END" || trimmed == "M END" {
                    didReachStructureEnd = true
                }
                index += 1
                continue
            }

            guard let fieldName = fieldName(from: lines[index]) else {
                index += 1
                continue
            }

            if dataFields[fieldName] == nil {
                dataFields[fieldName] = []
                fieldOrder.append(fieldName)
            }

            index += 1
            while index < lines.count {
                let valueLine = lines[index]
                let trimmedValue = valueLine.trimmingCharacters(in: .whitespacesAndNewlines)
                if trimmedValue == "$$$$" {
                    break
                }
                if trimmedValue.isEmpty {
                    index += 1
                    break
                }
                dataFields[fieldName, default: []].append(valueLine)
                index += 1
            }
        }

        return (dataFields, fieldOrder)
    }

    private static func fieldName(from headerLine: String) -> String? {
        let trimmed = headerLine.trimmingCharacters(in: .whitespacesAndNewlines)
        guard trimmed.hasPrefix(">") else { return nil }

        if let openAngle = trimmed.firstIndex(of: "<"),
           let closeAngle = trimmed[trimmed.index(after: openAngle)...].firstIndex(of: ">") {
            let fieldName = String(trimmed[trimmed.index(after: openAngle)..<closeAngle])
                .trimmingCharacters(in: .whitespacesAndNewlines)
            return fieldName.isEmpty ? nil : fieldName
        }

        let raw = trimmed.drop(while: { $0 == ">" || $0.isWhitespace })
        let fieldName = String(raw)
            .split(separator: "(", maxSplits: 1, omittingEmptySubsequences: true)
            .first
            .map { String($0).trimmingCharacters(in: .whitespacesAndNewlines) }

        guard let fieldName, !fieldName.isEmpty else { return nil }
        return fieldName
    }
}
