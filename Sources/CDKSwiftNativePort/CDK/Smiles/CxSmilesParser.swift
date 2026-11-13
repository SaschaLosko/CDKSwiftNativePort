import Foundation

public enum CDKCxSmilesParser {
    public struct SplitResult: Equatable {
        public let coreSmiles: String
        public let title: String?
        public let state: CDKCxSmilesState

        public init(coreSmiles: String, title: String?, state: CDKCxSmilesState) {
            self.coreSmiles = coreSmiles
            self.title = title
            self.state = state
        }
    }

    public static func split(_ input: String, enabled: Bool) throws -> SplitResult {
        let trimmed = input.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { throw ChemError.emptyInput }

        guard enabled else {
            return splitTitleField(from: trimmed)
        }

        let basic = splitTitleField(from: trimmed)
        guard let title = basic.title, title.hasPrefix("|") else {
            return basic
        }

        guard
            let secondPipe = CDKCxSmilesLayerDecoder.findClosingCxPipe(
                in: title, startingAt: title.index(after: title.startIndex))
        else {
            return basic
        }

        let cxBody = String(title[title.index(after: title.startIndex)..<secondPipe])
        guard let state = try? CDKCxSmilesLayerDecoder.decode(cxBody) else {
            return basic
        }

        let trailingStart = title.index(after: secondPipe)
        let trailing = String(title[trailingStart...]).trimmingCharacters(
            in: .whitespacesAndNewlines)
        return SplitResult(
            coreSmiles: basic.coreSmiles,
            title: trailing,
            state: state)
    }

    private static func splitTitleField(from input: String) -> SplitResult {
        let separators = CharacterSet.whitespacesAndNewlines
        guard let firstSeparatorRange = input.rangeOfCharacter(from: separators) else {
            return SplitResult(coreSmiles: input, title: nil, state: CDKCxSmilesState())
        }

        let core = String(input[..<firstSeparatorRange.lowerBound]).trimmingCharacters(
            in: separators)
        let trailing = String(input[firstSeparatorRange.upperBound...]).trimmingCharacters(
            in: separators)
        return SplitResult(
            coreSmiles: core,
            title: trailing.isEmpty ? nil : trailing,
            state: CDKCxSmilesState())
    }

    static func encodeStereoGroup(kind: String, number: Int) -> Int {
        let kindCode: Int
        switch kind {
        case "abs":
            kindCode = 1
        case "or":
            kindCode = 2
        case "and":
            kindCode = 3
        default:
            kindCode = 0
        }
        return (kindCode << 16) | max(0, number)
    }

    static func decodeStereoGroup(_ encoded: Int) -> (kind: String, number: Int)? {
        let kindCode = encoded >> 16
        let number = encoded & 0xFFFF
        switch kindCode {
        case 1:
            return ("abs", 0)
        case 2:
            return ("or", number)
        case 3:
            return ("and", number)
        default:
            return nil
        }
    }

}
