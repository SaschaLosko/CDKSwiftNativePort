import Foundation

enum CDKInChIKeyCodec {
    private struct Components {
        let major: String
        let minor: String
        let standardFlag: Character
        let protonFlag: Character
    }

    private enum KeyError: LocalizedError {
        case invalidPrefix
        case invalidVersion
        case invalidInchi
        case invalidStandardInchi
        case invalidProtonLayer

        var errorDescription: String? {
            switch self {
            case .invalidPrefix:
                return "Input is not a valid InChI string for InChIKey generation."
            case .invalidVersion:
                return "Only InChI version 1.x is supported for InChIKey generation."
            case .invalidInchi:
                return "Input InChI is not valid enough to derive an InChIKey."
            case .invalidStandardInchi:
                return "Standard InChI contains non-standard layers that are not valid for InChIKey generation."
            case .invalidProtonLayer:
                return "Unsupported protonation layer while deriving InChIKey."
            }
        }
    }

    static func makeKey(from inchi: String) throws -> String {
        let components = try split(inchi)
        let majorDigest = CDKSHA256.hash(bytes: Array(components.major.utf8))

        var minorBytes = Array(components.minor.utf8)
        if !minorBytes.isEmpty, minorBytes.count < 255 {
            minorBytes += minorBytes
        }
        let minorDigest = CDKSHA256.hash(bytes: minorBytes)

        let firstBlock =
            tripletString(for: triplet1(majorDigest)) +
            tripletString(for: triplet2(majorDigest)) +
            tripletString(for: triplet3(majorDigest)) +
            tripletString(for: triplet4(majorDigest)) +
            dubletString(for: dublet56to64(majorDigest))

        let secondBlock =
            tripletString(for: triplet1(minorDigest)) +
            tripletString(for: triplet2(minorDigest)) +
            dubletString(for: dublet28to36(minorDigest)) +
            String(components.standardFlag) +
            "A"

        return "\(firstBlock)-\(secondBlock)-\(components.protonFlag)"
    }

    private static func split(_ inchi: String) throws -> Components {
        let trimmed = inchi.trimmingCharacters(in: .whitespacesAndNewlines)
        let bytes = Array(trimmed.utf8)
        let prefix = Array("InChI=".utf8)

        guard bytes.count >= prefix.count + 3, bytes.starts(with: prefix) else {
            throw KeyError.invalidPrefix
        }

        guard bytes[prefix.count] == UInt8(ascii: "1") else {
            throw KeyError.invalidVersion
        }

        var slashIndex = prefix.count + 1
        var standardFlag: Character = "N"

        if slashIndex < bytes.count {
            if bytes[slashIndex] == UInt8(ascii: "S") {
                standardFlag = "S"
                slashIndex += 1
            } else if bytes[slashIndex] == UInt8(ascii: "B") {
                standardFlag = "B"
                slashIndex += 1
            }
        }

        guard slashIndex < bytes.count, bytes[slashIndex] == UInt8(ascii: "/") else {
            throw KeyError.invalidPrefix
        }

        guard slashIndex + 1 < bytes.count, isAllowedFirstPayloadByte(bytes[slashIndex + 1]) else {
            throw KeyError.invalidInchi
        }

        var protonStart: Int?
        var minorStart: Int?
        var scan = slashIndex + 1

        while scan < bytes.count - 1 {
            if bytes[scan] == UInt8(ascii: "/") {
                let layerKey = bytes[scan + 1]
                switch layerKey {
                case UInt8(ascii: "c"), UInt8(ascii: "h"), UInt8(ascii: "q"):
                    scan += 1
                    continue
                case UInt8(ascii: "p"):
                    protonStart = scan
                    scan += 1
                    continue
                case UInt8(ascii: "f") where standardFlag == "S":
                    throw KeyError.invalidStandardInchi
                case UInt8(ascii: "r") where standardFlag == "S":
                    throw KeyError.invalidStandardInchi
                default:
                    minorStart = scan
                    scan = bytes.count
                    continue
                }
            }
            scan += 1
        }

        let majorStart = slashIndex + 1
        let majorEnd = protonStart ?? minorStart ?? bytes.count
        let major = String(decoding: bytes[majorStart..<majorEnd], as: UTF8.self)

        let protonFlag: Character
        if let protonStart {
            let protonEnd = minorStart ?? bytes.count
            let protonLayer = String(decoding: bytes[protonStart..<protonEnd], as: UTF8.self)
            guard protonLayer.hasPrefix("/p"), protonLayer.count >= 3 else {
                throw KeyError.invalidProtonLayer
            }
            guard let protonCount = Int(protonLayer.dropFirst(2)), protonCount != 0 else {
                throw KeyError.invalidProtonLayer
            }
            protonFlag = protonationFlag(for: protonCount)
        } else {
            protonFlag = "N"
        }

        let minor = minorStart.map { String(decoding: bytes[$0..<bytes.count], as: UTF8.self) } ?? ""
        return Components(major: major, minor: minor, standardFlag: standardFlag, protonFlag: protonFlag)
    }

    private static func isAllowedFirstPayloadByte(_ byte: UInt8) -> Bool {
        isASCIIAlphaNumeric(byte) || byte == UInt8(ascii: "/") || byte == UInt8(ascii: "?")
    }

    private static func isASCIIAlphaNumeric(_ byte: UInt8) -> Bool {
        (byte >= UInt8(ascii: "0") && byte <= UInt8(ascii: "9")) ||
            (byte >= UInt8(ascii: "A") && byte <= UInt8(ascii: "Z")) ||
            (byte >= UInt8(ascii: "a") && byte <= UInt8(ascii: "z"))
    }

    private static func protonationFlag(for protonCount: Int) -> Character {
        if protonCount > 12 || protonCount < -12 {
            return "A"
        }

        if protonCount > 0 {
            let plusFlags = Array("OPQRSTUVWXYZ")
            return plusFlags[protonCount - 1]
        }

        let minusFlags = Array("MLKJIHGFEDCB")
        return minusFlags[(-protonCount) - 1]
    }

    private static func triplet1(_ digest: [UInt8]) -> Int {
        Int(digest[0]) | (Int(digest[1] & 0x3f) << 8)
    }

    private static func triplet2(_ digest: [UInt8]) -> Int {
        (Int(digest[1] & 0xc0) | (Int(digest[2]) << 8) | (Int(digest[3] & 0x0f) << 16)) >> 6
    }

    private static func triplet3(_ digest: [UInt8]) -> Int {
        (Int(digest[3] & 0xf0) | (Int(digest[4]) << 8) | (Int(digest[5] & 0x03) << 16)) >> 4
    }

    private static func triplet4(_ digest: [UInt8]) -> Int {
        (Int(digest[5] & 0xfc) | (Int(digest[6]) << 8)) >> 2
    }

    private static func dublet28to36(_ digest: [UInt8]) -> Int {
        (Int(digest[3] & 0xf0) | (Int(digest[4] & 0x1f) << 8)) >> 4
    }

    private static func dublet56to64(_ digest: [UInt8]) -> Int {
        Int(digest[7]) | (Int(digest[8] & 0x01) << 8)
    }

    private static func tripletString(for index: Int) -> String {
        precondition((0..<16384).contains(index))

        let lexicalIndex: Int
        if index >= 12_168 {
            lexicalIndex = index + 1_192
        } else if index >= 2_704 {
            lexicalIndex = index + 676
        } else {
            lexicalIndex = index
        }

        return base26String(for: lexicalIndex, length: 3)
    }

    private static func dubletString(for index: Int) -> String {
        precondition((0..<676).contains(index))
        return base26String(for: index, length: 2)
    }

    private static func base26String(for index: Int, length: Int) -> String {
        var value = index
        var bytes = Array(repeating: UInt8(ascii: "A"), count: length)

        for position in stride(from: length - 1, through: 0, by: -1) {
            bytes[position] = UInt8(ascii: "A") + UInt8(value % 26)
            value /= 26
        }

        return String(decoding: bytes, as: UTF8.self)
    }
}

private enum CDKSHA256 {
    private static let initialState: [UInt32] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
        0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
    ]

    private static let roundConstants: [UInt32] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
        0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
        0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
        0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
        0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
        0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
        0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
        0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
    ]

    static func hash(bytes: [UInt8]) -> [UInt8] {
        var message = bytes
        let bitLength = UInt64(message.count) * 8

        message.append(0x80)
        while message.count % 64 != 56 {
            message.append(0)
        }

        message.append(UInt8((bitLength >> 56) & 0xff))
        message.append(UInt8((bitLength >> 48) & 0xff))
        message.append(UInt8((bitLength >> 40) & 0xff))
        message.append(UInt8((bitLength >> 32) & 0xff))
        message.append(UInt8((bitLength >> 24) & 0xff))
        message.append(UInt8((bitLength >> 16) & 0xff))
        message.append(UInt8((bitLength >> 8) & 0xff))
        message.append(UInt8(bitLength & 0xff))

        var state = initialState
        var schedule = Array(repeating: UInt32(0), count: 64)

        for chunkStart in stride(from: 0, to: message.count, by: 64) {
            for wordIndex in 0..<16 {
                let offset = chunkStart + wordIndex * 4
                schedule[wordIndex] =
                    (UInt32(message[offset]) << 24) |
                    (UInt32(message[offset + 1]) << 16) |
                    (UInt32(message[offset + 2]) << 8) |
                    UInt32(message[offset + 3])
            }

            for wordIndex in 16..<64 {
                let s0 = rotateRight(schedule[wordIndex - 15], by: 7) ^
                    rotateRight(schedule[wordIndex - 15], by: 18) ^
                    (schedule[wordIndex - 15] >> 3)
                let s1 = rotateRight(schedule[wordIndex - 2], by: 17) ^
                    rotateRight(schedule[wordIndex - 2], by: 19) ^
                    (schedule[wordIndex - 2] >> 10)
                schedule[wordIndex] = schedule[wordIndex - 16] &+ s0 &+ schedule[wordIndex - 7] &+ s1
            }

            var a = state[0]
            var b = state[1]
            var c = state[2]
            var d = state[3]
            var e = state[4]
            var f = state[5]
            var g = state[6]
            var h = state[7]

            for index in 0..<64 {
                let sum1 = rotateRight(e, by: 6) ^ rotateRight(e, by: 11) ^ rotateRight(e, by: 25)
                let choose = (e & f) ^ ((~e) & g)
                let temp1 = h &+ sum1 &+ choose &+ roundConstants[index] &+ schedule[index]
                let sum0 = rotateRight(a, by: 2) ^ rotateRight(a, by: 13) ^ rotateRight(a, by: 22)
                let majority = (a & b) ^ (a & c) ^ (b & c)
                let temp2 = sum0 &+ majority

                h = g
                g = f
                f = e
                e = d &+ temp1
                d = c
                c = b
                b = a
                a = temp1 &+ temp2
            }

            state[0] &+= a
            state[1] &+= b
            state[2] &+= c
            state[3] &+= d
            state[4] &+= e
            state[5] &+= f
            state[6] &+= g
            state[7] &+= h
        }

        var digest: [UInt8] = []
        digest.reserveCapacity(32)
        for word in state {
            digest.append(UInt8((word >> 24) & 0xff))
            digest.append(UInt8((word >> 16) & 0xff))
            digest.append(UInt8((word >> 8) & 0xff))
            digest.append(UInt8(word & 0xff))
        }
        return digest
    }

    private static func rotateRight(_ value: UInt32, by amount: UInt32) -> UInt32 {
        (value >> amount) | (value << (32 - amount))
    }
}
