import Foundation

// CTfile reference mass numbers adapted from CDK 2.14 MDLIsotope.
// Copyright (C) 2026 John Mayfield; LGPL-2.1-or-later.
// https://github.com/cdk/cdk/pull/1350
// These are format offsets, not exact masses or necessarily major isotopes.
enum CDKMDLIsotope {
    static func referenceMass(for symbol: String) -> Int? {
        let normalized = symbol == "R" || symbol == "R#" ? "*" : CDKDescriptorSupport.canonicalElementSymbol(symbol)
        guard let element = QueryElements.number(normalized), reference.indices.contains(element) else { return nil }
        return reference[element]
    }

    private static let reference: [Int] = [
        0, 1, 4, 7, 9, 11, 12, 14, 16, 19, 20, 23, 24, 27, 28,
        31, 32, 35, 40, 39, 40, 45, 48, 51, 52, 55, 56, 59, 59, 64,
        65, 70, 73, 75, 79, 80, 84, 85, 88, 89, 91, 93, 96, 98, 101,
        103, 106, 108, 112, 115, 119, 122, 128, 127, 131, 133, 137, 139, 140, 141,
        144, 145, 150, 152, 157, 159, 162, 165, 167, 169, 173, 175, 178, 181, 184,
        186, 190, 192, 195, 197, 201, 204, 207, 209, 209, 210, 222, 223, 226, 227,
        232, 231, 238, 237, 244, 243, 247, 247, 251, 252, 257, 258, 259, 262, 261,
        262, 266, 264, 277, 268, 280, 280, 285, 286, 289, 289, 293, 294, 294,
    ]
}
