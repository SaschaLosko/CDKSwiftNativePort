import Foundation

// Derived from CDK 2.13 SquarePlanar/TrigonalBipyramidal/Octahedral conversions,
// Copyright 2017 John Mayfield. SPDX-License-Identifier: LGPL-2.1-or-later
extension QueryStereoSupport {
    static func coordinationEquivalent(
        _ query: CDKQueryAtomStereo, _ target: CDKQueryAtomStereo,
        mapping: [Int: Int]
    ) -> Bool {
        func octahedral(_ element: CDKQueryAtomStereo) -> CDKQueryAtomStereo {
            guard element.geometry == .squarePlanar else { return element }
            var converted = element
            converted.geometry = .octahedral
            converted.configuration = [1, 8, 4][element.configuration - 1]
            converted.carriers = [element.focus] + element.carriers + [element.focus]
            return converted
        }
        let expected = octahedral(query)
        var actual = octahedral(target)
        let mapped = Set(expected.carriers.compactMap { mapping[$0] })
        // Unspecified query ligands are wildcards at this focus, not requests to
        // remove target atoms. CDK filters target carriers to mapped query ligands.
        actual.carriers = actual.carriers.map { mapped.contains($0) ? $0 : actual.focus }
        if actual.geometry == expected.geometry { return equivalent(expected, actual, mapping: mapping) }
        guard [.octahedral, .trigonalBipyramidal].contains(expected.geometry),
            [.octahedral, .trigonalBipyramidal].contains(actual.geometry),
            actual.carriers.filter({ $0 != actual.focus }).count <= 3
        else { return false }
        let normalized = representatives(actual.geometry)[actual.configuration - 1].map { actual.carriers[$0] }
        for rotation in rotations(actual.geometry) {
            var carriers = rotation.map { normalized[$0] }
            guard carriers[1..<(carriers.count - 1)].filter({ $0 != actual.focus }).count <= 1 else { continue }
            if actual.geometry == .trigonalBipyramidal {
                carriers.insert(actual.focus, at: 1)
            } else if let implicit = [1, 2].first(where: { carriers[$0] == actual.focus }) {
                carriers.remove(at: implicit)
            } else {
                continue
            }
            let converted = CDKQueryAtomStereo(
                focus: actual.focus, carriers: carriers,
                geometry: expected.geometry, configuration: 1)
            if equivalent(expected, converted, mapping: mapping) { return true }
        }
        return false
    }
}
