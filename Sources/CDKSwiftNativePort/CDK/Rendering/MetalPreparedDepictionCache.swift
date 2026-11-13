#if canImport(CoreGraphics)
    import CoreGraphics
#endif
import Foundation

// The lock owns both entries and LRU order; expensive preparation stays outside it.
private struct MetalPreparedDepictionKey: Hashable {
    let molecule: Molecule
    let style: RenderStyle
}

struct MetalPreparedDepictionData {
    let depictionMolecule: Molecule
    let boundingBox: CGRect?
    let degreeByAtomID: [Int: Int]
    let aromaticRings: [[Int]]
    let aromaticBondIDs: Set<Int>
    let conjugatedRings: [[Int]]
}

enum MetalPreparedDepictionCache {
    private static let maxEntries = 16
    private static let lock = NSLock()
    nonisolated(unsafe) private static var entries: [MetalPreparedDepictionKey: MetalPreparedDepictionData] = [:]
    nonisolated(unsafe) private static var lruOrder: [MetalPreparedDepictionKey] = []

    static func preparedData(molecule: Molecule, style: RenderStyle) -> MetalPreparedDepictionData {
        let key = MetalPreparedDepictionKey(molecule: molecule, style: style)

        lock.lock()
        if let cached = entries[key] {
            touchLRU(key)
            lock.unlock()
            return cached
        }
        lock.unlock()

        let depictionMolecule = CDKDepictionPreprocessor.prepareForRendering(molecule: molecule, style: style)
        let aromaticRings = depictionMolecule.aromaticDisplayRings()
        let aromaticBondIDs = depictionMolecule.aromaticDisplayBondIDs()
        let conjugatedRings = depictionMolecule.simpleCycles(maxSize: 10).filter { ring in
            guard ring.count >= 5 else { return false }
            let ringBonds: [Bond] = (0..<ring.count).compactMap { idx in
                depictionMolecule.bond(between: ring[idx], and: ring[(idx + 1) % ring.count])
            }
            guard ringBonds.count == ring.count else { return false }
            let piCount = ringBonds.reduce(0) { partial, bond in
                partial + ((bond.order == .double || bond.order == .aromatic) ? 1 : 0)
            }
            return piCount >= 2 && ringBonds.allSatisfy { $0.order != .triple }
        }

        var degreeByAtomID: [Int: Int] = [:]
        for bond in depictionMolecule.bonds {
            degreeByAtomID[bond.a1, default: 0] += 1
            degreeByAtomID[bond.a2, default: 0] += 1
        }

        let prepared = MetalPreparedDepictionData(
            depictionMolecule: depictionMolecule,
            boundingBox: depictionMolecule.boundingBox(),
            degreeByAtomID: degreeByAtomID,
            aromaticRings: aromaticRings,
            aromaticBondIDs: aromaticBondIDs,
            conjugatedRings: conjugatedRings)

        lock.lock()
        if let cached = entries[key] {
            touchLRU(key)
            lock.unlock()
            return cached
        }
        entries[key] = prepared
        lruOrder.append(key)
        while lruOrder.count > maxEntries {
            let stale = lruOrder.removeFirst()
            entries.removeValue(forKey: stale)
        }
        lock.unlock()
        return prepared
    }

    private static func touchLRU(_ key: MetalPreparedDepictionKey) {
        guard let index = lruOrder.firstIndex(of: key) else { return }
        if index == lruOrder.count - 1 { return }
        let moved = lruOrder.remove(at: index)
        lruOrder.append(moved)
    }
}
