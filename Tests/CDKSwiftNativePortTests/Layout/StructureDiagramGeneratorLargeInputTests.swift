import Foundation
import XCTest

@testable import CDKSwiftNativePort

#if canImport(CoreGraphics)
    import CoreGraphics
#endif

/// Layout traversals must not recurse once per atom. Callers such as file previews
/// run layout on secondary threads, whose default stack is 512 KB; these tests run
/// the traversals on a 64 KB stack so a recursive form fails at test size.
final class StructureDiagramGeneratorLargeInputTests: XCTestCase {
    func testLongestChainSearchHandlesALongChainOnASmallStack() {
        let molecule = Self.chain(atomCount: 5_000)
        let length = runOnSmallStack {
            let graph = CDKMolecularGraph(molecule: molecule)
            return CDKStructureDiagramGenerator.longestUnplacedChain(
                anchor: 1, start: 2, component: Set(molecule.atoms.map(\.id)), graph: graph,
                molecule: molecule, ringAtoms: [], ringEdges: [], positions: [1: .zero]
            ).count
        }
        XCTAssertEqual(length, 5_000)
    }

    func testLongestChainSearchKeepsTheRecursiveTieBreak() {
        // A star: from the anchor, three branches of lengths 3, 3 and 2. The longest
        // path ends in the lexicographically smallest branch of the longest length.
        let bonds = [(1, 2), (2, 3), (3, 4), (2, 5), (5, 6), (2, 7)]
        let molecule = Self.molecule(atomCount: 7, bonds: bonds)
        let graph = CDKMolecularGraph(molecule: molecule)
        let path = CDKStructureDiagramGenerator.longestUnplacedChain(
            anchor: 1, start: 2, component: Set(1...7), graph: graph, molecule: molecule,
            ringAtoms: [], ringEdges: [], positions: [1: .zero])
        XCTAssertEqual(path, [1, 2, 3, 4])
    }

    func testRingSearchHandlesALongRingChainOnASmallStack() {
        // 1,500 cyclopropanes joined in a line: 4,500 atoms in one component.
        var bonds: [(Int, Int)] = []
        for ring in 0..<1_500 {
            let a = ring * 3 + 1
            bonds += [(a, a + 1), (a + 1, a + 2), (a + 2, a)]
            if ring > 0 { bonds.append((a - 1, a)) }
        }
        let molecule = Self.molecule(atomCount: 4_500, bonds: bonds)
        let rings = runOnSmallStack {
            CDKRingSearch.sssrLikeRings(in: molecule, graph: CDKMolecularGraph(molecule: molecule)).count
        }
        XCTAssertEqual(rings, 1_500)
    }

    func testRefinementLimitsBoundTheCostOfLargeComponents() {
        XCTAssertEqual(CDKStructureDiagramGenerator.Tuning.bondFlipComponentLimit, 256)
        XCTAssertEqual(
            CDKStructureDiagramGenerator.Tuning.relaxPairVisitBudget,
            CDKStructureDiagramGenerator.Tuning.relaxIterations * (512 * 511 / 2))
    }

    private static func chain(atomCount: Int) -> Molecule {
        molecule(atomCount: atomCount, bonds: (1..<atomCount).map { ($0, $0 + 1) })
    }

    private static func molecule(atomCount: Int, bonds: [(Int, Int)]) -> Molecule {
        Molecule(
            atoms: (1...atomCount).map { Atom(id: $0, element: "C", position: .zero) },
            bonds: bonds.enumerated().map {
                Bond(id: $0.offset + 1, a1: $0.element.0, a2: $0.element.1, order: .single)
            })
    }

    private final class ResultBox: @unchecked Sendable {
        var value = 0
    }

    /// The semaphore orders the worker's write before the caller's read.
    private func runOnSmallStack(_ work: @escaping @Sendable () -> Int) -> Int {
        let box = ResultBox()
        let finished = DispatchSemaphore(value: 0)
        let thread = Thread {
            box.value = work()
            finished.signal()
        }
        thread.stackSize = 64 * 1024
        thread.start()
        finished.wait()
        return box.value
    }
}
