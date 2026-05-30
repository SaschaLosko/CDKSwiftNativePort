import Foundation
import XCTest
#if canImport(CoreGraphics)
import CoreGraphics
#endif
@testable import CDKSwiftNativePort

final class StructureDiagramGeneratorTests: XCTestCase {
    private let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    // Regression for "C-O-C drawn as a straight line" complaints.
    func testAspirinEsterBridgeIsNotLinear() throws {
        let molecule = try parse("CC(=O)OC1=CC=CC=C1C(=O)O")

        let oxygen = try XCTUnwrap(molecule.atoms.first { atom in
            atom.element.uppercased() == "O"
                && molecule.neighbors(of: atom.id).count == 2
                && molecule.neighbors(of: atom.id).allSatisfy { neighbor in
                    molecule.atom(id: neighbor)?.element.uppercased() == "C"
                }
        })

        let neighbors = molecule.neighbors(of: oxygen.id)
        XCTAssertEqual(neighbors.count, 2)

        let p0 = try XCTUnwrap(molecule.atom(id: oxygen.id)?.position)
        let p1 = try XCTUnwrap(molecule.atom(id: neighbors[0])?.position)
        let p2 = try XCTUnwrap(molecule.atom(id: neighbors[1])?.position)

        let angle = angleDegrees(a: p1, center: p0, b: p2)
        XCTAssertGreaterThan(angle, 95)
        XCTAssertLessThan(angle, 170)
    }

    func testNaphthaleneAvoidsSevereNonBondedOverlap() throws {
        let molecule = try parse("c1cccc2ccccc12")

        let atomIDs = molecule.atoms.map(\.id)
        var minDistance = CGFloat.greatestFiniteMagnitude

        for i in 0..<atomIDs.count {
            for j in (i + 1)..<atomIDs.count {
                let a = atomIDs[i]
                let b = atomIDs[j]
                if molecule.bond(between: a, and: b) != nil { continue }

                let p1 = try XCTUnwrap(molecule.atom(id: a)?.position)
                let p2 = try XCTUnwrap(molecule.atom(id: b)?.position)
                minDistance = min(minDistance, p1.distance(to: p2))
            }
        }

        XCTAssertGreaterThan(minDistance, 0.45)
    }

    func testNorbornaneLayoutHasArea() throws {
        let molecule = try parse("C1CC2CCC1C2")
        let box = try XCTUnwrap(molecule.boundingBox())

        XCTAssertGreaterThan(box.width, 1.0)
        XCTAssertGreaterThan(box.height, 1.0)
    }

    func testCyclohexaneLayoutUsesCanonicalRegularRingOrientation() throws {
        let molecule = try parse("C1CCCCC1")

        let nearVerticalBondCount = molecule.bonds.reduce(into: 0) { count, bond in
            guard let p1 = molecule.atom(id: bond.a1)?.position,
                  let p2 = molecule.atom(id: bond.a2)?.position else {
                return
            }
            let dx = abs(p2.x - p1.x)
            let dy = abs(p2.y - p1.y)
            if dy > 0.5, dx < 0.15 {
                count += 1
            }
        }

        XCTAssertGreaterThanOrEqual(nearVerticalBondCount, 2,
                                    "Expected cyclohexane to seed with the conventional vertical-side hexagon orientation.")
    }

    func testExplicitHydrogensDoNotDistortIbuprofenLikeSideChainAngles() throws {
        let skeletal = try parse("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
        let explicitHydrogenExpanded = expandedWithExplicitHydrogens(from: skeletal)
        let laidOut = Depiction2DGenerator.generate(for: explicitHydrogenExpanded)

        let ringAtomIDs = Set(skeletal.simpleCycles(maxSize: 8).flatMap { $0 })
        let candidateAtomIDs = skeletal.atoms.compactMap { atom -> Int? in
            guard atom.element.uppercased() != "H",
                  !ringAtomIDs.contains(atom.id) else {
                return nil
            }
            let expected = heavyNeighborAngles(in: skeletal, around: atom.id)
            return expected.isEmpty ? nil : atom.id
        }

        XCTAssertFalse(candidateAtomIDs.isEmpty)

        for atomID in candidateAtomIDs {
            let expected = heavyNeighborAngles(in: skeletal, around: atomID)
            let actual = heavyNeighborAngles(in: laidOut, around: atomID)
            XCTAssertEqual(actual.count, expected.count, "Expected the same heavy-neighbor angle count for atom \(atomID).")
            for (expectedAngle, actualAngle) in zip(expected, actual) {
                XCTAssertEqual(actualAngle,
                               expectedAngle,
                               accuracy: 4.0,
                               "Explicit-hydrogen layout distorted the heavy-atom zig-zag at atom \(atomID).")
            }
        }
    }

    func testSingleAnchorRingAttachmentUsesExocyclicBisectorConvention() throws {
        let molecule = try parse("c1ccccc1NCC1OCOC1")
        let fiveMemberRing = try XCTUnwrap(molecule.simpleCycles(maxSize: 8).first { $0.count == 5 })
        let ringAtoms = Set(fiveMemberRing)
        let anchorID = try XCTUnwrap(fiveMemberRing.first { atomID in
            let neighbors = molecule.neighbors(of: atomID)
            let ringNeighborCount = neighbors.filter { ringAtoms.contains($0) }.count
            return ringNeighborCount == 2 && neighbors.contains { !ringAtoms.contains($0) }
        })

        let ringNeighbors = molecule.neighbors(of: anchorID).filter { ringAtoms.contains($0) }
        let externalNeighbor = try XCTUnwrap(molecule.neighbors(of: anchorID).first { !ringAtoms.contains($0) })

        XCTAssertEqual(ringNeighbors.count, 2)

        let center = try XCTUnwrap(molecule.atom(id: anchorID)?.position)
        let firstRingPoint = try XCTUnwrap(molecule.atom(id: ringNeighbors[0])?.position)
        let secondRingPoint = try XCTUnwrap(molecule.atom(id: ringNeighbors[1])?.position)
        let externalPoint = try XCTUnwrap(molecule.atom(id: externalNeighbor)?.position)

        let internalRingAngle = angleDegrees(a: firstRingPoint, center: center, b: secondRingPoint)
        let exocyclicAngles = [
            angleDegrees(a: firstRingPoint, center: center, b: externalPoint),
            angleDegrees(a: secondRingPoint, center: center, b: externalPoint)
        ].sorted()

        XCTAssertEqual(internalRingAngle, 108, accuracy: 8)
        XCTAssertGreaterThan(exocyclicAngles[0], 116)
        XCTAssertLessThan(exocyclicAngles[1], 136)
        XCTAssertEqual(exocyclicAngles[0],
                       exocyclicAngles[1],
                       accuracy: 6,
                       "Single-anchor ring placement should place the external bond on the ring-bond bisector, not collinear with one ring bond.")
    }

    func testTrigonalCarbonateOxygensUseOpenOneHundredTwentyDegreeLayout() throws {
        let molecule = try parse("O=C([O-])[O-]")
        let carbon = try XCTUnwrap(molecule.atoms.first { atom in
            atom.element.uppercased() == "C"
                && molecule.neighbors(of: atom.id).count == 3
                && molecule.neighbors(of: atom.id).allSatisfy { neighbor in
                    molecule.atom(id: neighbor)?.element.uppercased() == "O"
                }
        })

        let center = carbon.position
        let oxygenPoints = molecule.neighbors(of: carbon.id).compactMap { molecule.atom(id: $0)?.position }
        XCTAssertEqual(oxygenPoints.count, 3)

        var angles: [CGFloat] = []
        for leftIndex in 0..<oxygenPoints.count {
            for rightIndex in (leftIndex + 1)..<oxygenPoints.count {
                angles.append(angleDegrees(a: oxygenPoints[leftIndex],
                                           center: center,
                                           b: oxygenPoints[rightIndex]))
            }
        }

        for angle in angles {
            XCTAssertEqual(angle,
                           120,
                           accuracy: 8,
                           "Three-coordinate carbonate-like centers should depict as an open trigonal fan, not a collapsed or near-linear branch.")
        }
    }

    func testSubstitutedAlkeneParentChainKeepsConventionalZigZagLayout() throws {
        let molecule = Molecule(
            name: "5-ethyl-2,4-dimethyloct-2-ene",
            atoms: (1...12).map { Atom(id: $0, element: "C", position: .zero) },
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single),
                Bond(id: 2, a1: 2, a2: 3, order: .double),
                Bond(id: 3, a1: 3, a2: 4, order: .single),
                Bond(id: 4, a1: 4, a2: 5, order: .single),
                Bond(id: 5, a1: 5, a2: 6, order: .single),
                Bond(id: 6, a1: 6, a2: 7, order: .single),
                Bond(id: 7, a1: 7, a2: 8, order: .single),
                Bond(id: 8, a1: 2, a2: 9, order: .single),
                Bond(id: 9, a1: 4, a2: 10, order: .single),
                Bond(id: 10, a1: 5, a2: 11, order: .single),
                Bond(id: 11, a1: 11, a2: 12, order: .single),
            ]
        )

        let laidOut = Depiction2DGenerator.generate(for: molecule)
        let box = try XCTUnwrap(laidOut.boundingBox())
        XCTAssertGreaterThan(box.width,
                             box.height,
                             "A long substituted acyclic parent chain should be oriented as a readable horizontal zig-zag.")

        let parentChain = [1, 2, 3, 4, 5, 6, 7, 8]
        let terminalStart = try XCTUnwrap(laidOut.atom(id: parentChain[0])?.position)
        let terminalEnd = try XCTUnwrap(laidOut.atom(id: parentChain[parentChain.count - 1])?.position)
        XCTAssertGreaterThan(abs(terminalEnd.x - terminalStart.x),
                             abs(terminalEnd.y - terminalStart.y),
                             "The oct-2-ene parent chain should not be laid out vertically after the double bond.")

        for index in 1..<(parentChain.count - 1) {
            let previous = try XCTUnwrap(laidOut.atom(id: parentChain[index - 1])?.position)
            let center = try XCTUnwrap(laidOut.atom(id: parentChain[index])?.position)
            let next = try XCTUnwrap(laidOut.atom(id: parentChain[index + 1])?.position)
            let angle = angleDegrees(a: previous, center: center, b: next)
            XCTAssertGreaterThan(angle, 100, "Unexpectedly compressed parent-chain angle at atom \(parentChain[index]).")
            XCTAssertLessThan(angle, 140, "Unexpectedly linear parent-chain angle at atom \(parentChain[index]).")
        }
    }

    func testMarkushDefinitionsAreLaidOutBelowRootStructure() throws {
        let molecule = try parse("C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|")
        let components = connectedComponents(in: molecule)

        let rootComponents = components.filter { membership(of: $0, molecule: molecule) == nil }
        let r1Components = components.filter { membership(of: $0, molecule: molecule) == "R1" }

        XCTAssertEqual(rootComponents.count, 1)
        XCTAssertEqual(r1Components.count, 3)

        let rootMidY = try XCTUnwrap(rowBounds(for: rootComponents, molecule: molecule)?.midY)
        let r1MidY = try XCTUnwrap(rowBounds(for: r1Components, molecule: molecule)?.midY)
        XCTAssertGreaterThan(rootMidY, r1MidY)
    }

    func testMarkushRowsStayHorizontallyCentered() throws {
        let molecule = try parse("C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|")
        let components = connectedComponents(in: molecule)

        let rootComponents = components.filter { membership(of: $0, molecule: molecule) == nil }
        let r1Components = components.filter { membership(of: $0, molecule: molecule) == "R1" }

        let rootMidX = try XCTUnwrap(rowBounds(for: rootComponents, molecule: molecule)?.midX)
        let r1MidX = try XCTUnwrap(rowBounds(for: r1Components, molecule: molecule)?.midX)

        XCTAssertEqual(rootMidX, r1MidX, accuracy: 0.35)
    }

    func testMarkushRootOrientationMatchesCDKLayoutIntent() throws {
        let molecule = try parse("C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|")
        let components = connectedComponents(in: molecule)
        let rootComponent = try XCTUnwrap(components.first { membership(of: $0, molecule: molecule) == nil })
        let rootBounds = try XCTUnwrap(rowBounds(for: [rootComponent], molecule: molecule))

        let repeatAtomID = try XCTUnwrap(molecule.sgroups.first?.atomIDs.first)
        let repeatAtom = try XCTUnwrap(molecule.atom(id: repeatAtomID))
        let rGroupAtom = try XCTUnwrap(molecule.atoms.first { $0.rGroupMembership == nil && $0.symbolToDraw == "R1" })
        let nitrogen = try XCTUnwrap(molecule.atoms.first { $0.element.uppercased() == "N" && rootComponent.contains($0.id) })

        XCTAssertLessThan(repeatAtom.position.x, rootBounds.midX)
        XCTAssertGreaterThan(rGroupAtom.position.x, rootBounds.midX)
        XCTAssertLessThan(nitrogen.position.y, rootBounds.midY)
    }

    func testMultiAtomMarkushDefinitionsUseAttachmentDrivenOrientation() throws {
        let molecule = try parse("C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|")
        let r1Atoms = molecule.atoms.filter { $0.rGroupMembership == "R1" }
        let groups = Dictionary(grouping: r1Atoms, by: { $0.componentGroupID ?? -1 })

        let methoxy = try XCTUnwrap(groups.values.first { atoms in
            atoms.contains(where: { $0.element == "O" }) && atoms.contains(where: { $0.element == "C" })
        })
        let methoxyAttachment = try XCTUnwrap(methoxy.first { $0.attachmentPoint == 1 })
        let oxygen = try XCTUnwrap(methoxy.first { $0.element == "O" })
        let methyl = try XCTUnwrap(methoxy.first { $0.element == "C" })

        XCTAssertLessThan(methoxyAttachment.position.x, oxygen.position.x)
        XCTAssertGreaterThan(methoxyAttachment.position.y, oxygen.position.y)
        XCTAssertGreaterThan(methyl.position.x, oxygen.position.x)

        let nitrile = try XCTUnwrap(groups.values.first { atoms in
            atoms.contains(where: { $0.element == "N" }) && atoms.contains(where: { $0.element == "C" })
        })
        let nitrileAttachment = try XCTUnwrap(nitrile.first { $0.attachmentPoint == 1 })
        let nitrileNitrogen = try XCTUnwrap(nitrile.first { $0.element == "N" })

        XCTAssertLessThan(nitrileAttachment.position.x, nitrileNitrogen.position.x)
        XCTAssertGreaterThan(nitrileAttachment.position.y, nitrileNitrogen.position.y)
    }

    private func angleDegrees(a: CGPoint, center: CGPoint, b: CGPoint) -> CGFloat {
        let v1 = CGVector(dx: a.x - center.x, dy: a.y - center.y)
        let v2 = CGVector(dx: b.x - center.x, dy: b.y - center.y)
        let l1 = max(0.0001, hypot(v1.dx, v1.dy))
        let l2 = max(0.0001, hypot(v2.dx, v2.dy))
        let dot = (v1.dx * v2.dx + v1.dy * v2.dy) / (l1 * l2)
        let clamped = max(-1.0, min(1.0, dot))
        return acos(clamped) * 180.0 / .pi
    }

    private func parse(_ smiles: String) throws -> Molecule {
        try parser.parseSmiles(smiles)
    }

    private func expandedWithExplicitHydrogens(from molecule: Molecule) -> Molecule {
        var atoms = molecule.atoms.map { atom in
            var copied = atom
            copied.position = .zero
            copied.explicitHydrogenCount = nil
            return copied
        }
        var bonds = molecule.bonds
        var nextAtomID = (atoms.map(\.id).max() ?? 0) + 1
        var nextBondID = (bonds.map(\.id).max() ?? 0) + 1

        for atom in molecule.atoms where atom.element.uppercased() != "H" {
            let hydrogenCount = molecule.implicitHydrogenCount(for: atom.id)
            guard hydrogenCount > 0 else { continue }
            for _ in 0..<hydrogenCount {
                atoms.append(Atom(id: nextAtomID, element: "H", position: .zero))
                bonds.append(Bond(id: nextBondID, a1: atom.id, a2: nextAtomID, order: .single))
                nextAtomID += 1
                nextBondID += 1
            }
        }

        return Molecule(name: molecule.name, atoms: atoms, bonds: bonds)
    }

    private func heavyNeighborAngles(in molecule: Molecule, around atomID: Int) -> [CGFloat] {
        let heavyNeighborIDs = molecule.neighbors(of: atomID).filter { neighborID in
            molecule.atom(id: neighborID)?.element.uppercased() != "H"
        }
        guard heavyNeighborIDs.count >= 2,
              let center = molecule.atom(id: atomID)?.position else {
            return []
        }

        var angles: [CGFloat] = []
        for leftIndex in 0..<heavyNeighborIDs.count {
            for rightIndex in (leftIndex + 1)..<heavyNeighborIDs.count {
                guard let leftPoint = molecule.atom(id: heavyNeighborIDs[leftIndex])?.position,
                      let rightPoint = molecule.atom(id: heavyNeighborIDs[rightIndex])?.position else {
                    continue
                }
                angles.append(angleDegrees(a: leftPoint, center: center, b: rightPoint))
            }
        }
        return angles.sorted()
    }

    private func connectedComponents(in molecule: Molecule) -> [Set<Int>] {
        var adjacency: [Int: Set<Int>] = [:]
        for atom in molecule.atoms {
            _ = adjacency[atom.id, default: []]
        }
        for bond in molecule.bonds {
            adjacency[bond.a1, default: []].insert(bond.a2)
            adjacency[bond.a2, default: []].insert(bond.a1)
        }

        var seen = Set<Int>()
        var components: [Set<Int>] = []

        for atom in molecule.atoms where !seen.contains(atom.id) {
            var stack = [atom.id]
            var component = Set<Int>()
            seen.insert(atom.id)

            while let current = stack.popLast() {
                component.insert(current)
                for neighbor in adjacency[current, default: []] where !seen.contains(neighbor) {
                    seen.insert(neighbor)
                    stack.append(neighbor)
                }
            }

            components.append(component)
        }

        return components
    }

    private func membership(of component: Set<Int>, molecule: Molecule) -> String? {
        component.compactMap { molecule.atom(id: $0)?.rGroupMembership }.sorted().first
    }

    private func rowBounds(for components: [Set<Int>], molecule: Molecule) -> CGRect? {
        var result: CGRect?
        for component in components {
            let points = component.compactMap { molecule.atom(id: $0)?.position }
            guard let first = points.first else { continue }
            let box = points.dropFirst().reduce(CGRect(origin: first, size: .zero)) { partial, point in
                partial.union(CGRect(origin: point, size: .zero))
            }
            result = result?.union(box) ?? box
        }
        return result
    }
}
