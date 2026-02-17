import XCTest
import CoreGraphics
@testable import CDKPortExperimental

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
}
