import XCTest
@testable import CDKSwiftNativePort

final class Metal3DSceneBuilderTests: XCTestCase {
    func testBuildsAtomSpheresAndBondCylindersFrom3DCoordinates() throws {
        let molecule = Molecule(
            name: "Water",
            atoms: [
                Atom(id: 1, element: "O", position: CGPoint(x: 0.0, y: 0.0), zPosition: 0.10),
                Atom(id: 2, element: "H", position: CGPoint(x: 0.96, y: 0.0), zPosition: -0.12),
                Atom(id: 3, element: "H", position: CGPoint(x: -0.24, y: 0.93), zPosition: -0.10),
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single),
                Bond(id: 2, a1: 1, a2: 3, order: .single),
            ])

        let scene = CDKMetal3DSceneBuilder.build(molecule: molecule)

        XCTAssertTrue(scene.hasExplicit3DCoordinates)
        XCTAssertEqual(scene.atoms.count, 3)
        XCTAssertEqual(scene.bonds.count, 2)
        let firstAtom = try XCTUnwrap(scene.atoms.first)
        let firstBond = try XCTUnwrap(scene.bonds.first)
        let boundingBox = try XCTUnwrap(scene.boundingBox)
        XCTAssertEqual(firstAtom.center.z, 0.10, accuracy: 0.0001)
        XCTAssertEqual(firstBond.from.z, 0.10, accuracy: 0.0001)
        XCTAssertEqual(firstBond.to.z, -0.12, accuracy: 0.0001)
        XCTAssertEqual(boundingBox.min.x, -0.24, accuracy: 0.0001)
        XCTAssertEqual(boundingBox.max.y, 0.93, accuracy: 0.0001)
        XCTAssertEqual(boundingBox.min.z, -0.12, accuracy: 0.0001)
        XCTAssertEqual(boundingBox.max.z, 0.10, accuracy: 0.0001)
    }

    func testBuildsFlatSceneWhenZCoordinatesAreUnavailable() {
        let molecule = Molecule(
            name: "Ethane flat",
            atoms: [
                Atom(id: 1, element: "C", position: CGPoint(x: -0.7, y: 0.0)),
                Atom(id: 2, element: "C", position: CGPoint(x: 0.7, y: 0.0)),
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single),
            ])

        let scene = CDKMetal3DSceneBuilder.build(molecule: molecule)

        XCTAssertFalse(scene.hasExplicit3DCoordinates)
        XCTAssertEqual(scene.atoms.map(\.center.z), [0.0, 0.0])
        XCTAssertEqual(scene.boundingBox?.min.z, 0.0)
        XCTAssertEqual(scene.boundingBox?.max.z, 0.0)
    }

    func testZeroZCoordinatesAreNotTreatedAsExplicit3D() throws {
        let molecule = Molecule(
            name: "Cyclohexane zero z",
            atoms: (0..<6).map { index in
                let angle = Double(index) * Double.pi / 3.0
                return Atom(
                    id: index + 1,
                    element: "C",
                    position: CGPoint(x: cos(angle), y: sin(angle)),
                    zPosition: 0.0)
            },
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single),
                Bond(id: 2, a1: 2, a2: 3, order: .single),
                Bond(id: 3, a1: 3, a2: 4, order: .single),
                Bond(id: 4, a1: 4, a2: 5, order: .single),
                Bond(id: 5, a1: 5, a2: 6, order: .single),
                Bond(id: 6, a1: 6, a2: 1, order: .single),
            ])

        let scene = CDKMetal3DSceneBuilder.build(molecule: molecule)
        let boundingBox = try XCTUnwrap(scene.boundingBox)

        XCTAssertFalse(scene.hasExplicit3DCoordinates)
        XCTAssertGreaterThan(boundingBox.max.z - boundingBox.min.z, 0.50)
    }

    func testFlatStereoWedgeGeneratesOutOfPlaneDepth() throws {
        let molecule = Molecule(
            name: "Flat stereo",
            atoms: [
                Atom(id: 1, element: "C", position: CGPoint(x: 0.0, y: 0.0), zPosition: 0.0),
                Atom(id: 2, element: "C", position: CGPoint(x: 1.5, y: 0.0), zPosition: 0.0),
                Atom(id: 3, element: "C", position: CGPoint(x: -1.0, y: 1.0), zPosition: 0.0),
                Atom(id: 4, element: "O", position: CGPoint(x: -1.0, y: -1.0), zPosition: 0.0),
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single, stereo: .up),
                Bond(id: 2, a1: 1, a2: 3, order: .single),
                Bond(id: 3, a1: 1, a2: 4, order: .single),
            ])

        let scene = CDKMetal3DSceneBuilder.build(molecule: molecule)
        let center = try XCTUnwrap(scene.atoms.first { $0.id == 1 })
        let wedgeLigand = try XCTUnwrap(scene.atoms.first { $0.id == 2 })

        XCTAssertFalse(scene.hasExplicit3DCoordinates)
        XCTAssertGreaterThan(wedgeLigand.center.z, center.center.z)
        XCTAssertGreaterThan(abs(wedgeLigand.center.z - center.center.z), 0.50)
    }

    func testReversedFlatStereoWedgeGeneratesDepthOnLigandSide() throws {
        let molecule = Molecule(
            name: "Reversed flat stereo",
            atoms: [
                Atom(id: 1, element: "C", position: CGPoint(x: 0.0, y: 0.0), zPosition: 0.0),
                Atom(id: 2, element: "C", position: CGPoint(x: 1.5, y: 0.0), zPosition: 0.0),
                Atom(id: 3, element: "N", position: CGPoint(x: 3.0, y: 0.0), zPosition: 0.0),
            ],
            bonds: [
                Bond(id: 1, a1: 2, a2: 1, order: .single, stereo: .downReversed),
                Bond(id: 2, a1: 1, a2: 3, order: .single),
            ])

        let scene = CDKMetal3DSceneBuilder.build(molecule: molecule)
        let anchor = try XCTUnwrap(scene.atoms.first { $0.id == 1 })
        let ligand = try XCTUnwrap(scene.atoms.first { $0.id == 2 })

        XCTAssertLessThan(ligand.center.z, anchor.center.z)
        XCTAssertGreaterThan(abs(ligand.center.z - anchor.center.z), 0.50)
    }

    func testRendererModelControlsRadiiAndAtomColoring() throws {
        let molecule = Molecule(
            name: "Carbon monoxide",
            atoms: [
                Atom(id: 1, element: "C", position: CGPoint(x: 0.0, y: 0.0), zPosition: 0.0),
                Atom(id: 2, element: "O", position: CGPoint(x: 1.1, y: 0.0), zPosition: 0.0),
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .triple),
            ])
        let rendererModel = CDKRenderer3DModel(
            atomRadiusScale: 0.50,
            bondRadius: 0.14,
            atomColoringMode: .cdk2D,
            colorBondsByAtom: false)

        let scene = CDKMetal3DSceneBuilder.build(molecule: molecule, rendererModel: rendererModel)
        let firstBond = try XCTUnwrap(scene.bonds.first)

        XCTAssertEqual(firstBond.radius, 0.14, accuracy: 0.0001)
        XCTAssertEqual(firstBond.order, .triple)
        XCTAssertGreaterThan(scene.atoms[0].radius, 0.30)
        XCTAssertGreaterThan(scene.atoms[1].color.red, 0.80)
        XCTAssertLessThan(scene.atoms[1].color.green, 0.30)
    }
}
