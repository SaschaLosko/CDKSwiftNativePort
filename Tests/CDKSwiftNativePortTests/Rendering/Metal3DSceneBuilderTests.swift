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
