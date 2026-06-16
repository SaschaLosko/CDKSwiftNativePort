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

    func testRendererModelCanBuildSpaceFillingScene() throws {
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
            atomRadiusScale: 0.01,
            bondRadius: 0.20,
            minimumAtomRadius: 0.01,
            maximumAtomRadius: 0.02,
            atomColoringMode: .cdk2D,
            atomColorPalette: .cpk,
            representationMode: .spaceFilling)

        let scene = CDKMetal3DSceneBuilder.build(molecule: molecule, rendererModel: rendererModel)
        let carbon = try XCTUnwrap(scene.atoms.first { $0.id == 1 })
        let oxygen = try XCTUnwrap(scene.atoms.first { $0.id == 2 })

        XCTAssertEqual(CDK3DRepresentationMode.allCases, [.ballAndStick, .stick, .spaceFilling])
        XCTAssertEqual(CDK3DRepresentationMode.spaceFilling.displayName, "Space Filling")
        XCTAssertEqual(scene.bonds.count, 0)
        XCTAssertEqual(carbon.radius, 1.70, accuracy: 0.0001)
        XCTAssertEqual(oxygen.radius, 1.55, accuracy: 0.0001)
        assertColor(carbon.color, hex: 0xC8C8C8)
        assertColor(oxygen.color, hex: 0xF00000)
    }

    func testRendererModelCanBuildStickScene() throws {
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
            bondRadius: 0.20,
            atomColoringMode: .cdk2D,
            atomColorPalette: .cpk,
            representationMode: .stick)

        let scene = CDKMetal3DSceneBuilder.build(molecule: molecule, rendererModel: rendererModel)
        let bond = try XCTUnwrap(scene.bonds.first)

        XCTAssertEqual(CDK3DRepresentationMode.allCases, [.ballAndStick, .stick, .spaceFilling])
        XCTAssertEqual(CDK3DRepresentationMode.stick.displayName, "Stick")
        XCTAssertEqual(scene.atoms.count, 0)
        XCTAssertEqual(scene.bonds.count, 1)
        XCTAssertEqual(bond.radius, 0.27, accuracy: 0.0001)
        XCTAssertNotNil(scene.boundingBox)
    }

    func testDefaultRendererModelUsesBallAndStickRepresentation() {
        XCTAssertEqual(CDKRenderer3DModel().representationMode, .ballAndStick)
    }

    func testRendererModelUsesSelected3DAtomColorPalette() throws {
        let molecule = Molecule(
            name: "Chloromethane",
            atoms: [
                Atom(id: 1, element: "C", position: CGPoint(x: 0.0, y: 0.0), zPosition: 0.0),
                Atom(id: 2, element: "Cl", position: CGPoint(x: 1.5, y: 0.0), zPosition: 0.0),
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single),
            ])

        let jmolScene = CDKMetal3DSceneBuilder.build(
            molecule: molecule,
            rendererModel: CDKRenderer3DModel(
                atomColoringMode: .cdk2D,
                atomColorPalette: .jmol,
                colorBondsByAtom: true))
        let cpkScene = CDKMetal3DSceneBuilder.build(
            molecule: molecule,
            rendererModel: CDKRenderer3DModel(
                atomColoringMode: .cdk2D,
                atomColorPalette: .cpk,
                colorBondsByAtom: true))

        let jmolCarbon = try XCTUnwrap(jmolScene.atoms.first { $0.id == 1 })
        let jmolChlorine = try XCTUnwrap(jmolScene.atoms.first { $0.id == 2 })
        let cpkCarbon = try XCTUnwrap(cpkScene.atoms.first { $0.id == 1 })
        let cpkChlorine = try XCTUnwrap(cpkScene.atoms.first { $0.id == 2 })
        let jmolBond = try XCTUnwrap(jmolScene.bonds.first)
        let cpkBond = try XCTUnwrap(cpkScene.bonds.first)

        assertColor(jmolCarbon.color, hex: 0x909090)
        assertColor(jmolChlorine.color, hex: 0x1FF01F)
        assertColor(cpkCarbon.color, hex: 0xC8C8C8)
        assertColor(cpkChlorine.color, hex: 0x00FF00)
        XCTAssertNotEqual(jmolCarbon.color, cpkCarbon.color)
        XCTAssertNotEqual(jmolChlorine.color, cpkChlorine.color)
        XCTAssertNotEqual(jmolBond.color, cpkBond.color)
    }

    func testRendererModelUsesAdditional3DAtomColorPalettes() throws {
        let molecule = Molecule(
            name: "Carbon monoxide",
            atoms: [
                Atom(id: 1, element: "C", position: CGPoint(x: 0.0, y: 0.0), zPosition: 0.0),
                Atom(id: 2, element: "O", position: CGPoint(x: 1.1, y: 0.0), zPosition: 0.0),
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .triple),
            ])

        XCTAssertEqual(CDK3DAtomColorPalette.allCases.count, 14)
        XCTAssertEqual(CDK3DAtomColorPalette.okabeIto.displayName, "Okabe-Ito")
        XCTAssertEqual(CDK3DAtomColorPalette.colorBrewerSet2.displayName, "ColorBrewer Set2")

        for palette in CDK3DAtomColorPalette.allCases {
            let scene = CDKMetal3DSceneBuilder.build(
                molecule: molecule,
                rendererModel: CDKRenderer3DModel(
                    atomColoringMode: .cdk2D,
                    atomColorPalette: palette,
                    colorBondsByAtom: true))
            XCTAssertEqual(scene.atoms.count, 2)
            XCTAssertEqual(scene.bonds.count, 1)
            XCTAssertFalse(palette.displayName.isEmpty)
        }

        assertCarbonColor(for: molecule, palette: .okabeIto, hex: 0xD55E00)
        assertCarbonColor(for: molecule, palette: .viridis, hex: 0x1F9E89)
        assertCarbonColor(for: molecule, palette: .cividis, hex: 0x8A8678)
        assertCarbonColor(for: molecule, palette: .cartoVivid, hex: 0x24796C)
        assertCarbonColor(for: molecule, palette: .matplotlibTab20, hex: 0x98DF8A)
    }

    private func assertCarbonColor(
        for molecule: Molecule,
        palette: CDK3DAtomColorPalette,
        hex: UInt32,
        file: StaticString = #filePath,
        line: UInt = #line
    ) {
        let scene = CDKMetal3DSceneBuilder.build(
            molecule: molecule,
            rendererModel: CDKRenderer3DModel(atomColoringMode: .cdk2D, atomColorPalette: palette))
        guard let carbon = scene.atoms.first(where: { $0.element == "C" }) else {
            XCTFail("Missing carbon atom", file: file, line: line)
            return
        }
        assertColor(carbon.color, hex: hex, file: file, line: line)
    }

    private func assertColor(
        _ color: CDKRenderColor,
        hex: UInt32,
        file: StaticString = #filePath,
        line: UInt = #line
    ) {
        XCTAssertEqual(
            color.red,
            CGFloat((hex >> 16) & 0xFF) / 255.0,
            accuracy: 0.0001,
            file: file,
            line: line)
        XCTAssertEqual(
            color.green,
            CGFloat((hex >> 8) & 0xFF) / 255.0,
            accuracy: 0.0001,
            file: file,
            line: line)
        XCTAssertEqual(
            color.blue,
            CGFloat(hex & 0xFF) / 255.0,
            accuracy: 0.0001,
            file: file,
            line: line)
        XCTAssertEqual(color.alpha, 1.0, accuracy: 0.0001, file: file, line: line)
    }
}
