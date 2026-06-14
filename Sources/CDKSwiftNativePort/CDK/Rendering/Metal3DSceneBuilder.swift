#if canImport(CoreGraphics)
import CoreGraphics
#endif
import Foundation

public struct CDKPoint3D: Hashable, Codable, Sendable {
    public var x: Double
    public var y: Double
    public var z: Double

    public init(x: Double, y: Double, z: Double) {
        self.x = x
        self.y = y
        self.z = z
    }
}

public struct CDK3DBoundingBox: Hashable, Sendable {
    public let min: CDKPoint3D
    public let max: CDKPoint3D

    public var center: CDKPoint3D {
        CDKPoint3D(
            x: (min.x + max.x) * 0.5,
            y: (min.y + max.y) * 0.5,
            z: (min.z + max.z) * 0.5)
    }

    public var size: CDKPoint3D {
        CDKPoint3D(
            x: max.x - min.x,
            y: max.y - min.y,
            z: max.z - min.z)
    }

    public var maximumExtent: Double {
        Swift.max(size.x, Swift.max(size.y, size.z))
    }

    public init(min: CDKPoint3D, max: CDKPoint3D) {
        self.min = min
        self.max = max
    }
}

public struct CDKRenderer3DModel: Hashable, Sendable {
    public var atomRadiusScale: Double
    public var bondRadius: Double
    public var minimumAtomRadius: Double
    public var maximumAtomRadius: Double
    public var atomColoringMode: CDKAtomColoringMode
    public var colorBondsByAtom: Bool

    public init(
        atomRadiusScale: Double = 0.34,
        bondRadius: Double = 0.10,
        minimumAtomRadius: Double = 0.22,
        maximumAtomRadius: Double = 0.58,
        atomColoringMode: CDKAtomColoringMode = .cdk2D,
        colorBondsByAtom: Bool = true
    ) {
        self.atomRadiusScale = atomRadiusScale
        self.bondRadius = bondRadius
        self.minimumAtomRadius = minimumAtomRadius
        self.maximumAtomRadius = maximumAtomRadius
        self.atomColoringMode = atomColoringMode
        self.colorBondsByAtom = colorBondsByAtom
    }
}

public struct CDKMetal3DScene: Hashable, Sendable {
    public struct AtomSphere: Hashable, Sendable {
        public let id: Int
        public let element: String
        public let center: CDKPoint3D
        public let radius: Double
        public let color: CDKRenderColor

        public init(id: Int, element: String, center: CDKPoint3D, radius: Double, color: CDKRenderColor) {
            self.id = id
            self.element = element
            self.center = center
            self.radius = radius
            self.color = color
        }
    }

    public struct BondCylinder: Hashable, Sendable {
        public let id: Int
        public let fromAtomID: Int
        public let toAtomID: Int
        public let from: CDKPoint3D
        public let to: CDKPoint3D
        public let radius: Double
        public let order: BondOrder
        public let color: CDKRenderColor

        public init(
            id: Int,
            fromAtomID: Int,
            toAtomID: Int,
            from: CDKPoint3D,
            to: CDKPoint3D,
            radius: Double,
            order: BondOrder,
            color: CDKRenderColor
        ) {
            self.id = id
            self.fromAtomID = fromAtomID
            self.toAtomID = toAtomID
            self.from = from
            self.to = to
            self.radius = radius
            self.order = order
            self.color = color
        }
    }

    public let atoms: [AtomSphere]
    public let bonds: [BondCylinder]
    public let boundingBox: CDK3DBoundingBox?
    public let hasExplicit3DCoordinates: Bool

    public init(
        atoms: [AtomSphere],
        bonds: [BondCylinder],
        boundingBox: CDK3DBoundingBox?,
        hasExplicit3DCoordinates: Bool
    ) {
        self.atoms = atoms
        self.bonds = bonds
        self.boundingBox = boundingBox
        self.hasExplicit3DCoordinates = hasExplicit3DCoordinates
    }
}

public enum CDKMetal3DSceneBuilder {
    public static func build(
        molecule: Molecule,
        rendererModel: CDKRenderer3DModel = CDKRenderer3DModel()
    ) -> CDKMetal3DScene {
        let hasExplicit3DCoordinates = CDKModelBuilder3D.hasMeaningful3DCoordinates(molecule)
        let sceneMolecule = hasExplicit3DCoordinates
            ? molecule
            : CDKModelBuilder3D.generate3DCoordinates(molecule: molecule)
        let style = RenderStyle(
            atomColoringMode: rendererModel.atomColoringMode,
            colorBondsByAtom: rendererModel.colorBondsByAtom)
        let atomPoints = Dictionary(uniqueKeysWithValues: sceneMolecule.atoms.map { atom in
            (atom.id, point3D(for: atom))
        })

        let atoms = sceneMolecule.atoms.map { atom in
            CDKMetal3DScene.AtomSphere(
                id: atom.id,
                element: atom.element,
                center: atomPoints[atom.id] ?? point3D(for: atom),
                radius: radius(for: atom, rendererModel: rendererModel),
                color: CDKRenderingStyleResolver.atomColor(for: atom, style: style))
        }

        let bonds = sceneMolecule.bonds.compactMap { bond -> CDKMetal3DScene.BondCylinder? in
            guard let from = atomPoints[bond.a1],
                  let to = atomPoints[bond.a2] else {
                return nil
            }
            return CDKMetal3DScene.BondCylinder(
                id: bond.id,
                fromAtomID: bond.a1,
                toAtomID: bond.a2,
                from: from,
                to: to,
                radius: max(0.025, rendererModel.bondRadius),
                order: bond.order,
                color: CDKRenderingStyleResolver.bondColor(for: bond, molecule: sceneMolecule, style: style))
        }

        return CDKMetal3DScene(
            atoms: atoms,
            bonds: bonds,
            boundingBox: boundingBox(points: Array(atomPoints.values)),
            hasExplicit3DCoordinates: hasExplicit3DCoordinates)
    }

    private static func point3D(for atom: Atom) -> CDKPoint3D {
        CDKPoint3D(
            x: Double(atom.position.x),
            y: Double(atom.position.y),
            z: atom.zPosition ?? 0.0)
    }

    private static func radius(for atom: Atom, rendererModel: CDKRenderer3DModel) -> Double {
        let normalized = atom.element.trimmingCharacters(in: .whitespacesAndNewlines).uppercased()
        let covalentRadius: Double =
            switch normalized {
            case "H": 0.31
            case "B": 0.84
            case "C": 0.76
            case "N": 0.71
            case "O": 0.66
            case "F": 0.57
            case "P": 1.07
            case "S": 1.05
            case "CL": 1.02
            case "BR": 1.20
            case "I": 1.39
            default: 0.80
            }
        let scaled = covalentRadius * rendererModel.atomRadiusScale
        return min(rendererModel.maximumAtomRadius, max(rendererModel.minimumAtomRadius, scaled))
    }

    private static func boundingBox(points: [CDKPoint3D]) -> CDK3DBoundingBox? {
        guard var minX = points.first?.x,
              var minY = points.first?.y,
              var minZ = points.first?.z,
              var maxX = points.first?.x,
              var maxY = points.first?.y,
              var maxZ = points.first?.z else {
            return nil
        }
        for point in points.dropFirst() {
            minX = min(minX, point.x)
            minY = min(minY, point.y)
            minZ = min(minZ, point.z)
            maxX = max(maxX, point.x)
            maxY = max(maxY, point.y)
            maxZ = max(maxZ, point.z)
        }
        return CDK3DBoundingBox(
            min: CDKPoint3D(x: minX, y: minY, z: minZ),
            max: CDKPoint3D(x: maxX, y: maxY, z: maxZ))
    }
}
