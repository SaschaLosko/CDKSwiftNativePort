// swift-tools-version: 5.9
import PackageDescription

let package = Package(
    name: "CDKPortExperimental",
    platforms: [
        .macOS(.v14)
    ],
    products: [
        .library(
            name: "CDKPortExperimental",
            targets: ["CDKPortExperimental"]
        )
    ],
    targets: [
        .target(
            name: "CDKPortExperimental"
        ),
        .testTarget(
            name: "CDKPortExperimentalTests",
            dependencies: ["CDKPortExperimental"],
            exclude: ["Smiles/port_metadata.json"]
        )
    ]
)
