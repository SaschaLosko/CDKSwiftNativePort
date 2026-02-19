// swift-tools-version: 5.9
import PackageDescription

let package = Package(
    name: "CDKSwiftNativePort",
    platforms: [
        .macOS(.v14)
    ],
    products: [
        .library(
            name: "CDKSwiftNativePort",
            targets: ["CDKSwiftNativePort"]
        )
    ],
    targets: [
        .target(
            name: "CDKSwiftNativePort"
        ),
        .testTarget(
            name: "CDKSwiftNativePortTests",
            dependencies: ["CDKSwiftNativePort"],
            exclude: [
                "Smiles/port_metadata.json",
                "MDL/port_metadata.json",
                "InChI/port_metadata.json"
            ]
        )
    ]
)
