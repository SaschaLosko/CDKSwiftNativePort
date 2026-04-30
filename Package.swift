// swift-tools-version: 6.0
import PackageDescription

let package = Package(
    name: "CDKSwiftNativePort",
    platforms: [
        .macOS(.v15),
        .iOS(.v13)
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
                "CML/port_metadata.json",
                "Smiles/port_metadata.json",
                "MDL/port_metadata.json",
                "InChI/port_metadata.json",
                "InChI/OfficialReference/official_reference_cases.json",
                "InChI/OfficialReference/known_gap_inventory.json",
                "Reaction/UpstreamReference",
                "Reaction/port_metadata.json",
                "RInChI/UpstreamReference"
            ]
        )
    ],
    swiftLanguageModes: [.v6]
)
