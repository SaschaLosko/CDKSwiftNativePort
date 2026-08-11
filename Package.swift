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
            type: .dynamic,
            targets: ["CDKSwiftNativePort"]
        )
    ],
    targets: [
        .target(
            name: "IUPACInChI",
            path: "Sources/IUPACInChI",
            exclude: [
                "External-contributors-IUPAC-InChI.txt",
                "LICENSE-IUPAC-InChI.txt"
            ],
            publicHeadersPath: "include",
            cSettings: [
                .define("COMPILE_ANSI_ONLY"),
                .define("TARGET_API_LIB"),
                .headerSearchPath("include"),
                .headerSearchPath("official/INCHI-1-SRC/INCHI_BASE/src"),
                .headerSearchPath("official/INCHI-1-SRC/INCHI_API/libinchi/src"),
                .headerSearchPath("official/INCHI-1-SRC/INCHI_API/libinchi/src/ixa")
            ],
            linkerSettings: [
                .linkedLibrary("m", .when(platforms: [.linux]))
            ]
        ),
        .target(
            name: "CDKSwiftNativePort",
            dependencies: ["IUPACInChI"]
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
