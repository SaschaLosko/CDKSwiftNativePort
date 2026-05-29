# Linux Guide

This guide shows how to use `CDKSwiftNativePort` as a standalone Swift package
on Linux.

The package is headless. It does not require AtomLens, Xcode, `AppKit`,
`UIKit`, `SwiftUI`, Quick Look, or Spotlight infrastructure.

## Verified Environment

The package has been verified on:

- Ubuntu 26.04 LTS
- Swift 6.3.2

The same public API is intended to work in other Swift-on-Linux environments as
long as the package can be built with Swift Package Manager.

## What Works on Linux

The supported Linux path includes:

- chemistry models and Sgroup metadata
- SMILES, CXSMILES, reaction SMILES, MDL, SDF, CML, RXN, RDF, MOL2, PDB, XYZ,
  and InChI import/export
- 2D coordinate generation
- SVG generation
- descriptor and identifier services
- scene generation through `CDKMetalDepictionSceneBuilder` and
  `CDKMetalReactionDepictionSceneBuilder`

What stays outside the package:

- UI frameworks
- renderer/window ownership
- command-line argument parsing
- app and extension packaging

## Create a Linux Executable Package

Create a new executable package:

```bash
mkdir ChemCLI
cd ChemCLI
swift package init --type executable
```

Update `Package.swift`:

```swift
// swift-tools-version: 6.0
import PackageDescription

let package = Package(
    name: "ChemCLI",
    dependencies: [
        .package(url: "https://github.com/SaschaLosko/CDKSwiftNativePort.git", from: "1.4.0")
    ],
    targets: [
        .executableTarget(
            name: "ChemCLI",
            dependencies: [
                .product(name: "CDKSwiftNativePort", package: "CDKSwiftNativePort")
            ]
        )
    ]
)
```

## Example 1: SMILES to SVG CLI

Replace `Sources/main.swift` with:

```swift
import Foundation
import CDKSwiftNativePort

let args = CommandLine.arguments
guard args.count == 2 else {
    fputs("usage: ChemCLI '<SMILES>'\n", stderr)
    exit(64)
}

let molecules = try CDKFileImporter.readMolecules(text: args[1], fileExtension: "smi")
guard let molecule = molecules.first else {
    fputs("No molecule could be parsed.\n", stderr)
    exit(65)
}

let laidOut = Depiction2DGenerator.generate(for: molecule)
let identifiers = CDKMoleculeIdentifierService.compute(for: laidOut)
let properties = CDKMoleculePropertyService.compute(for: laidOut)

let svg = CDKDepictionGenerator.toSVG(
    molecule: laidOut,
    style: RenderStyle(),
    canvasSize: CGSize(width: 960, height: 640),
    includeBackground: true
)

let outputURL = URL(fileURLWithPath: "molecule.svg")
try svg.write(to: outputURL, atomically: true, encoding: .utf8)

print("SMILES:", identifiers.smiles)
print("InChIKey:", identifiers.inchiKey)
print("Formula:", properties.formula)
print("Wrote:", outputURL.path)
```

Build and run it:

```bash
swift build
swift run ChemCLI "CC(=O)OC1=CC=CC=C1C(=O)O"
```

That produces `molecule.svg` in the current working directory.

## Example 2: Read a File and Export V3000

This example reads any supported molecule file, generates 2D coordinates if
needed, and prints a V3000 molfile.

`Sources/main.swift`:

```swift
import Foundation
import CDKSwiftNativePort

let args = CommandLine.arguments
guard args.count == 2 else {
    fputs("usage: ChemCLI <input-file>\n", stderr)
    exit(64)
}

let inputURL = URL(fileURLWithPath: args[1])
let molecules = try CDKFileImporter.readMolecules(from: inputURL)
guard let molecule = molecules.first else {
    fputs("No molecule found in input file.\n", stderr)
    exit(65)
}

let laidOut = Depiction2DGenerator.generate(for: molecule)
let molfile = try CDKFileExporter.write(molecule: laidOut, as: .molV3000)

print(molfile)
```

Run it:

```bash
swift run ChemCLI ./sample.sdf
```

This pattern works well for:

- batch conversion tools
- CI validation utilities
- server-side chemistry preprocessing

## Example 3: Build a Renderer-Neutral Scene on Linux

If your Linux host owns its own drawing backend, you can ask the package to
build scene geometry and labels without taking ownership of the renderer:

```swift
import Foundation
import CDKSwiftNativePort

let molecule = try CDKSmilesParserFactory.shared
    .newSmilesParser(flavor: .cdkDefault)
    .parseSmiles("C1=CC=CC=C1")

let laidOut = Depiction2DGenerator.generate(for: molecule)
let scene = CDKMetalDepictionSceneBuilder.build(
    molecule: laidOut,
    style: RenderStyle(showCarbons: false, showImplicitHydrogens: true),
    canvasRect: CGRect(x: 0, y: 0, width: 800, height: 600),
    zoom: 1.0,
    pan: .zero
)

print("bond segments:", scene.bondSegments.count)
print("labels:", scene.labels.count)
print("background boxes:", scene.backgroundBoxes.count)
```

The package computes chemistry layout and scene data. Your Linux host decides
how to rasterize or display it.

## Running the Package Tests on Linux

From the package root:

```bash
swift test
```

For release preparation, run the package suite on Linux as well as on Apple
platforms so the public examples and package boundary stay in sync.
