# CDKSwiftNativePort

`CDKSwiftNativePort` is a Swift-native chemistry package that keeps
CDK-derived parsing, layout, depiction, identifiers, descriptors, and file I/O
inside a reusable package boundary.

GitHub repository: <https://github.com/SaschaLosko/CDKSwiftNativePort>

## Overview

This package exists for two reasons:

- keep all CDK-derived chemistry logic out of the host app
- expose a Swift-first chemistry API that can be reused from Apple apps,
  command-line tools, and Swift on Linux

The first-party consumer is AtomLens, but the package itself has no dependency
on AtomLens modules, Quick Look targets, Spotlight targets, bundle identifiers,
or app session logic.

## Feature Set

### Core chemistry model

- `Molecule`, `Atom`, `Bond`, `MoleculeSgroup`
- `CDKReaction`, `CDKReactionParticipant`, `CDKReactionDirection`
- stereochemistry, atom maps, aliases, attachment points
- query atoms and query bonds
- SDF-style data fields and Sgroup metadata

### SMILES and CXSMILES

- SMILES parsing and generation
- reaction SMILES parsing and generation
- native reaction hierarchy types for sets, lists, schemes, and lossless list entries
- native reaction manipulator utilities for counts, reversal, inline conversion,
  mapped-object lookup, and set-level queries
- CXSMILES parsing for:
  - atom labels and atom values
  - 2D and 3D coordinates
  - radicals
  - enhanced stereo metadata
  - fragment grouping
  - highlight layers (`ha:` / `hb:`)
  - ligand ordering
  - polymer, data, generic, and positional-variation Sgroups
  - Markush `RG:` and link-node `LN:` handling

### Layout and depiction

- CDK-derived 2D coordinate generation via `Depiction2DGenerator`
- SVG depiction via `CDKDepictionGenerator`
- renderer-neutral scene generation for host Metal/Core Graphics pipelines:
  - `CDKMetalDepictionSceneBuilder`
  - `CDKMetalReactionDepictionSceneBuilder`
- depiction coverage includes:
  - aromatic circle mode
  - Markush legends and link-node annotations
  - non-Markush Sgroup brackets
  - query atom and query bond visuals
  - highlight styles: none, colored, glow, glow with white edge

### File formats

Import and export coverage currently includes:

- MDL Molfile V2000
- MDL Molfile V3000
- SDF
- SMILES and CXSMILES
- reaction SMILES
- RInChI text files
- InChI
- MOL2
- PDB
- XYZ
- CML molecules
- CML reactions
- RXN V2000 and V3000
- RDF
- SVG depiction export

### Identifiers and descriptors

- SMILES, isomeric SMILES, InChI, InChIKey
- RInChI, RAuxInfo, Long-RInChIKey, Short-RInChIKey, Web-RInChIKey
- formula, molecular weight, exact mass
- heavy atom count, donor/acceptor counts, ring count, rotatable bonds
- XLogP, Mannhold LogP, TPSA, van der Waals volume
- Rule-of-Five summary

## Encapsulation and App Boundary

All CDK-derived implementation lives under `Sources/CDKSwiftNativePort`.

The package owns:

- molecule and reaction models
- reaction hierarchy models (`CDKReactionSet`, `CDKReactionList`, `CDKReactionScheme`) and entry enums
- reaction manipulator helpers (`CDKReactionManipulator`, `CDKReactionSetManipulator`, `CDKReactionSchemeManipulator`)
- readers, writers, parsers, generators
- layout and depiction algorithms
- scene building for host renderers
- identifier and descriptor services

The host app or tool owns:

- window and document behavior
- UI state and toolbar controls
- render loops and view ownership
- Quick Look and Spotlight extension targets
- bundle identifiers, entitlements, and app metadata
- command-line argument parsing and process lifecycle

This boundary is enforced in tests:

- package sources are scanned for forbidden host-coupling imports
- `Package.swift` is checked for AtomLens target references
- the monorepo host chemistry layer is checked to ensure it only contains the
  thin alias adapter
- the AtomLens Xcode project is checked to ensure it links the package product
  rather than compiling package implementation files directly

## Platform and Release Status

- Swift tools: `5.9`
- Apple deployment targets declared in `Package.swift`:
  - `macOS 14+`
  - `iOS 13+`
- Linux:
  - builds and runs as a Swift package under Swift on Linux
  - verified on Ubuntu 24.04 with Swift 6.3

The package is headless. It does not depend on `AppKit`, `UIKit`, `SwiftUI`,
Quick Look, Spotlight, or AtomLens internals. That makes it suitable for:

- macOS and iOS apps
- helper tools and import/export utilities
- Quick Look and Spotlight-style host integrations on Apple platforms
- Linux command-line tools, batch converters, and service-side workflows

## InChI Status

The package includes a native Swift InChI implementation and an official
reference parity harness. The runtime stays Swift native, while the vendored
official InChI corpus is used as the correctness oracle for the
reference-grade work tracked in
[`Documentation/INCHI_REFERENCE_GRADE.md`](Documentation/INCHI_REFERENCE_GRADE.md).

## Reaction Parity Status

Reaction hierarchy and reaction-IO parity now have a committed upstream-derived
reference corpus as well. The vendored reaction fixtures now cover CML
hierarchy cases, RXN V2000 and V3000 resources, strict/relaxed RXN counts-line
behavior, RDF parsing, and reaction SMILES behavior, with zero committed gaps
across the current `22`-case corpus. Native reaction manipulator and
scheme-manipulator coverage is tracked separately through executable upstream
port metadata and direct Swift regression tests. Inventory/gate commands are
documented in
[`Documentation/REACTION_PARITY.md`](Documentation/REACTION_PARITY.md).

## RInChI Status

The package now also includes a native Swift `RInChI` surface for reactions,
including `RInChI` generation, `RAuxInfo`, Long/Short/Web keys,
decomposition, and `RInChI`-to-reaction reconstruction. Vendored upstream CDK
`storage/rinchi` fixtures are exercised through a strict executable gate
documented in [`Documentation/RINCHI_PARITY.md`](Documentation/RINCHI_PARITY.md).

## Installation

### Xcode

Add the package dependency:

```text
https://github.com/SaschaLosko/CDKSwiftNativePort.git
```

### Swift Package Manager

Add the dependency in `Package.swift`:

```swift
dependencies: [
    .package(url: "https://github.com/SaschaLosko/CDKSwiftNativePort.git", from: "1.3.0")
]
```

Then add the product to your target:

```swift
.target(
    name: "MyTarget",
    dependencies: [
        .product(name: "CDKSwiftNativePort", package: "CDKSwiftNativePort")
    ]
)
```

### SwiftPM Executable on Linux

For a Linux command-line tool, use an executable target:

```swift
// Package.swift
// swift-tools-version: 6.0
import PackageDescription

let package = Package(
    name: "ChemCLI",
    dependencies: [
        .package(url: "https://github.com/SaschaLosko/CDKSwiftNativePort.git", from: "1.3.0")
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

## Quick Start

### Parse, Layout, Identify, and Export SVG

This example works in a Swift package on macOS and Linux.

```swift
import Foundation
import CDKSwiftNativePort

let molecules = try CDKFileImporter.readMolecules(
    text: "CC(=O)OC1=CC=CC=C1C(=O)O",
    fileExtension: "smi"
)

guard let molecule = molecules.first else {
    throw NSError(domain: "CDKSwiftNativePort", code: 1, userInfo: [
        NSLocalizedDescriptionKey: "No molecule could be parsed."
    ])
}

let laidOut = Depiction2DGenerator.generate(for: molecule)
let identifiers = CDKMoleculeIdentifierService.compute(for: laidOut)
let properties = CDKMoleculePropertyService.compute(for: laidOut)

let svg = CDKDepictionGenerator.toSVG(
    molecule: laidOut,
    style: RenderStyle(),
    canvasSize: CGSize(width: 1280, height: 840),
    includeBackground: true
)

print(identifiers.smiles)
print(identifiers.inchiKey)
print(properties.formula)
print(svg.prefix(120))
```

### Build a Linux CLI That Writes an SVG File

`Sources/ChemCLI/main.swift`:

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
print("Wrote:", outputURL.path)
```

Run it with SwiftPM:

```bash
swift build
swift run ChemCLI "CC(=O)OC1=CC=CC=C1C(=O)O"
```

For more detailed Linux examples, including file import/export workflows, see
[`Documentation/LINUX.md`](Documentation/LINUX.md).

### Parse a Markush CXSMILES and Build a Scene

```swift
import Foundation
import CDKSwiftNativePort

let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
let molecule = try parser.parseSmiles(
    "C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|"
)

let laidOut = Depiction2DGenerator.generate(for: molecule)
let scene = CDKMetalDepictionSceneBuilder.build(
    molecule: laidOut,
    style: RenderStyle(showCarbons: false, showImplicitHydrogens: true),
    canvasRect: CGRect(x: 0, y: 0, width: 960, height: 720),
    zoom: 1.0,
    pan: .zero
)

print(scene.backgroundBoxes.count)
print(scene.labels.count)
```

### Export a Highlighted Molecule as V3000

```swift
import CDKSwiftNativePort

let molecule = try CDKSmilesParserFactory.shared
    .newSmilesParser(flavor: .cdkDefault)
    .parseSmiles("C1NCNC1 |ha:0,1,3,hb:2,4|")

let molV3000 = try CDKFileExporter.write(
    molecule: molecule,
    as: .molV3000
)

print(molV3000)
```

### Read and Export a Reaction as CML

```swift
import CDKSwiftNativePort

let reaction = try CDKFileImporter.readReaction(
    text: "CCO>>CC=O",
    fileExtension: "rsmi"
)

let cml = try CDKFileExporter.write(reaction: reaction, as: .cml)
print(cml.prefix(120))
```

### Read and Export a Reaction as RInChI

```swift
import CDKSwiftNativePort

let reaction = try CDKFileImporter.readReaction(
    text: "RInChI=1.00.1S/C7H13BrN2O2/c1-3-7(8,4-2)5(11)10-6(9)12<>C7H14N2O2/c1-3-5(4-2)6(10)9-7(8)11/d+",
    fileExtension: "rinchi"
)

let rinchi = try CDKFileExporter.write(reaction: reaction, as: .rinchi)
print(rinchi)
```

### Preserve a CML Reaction Scheme Hierarchy

```swift
import CDKSwiftNativePort

let hierarchy = try CDKCMLReactionReader.readHierarchy(text: cmlText)

if case .scheme(let scheme) = hierarchy {
    print(scheme.flattenedReactions.count)
}
```

Step lists now preserve nested scheme/list entries as well, so branching and
step grouping survive `readHierarchy` -> `write(_:)` -> `readHierarchy`
round-trips without being collapsed to flat reaction arrays.

## Documentation

- API reference: [`Documentation/API.md`](Documentation/API.md)
- Integration guide: [`Documentation/INTEGRATION.md`](Documentation/INTEGRATION.md)
- Linux guide: [`Documentation/LINUX.md`](Documentation/LINUX.md)
- Architecture and boundary contract: [`Documentation/ARCHITECTURE.md`](Documentation/ARCHITECTURE.md)
- Comparison with original CDK: [`Documentation/CDK_COMPARISON.md`](Documentation/CDK_COMPARISON.md)
- macOS-specific notes: [`Documentation/MACOS.md`](Documentation/MACOS.md)
- Publishing guide: [`PUBLISHING.md`](PUBLISHING.md)
- Attribution and notice: [`NOTICE.md`](NOTICE.md)
- Changelog: [`CHANGELOG.md`](CHANGELOG.md)

## Release Quality Gate

Package-only:

```bash
swift test
```

Before publishing a GitHub release, run the same package suite on at least one
Linux environment as well. The package has been verified on Ubuntu 24.04 with
Swift 6.3.

If you are validating the first-party monorepo as well, also build AtomLens
against the package product:

```bash
xcodebuild -project /path/to/AtomLens.xcodeproj -scheme AtomLens -configuration Debug -sdk macosx build
```

## Upstream Reference and Attribution

`CDKSwiftNativePort` includes code derived from and inspired by the Chemistry
Development Kit (CDK). The current parity target for the supported workflows in
this package is CDK `2.12`.

Upstream project:

- <https://github.com/cdk/cdk>

Reference citations:

- Willighagen et al. (2017), *The Chemistry Development Kit (CDK) v2.0: atom
  typing, depiction, molecular formulas, and substructure searching*.
  DOI: [10.1186/s13321-017-0220-4](https://doi.org/10.1186/s13321-017-0220-4)
- May and Steinbeck (2014), *Efficient ring perception for the Chemistry
  Development Kit*.
  DOI: [10.1186/1758-2946-6-3](https://doi.org/10.1186/1758-2946-6-3)
- Steinbeck et al. (2006), *Recent Developments of the Chemistry Development
  Kit (CDK) - An Open-Source Java Library for Chemo- and Bioinformatics*.
  DOI: [10.2174/138161206777585274](https://doi.org/10.2174/138161206777585274)
- Steinbeck et al. (2003), *The Chemistry Development Kit (CDK): An Open-Source
  Java Library for Chemo- and Bioinformatics*.
  DOI: [10.1021/ci025584y](https://doi.org/10.1021/ci025584y)

For the full upstream attribution and redistribution notice, see
[`NOTICE.md`](NOTICE.md).
