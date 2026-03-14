# CDKSwiftNativePort

`CDKSwiftNativePort` is a Swift-native chemistry package that keeps
CDK-derived parsing, layout, depiction, identifiers, descriptors, and file IO
inside a reusable package boundary.

Repository: <https://github.com/SaschaLosko/CDKSwiftNativePort>

## Overview

This package exists for two reasons:

- keep all CDK-derived chemistry logic out of the host app
- expose a Swift-first API for native macOS chemistry workflows

The first-party consumer is AtomLens, but the package itself has no dependency
on AtomLens modules, Quick Look targets, Spotlight targets, bundle identifiers,
or app session logic.

## Feature Set

### Core chemistry model

- `Molecule`, `Atom`, `Bond`, `MoleculeSgroup`
- stereochemistry, atom maps, aliases, attachment points
- query atoms and query bonds
- SDF-style data fields and Sgroup metadata

### SMILES and CXSMILES

- SMILES parsing and generation
- reaction SMILES parsing and generation
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
- InChI
- MOL2
- PDB
- XYZ
- CML molecules
- CML reactions
- RXN
- RDF
- SVG depiction export

### Identifiers and descriptors

- SMILES, isomeric SMILES, InChI, InChIKey
- formula, molecular weight, exact mass
- heavy atom count, donor/acceptor counts, ring count, rotatable bonds
- XLogP, Mannhold LogP, TPSA, van der Waals volume
- Rule-of-Five summary

## Encapsulation and App Boundary

All CDK-derived implementation lives under `Sources/CDKSwiftNativePort`.

The package owns:

- molecule and reaction models
- readers, writers, parsers, generators
- layout and depiction algorithms
- scene building for host renderers
- identifier and descriptor services

The host app owns:

- window and document behavior
- UI state and toolbar controls
- `MTKView` / Core Graphics drawing loops
- Quick Look and Spotlight extension targets
- bundle identifiers, entitlements, and App Store metadata

This boundary is enforced in tests:

- package sources are scanned for forbidden host-coupling imports
- `Package.swift` is checked for AtomLens target references
- the monorepo host chemistry layer is checked to ensure it only contains the
  thin alias adapter
- the AtomLens Xcode project is checked to ensure it links the package product
  rather than compiling package implementation files directly

## Platform

- Swift tools: `5.9`
- declared platform: `macOS 14+`

The chemistry core is headless and does not depend on `AppKit` or `SwiftUI`.
The package uses Foundation/CoreGraphics/CoreText-level APIs and is well suited
to app, Quick Look, and metadata-importer targets that want to share the same
chemistry core.

## Installation

### Xcode

Add the package dependency:

```text
https://github.com/SaschaLosko/CDKSwiftNativePort.git
```

### `Package.swift`

```swift
dependencies: [
    .package(url: "https://github.com/SaschaLosko/CDKSwiftNativePort.git", from: "1.2.0")
]
```

## Quick Start

### Parse, layout, identify, and export SVG

```swift
import CDKSwiftNativePort
import CoreGraphics

let molecules = try CDKFileImporter.readMolecules(
    text: "CC(=O)OC1=CC=CC=C1C(=O)O",
    fileExtension: "smi"
)

guard let molecule = molecules.first else {
    throw ChemError.emptyInput
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

### Parse a Markush CXSMILES and build a scene

```swift
import CDKSwiftNativePort
import CoreGraphics

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

### Export a highlighted molecule as V3000

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

### Read and export a reaction as CML

```swift
import CDKSwiftNativePort

let reaction = try CDKFileImporter.readReaction(
    text: "CCO>>CC=O",
    fileExtension: "rsmi"
)

let cml = try CDKFileExporter.write(reaction: reaction, as: .cml)
print(cml.prefix(120))
```

## Documentation

- API reference: [`Documentation/API.md`](Documentation/API.md)
- Integration guide: [`Documentation/INTEGRATION.md`](Documentation/INTEGRATION.md)
- Architecture and boundary contract: [`Documentation/ARCHITECTURE.md`](Documentation/ARCHITECTURE.md)
- Comparison with original CDK: [`Documentation/CDK_COMPARISON.md`](Documentation/CDK_COMPARISON.md)
- macOS-specific notes: [`Documentation/MACOS.md`](Documentation/MACOS.md)
- Publishing guide: [`PUBLISHING.md`](PUBLISHING.md)
- Attribution and notice: [`NOTICE.md`](NOTICE.md)
- Changelog: [`CHANGELOG.md`](CHANGELOG.md)

## Quality Gate

Package-only:

```bash
swift test
```

Monorepo verification before shipping AtomLens:

```bash
xcodebuild -project /path/to/AtomLens.xcodeproj -scheme AtomLens -configuration Debug -sdk macosx build
```

## Upstream Reference

The current parity target for the supported workflows in this package is
CDK `2.12`.

This is not a full port of every CDK Java module. See
[`Documentation/CDK_COMPARISON.md`](Documentation/CDK_COMPARISON.md) for the
current parity surface and known scope limits.

## License and Attribution

This repository contains CDK-derived work.

- upstream CDK: <https://github.com/cdk/cdk>
- reference parity target: CDK `2.12`
- license: `LGPL-2.1-or-later`
- attribution notes: [`NOTICE.md`](NOTICE.md)
