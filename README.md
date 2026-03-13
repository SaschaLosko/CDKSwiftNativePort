# CDKSwiftNativePort

`CDKSwiftNativePort` is a Swift-native chemistry package that encapsulates
CDK-derived parsing, layout, depiction, identifier, IO, and selected descriptor
work behind a reusable library boundary.

Repository: <https://github.com/SaschaLosko/CDKSwiftNativePort>

## What This Package Is For

- Keep all CDK-derived chemistry implementation inside one standalone package.
- Let host apps link a package product instead of embedding chemistry-port source.
- Preserve a clean host boundary: package code has no dependency on AtomLens app
  modules, app extension targets, app bundle identifiers, or app session logic.
- Provide a native Swift API for the workflows most often needed in macOS
  chemistry apps: import, layout, depiction, identifiers, properties, and export.

The first-party consumer is AtomLens, but the package is intentionally authored
as a host-agnostic library.

## Feature Highlights

- Core chemistry model:
  - `Molecule`, `Atom`, `Bond`, `MoleculeSgroup`
  - stereochemistry, query atoms/bonds, atom maps, attachment points
  - data fields for SDF-style metadata
- SMILES and CXSMILES:
  - CDK-style SMILES parser/generator
  - reaction SMILES support
  - CXSMILES state parsing
  - CDK 2.12-style Markush / R-group support, including `RG:` and `LN:` layers
- Layout and depiction:
  - `Depiction2DGenerator` for CDK-style 2D coordinates
  - SVG export with `CDKDepictionGenerator`
  - SwiftUI `GraphicsContext` rendering with `CDKStandardGenerator`
  - renderer-agnostic Metal scene primitives with `CDKMetalDepictionSceneBuilder`
    and `CDKMetalReactionDepictionSceneBuilder`
- File IO:
  - MDL Molfile / SDF / RXN / RDF
  - SMILES / reaction SMILES
  - InChI
  - MOL2, PDB, XYZ, CML
  - unified import/export facades (`CDKFileImporter`, `CDKFileExporter`)
- Identifiers and descriptors:
  - SMILES, isomeric SMILES, InChI, InChIKey
  - formula, molecular weight, exact mass
  - H-bond donors/acceptors, heavy atoms, rings, rotatable bonds
  - XLogP, TPSA, van der Waals volume, Rule-of-Five summary

## Encapsulation and Host Boundary

All CDK-derived implementation lives under `Sources/CDKSwiftNativePort`.

What stays inside the package:

- parsers, generators, readers, writers
- CDK-derived 2D layout and depiction logic
- descriptor and identifier implementations
- scene-building logic for host renderers

What stays outside the package:

- Spotlight importers and indexers
- Quick Look preview / thumbnail extensions
- `MTKView` ownership, AppKit / SwiftUI app wiring, window state, entitlements
- App Store metadata, bundle IDs, sandbox policy, and host-session logic

This is enforced in tests:

- package sources are scanned for forbidden AtomLens and app-extension coupling
- the package manifest is checked for host-target references
- the monorepo host chemistry layer is checked to ensure it only contains a thin
  alias adapter
- the host Xcode project is checked to ensure it links the package product rather
  than compiling package source files directly

## Platform and Toolchain

- Swift tools: `5.9`
- Declared platform: `macOS 14+`

The package is currently optimized for native macOS hosts. Core chemistry APIs
are Swift/Apple-framework based, while the rendering surface includes a SwiftUI
renderer and Metal scene builders for host-owned rendering pipelines.

## Installation

### Xcode

Add the package dependency:

```text
https://github.com/SaschaLosko/CDKSwiftNativePort.git
```

### `Package.swift`

```swift
dependencies: [
    .package(url: "https://github.com/SaschaLosko/CDKSwiftNativePort.git", from: "1.1.1")
]
```

If you need functionality that is newer than the latest tag, pin a branch or
revision until the next release is cut.

## Quick Start

### Parse, Layout, Identify, and Export

```swift
import CDKSwiftNativePort

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
let svg = CDKDepictionGenerator.toSVG(molecule: laidOut)

print(identifiers.smiles)
print(identifiers.inchiKey)
print(properties.formula)
print(svg.prefix(120))
```

### Parse a Markush CXSMILES

```swift
import CDKSwiftNativePort
import CoreGraphics

let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
let markush = try parser.parseSmiles(
    "C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|"
)

let laidOut = Depiction2DGenerator.generate(for: markush)
let scene = CDKMetalDepictionSceneBuilder.build(
    molecule: laidOut,
    style: RenderStyle(showCarbons: false, showImplicitHydrogens: true),
    canvasRect: CGRect(x: 0, y: 0, width: 960, height: 720),
    zoom: 1.0,
    pan: .zero
)

print(scene.backgroundBoxes.count) // Markush legend boxes
```

### Build a Reaction Scene for a Host Renderer

```swift
import CDKSwiftNativePort
import CoreGraphics

let reaction = try CDKFileImporter.readReaction(
    text: "CCO>>CC=O",
    fileExtension: "rsmi"
)

let scene = CDKMetalReactionDepictionSceneBuilder.build(
    reaction: reaction,
    style: RenderStyle(),
    canvasRect: CGRect(x: 0, y: 0, width: 1200, height: 680),
    zoom: 1.0,
    pan: .zero
)

print(scene.bondSegments.count)
print(scene.labels.count)
```

## Documentation

- API reference: [`Documentation/API.md`](Documentation/API.md)
- Integration guide: [`Documentation/INTEGRATION.md`](Documentation/INTEGRATION.md)
- Architecture and boundary contract: [`Documentation/ARCHITECTURE.md`](Documentation/ARCHITECTURE.md)
- Comparison with original CDK: [`Documentation/CDK_COMPARISON.md`](Documentation/CDK_COMPARISON.md)
- macOS-specific notes: [`Documentation/MACOS.md`](Documentation/MACOS.md)
- Contributing: [`CONTRIBUTING.md`](CONTRIBUTING.md)
- Publishing: [`PUBLISHING.md`](PUBLISHING.md)
- Changelog: [`CHANGELOG.md`](CHANGELOG.md)

## Quality Gate

```bash
swift test
```

## Upstream Reference and Scope

This package aims at practical parity for the workflows currently implemented in
the Swift port, not one-to-one coverage of all Java CDK modules. The current
reference target for parity work is CDK `2.12`.

See [`Documentation/CDK_COMPARISON.md`](Documentation/CDK_COMPARISON.md) for a
feature-by-feature comparison and known scope limits.

## Attribution and License

This package contains CDK-derived work.

- Upstream CDK: <https://github.com/cdk/cdk>
- Reference parity target: CDK `2.12`
- Package license: `LGPL-2.1-or-later` (`LICENSE`)
- Attribution notes: [`NOTICE.md`](NOTICE.md)

Reference citations:

- Willighagen et al. (2017), The Chemistry Development Kit (CDK) v2.0:
  atom typing, depiction, molecular formulas, and substructure searching.
  DOI: [10.1186/s13321-017-0220-4](https://doi.org/10.1186/s13321-017-0220-4)
- May and Steinbeck (2014), Efficient ring perception for the Chemistry
  Development Kit. DOI: [10.1186/1758-2946-6-3](https://doi.org/10.1186/1758-2946-6-3)
- Steinbeck et al. (2006), Recent Developments of the Chemistry Development
  Kit (CDK) - An Open-Source Java Library for Chemo- and Bioinformatics.
  DOI: [10.2174/138161206777585274](https://doi.org/10.2174/138161206777585274)
- Steinbeck et al. (2003), The Chemistry Development Kit (CDK): An
  Open-Source Java Library for Chemo- and Bioinformatics.
  DOI: [10.1021/ci025584y](https://doi.org/10.1021/ci025584y)
