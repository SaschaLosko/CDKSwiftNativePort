# CDKSwiftNativePort

`CDKSwiftNativePort` is a Swift-native chemistry package that encapsulates CDK-derived functionality behind a reusable API for macOS hosts.

Repository: <https://github.com/SaschaLosko/CDKSwiftNativePort>

## Why This Package Exists

- Keep CDK-derived chemistry logic in one isolated, reusable package.
- Let host apps (for example AtomLens) consume a stable API instead of embedding chemistry-port internals.
- Preserve a clean boundary: package code has no dependency on host app modules, UI targets, Spotlight, Quick Look, or app-specific identifiers.

## Feature Overview

- Core model: `Molecule`, `Atom`, `Bond`, stereochemistry and query metadata.
- Layout/depiction:
  - `Depiction2DGenerator` (CDK-style 2D layout)
  - `CDKDepictionGenerator` (SVG output)
  - `CDKStandardGenerator` (SwiftUI `GraphicsContext` drawing)
  - `CDKMetalDepictionSceneBuilder` and `CDKMetalReactionDepictionSceneBuilder` (renderer-agnostic scene primitives for Metal hosts)
- Parsing/import:
  - MDL Molfile/SDF, SMILES and reaction SMILES, InChI, MOL2, PDB, XYZ, CML, RXN, RDF
  - unified `CDKFileImporter`
- Export:
  - Molfile, SDF, SMILES/isomeric SMILES, InChI, MOL2, PDB, XYZ, CML, RXN, RDF, SVG
  - unified `CDKFileExporter`
- Identifier/property services:
  - SMILES, isomeric SMILES, InChI, InChIKey (`CDKMoleculeIdentifierService`)
  - molecular formula/mass, H-bond counts, rotatable bonds, rings, Rule-of-Five, XLogP (`CDKMoleculePropertyService`)

## Platform and Toolchain

- Swift tools: `5.9`
- Platform: `macOS 14+`

## Installation

### Xcode

Add package dependency:

```text
https://github.com/SaschaLosko/CDKSwiftNativePort.git
```

### `Package.swift`

```swift
.package(url: "https://github.com/SaschaLosko/CDKSwiftNativePort.git", from: "1.0.0")
```

## Quick Start

```swift
import CDKSwiftNativePort

let molecules = try CDKFileImporter.readMolecules(
    text: "CC(=O)OC1=CC=CC=C1C(=O)O",
    fileExtension: "smi"
)
guard let parsed = molecules.first else { fatalError("No molecule parsed") }

let laidOut = Depiction2DGenerator.generate(for: parsed)
let ids = CDKMoleculeIdentifierService.compute(for: laidOut)
let props = CDKMoleculePropertyService.compute(for: laidOut)
let svg = CDKDepictionGenerator.toSVG(molecule: laidOut)

print(ids.inchiKey)
print(props.formula)
print(svg.prefix(80))
```

### Reaction + Metal Scene (Host Rendering)

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
    canvasRect: CGRect(x: 0, y: 0, width: 960, height: 540),
    zoom: 1.0,
    pan: .zero
)

print(scene.bondSegments.count)
```

## Boundary Contract

`CDKSwiftNativePort` intentionally excludes host-app integration concerns:

- no Spotlight indexer implementation
- no Quick Look provider implementation
- no app bundle identifiers/entitlements/window/session logic

Boundary checks are test-enforced (`PackageBoundaryTests`).

## Documentation

- API reference: `Documentation/API.md`
- Integration guide: `Documentation/INTEGRATION.md`
- Architecture: `Documentation/ARCHITECTURE.md`
- CDK comparison: `Documentation/CDK_COMPARISON.md`
- macOS notes: `Documentation/MACOS.md`
- Contributing: `CONTRIBUTING.md`
- Publishing: `PUBLISHING.md`
- Changelog: `CHANGELOG.md`

## CDK Comparison

This package targets practical parity for the chemistry workflows used by AtomLens, not full one-to-one coverage of all Java CDK modules.

- See `Documentation/CDK_COMPARISON.md` for a detailed comparison and scope notes.

## macOS Notes

- See `Documentation/MACOS.md` for sandbox/security-scoped IO behavior, rendering integration patterns, and package/host separation guidance.

## Quality Gate

```bash
swift test
```

## Attribution and License

This package contains CDK-derived work.

- Upstream CDK: <https://github.com/cdk/cdk>
- Port parity target: CDK `2.11`
- Package license: `LGPL-2.1-or-later` (`LICENSE`)
- Attribution notes: `NOTICE.md`
