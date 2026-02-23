# CDKSwiftNativePort

`CDKSwiftNativePort` is a Swift-native chemistry package that encapsulates CDK-derived logic behind a reusable API surface for macOS applications.

Repository: `https://github.com/SaschaLosko/CDKSwiftNativePort`

## Current Feature Set

- Core chemical graph model: `Molecule`, `Atom`, `Bond` (+ query/stereo/valence helpers).
- CDK-like 2D layout pipeline: `Depiction2DGenerator` (`StructureDiagramGenerator` port).
- Depiction/export pipeline:
  - SVG output (`CDKDepictionGenerator`)
  - SwiftUI `GraphicsContext` renderer (`CDKStandardGenerator`)
  - Metal-ready scene model (`CDKMetalDepictionSceneBuilder`)
  - label clipping with glyph-outline-aware obstacles.
- SMILES family:
  - SMILES parser/generator
  - CXSMILES split/state handling
  - reaction SMILES parsing (`CDKReaction`).
- InChI family:
  - InChI -> structure
  - structure -> InChI / InChIKey
  - native Swift generation/parsing path.
- IO coverage:
  - unified importer/exporter facades
  - format-specific readers/writers (MDL V2000/V3000, SDF, MOL2, PDB, XYZ, CML, RXN/RDF, SMILES, InChI).
- Descriptor/property services:
  - molecular formula/mass descriptors
  - Rule-of-Five aggregation
  - XLogP and Mannhold LogP.

## Supported Formats

Input (read):
- MDL Molfile / SDFile (`.mol`, `.sdf`, `.sd`)
- SMILES (`.smi`, `.smiles`, `.ism`, `.can`)
- InChI (`.inchi`, `.ich`)
- MOL2 (`.mol2`)
- PDB (`.pdb`, `.ent`)
- XYZ (`.xyz`)
- CML (`.cml`)
- RXN (`.rxn`)
- RDF (`.rdf`)
- plain text auto-detection for supported chemistry payloads.

Output (write):
- MDL Molfile / SDFile
- SMILES / isomeric SMILES
- InChI
- MOL2 / PDB / XYZ / CML
- RXN / RDF
- SVG depiction

## Package Boundary Contract

`CDKSwiftNativePort` intentionally excludes host-app integration concerns:
- no Spotlight indexer implementation
- no Quick Look provider implementation
- no app bundle identifiers, entitlements, or app window/session management

The package is consumed by host apps; app-specific integration remains outside this repository.

## Platform and Toolchain

- Swift tools: `5.9`
- Platform: `macOS 14+`

## Installation (SwiftPM)

In Xcode, add:

```text
https://github.com/SaschaLosko/CDKSwiftNativePort.git
```

or in `Package.swift`:

```swift
.package(url: "https://github.com/SaschaLosko/CDKSwiftNativePort.git", from: "1.0.0")
```

## Quick Start

```swift
import CDKSwiftNativePort

let molecules = try CDKSMILESReader.read(text: "CC(=O)OC1=CC=CC=C1C(=O)O")
guard let molecule = molecules.first else { fatalError("No molecule parsed") }

let laidOut = Depiction2DGenerator.generate(for: molecule)
let ids = CDKMoleculeIdentifierService.compute(for: laidOut)
let props = CDKMoleculePropertyService.compute(for: laidOut)
let svg = CDKDepictionGenerator.toSVG(molecule: laidOut)

print(ids.smiles, ids.inchi, props.molecularWeight, svg.count)
```

## Documentation

- API reference: `Documentation/API.md`
- Architecture: `Documentation/ARCHITECTURE.md`
- Contributing: `CONTRIBUTING.md`
- Security policy: `SECURITY.md`
- Publishing process: `PUBLISHING.md`
- Changelog: `CHANGELOG.md`

## Quality Gates

```bash
swift test
```

The test suite includes package-boundary guards that fail if app-specific coupling markers leak into package sources.

## CDK Attribution and License

This package contains CDK-derived work.

- Upstream CDK: `https://github.com/cdk/cdk`
- Port parity target: CDK `2.11`
- Package license: `LGPL-2.1-or-later` (`LICENSE`)
- Attribution notes: `NOTICE.md`
