# CDKSwiftNativePort

`CDKSwiftNativePort` is a Swift-native package that encapsulates CDK-derived chemistry functionality and keeps it reusable across host applications.

Repository: `https://github.com/SaschaLosko/CDKSwiftNativePort`

## What It Provides

- Core molecular graph model (`Molecule`, `Atom`, `Bond`) and graph/ring utilities.
- CDK-style 2D coordinate generation and depiction pipelines.
- SMILES/CXSMILES/reaction SMILES parsing and generation.
- InChI parsing + native Swift InChI/InChIKey generation path.
- Unified multi-format import/export APIs.
- Descriptor/property services including XLogP and Rule-of-Five support.

## Supported Formats

Input (read):
- MDL Molfile / SDFile (`.mol`, `.sdf`, `.sd`)
- SMILES (`.smi`, `.smiles`)
- InChI (`.inchi`, `.ich`)
- MOL2 (`.mol2`)
- PDB (`.pdb`)
- XYZ (`.xyz`)
- CML (`.cml`, `.xml`)
- RXN (`.rxn`)
- RDF (`.rdf`)

Output (write):
- MDL Molfile / SDFile
- SMILES / ISO SMILES
- InChI
- MOL2, PDB, XYZ, CML
- RXN, RDF
- SVG depiction export

## Package Boundary Contract

The package intentionally does **not** depend on host-app integration layers, including:
- Spotlight APIs (`CoreSpotlight`)
- Quick Look APIs (`QuickLook`, `QuickLookThumbnailing`)
- App-level identifiers, entitlements, or sandbox/bookmark management
- Host UI state and window lifecycle logic

These concerns stay in the consuming app.

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

The test suite includes package-boundary checks that fail if app-specific coupling leaks into package sources.

## CDK Attribution and License

This package contains CDK-derived work.

- Upstream CDK: `https://github.com/cdk/cdk`
- Port parity target: CDK `2.11`
- Package license: `LGPL-2.1-or-later` (`LICENSE`)
- Attribution notes: `NOTICE.md`
