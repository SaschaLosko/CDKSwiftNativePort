# CDKSwiftNativePort

`CDKSwiftNativePort` is a Swift Package that encapsulates CDK-derived chemistry functionality in native Swift.

It provides:
- molecule data model and graph utilities
- import/export for multiple chemistry formats
- SMILES and InChI parsing/generation paths
- CDK-style 2D layout + depiction scene generation
- molecular property and descriptor calculation (including XLogP and Rule-of-Five helpers)

## Design Goals

- Keep CDK-derived chemistry logic in one reusable package.
- Keep package boundaries clean: no dependency on app-level code (Spotlight, Quick Look, app bundle identifiers, window/UI app state).
- Preserve CDK-inspired organization while staying idiomatic in Swift.

## Platform and Toolchain

- Swift tools: `5.9`
- Platform: `macOS 14+` (as declared in `Package.swift`)

## Installation (SwiftPM)

In Xcode, add a package dependency:

```text
https://github.com/<your-org-or-user>/CDKSwiftNativePort.git
```

or in `Package.swift`:

```swift
.package(url: "https://github.com/<your-org-or-user>/CDKSwiftNativePort.git", from: "1.0.0")
```

## Quick Start

```swift
import CDKSwiftNativePort

let input = "CC(=O)OC1=CC=CC=C1C(=O)O"
let molecules = try CDKSMILESReader.read(text: input)
guard let first = molecules.first else { fatalError("No molecule") }

// 2D layout
let laidOut = Depiction2DGenerator.generate(for: first)

// Identifiers
let ids = CDKMoleculeIdentifierService.compute(for: laidOut)
print(ids.smiles)
print(ids.inchi)

// Properties
let props = CDKMoleculePropertyService.compute(for: laidOut)
print(props.formula, props.molecularWeight, props.ruleOfFive.statusText)

// Export SVG
let svg = CDKDepictionGenerator.toSVG(molecule: laidOut)
```

## Import and Export APIs

- Generic importer: `CDKFileImporter`
- Generic exporter: `CDKFileExporter`

Supported families include:
- MDL Molfile / SDFile (`.mol`, `.sdf`, `.sd`)
- SMILES / isomeric SMILES
- InChI
- MOL2, PDB, XYZ, CML
- RXN / RDF
- SVG depiction export

## API Documentation

- High-level API reference: `Documentation/API.md`
- Architecture notes: `Documentation/ARCHITECTURE.md`
- Contribution guide: `CONTRIBUTING.md`
- Security policy: `SECURITY.md`
- Publishing/release process: `PUBLISHING.md`
- Changelog: `CHANGELOG.md`

## Quality Gates

- Run tests:

```bash
swift test
```

- The suite includes package-boundary checks that prevent app-specific coupling from leaking into package source code.

## Repository Standards

This repository includes:
- CI workflow (`.github/workflows/ci.yml`)
- Pull request template (`.github/pull_request_template.md`)
- Issue templates (`.github/ISSUE_TEMPLATE/`)
- Changelog-based release notes (`CHANGELOG.md`)

## CDK Attribution and License

This package contains CDK-derived work.

- Upstream project: `https://github.com/cdk/cdk`
- Reference parity target used in this porting effort: CDK `2.11`

Licensing:
- Package license: `LGPL-2.1-or-later` (see `LICENSE`)
- Additional attribution notes: `NOTICE.md`
