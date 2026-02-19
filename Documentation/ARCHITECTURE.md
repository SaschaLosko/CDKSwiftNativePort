# Architecture Overview

## Top-level Layout

- `Sources/CDKSwiftNativePort/Molecule.swift`
  - core domain model and shared error types.
- `Sources/CDKSwiftNativePort/CDK/Common`
  - service facades and shared utilities.
- `Sources/CDKSwiftNativePort/CDK/Smiles`
  - SMILES/CXSMILES/reaction parse + generation.
- `Sources/CDKSwiftNativePort/CDK/InChI`
  - InChI parse/generate services.
- `Sources/CDKSwiftNativePort/CDK/Layout`
  - structure diagram generation and placement heuristics.
- `Sources/CDKSwiftNativePort/CDK/Rendering`
  - depiction, style resolution, clipping, scene generation.
- `Sources/CDKSwiftNativePort/CDK/IO`
  - format-specific readers/writers and unified importer/exporter APIs.
- `Sources/CDKSwiftNativePort/CDK/QSAR`
  - descriptors and molecular property calculations.

## Boundary Principles

`CDKSwiftNativePort` is intentionally independent from any host app.

Not allowed in package source:
- app bundle identifiers/names
- Spotlight/Quick Look integration logic
- app window/session state management

App integrations belong in host projects and should consume this package through public APIs.

## Testing Strategy

- Unit tests grouped by chemistry area (`Smiles`, `InChI`, `MDL`, `Rendering`, `Layout`, `QSAR`).
- Parity metadata tests track links to upstream CDK tests.
- Round-trip tests validate IO reader/writer consistency.
- Boundary guard tests detect app-level coupling markers in source files.
