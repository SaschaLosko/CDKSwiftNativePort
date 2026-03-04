# Changelog

All notable changes to this project are documented in this file.

The format is based on Keep a Changelog and this project follows Semantic Versioning.

## [Unreleased]

## [1.1.0] - 2026-03-04

### Added

- Reaction model value types (`CDKReactionRole`, `CDKReactionParticipant`) and participant-preserving reaction parsing across SMILES/RXN/RDF workflows.
- Metal reaction depiction scene builder and participant hit-testing API for host-side interaction (`CDKMetalReactionDepictionSceneBuilder`).
- Extended importer coverage and tests for reaction parsing edge cases, stoichiometry handling, and RXN/RDF participant metadata.
- Additional rendering test coverage for aromatic display and reaction visualization consistency.
- New documentation set for GitHub publishing:
  - `Documentation/INTEGRATION.md`
  - `Documentation/CDK_COMPARISON.md`
  - `Documentation/MACOS.md`

### Changed

- Updated public documentation (`README.md`, `Documentation/API.md`, `Documentation/ARCHITECTURE.md`) to reflect the current API surface and package boundary contract.
- Strengthened package boundary guard tests to catch host-app coupling markers and accidental workspace-level chemistry leakage.
- Refined rendering/layout internals and label clipping behavior used by SVG and Metal scene outputs.
- Updated installation and publishing references to the canonical GitHub repository URL.

### Notes

- Package boundary policy remains enforced by tests (`PackageBoundaryTests`) to keep CDK-derived code app-agnostic.

## [1.0.0] - 2026-02-19

### Added

- Unified import/export dispatch APIs (`CDKFileImporter`, `CDKFileExporter`) across MDL, SMILES, InChI, MOL2, PDB, XYZ, CML, RXN, RDF, SVG.
- CDK-derived layout/depiction pipeline, including scene generation and style controls.
- Native Swift InChI parse/generate path and InChIKey generation service.
- Molecular property service and descriptor set including XLogP and Rule-of-Five helpers.
- Expanded CDK parity-oriented test suite and metadata tests.
- Package boundary guard tests to prevent app-level coupling from entering package sources.
- Release documentation set (`Documentation/API.md`, `CONTRIBUTING.md`, `SECURITY.md`).

### Changed

- Removed package source coupling to app-specific identifiers and app-branded defaults.
- Replaced app-scoped content-type identifiers in package IO metadata with package-neutral/common identifiers.
- Cleaned renderer dependency surface by removing `AppKit` import from package rendering implementation.

### Notes

- This package contains CDK-derived functionality and remains licensed under LGPL-2.1-or-later.
