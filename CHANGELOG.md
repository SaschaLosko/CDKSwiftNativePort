# Changelog

All notable changes to this project are documented in this file.

The format is based on Keep a Changelog and this project follows Semantic Versioning.

## [Unreleased]

## [1.4.7] - 2026-06-09

### Fixed

- Aligned reaction agent depiction with Java CDK behavior by placing agents in
  the lane directly above the reaction arrow instead of a separate canvas-wide
  top band.

## [1.4.6] - 2026-06-09

### Fixed

- Preserved CML atom parity when generating native InChI identifiers through the
  official libinchi bridge.

## [1.4.5] - 2026-06-01

### Fixed

- Applied a shared structural scale across reaction participants so reactants,
  agents, and products preserve consistent bond lengths and aromatic ring sizes
  within one reaction depiction.

## [1.4.4] - 2026-05-31

### Fixed

- Converted large vendored InChI line buffers from tentative common symbols to
  explicit zero-initialized storage so Apple linkers no longer emit
  `__DATA,__common` alignment warnings in downstream app targets.

## [1.4.3] - 2026-05-31

### Changed

- Removed the package-owned Swift InChI generator fallback and Swift InChIKey
  codec from the runtime path. InChI and component InChIKey generation now go
  through the vendored official IUPAC libinchi bridge.
- Routed RInChI long-key component InChIKey derivation through the official
  libinchi bridge while keeping the RInChI-specific short/web hash algorithm in
  native Swift.

### Fixed

- Preserved two-letter element symbols that are not present in the package's
  descriptor mass table when passing structures to libinchi, fixing reference
  isotope cases such as dysprosium.

## [1.4.2] - 2026-05-30

### Added

- Vendored the official IUPAC InChI `v1.07.5` C sources as a native
  SwiftPM C target (`IUPACInChI`) and routed default InChI/InChIKey
  generation through that in-process library.

### Changed

- Tightened native InChI generation for simple alcohol and carboxylic-acid
  cases so ethanol and acetic acid now produce exact official InChI/InChIKey
  smoke-test outputs.

## [1.4.1] - 2026-05-30

### Fixed

- Fixed acyclic alkene 2D layout so substituted open-chain alkenes keep a more
  conventional zig-zag depiction.

## [1.4.0] - 2026-05-29

### Added

- Native reaction hierarchy parity coverage for reaction sets, reaction lists,
  reaction schemes, CML hierarchy round-trips, RXN V2000/V3000, RDF, and
  reaction SMILES workflows.
- Native `RInChI` support for generation, decomposition,
  `RInChI`-to-reaction reconstruction, `RAuxInfo`, and long/short/web
  `RInChIKey` derivation.
- Official InChI reference parity fixtures and an opt-in strict release gate
  for the maintained native Swift InChI surface.
- Native InChIKey derivation aligned with the official major/minor block split.

### Changed

- Adopted Swift tools 6.0, Swift 6 language mode, and `macOS 15+` as the
  declared Apple deployment baseline.
- Improved molecule and reaction depiction geometry for stereo wedges, ring
  stereo, explicit hydrogens, hidden-carbon defaults, reaction packing, and
  trigonal carbonate branches.
- Refined CML, reaction, and InChI round-tripping so source metadata and
  hierarchy information survive maintained import/export workflows more
  consistently.
- Kept generated reaction gap reports out of source control while preserving
  committed parity inventories used by the test suite.

## [1.3.0] - 2026-03-28

### Added

- Linux Swift Package Manager support for the headless chemistry core, including
  cross-platform geometry compatibility, conditional platform imports, and
  Linux verification coverage.
- RGFile import/export support with round-trip-focused tests.
- A dedicated Linux integration guide with runnable SwiftPM executable examples.

### Changed

- Refined rendering and label-clipping fallbacks so the depiction and
  scene-generation paths behave consistently on Linux as well as Apple
  platforms.
- Updated the GitHub-facing documentation set to present the package as a
  standalone cross-platform release, including Linux installation guidance and
  upstream CDK attribution/citations.
- Strengthened package-boundary checks to keep package-source scanning robust on
  Linux and to preserve the no-AtomLens-coupling contract.

## [1.2.0] - 2026-03-14

### Added

- Expanded CDK 2.12-style CXSMILES support, including highlights, richer
  Sgroups, ligand ordering, and Markush/link-node handling across parsing,
  layout, SVG, and scene-generation paths.
- V2000/V3000/SDF parity work for import/export, including round-trip-focused
  coverage and V3000 `HILITE` support.
- CML molecule and reaction import/export support aligned with the current
  package model and tested for round-trip workflows.
- Rendering support for non-Markush Sgroups, query atoms/bonds, and additional
  highlight styles backed by the newer format support.
- A package-boundary test that verifies the host Xcode project links the
  `CDKSwiftNativePort` package product instead of compiling package source files
  directly.

### Changed

- Refreshed the GitHub-facing documentation set (`README.md`,
  `Documentation/API.md`, `Documentation/INTEGRATION.md`,
  `Documentation/ARCHITECTURE.md`, `Documentation/CDK_COMPARISON.md`,
  `Documentation/MACOS.md`, `PUBLISHING.md`) to reflect the current public API
  surface, host-boundary contract, and macOS integration model.
- Updated attribution notes to reflect CDK `2.12` as the current parity target
  for supported workflows.

## [1.1.1] - 2026-03-09

### Changed

- Updated installation snippets in `README.md` and `Documentation/INTEGRATION.md` to reference the current package release line.
- Refined `PUBLISHING.md` to use a reusable `VERSION` variable in release commands and to document `1.1.1` as the current consumer pin target.

### Notes

- This is a documentation-only follow-up release to keep published guidance aligned with the latest tag.

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
