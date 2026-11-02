# Changelog

All notable changes to this project are documented in this file.

The format is based on Keep a Changelog and this project follows Semantic Versioning.

## [Unreleased]

## [1.4.38] - 2026-11-02

### Changed

- Separate atom, bond, shared chemistry types and depiction entry points from the
  molecule model. Split structure-diagram generation into graph analysis and
  focused layout stages without changing public APIs or chemistry behavior.
- Normalize maintained Swift source and tests with the checked-in formatter;
  full strict Swift lint passes. Preserve the official native sources unchanged.

### Added

- Add bounded finite family overlap, containment, equivalence and disjointness,
  with shared budgets, witnesses/counterexamples and original atom provenance.
- Add reviewed finite definition sidecars, nested choices, explicit exclusions,
  scaffold OR/ABS stereo, neutral imine stereo and explicit dummy position variation.
- Retain the 103 M2S annotations and qualify five finite interpretations; unknown
  source conditions, mixtures and ambiguous attachment semantics remain unresolved.
- Check expanded and retained atom limits before allocation and preserve omission
  and repeated-atom provenance.

- Add finite `CDKMarkushMatcher` membership with match/no-match/unresolved results,
  concrete witnesses, atom mappings, traversal completeness and cooperative bounds.
- Support numeric variables beyond R32 and explicit named-pseudoatom resolution;
  retain the R0-R32 boundary for RGfile output.
- Add focused tetrahedral/alkene, REST_H, repeat, variable-label and malformed-input
  regressions. Correct CX fragment ligand-ID remapping and preserve atom properties.
- Make package checks and documentation fully standalone, remove host-checkout
  assumptions, and derive reference-import output paths from the package location.
- Remove the development Java reference harness and Maven configuration; retain
  attributed reference data without a Java source or runtime dependency.

- Added additional 3D atom color palettes for host renderers: Okabe-Ito,
  Viridis, Cividis, Magma, Inferno, Plasma, ColorBrewer Set2/Dark2,
  CARTO Safe/Vivid, and Matplotlib Tab10/Tab20.
- Added `CDK3DRepresentationMode` with a space-filling/CPK representation that
  uses van der Waals atom radii and omits bond cylinders.

## [1.4.37] - 2026-11-01

- Use project-neutral terminology in the Markush regression documentation.
- Documentation-only update; runtime code and tests are unchanged from 1.4.36.

## [1.4.36] - 2026-11-01

- Add bounded scaffold repeat/link-node expansion combined with R-group enumeration,
  including streaming, cancellation, output/atom limits and explicit unresolved errors.
- Correct CX LN neighbor-index interpretation; retain the earlier colon/bond-index variant.
- Parse CX LOG, record rule provenance, and preserve alternatives, attachments and
  conditions in supported flat-family CXSMILES output. Reject unsafe exports.
- Preserve named wildcard provenance; all 103 M2S regression backbones are classified
  as variable, including 12 cases previously exposed as specified compounds.
- Verify nine distinct molecules from the original ring/repeat/R1 example and all
  twelve identities in the assessment's RGfile-to-CX round trip.
- Full macOS/Linux suites: 613 tests each, two existing skips, zero failures.

## [1.4.35] - 2026-11-01

### Added

- Native CDK 2.13 R-group configuration enumeration with occurrence rules,
  positional assignments, dependencies, ordered attachments, REST_H and hydrogen
  handling, deterministic streaming, cancellation and configurable output limits.
- Java-derived regression fixtures and independently generated graph/REST_H
  histogram references covering 17,866 configurations. Repeat units and nested
  definitions report unsupported features instead of losing their semantics.

### Fixed

- Avoid overflow and invalid ranges in R-group occurrence matching.
- Preserve atom properties and stereo ligand references when flattening R-group
  definitions for round trips.

## [1.4.34] - 2026-10-11

### Fixed

- Reconstruct neutral sulfur oxo bonds from InChI connectivity and hydrogen
  assignments using compatible sulfur valences. Sulfuric/sulfurous acid and
  sulfoxides retain their formulas through CML serialization; sulfides retain
  divalent sulfur. The bundled RInChI esterification agent remains H2O4S.
- Regression checks include uncached native InChI regeneration and the full
  RInChI → CML reaction round trip.

## [1.4.32] - 2026-09-14

### Fixed

- Preserve reaction and participant display names in CML even when the name is
  also used as the generated identifier.
- Keep distinct participants with colliding identifiers separate while retaining
  references to identical molecules shared across reaction lists and schemes.
- Verified with 598 tests on macOS and Linux (two existing skips on each).

## [1.4.31] - 2026-09-14

### Fixed

- Reject reaction-level CXSMILES annotations and unprojected molecular variation,
  link-node, R-group, and Sgroup data instead of silently omitting their chemistry
  from ChemDraw exports.

## [1.4.30] - 2026-09-14

### Added

- Native Swift CDX and CDXML export for 2D molecules and reactions, with binary
  Data APIs, Unicode titles, stereo, multi-page layout, and reaction atom maps.
- Explicit unsupported-feature errors and independent Indigo read-back tooling;
  no additional runtime dependency.

### Fixed

- Preserve SMILES neighbor encounter order, including ring-closure positions and
  virtual hydrogen ligands, for parity-correct ChemDraw wedges. Mark generated
  SMILES coordinates so unspecified double bonds remain unspecified on export.

## [1.4.29] - 2026-09-05

### Fixed

- Replaced the VABC descriptor's exhaustive simple-cycle fallback with
  shortest-path-tree cycle candidates and independent minimum-basis selection.
  Dense or fused ring systems no longer enumerate every cycle up to 24 atoms,
  avoiding runaway CPU and memory use during molecular property calculation.
- Verified basis rank, total length, connectivity, macrocycles, and explicit
  hydrogen exclusion, including all 1,024 five-atom graphs and a 70-atom
  disconnected stress graph. Existing VABC reference values remain unchanged.

## [1.4.27] - 2026-08-22

### Added

- Promoted the maintained upstream parity baseline to CDK 2.13, including
  CXSMILES coordinate bonds, MDL V2000/V3000 reacting-center metadata,
  R0/partial R-group handling, small-ring queries, fixed-H AuxInfo numbering,
  and coordinate-bonded sandwich-ring layout.
- Added CDK 2.13 regression coverage for adjacent stereocenters, double-bond
  hydrogen suppression, degenerate double-bond geometry, and heavy-atom counts
  that exclude pseudo atoms.

## [1.4.26] - 2026-08-15

### Fixed

- Preserved explicit MDL V2000 atom valence fields even when they equal the
  element's default valence, preventing hydrogen changes after metal
  disconnection during Standard and Fixed-H InChI generation.

## [1.4.23] - 2026-08-11

### Changed

- Declared `CDKSwiftNativePort` as an explicitly dynamic SwiftPM library product
  so host applications do not absorb the package implementation into their
  executables.

## [1.4.12] - 2026-06-14

### Added

- Added a CDK-style 3D renderer model and Metal-ready 3D scene extraction for
  molecule atom spheres and bond cylinders.

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
  Linux and to preserve the standalone-package contract.

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
