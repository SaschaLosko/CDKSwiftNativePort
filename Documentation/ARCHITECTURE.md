# Architecture Overview

## High-Level Goal

`CDKSwiftNativePort` exists to keep CDK-derived chemistry functionality inside a
standalone Swift package while exposing a stable, host-consumable API.

The package owns:

- chemistry model types
- parsing and generation
- layout and depiction logic
- readers and writers
- identifiers and descriptors

The host app owns:

- UI composition
- document and window lifecycle
- Metal / Core Graphics drawing loops
- Spotlight and Quick Look targets
- sandbox policy and App Store-facing app behavior

## Top-Level Source Layout

- `Sources/CDKSwiftNativePort/Molecule.swift`
  - domain model, graph helpers, shared errors, layout entry point
- `Sources/CDKSwiftNativePort/CDK/Common`
  - identifier and property service facades
- `Sources/CDKSwiftNativePort/CDK/Smiles`
  - SMILES / CXSMILES / reaction parser and generator stack
- `Sources/CDKSwiftNativePort/CDK/InChI`
  - InChI parse / generate services
- `Sources/CDKSwiftNativePort/CDK/Layout`
  - CDK-derived structure diagram generation
- `Sources/CDKSwiftNativePort/CDK/Rendering`
  - depiction preprocessing, SVG generation, SwiftUI rendering, Metal scene generation
- `Sources/CDKSwiftNativePort/CDK/IO`
  - format readers / writers and unified importer / exporter facades
- `Sources/CDKSwiftNativePort/CDK/QSAR`
  - descriptor implementations and molecular property calculations

## Public API Design

The package intentionally favors a small number of host-facing facade types:

- `CDKFileImporter`
- `CDKFileExporter`
- `Depiction2DGenerator`
- `CDKDepictionGenerator`
- `CDKMetalDepictionSceneBuilder`
- `CDKMetalReactionDepictionSceneBuilder`
- `CDKMoleculeIdentifierService`
- `CDKMoleculePropertyService`

This keeps host code straightforward while allowing the internal port to evolve
without exposing every implementation detail.

## Boundary Contract

### Not allowed in package sources

- app module imports such as `AtomLens`
- Quick Look or Spotlight implementation code
- host bundle identifiers or entitlements
- window, session, document, or app-restoration logic

### Allowed in host apps

- importing `CDKSwiftNativePort`
- turning package scene data into host-specific rendering
- building metadata importers, Quick Look providers, or app-document logic on
  top of the package

## Monorepo Encapsulation Model

In the first-party monorepo:

- all CDK-derived implementation lives in `CDKSwiftNativePort`
- the AtomLens host keeps only a thin compatibility layer in
  `AtomLens/Chemistry/CDKPortAliases.swift`
- the AtomLens Xcode project links the package product through a local Swift
  package reference instead of compiling package sources directly

That arrangement is important for reuse, package publishing, and App Store
review clarity.

## Boundary Enforcement

- `PackageBoundaryTests.testSourcesContainNoAppLevelCouplingMarkers`
  - scans package sources for app-coupling markers
- `PackageBoundaryTests.testPackageManifestDoesNotDependOnAtomLensTargets`
  - verifies `Package.swift` does not reference AtomLens targets or local app
    paths
- `PackageBoundaryTests.testWorkspaceChemistryLayerContainsOnlyAliasAdapter`
  - monorepo guard: ensures host `AtomLens/Chemistry` does not accumulate
    CDK-derived implementation files
- `PackageBoundaryTests.testHostXcodeProjectUsesPackageProductInsteadOfPackageSourceFiles`
  - verifies the host Xcode project links the local package product and does not
    compile package source files directly

## Testing Strategy

- unit tests grouped by feature area (`Smiles`, `InChI`, `MDL`, `Rendering`,
  `Layout`, `QSAR`, `IO`)
- parity metadata tests that keep links back to upstream CDK tests where
  applicable
- IO round-trip tests for reader / writer consistency
- rendering tests for SVG and scene-builder behavior
- boundary tests that enforce package / host separation
