# Architecture Overview

## Top-Level Layout

- `Sources/CDKSwiftNativePort/Molecule.swift`
  - domain model, graph helpers, shared errors, layout entry point.
- `Sources/CDKSwiftNativePort/CDK/Common`
  - identifier and property service facades.
- `Sources/CDKSwiftNativePort/CDK/Smiles`
  - SMILES/CXSMILES/reaction parser + generator stack.
- `Sources/CDKSwiftNativePort/CDK/InChI`
  - InChI parse/generate services.
- `Sources/CDKSwiftNativePort/CDK/Layout`
  - structure diagram generation.
- `Sources/CDKSwiftNativePort/CDK/Rendering`
  - depiction preprocessing, style resolution, SVG and scene generation.
- `Sources/CDKSwiftNativePort/CDK/IO`
  - format readers/writers and unified importer/exporter facades.
- `Sources/CDKSwiftNativePort/CDK/QSAR`
  - descriptors and molecular property calculations.

## Design Goals

- Encapsulate CDK-derived chemistry logic in one package.
- Provide stable, host-consumable Swift APIs.
- Keep host application wiring out of chemistry-port internals.

## Boundary Contract

Not allowed in package sources:

- app module imports (for example `AtomLens`)
- Spotlight / Quick Look integration code
- host bundle IDs, entitlements, or host window/session logic

Allowed in host apps:

- importing `CDKSwiftNativePort`
- rendering package scene models in host-specific UI code
- host-specific metadata/indexing/extension layers

## Boundary Enforcement

- `PackageBoundaryTests.testSourcesContainNoAppLevelCouplingMarkers`
  - scans package sources for app-coupling markers.
- `PackageBoundaryTests.testPackageManifestDoesNotDependOnAtomLensTargets`
  - verifies `Package.swift` does not reference AtomLens targets or local app paths.
- `PackageBoundaryTests.testWorkspaceChemistryLayerContainsOnlyAliasAdapter`
  - monorepo guard: ensures host `AtomLens/Chemistry` does not accumulate CDK-derived implementation files.

## Testing Strategy

- Unit tests grouped by feature area (`Smiles`, `InChI`, `MDL`, `Rendering`, `Layout`, `QSAR`, `IO`).
- Port metadata tests keep parity traceability to upstream CDK tests.
- IO round-trip tests validate reader/writer consistency.
- Boundary tests enforce package/host separation.
