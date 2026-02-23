# API Reference

This document summarizes the current public API surface of `CDKSwiftNativePort`.

## 1) Core Types

Defined in `Sources/CDKSwiftNativePort/Molecule.swift`.

- `Molecule`
- `Atom`
- `Bond`
- `BondOrder`
- `BondStereo`
- `AtomChirality`
- `AtomQueryType`
- `BondQueryType`
- `ChemFormat`
- `ChemError`
- `Depiction2DGenerator`

Key `Molecule` helpers include:
- connectivity/lookup: `bonds(forAtom:)`, `bond(between:and:)`, `neighbors(of:)`
- ring/cycle helpers: `simpleCycles(maxSize:)`, `aromaticDisplayRings()`, `aromaticDisplayBondIDs()`
- depiction-related chemistry helpers: `implicitHydrogenCount(for:)`, `assignWedgeHashFromChiralCenters()`
- geometry: `boundingBox()`

## 2) Layout + Rendering

### Layout
- `Depiction2DGenerator.generate(for:) -> Molecule`
  - CDK-style 2D placement entry point.

### Rendering outputs
- `CDKDepictionGenerator.toSVG(...) -> String`
  - vector depiction export.
- `CDKMetalDepictionSceneBuilder.build(...) -> CDKMetalDepictionScene`
  - renderer-agnostic primitive scene suitable for Metal hosts.
- `CDKStandardGenerator.draw(molecule:in:style:context:)`
  - SwiftUI `GraphicsContext` drawing path.

### Rendering style and color APIs
- `RenderStyle`
- `CDKAtomColoringMode`
- `CDKAromaticDisplayMode`
- `CDKRenderColor`
- `CDKRenderingStyleResolver`
- `CDKMetalDepictionScene` and nested:
  - `CDKMetalDepictionScene.LineSegment`
  - `CDKMetalDepictionScene.AtomLabel`

## 3) SMILES / CXSMILES / Reaction SMILES

### SMILES parser/generator
- `CDKSmilesParser`
  - `parseSmiles(_:)`
  - `parseCoreSmiles(_:)`
  - `parseReactionSmiles(_:)`
- `CDKSmilesParserFactory.shared.newSmilesParser(flavor:)`
- `CDKSmilesGenerator`
  - `create(_:)`
- `CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor:)`
- `CDKSmiFlavor`

### CXSMILES support
- `CDKCxSmilesState`
- `CDKCxSmilesParser`
  - `split(_:enabled:)`
  - `applyAtomLabels(to:state:)`
  - `SplitResult`

### Reaction container
- `CDKReaction`
  - `reactants`, `agents`, `products`, `cxState`

## 4) InChI APIs

- `CDKInChIGeneratorFactory.shared`
  - `getInChIToStructure(_:)`
  - `getInChIGenerator(_:)`
- `CDKInChIToStructure`
  - `getStatus()`
  - `getMessage()`
  - `getAtomContainer()`
- `CDKInChIGenerator`
  - `getStatus()`
  - `getMessage()`
  - `getInchi()`
  - `getInchiKey()`
- `CDKInChIStatus`

## 5) Import APIs

### Unified importer facade
- `CDKFileImporterFormat`
- `CDKFileImporter`
  - `formats`
  - `supportedFileExtensions`
  - `supportedUTIIdentifiers`
  - `supports(fileExtension:)`
  - `preferredInputFormat(forFileExtension:text:)`
  - `readMolecules(from:)`
  - `readMolecules(text:fileExtension:)`

### Format-specific readers
- MDL/SDF: `CDKMDLReader`, `CDKMDLV2000Reader`, `CDKIteratingSDFReader`
- V3000: `CDKMDLV3000Reader`, `CDKMDLV3000Scaffold`
- SMILES: `CDKSMILESReader`
- InChI: `CDKInChIReader`
- MOL2: `CDKMol2Reader`
- PDB: `CDKPDBReader`
- XYZ: `CDKXYZReader`
- CML: `CDKCMLReader`
- RXN/RDF: `CDKRXNReader`, `CDKRDFReader`

## 6) Export APIs

### Unified exporter facade
- `CDKFileExportFormat`
- `CDKFileExporterFormat`
- `CDKFileExportOptions`
- `CDKFileExporter`
  - `formats`
  - `supportedFileExtensions`
  - `supportedUTIIdentifiers`
  - `format(forFileExtension:)`
  - `format(forUTIIdentifier:)`
  - `write(molecule:as:options:)`
  - `write(molecules:as:options:)`
  - `write(molecule:to:as:options:)`
  - `write(molecules:to:as:options:)`

### Format-specific writers
- MDL/SDF: `CDKMDLV2000Writer`, `CDKSDFWriter`
- SMILES: `CDKSMILESWriter`
- InChI: `CDKInChIWriter`
- MOL2: `CDKMol2Writer`
- PDB: `CDKPDBWriter`
- XYZ: `CDKXYZWriter`
- CML: `CDKCMLWriter`
- RXN/RDF: `CDKRXNWriter`, `CDKRDFWriter`

## 7) Identifier + Property Facades

### Identifier service
- `CDKMoleculeIdentifiers`
- `CDKMoleculeIdentifierService`
  - `compute(for:smilesFlavor:isoSmilesFlavor:)`
  - `unavailableText(from:)`

### Property service
- `CDKMolecularProperties`
- `CDKRuleOfFiveResult`
- `CDKRuleOfFiveDescriptor`
  - `evaluate(for:xlogP:)`
  - `evaluate(molecularWeight:hBondDonorCount:hBondAcceptorCount:xlogP:)`
- `CDKMoleculePropertyService`
  - `compute(for:xlogP:)`

### Descriptor APIs
- `CDKMolecularFormulaDescriptor`
- `CDKMolecularWeightDescriptor`
- `CDKExactMassDescriptor`
- `CDKHeavyAtomCountDescriptor`
- `CDKHBondDonorCountDescriptor`
- `CDKHBondAcceptorCountDescriptor`
- `CDKRotatableBondsCountDescriptor`
- `CDKRingCountDescriptor`
- `CDKXLogPDescriptor`
- `CDKMannholdLogPDescriptor`

## 8) Errors

Most parsing and IO APIs throw `ChemError`:
- `.emptyInput`
- `.unsupported(String)`
- `.parseFailed(String)`

## 9) Host-App Boundary

This package intentionally excludes app integration layers:
- Spotlight extension/indexing logic
- Quick Look thumbnail/preview providers
- app-level identifiers/entitlements/window/session logic

`CDKSwiftNativePort` can be used from sandboxed and non-sandboxed hosts, but host-level policies and extension wiring remain the app’s responsibility.

## 10) Stability Notes

- The API families above are the intended integration surface for host apps.
- Internal implementation details under `Sources/CDKSwiftNativePort/CDK/...` may evolve as parity work continues.
- See `CHANGELOG.md` for compatibility notes.
