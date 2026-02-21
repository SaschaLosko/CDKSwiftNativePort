# API Reference (High-Level)

This document summarizes the primary public API surface of `CDKSwiftNativePort`.

## 1) Core Model

Defined in `Sources/CDKSwiftNativePort/Molecule.swift`.

- `Molecule`
- `Atom`
- `Bond`
- `BondOrder`
- `BondStereo`
- `AtomChirality`
- `ChemFormat`
- `ChemError`

`Molecule` exposes graph helpers (neighbors, bond lookup), ring/cycle utilities, implicit-hydrogen helpers, and geometry helpers such as `boundingBox()`.

## 2) Layout and Depiction

- `Depiction2DGenerator.generate(for:)`
  - CDK-style 2D coordinate generation.
- `CDKDepictionGenerator.toSVG(...)`
  - SVG depiction output.
- `CDKMetalDepictionSceneBuilder.build(...)`
  - render-agnostic scene generation for Metal renderers in host apps.
- `CDKStandardGenerator.draw(...)`
  - SwiftUI `GraphicsContext` rendering path.

Style controls:
- `RenderStyle`
- `CDKAtomColoringMode`
- `CDKAromaticDisplayMode`

## 3) Import APIs

Generic dispatch:
- `CDKFileImporter.readMolecules(from:)`
- `CDKFileImporter.readMolecules(text:fileExtension:)`
- `CDKFileImporter.preferredInputFormat(...)`
- `CDKFileImporter.formats`

Format-specific readers:
- MDL: `CDKMDLReader`, `CDKMDLV2000Reader`, `CDKMDLV3000Reader`, `CDKIteratingSDFReader`
- SMILES: `CDKSMILESReader`
- InChI: `CDKInChIReader`
- MOL2: `CDKMol2Reader`
- PDB: `CDKPDBReader`
- XYZ: `CDKXYZReader`
- CML: `CDKCMLReader`
- RXN/RDF: `CDKRXNReader`, `CDKRDFReader`

## 4) Export APIs

Generic dispatch:
- `CDKFileExporter.write(...)`
- `CDKFileExporter.write(..., to:as:options:)`
- `CDKFileExporter.formats`
- `CDKFileExportFormat`
- `CDKFileExportOptions`

Format-specific writers:
- MDL/SDF: `CDKMDLV2000Writer`, `CDKSDFWriter`
- SMILES: `CDKSMILESWriter`
- InChI: `CDKInChIWriter`
- MOL2: `CDKMol2Writer`
- PDB: `CDKPDBWriter`
- XYZ: `CDKXYZWriter`
- CML: `CDKCMLWriter`
- RXN/RDF: `CDKRXNWriter`, `CDKRDFWriter`

## 5) SMILES APIs

Parsing:
- `CDKSmilesParser`
- `CDKSmilesParserFactory`
- `CDKCxSmilesParser`
- `CDKSmilesReactionParser`

Generation:
- `CDKSmilesGenerator`
- `CDKSmilesGeneratorFactory`
- `CDKSmiFlavor`

## 6) InChI APIs

- `CDKInChIToStructure`
- `CDKInChIGenerator`
- `CDKInChIGeneratorFactory`
- `CDKInChIStatus`

Identifier facade:
- `CDKMoleculeIdentifierService.compute(for:)`
  - returns `CDKMoleculeIdentifiers` (`smiles`, `isoSmiles`, `inchi`, `inchiKey`).

## 7) Descriptors and Property Services

Property facade:
- `CDKMoleculePropertyService.compute(for:)`
  - returns `CDKMolecularProperties`.

Rule-of-Five:
- `CDKRuleOfFiveDescriptor`
- `CDKRuleOfFiveResult`

Descriptor set includes:
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

## 8) Error Handling

Parsing and writing APIs throw `ChemError`:
- `.emptyInput`
- `.unsupported(String)`
- `.parseFailed(String)`

Host apps should use `do/try/catch` and expose `LocalizedError.errorDescription`.

## 9) Host App Boundary

The package intentionally excludes host-application integration concerns:
- no Spotlight indexer code
- no Quick Look providers
- no security-scoped bookmark/sandbox file-access policy code
- no app bundle identifiers or app-specific metadata keys

These are expected to live in the consuming application.

## 10) Stability Notes

- The API families listed above are the intended integration surface.
- Internal implementation details under `Sources/CDKSwiftNativePort/CDK/...` may evolve as parity work continues.
- See `CHANGELOG.md` for compatibility notes and release history.
