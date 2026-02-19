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

Useful methods on `Molecule` include:
- graph neighborhood/bond lookup helpers
- ring/cycle helpers
- implicit hydrogen helpers
- geometry helpers (`boundingBox`, normalization-friendly traversal utilities)

## 2) Layout and Depiction

- `Depiction2DGenerator.generate(for:)`
  - applies CDK-style 2D coordinate generation.
- `CDKDepictionGenerator.toSVG(...)`
  - produces SVG depictions.
- `CDKMetalDepictionSceneBuilder.build(...)`
  - builds a render-agnostic scene representation suitable for Metal-style renderers.
- `CDKStandardGenerator.draw(...)`
  - SwiftUI `GraphicsContext` drawing path.

Style controls:
- `RenderStyle`
- `CDKAtomColoringMode`
- `CDKAromaticDisplayMode`

## 3) Import APIs

Generic dispatch:
- `CDKFileImporter.readMolecules(from:)`
- `CDKFileImporter.readMolecules(text:fileExtension:)`
- `CDKFileImporter.preferredInputFormat(...)`

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

- Parsing:
  - `CDKSmilesParser`
  - `CDKSmilesParserFactory`
  - `CDKCxSmilesParser`
  - `CDKSmilesReactionParser`
- Generation:
  - `CDKSmilesGenerator`
  - `CDKSmilesGeneratorFactory`
  - `CDKSmiFlavor`

## 6) InChI APIs

- `CDKInChIToStructure`
- `CDKInChIGenerator`
- `CDKInChIGeneratorFactory`
- `CDKInChIStatus`

Facade:
- `CDKMoleculeIdentifierService.compute(for:)`
  - returns `CDKMoleculeIdentifiers` (SMILES, ISO SMILES, InChI, InChIKey).

## 7) Descriptors and Property Services

Facade:
- `CDKMoleculePropertyService.compute(for:)`
  - returns `CDKMolecularProperties`.

Rule-of-Five helpers:
- `CDKRuleOfFiveDescriptor`
- `CDKRuleOfFiveResult`

Descriptors:
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

Most parsing and writing APIs throw `ChemError`:
- `.emptyInput`
- `.unsupported(String)`
- `.parseFailed(String)`

Call sites should use `do / try / catch` and surface `LocalizedError.errorDescription` where appropriate.

## 9) Stability Notes

- Public symbols listed above are intended as the stable integration surface for application use.
- Internals under the `CDK/...` tree may evolve as parity with upstream CDK improves.
- For release history and compatibility notes, see `CHANGELOG.md`.
