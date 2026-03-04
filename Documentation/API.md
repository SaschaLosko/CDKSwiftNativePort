# API Reference

This document summarizes the public integration surface of `CDKSwiftNativePort`.

## 1. Core Domain Model

Defined primarily in `Sources/CDKSwiftNativePort/Molecule.swift`.

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

Key `Molecule` helpers:

- graph and lookup: `bonds(forAtom:)`, `bond(between:and:)`, `neighbors(of:)`, `atom(id:)`
- ring/aromatic support: `simpleCycles(maxSize:)`, `aromaticDisplayRings()`, `aromaticDisplayBondIDs()`
- depiction chemistry helpers: `implicitHydrogenCount(for:)`, `assignWedgeHashFromChiralCenters()`
- geometry: `boundingBox()`

## 2. Layout and Rendering APIs

### Layout

- `Depiction2DGenerator.generate(for:) -> Molecule`

### Rendering outputs

- `CDKDepictionGenerator.toSVG(molecule:style:canvasSize:includeBackground:) -> String`
- `CDKStandardGenerator.draw(molecule:in:style:context:)`
- `CDKMetalDepictionSceneBuilder.build(molecule:style:canvasRect:zoom:pan:rotationDegrees:...) -> CDKMetalDepictionScene`
- `CDKMetalReactionDepictionSceneBuilder.build(reaction:style:canvasRect:zoom:pan:rotationDegrees:...) -> CDKMetalDepictionScene`
- `CDKMetalReactionDepictionSceneBuilder.participant(at:in:style:canvasRect:zoom:pan:rotationDegrees:) -> CDKReactionParticipantSelection?`

### Rendering configuration and scene types

- `RenderStyle`
- `CDKAtomColoringMode`
- `CDKAromaticDisplayMode`
- `CDKRenderColor`
- `CDKRenderingStyleResolver`
- `CDKMetalDepictionScene`
  - `CDKMetalDepictionScene.LineSegment`
  - `CDKMetalDepictionScene.AtomLabel`

## 3. SMILES and Reaction APIs

### Parsers and generators

- `CDKSmilesParser`
  - `parseSmiles(_:)`
  - `parseCoreSmiles(_:)`
  - `parseReactionSmiles(_:)`
- `CDKSmilesParserFactory.shared.newSmilesParser(flavor:)`
- `CDKSmilesGenerator`
  - `create(_:)`
- `CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor:)`
- `CDKSmiFlavor`

### CXSMILES

- `CDKCxSmilesState`
- `CDKCxSmilesParser`
  - `split(_:enabled:)`
  - `applyAtomLabels(to:state:)`
  - `SplitResult`

### Reaction model

- `CDKReaction`
- `CDKReactionRole`
- `CDKReactionParticipant`
- `CDKReactionParticipantSelection`

## 4. InChI APIs

- `CDKInChIGeneratorFactory.shared`
  - `getInChIToStructure(_:) -> CDKInChIToStructure`
  - `getInChIGenerator(_:) -> CDKInChIGenerator`
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

## 5. Import APIs

### Unified importer

- `CDKFileImporterFormat`
- `CDKFileImporter`
  - `formats`
  - `supportedFileExtensions`
  - `supportedUTIIdentifiers`
  - `supports(fileExtension:)`
  - `preferredInputFormat(forFileExtension:text:)`
  - `readMolecules(from:)`
  - `readMolecules(text:fileExtension:)`
  - `readReaction(from:)`
  - `readReaction(text:fileExtension:)`

### Format readers

- MDL/SDF: `CDKMDLReader`, `CDKMDLV2000Reader`, `CDKMDLV3000Reader`, `CDKMDLV3000Scaffold`, `CDKIteratingSDFReader`
- SMILES: `CDKSMILESReader`
- InChI: `CDKInChIReader`
- MOL2: `CDKMol2Reader`
- PDB: `CDKPDBReader`
- XYZ: `CDKXYZReader`
- CML: `CDKCMLReader`
- RXN/RDF: `CDKRXNReader`, `CDKRDFReader`

## 6. Export APIs

### Unified exporter

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

### Format writers

- MDL/SDF: `CDKMDLV2000Writer`, `CDKSDFWriter`
- SMILES: `CDKSMILESWriter`
- InChI: `CDKInChIWriter`
- MOL2: `CDKMol2Writer`
- PDB: `CDKPDBWriter`
- XYZ: `CDKXYZWriter`
- CML: `CDKCMLWriter`
- RXN/RDF: `CDKRXNWriter`, `CDKRDFWriter`

## 7. Identifier and Property APIs

### Identifiers

- `CDKMoleculeIdentifiers`
- `CDKMoleculeIdentifierService`
  - `compute(for:smilesFlavor:isoSmilesFlavor:)`
  - `unavailableText(from:)`

### Properties and Rule-of-Five

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

## 8. Error Model

Most parsing and IO APIs throw `ChemError`:

- `.emptyInput`
- `.unsupported(String)`
- `.parseFailed(String)`

## 9. Minimal Integration Example

```swift
import CDKSwiftNativePort

let molecules = try CDKFileImporter.readMolecules(text: "CCO", fileExtension: "smi")
guard let molecule = molecules.first else { throw ChemError.emptyInput }

let laidOut = Depiction2DGenerator.generate(for: molecule)
let ids = CDKMoleculeIdentifierService.compute(for: laidOut)
let properties = CDKMoleculePropertyService.compute(for: laidOut)
let svg = CDKDepictionGenerator.toSVG(molecule: laidOut)

print(ids.smiles, properties.formula, svg.count)
```

## 10. Package Boundary

`CDKSwiftNativePort` excludes host-app integration code (Spotlight, Quick Look, app window/session logic, app identifiers). Host apps should consume this package only through the public APIs listed above.
