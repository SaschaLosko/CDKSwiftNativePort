# API Reference

This document summarizes the public integration surface of
`CDKSwiftNativePort`. It is organized around the way host applications
typically consume the package: model, layout, rendering, identifiers,
descriptors, and IO.

## 1. Core Domain Model

Defined primarily in [`Sources/CDKSwiftNativePort/Molecule.swift`](../Sources/CDKSwiftNativePort/Molecule.swift).

### Primary value types

- `Molecule`
- `Atom`
- `Bond`
- `MoleculeSgroup`
- `MoleculeSgroupBracket`

### Supporting enums and value types

- `ChemFormat`
- `BondOrder`
- `BondStereo`
- `BondTopology`
- `BondQueryType`
- `AtomChirality`
- `AtomQueryType`
- `CxCoordinate`
- `CxRadicalType`
- `ChemError`

### Common `Molecule` helpers

- graph queries:
  - `atom(id:)`
  - `bonds(forAtom:)`
  - `bond(between:and:)`
  - `neighbors(of:)`
- chemistry/layout helpers:
  - `implicitHydrogenCount(for:)`
  - `simpleCycles(maxSize:)`
  - `aromaticDisplayRings()`
  - `aromaticDisplayBondIDs()`
  - `boundingBox()`
- metadata:
  - `dataFields`
  - `orderedDataFieldNames`
  - `dataFieldValues(named:)`
  - `appendDataFieldValue(_:named:)`

## 2. Layout and Rendering

### 2.1 Layout

- `Depiction2DGenerator.generate(for:) -> Molecule`

Use this to generate or normalize 2D coordinates before depiction or export.

### 2.2 Rendering configuration

- `RenderStyle`
- `CDKHighlightStyle`
- `CDKAtomColoringMode`
- `CDKAromaticDisplayMode`
- `CDKRenderColor`
- `CDKRenderingStyleResolver`

Important `RenderStyle` controls include:

- carbon visibility
- implicit/explicit hydrogen visibility
- atom IDs and atom-map numbers
- highlight mode
- aromatic display mode
- bond width, font size, and padding

### 2.3 SVG depiction

- `CDKDepictionGenerator.toSVG(molecule:style:canvasSize:includeBackground:) -> String`

The SVG path includes:

- normal molecule depiction
- query atom/query bond visuals
- Markush legends and link-node annotations
- non-Markush Sgroup bracket annotations
- aromatic circle mode
- highlight styles

### 2.4 Renderer-neutral scene generation

- `CDKMetalDepictionSceneBuilder.build(molecule:style:canvasRect:zoom:pan:rotationDegrees:minimumLabelFontSize:includeAromaticCarbonLabelsWhenCarbonsHidden:includeTerminalCarbonLabelsWhenCarbonsHidden:) -> CDKMetalDepictionScene`
- `CDKMetalReactionDepictionSceneBuilder.build(reaction:style:canvasRect:zoom:pan:rotationDegrees:highlightedParticipant:withOuterGlowHighlight:) -> CDKMetalDepictionScene`
- `CDKMetalReactionDepictionSceneBuilder.participant(at:in:style:canvasRect:zoom:pan:rotationDegrees:) -> CDKReactionParticipantSelection?`

These scene builders do not require `MetalKit`. They return geometry and label
data for host-owned rendering pipelines.

### 2.5 Scene data types

- `CDKMetalDepictionScene`
  - `BackgroundBox`
  - `LineSegment`
  - `AtomLabel`
- `CDKReactionParticipantSelection`

### 2.6 Historical compatibility symbol

- `CDKStandardGenerator`

This symbol is currently retained as a compatibility name only. The active
public rendering APIs are `CDKDepictionGenerator` and the scene builders above.

## 3. SMILES, CXSMILES, and Reactions

### 3.1 SMILES parser/generator

- `CDKSmilesParser`
  - `parseSmiles(_:)`
  - `parseCoreSmiles(_:)`
  - `parseReactionSmiles(_:)`
- `CDKSmilesParserFactory.shared.newSmilesParser(flavor:)`
- `CDKSmilesGenerator`
  - `create(_ molecule:)`
  - `create(_ reaction:)`
- `CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor:)`
- `CDKSmiFlavor`

### 3.2 CXSMILES support

- `CDKCxSmilesState`
- `CDKCxSmilesParser`
  - `split(_:enabled:)`
  - `apply(to:state:)`
  - `applyAtomLabels(to:state:)`
  - `applyAtomValues(to:state:)`
  - `applyCoordinates(to:state:)`
  - `applyRadicals(to:state:)`
  - `applyStereoMetadata(to:state:)`
  - `applyLigandOrdering(to:state:)`
  - `applyHighlights(to:state:)`
  - `applySgroups(to:state:)`
  - `applyLinkNodes(to:state:)`

Current CX support includes:

- atom labels and atom values
- coordinates
- radicals
- enhanced stereo metadata
- fragment grouping
- ligand ordering
- highlight layers
- Sgroups
- link nodes and Markush `RG:` data

### 3.3 Reaction model

- `CDKReactionRole`
- `CDKReactionParticipant`
- `CDKReaction`
- `CDKReactionParticipantSelection`

## 4. InChI

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

## 5. Unified Import APIs

### 5.1 Import facade

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
  - `looksLikeReaction(text:fileExtension:)`

### 5.2 Supported molecule readers

- MDL:
  - `CDKMDLReader`
  - `CDKMDLV2000Reader`
  - `CDKMDLV3000Reader`
  - `CDKMDLV3000Scaffold`
  - `CDKIteratingSDFReader`
- SMILES:
  - `CDKSMILESReader`
- InChI:
  - `CDKInChIReader`
- MOL2 / PDB / XYZ:
  - `CDKMol2Reader`
  - `CDKPDBReader`
  - `CDKXYZReader`
- CML:
  - `CDKCMLReader`
- reaction readers:
  - `CDKCMLReactionReader`
  - `CDKRXNReader`
  - `CDKRDFReader`

### 5.3 File access helper

- `CDKFileAccess`
  - `withReadableURL(at:_:)`
  - `readData(from:)`
  - `decodeText(from:)`

This is the package-level, sandbox-aware file access helper used by importer
and extension-style host integrations.

## 6. Unified Export APIs

### 6.1 Export facade

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
  - `write(reaction:as:options:)`
  - `write(reactions:as:options:)`
  - `write(molecule:to:as:options:)`
  - `write(molecules:to:as:options:)`

`CDKFileExportFormat` currently includes:

- `.mol`
- `.molV3000`
- `.sdf`
- `.smiles`
- `.isomericSmiles`
- `.inchi`
- `.mol2`
- `.pdb`
- `.xyz`
- `.cml`
- `.rxn`
- `.rdf`
- `.svg`

### 6.2 Format-specific writers

- MDL / SDF:
  - `CDKMDLV2000Writer`
  - `CDKMDLV3000Writer`
  - `CDKSDFWriter`
  - `CDKSDFWriterOptions`
- SMILES / InChI:
  - `CDKSMILESWriter`
  - `CDKInChIWriter`
- MOL2 / PDB / XYZ / CML:
  - `CDKMol2Writer`
  - `CDKPDBWriter`
  - `CDKXYZWriter`
  - `CDKCMLWriter`
- reaction writers:
  - `CDKCMLReactionWriter`
  - `CDKRXNWriter`
  - `CDKRDFWriter`

## 7. Identifier and Property Services

### 7.1 Identifier service

- `CDKMoleculeIdentifiers`
- `CDKMoleculeIdentifierService`
  - `compute(for:smilesFlavor:isoSmilesFlavor:)`
  - `unavailableText(from:)`

### 7.2 Property service

- `CDKMolecularProperties`
- `CDKRuleOfFiveResult`
- `CDKRuleOfFiveDescriptor`
- `CDKMoleculePropertyService`

### 7.3 Low-level descriptors

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
- `CDKTPSADescriptor`
- `CDKVABCDescriptor`

## 8. Recommended Entry Points

For most host applications, the smallest useful public surface is:

- `CDKFileImporter`
- `Depiction2DGenerator`
- `CDKDepictionGenerator`
- `CDKMetalDepictionSceneBuilder`
- `CDKMetalReactionDepictionSceneBuilder`
- `CDKMoleculeIdentifierService`
- `CDKMoleculePropertyService`
- `CDKFileExporter`

That set covers the common parse, layout, inspect, depict, and export loop
without forcing host code to bind to every lower-level reader or writer.
