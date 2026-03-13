# API Reference

This document summarizes the public integration surface of
`CDKSwiftNativePort`. It is organized by the way host applications typically
consume the package: model, import, layout, rendering, identifiers, properties,
and export.

## 1. Core Domain Model

Defined primarily in `Sources/CDKSwiftNativePort/Molecule.swift`.

### Primary value types

- `Molecule`
- `Atom`
- `Bond`
- `MoleculeSgroup`

### Supporting enums

- `ChemFormat`
- `BondOrder`
- `BondStereo`
- `AtomChirality`
- `AtomQueryType`
- `BondQueryType`
- `ChemError`

### Important `Molecule` helpers

- Graph queries:
  - `atom(id:)`
  - `bonds(forAtom:)`
  - `bond(between:and:)`
  - `neighbors(of:)`
- Layout/render support:
  - `boundingBox()`
  - `implicitHydrogenCount(for:)`
  - `assignWedgeHashFromChiralCenters()`
- Ring / aromatic helpers:
  - `simpleCycles(maxSize:)`
  - `aromaticDisplayRings()`
  - `aromaticDisplayBondIDs()`
- Data fields:
  - `dataFields`
  - `orderedDataFieldNames`
  - `appendDataFieldValue(_:named:)`

## 2. Layout and Rendering

### 2.1 Layout

- `Depiction2DGenerator.generate(for:) -> Molecule`

This is the main entry point for generating or refining 2D coordinates before
rendering or export.

### 2.2 Rendering configuration

- `RenderStyle`
- `CDKAtomColoringMode`
- `CDKAromaticDisplayMode`
- `CDKRenderColor`
- `CDKRenderingStyleResolver`

### 2.3 SVG output

- `CDKDepictionGenerator.toSVG(molecule:style:canvasSize:includeBackground:) -> String`

Use this for GitHub-ready previews, web embedding, snapshots, or export flows.

### 2.4 SwiftUI `GraphicsContext` rendering

- `CDKStandardGenerator.draw(molecule:in:style:context:)`

This is the package's native SwiftUI renderer.

### 2.5 Renderer-agnostic Metal scene generation

- `CDKMetalDepictionSceneBuilder.build(molecule:style:canvasRect:zoom:pan:rotationDegrees:minimumLabelFontSize:includeAromaticCarbonLabelsWhenCarbonsHidden:includeTerminalCarbonLabelsWhenCarbonsHidden:) -> CDKMetalDepictionScene`
- `CDKMetalReactionDepictionSceneBuilder.build(reaction:style:canvasRect:zoom:pan:rotationDegrees:highlightedParticipant:withOuterGlowHighlight:) -> CDKMetalDepictionScene`
- `CDKMetalReactionDepictionSceneBuilder.participant(at:in:style:canvasRect:zoom:pan:rotationDegrees:) -> CDKReactionParticipantSelection?`

### 2.6 Scene data types

- `CDKMetalDepictionScene`
  - `CDKMetalDepictionScene.BackgroundBox`
  - `CDKMetalDepictionScene.LineSegment`
  - `CDKMetalDepictionScene.AtomLabel`
- `CDKReactionParticipantSelection`

These types are intentionally rendering-backend neutral. The package computes
geometry, labels, and annotation placement; the host app owns drawing.

## 3. SMILES, CXSMILES, and Reactions

### 3.1 SMILES parser/generator

- `CDKSmilesParser`
  - `parseSmiles(_:)`
  - `parseCoreSmiles(_:)`
  - `parseReactionSmiles(_:)`
- `CDKSmilesParserFactory.shared.newSmilesParser(flavor:)`
- `CDKSmilesGenerator`
  - `create(_:)`
- `CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor:)`
- `CDKSmiFlavor`

`CDKSmiFlavor` includes:

- `.useAromaticSymbols`
- `.isomeric`
- `.strict`
- `.cxsmiles`
- `.cdkDefault`

### 3.2 CXSMILES state and Markush support

- `CDKCxSmilesState`
- `CDKCxSmilesParser`
  - `split(_:enabled:)`
  - `applyAtomLabels(to:state:)`
  - `SplitResult`

Current CDK-derived CX coverage includes atom labels, CX suffix handling, link
nodes, and CDK 2.12-style Markush / R-group parsing used in the supported
depiction path.

### 3.3 Reaction model

- `CDKReactionRole`
- `CDKReactionParticipant`
- `CDKReaction`
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

### 5.2 Format-specific readers

- MDL:
  - `CDKMDLReader`
  - `CDKMDLV2000Reader`
  - `CDKMDLV3000Reader`
  - `CDKMDLV3000Scaffold`
  - `CDKIteratingSDFReader`
- SMILES / reaction SMILES:
  - `CDKSMILESReader`
- InChI:
  - `CDKInChIReader`
- MOL2 / PDB / XYZ / CML:
  - `CDKMol2Reader`
  - `CDKPDBReader`
  - `CDKXYZReader`
  - `CDKCMLReader`
- RXN / RDF:
  - `CDKRXNReader`
  - `CDKRDFReader`

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
  - `write(molecule:to:as:options:)`
  - `write(molecules:to:as:options:)`

### 6.2 Format-specific writers

- MDL / SDF:
  - `CDKMDLV2000Writer`
  - `CDKSDFWriter`
- SMILES / InChI:
  - `CDKSMILESWriter`
  - `CDKInChIWriter`
- MOL2 / PDB / XYZ / CML:
  - `CDKMol2Writer`
  - `CDKPDBWriter`
  - `CDKXYZWriter`
  - `CDKCMLWriter`
- Reaction writers:
  - `CDKRXNWriter`
  - `CDKRDFWriter`

For full reaction export, call `CDKRXNWriter.write(reactants:products:agents:)`
or `CDKRDFWriter.write(...)` directly.

## 7. Identifiers and Property Services

### 7.1 Identifier service

- `CDKMoleculeIdentifiers`
- `CDKMoleculeIdentifierService`
  - `compute(for:smilesFlavor:isoSmilesFlavor:)`
  - `unavailableText(from:)`

### 7.2 Property summary service

- `CDKMolecularProperties`
- `CDKRuleOfFiveResult`
- `CDKRuleOfFiveDescriptor`
  - `evaluate(for:xlogP:)`
  - `evaluate(molecularWeight:hBondDonorCount:hBondAcceptorCount:xlogP:)`
- `CDKMoleculePropertyService`
  - `compute(for:xlogP:)`

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

## 8. Common Usage Patterns

### Parse and layout a molecule

```swift
let molecules = try CDKFileImporter.readMolecules(text: "CCO", fileExtension: "smi")
guard let molecule = molecules.first else { throw ChemError.emptyInput }
let laidOut = Depiction2DGenerator.generate(for: molecule)
```

### Compute identifiers and properties

```swift
let ids = CDKMoleculeIdentifierService.compute(for: laidOut)
let properties = CDKMoleculePropertyService.compute(for: laidOut)
```

### Export SVG

```swift
import CoreGraphics

let svg = CDKDepictionGenerator.toSVG(
    molecule: laidOut,
    style: RenderStyle(),
    canvasSize: CGSize(width: 1200, height: 840),
    includeBackground: true
)
```

### Parse a Markush CXSMILES

```swift
let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
let molecule = try parser.parseSmiles(
    "C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|"
)
```

### Build a Metal-friendly scene

```swift
import CoreGraphics

let scene = CDKMetalDepictionSceneBuilder.build(
    molecule: laidOut,
    style: RenderStyle(),
    canvasRect: CGRect(x: 0, y: 0, width: 960, height: 720),
    zoom: 1.0,
    pan: .zero
)
```

## 9. Error Model

Most parsing and IO APIs throw `ChemError`:

- `.emptyInput`
- `.unsupported(String)`
- `.parseFailed(String)`

## 10. Boundary Contract

`CDKSwiftNativePort` excludes host-app implementation such as Spotlight, Quick
Look, app window/session logic, and app-branded identifiers. Host apps should
consume the package through the public APIs listed above and keep platform- or
product-specific wiring outside the package.
