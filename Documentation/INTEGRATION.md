# Integration Guide

This guide shows how to integrate `CDKSwiftNativePort` into a macOS host and
use the current public APIs end-to-end.

## 1. Add the Package

Use Swift Package Manager with:

```text
https://github.com/SaschaLosko/CDKSwiftNativePort.git
```

or in `Package.swift`:

```swift
dependencies: [
    .package(url: "https://github.com/SaschaLosko/CDKSwiftNativePort.git", from: "1.1.1")
]
```

Then import:

```swift
import CDKSwiftNativePort
```

## 2. Choose an Integration Style

Most hosts use one of these patterns:

- Import / inspect / export only
- Parse -> layout -> SVG export
- Parse -> layout -> build `CDKMetalDepictionScene` -> render in a host-owned
  Metal or Core Graphics layer
- Parse reactions and use `CDKMetalReactionDepictionSceneBuilder` for reaction UI

## 3. Read Molecules from Text

```swift
let molecules = try CDKFileImporter.readMolecules(
    text: "CC(=O)OC1=CC=CC=C1C(=O)O",
    fileExtension: "smi"
)

guard let molecule = molecules.first else {
    throw ChemError.emptyInput
}
```

`CDKFileImporter` also supports text autodetection for SDF, RXN, RDF, InChI,
SMILES, reaction SMILES, MOL2, CML, PDB, and XYZ payloads.

## 4. Read from Sandboxed File URLs

```swift
let fileURL = URL(fileURLWithPath: "/path/to/molecule.sdf")
let molecules = try CDKFileImporter.readMolecules(from: fileURL)
```

`readMolecules(from:)` and `readReaction(from:)` start and stop security-scoped
resource access for file URLs when needed, which makes them suitable for App
Sandbox hosts using user-selected files or security-scoped bookmarks.

## 5. Generate 2D Coordinates

```swift
let laidOut = Depiction2DGenerator.generate(for: molecule)
```

Layout is recommended before SVG export or scene generation unless your input
already carries usable 2D coordinates.

## 6. Compute Identifiers and Properties

```swift
let identifiers = CDKMoleculeIdentifierService.compute(for: laidOut)
let properties = CDKMoleculePropertyService.compute(for: laidOut)

print(identifiers.smiles)
print(identifiers.isoSmiles)
print(identifiers.inchi)
print(identifiers.inchiKey)

print(properties.formula)
print(properties.molecularWeight)
print(properties.tpsa)
print(properties.ruleOfFive.statusText)
```

For direct descriptor access:

```swift
let xlogP = CDKXLogPDescriptor.calculate(for: laidOut)
let tpsa = CDKTPSADescriptor.calculate(for: laidOut, checkAromaticity: true)
let vabc = CDKVABCDescriptor.calculate(for: laidOut)
```

## 7. Generate SVG

```swift
import CoreGraphics

var options = CDKFileExportOptions()
options.svgCanvasSize = CGSize(width: 1280, height: 840)
options.svgIncludeBackground = true

let svg = try CDKFileExporter.write(molecule: laidOut, as: .svg, options: options)
```

You can also call the lower-level generator directly:

```swift
import CoreGraphics

let svg = CDKDepictionGenerator.toSVG(
    molecule: laidOut,
    style: RenderStyle(),
    canvasSize: CGSize(width: 1280, height: 840),
    includeBackground: true
)
```

## 8. Build a Scene for a Host Renderer

```swift
import CoreGraphics

let scene = CDKMetalDepictionSceneBuilder.build(
    molecule: laidOut,
    style: RenderStyle(),
    canvasRect: CGRect(x: 0, y: 0, width: 960, height: 720),
    zoom: 1.0,
    pan: .zero,
    rotationDegrees: 0
)

print(scene.backgroundBoxes.count)
print(scene.bondSegments.count)
print(scene.labels.count)
```

Important: the package computes scene geometry only. Your app owns the actual
Metal, Core Graphics, or raster drawing implementation.

## 9. Parse a Markush CXSMILES

```swift
import CoreGraphics

let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
let markush = try parser.parseSmiles(
    "C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|"
)

let laidOutMarkush = Depiction2DGenerator.generate(for: markush)
let markushScene = CDKMetalDepictionSceneBuilder.build(
    molecule: laidOutMarkush,
    style: RenderStyle(showCarbons: false, showImplicitHydrogens: true),
    canvasRect: CGRect(x: 0, y: 0, width: 960, height: 720),
    zoom: 1.0,
    pan: .zero
)
```

Supported Markush-related behavior includes:

- CXSMILES `RG:` parsing
- CXSMILES `LN:` repeat/link-node parsing for the supported path
- R-group legend boxes in SVG / SwiftUI / Metal scene outputs
- fixed-legend Metal scene behavior for interactive host rotation

## 10. Parse and Render Reactions

```swift
import CoreGraphics

let reaction = try CDKFileImporter.readReaction(
    text: "CCO>>CC=O",
    fileExtension: "rsmi"
)

let reactionScene = CDKMetalReactionDepictionSceneBuilder.build(
    reaction: reaction,
    style: RenderStyle(),
    canvasRect: CGRect(x: 0, y: 0, width: 1200, height: 680),
    zoom: 1.0,
    pan: .zero
)
```

Optional participant hit-testing:

```swift
let hit = CDKMetalReactionDepictionSceneBuilder.participant(
    at: CGPoint(x: 320, y: 240),
    in: reaction,
    style: RenderStyle(),
    canvasRect: CGRect(x: 0, y: 0, width: 1200, height: 680),
    zoom: 1.0,
    pan: .zero
)
```

## 11. Export to Common Formats

```swift
let smiles = try CDKFileExporter.write(molecule: laidOut, as: .smiles)
let isoSmiles = try CDKFileExporter.write(molecule: laidOut, as: .isomericSmiles)
let inchi = try CDKFileExporter.write(molecule: laidOut, as: .inchi)
let molfile = try CDKFileExporter.write(molecule: laidOut, as: .mol)
```

For full reaction export:

```swift
let rxn = try CDKRXNWriter.write(
    reactants: reaction.reactants,
    products: reaction.products,
    agents: reaction.agents,
    reactionName: "Oxidation"
)
```

## 12. Error Handling Pattern

```swift
do {
    _ = try CDKFileImporter.readMolecules(text: "", fileExtension: "smi")
} catch let error as ChemError {
    switch error {
    case .emptyInput:
        print("Input is empty")
    case .unsupported(let details):
        print("Unsupported format: \(details)")
    case .parseFailed(let details):
        print("Parse failed: \(details)")
    }
} catch {
    print("Unexpected error: \(error)")
}
```

## 13. Host Boundary Guidance

Keep host-specific concerns outside this package:

- Spotlight metadata and index extensions
- Quick Look preview / thumbnail extension wiring
- `MTKView` ownership and renderer lifecycle
- app entitlements, bundle IDs, document/window/session logic

Use `CDKSwiftNativePort` as the chemistry core and build host features around
its public APIs.
