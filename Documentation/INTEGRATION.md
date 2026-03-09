# Integration Guide

This guide shows how to integrate `CDKSwiftNativePort` into a macOS app and use the main APIs end-to-end.

## 1. Add the Package

Use Swift Package Manager with:

```text
https://github.com/SaschaLosko/CDKSwiftNativePort.git
```

or in `Package.swift`:

```swift
.package(url: "https://github.com/SaschaLosko/CDKSwiftNativePort.git", from: "1.1.1")
```

Then import:

```swift
import CDKSwiftNativePort
```

## 2. Parse Molecules from Text

```swift
let molecules = try CDKFileImporter.readMolecules(
    text: "CC(=O)OC1=CC=CC=C1C(=O)O",
    fileExtension: "smi"
)

guard let molecule = molecules.first else {
    throw ChemError.emptyInput
}
```

## 3. Parse from File URLs

```swift
let fileURL = URL(fileURLWithPath: "/path/to/molecule.sdf")
let molecules = try CDKFileImporter.readMolecules(from: fileURL)
```

`CDKFileImporter.readMolecules(from:)` and `readReaction(from:)` handle security-scoped resource access for file URLs when the host is sandboxed.

## 4. Generate Layout, IDs, and Properties

```swift
let laidOut = Depiction2DGenerator.generate(for: molecule)
let identifiers = CDKMoleculeIdentifierService.compute(for: laidOut)
let properties = CDKMoleculePropertyService.compute(for: laidOut)

print(identifiers.smiles)
print(identifiers.inchiKey)
print(properties.formula)
print(properties.ruleOfFive.statusText)
```

## 5. Generate SVG Output

```swift
var options = CDKFileExportOptions()
options.svgCanvasSize = CGSize(width: 1280, height: 840)
options.svgIncludeBackground = true

let svg = try CDKFileExporter.write(molecule: laidOut, as: .svg, options: options)
```

## 6. Build Metal-Friendly Scenes

```swift
import CoreGraphics

let scene = CDKMetalDepictionSceneBuilder.build(
    molecule: laidOut,
    style: RenderStyle(),
    canvasRect: CGRect(x: 0, y: 0, width: 900, height: 620),
    zoom: 1.0,
    pan: .zero
)

print(scene.bondSegments.count, scene.labels.count)
```

## 7. Parse and Render Reactions

```swift
import CoreGraphics

let reaction = try CDKFileImporter.readReaction(text: "CCO>>CC=O", fileExtension: "rsmi")

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

## 8. Export to Common Formats

```swift
let smiles = try CDKFileExporter.write(molecule: laidOut, as: .smiles)
let inchi = try CDKFileExporter.write(molecule: laidOut, as: .inchi)
let molfile = try CDKFileExporter.write(molecule: laidOut, as: .mol)

print(smiles)
print(inchi)
print(molfile.prefix(80))
```

## 9. Error Handling Pattern

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

## 10. Host Boundary Guidance

Keep host-specific concerns outside this package:

- Spotlight metadata and index extensions
- Quick Look preview/thumbnail extension wiring
- app entitlements, bundle IDs, and app-window/session logic

Use `CDKSwiftNativePort` as the chemistry core and build host features around it.
