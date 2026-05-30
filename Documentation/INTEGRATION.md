# Integration Guide

This guide shows how to integrate `CDKSwiftNativePort` into a native app,
helper tool, or Linux command-line host without pulling chemistry logic into
the host target itself.

## 1. Add the Package

Use Swift Package Manager with:

```text
https://github.com/SaschaLosko/CDKSwiftNativePort.git
```

or in `Package.swift`:

```swift
dependencies: [
    .package(url: "https://github.com/SaschaLosko/CDKSwiftNativePort.git", from: "1.4.2")
]
```

Then import it:

```swift
import CDKSwiftNativePort
```

For a full Swift-on-Linux walkthrough, including runnable executable-package
examples, see [`LINUX.md`](LINUX.md).

## 2. Pick an Integration Style

Most hosts fit one of these patterns:

- import, inspect, and export only
- parse -> layout -> SVG
- parse -> layout -> build `CDKMetalDepictionScene` -> render in host-owned
  Metal/Core Graphics
- parse reactions and use `CDKMetalReactionDepictionSceneBuilder`
- build Linux command-line tools or batch converters
- build Quick Look / Spotlight / importer targets on top of the package APIs

## 3. Read Molecules From Text

```swift
import Foundation

let molecules = try CDKFileImporter.readMolecules(
    text: "CC(=O)OC1=CC=CC=C1C(=O)O",
    fileExtension: "smi"
)

guard let molecule = molecules.first else {
    throw ChemError.emptyInput
}
```

`CDKFileImporter` supports text autodetection for common payloads including:

- SMILES and CXSMILES
- reaction SMILES
- MDL Molfile and SDF
- RXN and RDF
- InChI
- MOL2, PDB, XYZ
- CML molecules and CML reactions

## 4. Read From Sandboxed File URLs

```swift
import Foundation

let fileURL = URL(fileURLWithPath: "/path/to/molecule.sdf")
let molecules = try CDKFileImporter.readMolecules(from: fileURL)
```

The package uses `CDKFileAccess` internally, which coordinates file access and
works with security-scoped URLs in App Sandbox hosts.

On Linux and other non-sandboxed hosts, the same API uses standard Foundation
file access without introducing host-app dependencies.

## 5. Generate 2D Coordinates

```swift
import Foundation

let laidOut = Depiction2DGenerator.generate(for: molecule)
```

Run layout before depiction unless the input already provides good 2D
coordinates that you want to preserve.

## 6. Compute Identifiers and Properties

```swift
import Foundation

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

Low-level descriptors remain available directly:

```swift
let xlogP = CDKXLogPDescriptor.calculate(for: laidOut)
let tpsa = CDKTPSADescriptor.calculate(for: laidOut, checkAromaticity: true)
let vabc = CDKVABCDescriptor.calculate(for: laidOut)
```

## 7. Generate SVG

Using the export facade:

```swift
import Foundation

var options = CDKFileExportOptions()
options.svgCanvasSize = CGSize(width: 1280, height: 840)
options.svgIncludeBackground = true

let svg = try CDKFileExporter.write(molecule: laidOut, as: .svg, options: options)
```

Using the lower-level depiction generator directly:

```swift
import Foundation

let svg = CDKDepictionGenerator.toSVG(
    molecule: laidOut,
    style: RenderStyle(),
    canvasSize: CGSize(width: 1280, height: 840),
    includeBackground: true
)
```

## 8. Build a Scene for a Host Renderer

```swift
import Foundation

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

Important: the package computes scene geometry only. The host owns the drawing
backend and event handling.

## 9. Configure Rendering

```swift
let style = RenderStyle(
    showCarbons: false,
    showImplicitHydrogens: true,
    atomColoringMode: .cdk2D,
    highlightStyle: .outerGlowWhiteEdge,
    colorBondsByAtom: true,
    aromaticDisplayMode: .circle,
    bondWidth: 2.2,
    fontSize: 18,
    padding: 24
)
```

Useful rendering features include:

- aromatic circle mode
- colored or glow highlights
- query atom and query bond depiction
- Markush legends and link-node annotations
- non-Markush Sgroup brackets and labels

## 10. Parse a Markush CXSMILES

```swift
import Foundation

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

Current CDK 2.12-style CX coverage in the supported path includes:

- atom labels and atom values
- coordinates
- radicals
- enhanced stereo state
- highlight layers
- ligand ordering
- polymer/data/generic/positional-variation Sgroups
- Markush `RG:` and `LN:`

## 11. Parse and Render Reactions

```swift
import Foundation

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
import Foundation

let hit = CDKMetalReactionDepictionSceneBuilder.participant(
    at: CGPoint(x: 320, y: 240),
    in: reaction,
    style: RenderStyle(),
    canvasRect: CGRect(x: 0, y: 0, width: 1200, height: 680),
    zoom: 1.0,
    pan: .zero
)
```

## 12. Export Common Formats

### Molecules

```swift
let smiles = try CDKFileExporter.write(molecule: laidOut, as: .smiles)
let isoSmiles = try CDKFileExporter.write(molecule: laidOut, as: .isomericSmiles)
let inchi = try CDKFileExporter.write(molecule: laidOut, as: .inchi)
let molV2000 = try CDKFileExporter.write(molecule: laidOut, as: .mol)
let molV3000 = try CDKFileExporter.write(molecule: laidOut, as: .molV3000)
let cml = try CDKFileExporter.write(molecule: laidOut, as: .cml)
```

### SDF with options

```swift
var options = CDKFileExportOptions()
options.sdfOptions = CDKSDFWriterOptions(alwaysV3000: true, programName: "MyHost")

let sdf = try CDKFileExporter.write(molecules: [laidOut], as: .sdf, options: options)
```

### Reactions

```swift
let rsmi = try CDKFileExporter.write(reaction: reaction, as: .smiles)
let cmlReaction = try CDKFileExporter.write(reaction: reaction, as: .cml)
let rxn = try CDKFileExporter.write(reaction: reaction, as: .rxn)
let rdf = try CDKFileExporter.write(reaction: reaction, as: .rdf)
```

## 13. Use the Package in Extensions

The package is suitable for:

- Quick Look preview or thumbnail targets
- Spotlight metadata importers
- helper tools and batch converters
- Linux CLI and service-side workflows

Those targets should import `CDKSwiftNativePort` directly instead of copying
chemistry code into the extension target.

## 14. Host Boundary Guidance

Keep host-specific concerns outside the package:

- Spotlight indexing policy
- Quick Look extension lifecycle
- `MTKView` ownership and render loops
- app entitlements, bundle IDs, document logic, and window/session state

Use `CDKSwiftNativePort` as the chemistry core and build host features around
its public APIs.
