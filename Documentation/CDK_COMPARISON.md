# Comparison With Original CDK

This document compares `CDKSwiftNativePort` with upstream Java CDK.

Upstream reference: <https://github.com/cdk/cdk>

## Scope Summary

- **Original CDK**: broad Java chemistry toolkit covering many domains and algorithms.
- **CDKSwiftNativePort**: Swift-native port focused on the feature set required by AtomLens workflows (parsing, layout/depiction, identifiers, and core molecular properties).

`CDKSwiftNativePort` targets practical parity for supported workflows, not complete API/module parity with every CDK package.

## Side-by-Side Overview

| Area | Original CDK (Java) | CDKSwiftNativePort (Swift) |
|---|---|---|
| Runtime | JVM | Native Swift (no JVM dependency) |
| Primary model | `IAtomContainer` ecosystem | `Molecule` / `Atom` / `Bond` value types |
| API style | Java OO APIs and interfaces | Swift value-oriented APIs + facades |
| SMILES | Mature parser/generator stack | CDK-style parser/generator + CXSMILES state support |
| InChI | CDK wrappers and integrations | Native Swift-style InChI parse/generate facades |
| File IO | Many readers/writers and modules | Unified importer/exporter + supported format readers/writers |
| 2D layout | CDK structure diagram APIs | `Depiction2DGenerator` (CDK-style layout) |
| Rendering | CDK depiction components | SVG export + SwiftUI `GraphicsContext` + Metal scene builders |
| Descriptor breadth | Very broad QSAR/descriptors | Focused descriptor set used by host app workflows |
| Host integration | N/A (library only) | Explicitly excludes host app wiring (Spotlight/Quick Look/etc.) |

## Supported Feature Families in This Port

- Molecule graph model, bond order/stereo, aromatic display helpers.
- SMILES parsing/generation, reaction SMILES parsing.
- InChI import/export APIs and identifier services.
- MDL/SDF, MOL2, PDB, XYZ, CML, RXN, RDF reader/writer coverage for current workflows.
- CDK-style 2D layout and depiction exports.
- Property services and selected descriptor calculations.

## Intentional Differences

- API naming and ergonomics are Swift-first.
- Rendering outputs are tailored to native Apple host stacks (SVG, SwiftUI context, Metal scene primitives).
- The package keeps a strict boundary from host application targets and extension code.
- Not all CDK modules are ported; unsupported areas should be treated as out-of-scope until explicitly added.

## Parity Tracking

Parity tests and metadata references live in:

- `Tests/CDKSwiftNativePortTests/**/port_metadata.json`

These metadata files document links back to upstream CDK tests where applicable.
