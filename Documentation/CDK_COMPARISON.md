# Comparison With Original CDK

This document compares `CDKSwiftNativePort` with upstream Java CDK.

Upstream reference: <https://github.com/cdk/cdk>

## Scope Summary

- **Original CDK**: a broad Java chemistry toolkit with extensive modules for
  IO, structure handling, descriptors, rendering, atom typing, fingerprints,
  substructure search, and more.
- **CDKSwiftNativePort**: a Swift-native package that ports and re-expresses the
  subset of CDK functionality needed for current native-host workflows.

The current parity target for this port is CDK `2.12` for the areas covered by
the package today. That does not imply full port parity for every CDK module.

## Side-by-Side Overview

| Area | Original CDK (Java) | CDKSwiftNativePort (Swift) |
|---|---|---|
| Runtime | JVM | Native Swift on Apple platforms |
| Core model | `IAtom`, `IBond`, `IAtomContainer`, interface-heavy graph model | `Atom`, `Bond`, `Molecule` value types |
| API style | Java OO interfaces and factories | Swift value types plus focused facade APIs |
| SMILES | Mature parser/generator stack | CDK-style parser/generator, reaction SMILES, CXSMILES state handling |
| CXSMILES Markush | CDK 2.12 supports `RG:` / link-node depiction path | Ported support for the currently implemented Markush / R-group path |
| InChI | CDK wrapper APIs | Native Swift-facing parse/generate facades |
| File IO | Broad reader/writer module set | Unified importer/exporter plus focused reader/writer coverage |
| 2D layout | Structure diagram generator | `Depiction2DGenerator` backed by the Swift port |
| Rendering | CDK depiction infrastructure | SVG output, SwiftUI rendering, Metal scene primitives |
| Descriptors | Very broad | Focused descriptor set used by current hosts |
| Host integration | N/A | Explicitly excludes host-app wiring |

## Areas With Strong Current Coverage

- Molecule graph model and stereochemistry
- SMILES parsing and generation
- reaction SMILES parsing
- CXSMILES parsing used by the supported depiction path
- InChI parse/generate integration
- MDL Molfile / SDF / RXN / RDF support
- MOL2, PDB, XYZ, and CML IO
- CDK-derived 2D layout
- SVG / SwiftUI / Metal-scene depiction
- identifier and selected descriptor services

## CDK 2.12-Specific Highlights in This Port

The current branch includes CDK 2.12-derived work for the supported Markush
path, including:

- CXSMILES `RG:` R-group parsing
- link-node / repeat-unit handling used by the current depiction flow
- Markush legend layout and depiction
- host-friendly Metal scene behavior for rotating scaffolds while keeping the
  Markush legend fixed

## Intentional Differences

- The API surface is Swift-first rather than a Java API mirror.
- The package favors a small number of host-facing facade types
  (`CDKFileImporter`, `CDKFileExporter`, `CDKMoleculeIdentifierService`,
  `CDKMoleculePropertyService`) over exposing every internal port layer.
- Rendering is expressed in outputs that fit native Apple hosts:
  SVG strings, SwiftUI `GraphicsContext`, and renderer-neutral scene primitives.
- Host integration concerns such as Spotlight, Quick Look, document windows,
  and bundle metadata stay outside the package.

## Not Yet a Full CDK Port

The following CDK areas are not exposed as package-wide parity targets today:

- generalized SMARTS and substructure search modules
- fingerprint / similarity modules
- force-field / 3D generation workflows
- full atom typing / matcher breadth
- the broader long tail of CDK QSAR descriptors
- complete parity for every CXSMILES extension outside the supported paths

These areas should be treated as out of scope unless explicitly added.

## Parity Tracking

Parity-oriented tests and metadata references live in:

- `Tests/CDKSwiftNativePortTests/**/port_metadata.json`

These metadata files document links back to upstream CDK tests where applicable
and help keep the supported feature surface traceable to the original toolkit.
