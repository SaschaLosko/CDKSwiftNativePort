# Comparison With Original CDK

This document compares `CDKSwiftNativePort` with upstream Java CDK.

Upstream reference: <https://github.com/cdk/cdk>

## Scope Summary

- **Original CDK** is a broad Java chemistry toolkit with a large module set.
- **CDKSwiftNativePort** is a Swift-native port that focuses on the workflows
  needed by native Apple hosts: parsing, layout, depiction, identifiers,
  descriptors, and major file formats.

The current parity target for the supported surface in this package is CDK
`2.12`. That does not mean every CDK module is ported.

## Side-by-Side Overview

| Area | Original CDK (Java) | CDKSwiftNativePort (Swift) |
|---|---|---|
| Runtime | JVM | Native Swift |
| Core model | interface-heavy `IAtomContainer` ecosystem | Swift value types (`Molecule`, `Atom`, `Bond`) |
| API style | Java OO interfaces and factories | Swift-first facades plus lower-level readers/writers |
| SMILES | mature parser/generator stack | SMILES, reaction SMILES, CXSMILES |
| MDL | V2000, V3000, SDF, RXN, RDF | same major formats with round-trip coverage |
| CML | broad molecule/reaction support | molecule and reaction CML support |
| Layout | structure diagram generator | `Depiction2DGenerator` |
| Rendering | CDK depiction infrastructure | SVG plus renderer-neutral scene generation |
| Descriptors | extensive descriptor set | focused descriptor/property set used by hosts |
| Host integration | not applicable | explicitly excludes app/extension implementation |

## Areas With Strong Current Coverage

- molecule and reaction model types
- SMILES and reaction SMILES parsing/generation
- CXSMILES layers used by CDK 2.12-backed host workflows:
  - atom labels and atom values
  - coordinates
  - radicals
  - enhanced stereo metadata
  - fragment grouping
  - `ha:` / `hb:` highlights
  - ligand ordering
  - Markush `RG:` and link-node `LN:`
  - polymer/data/generic/positional-variation Sgroups
- MDL Molfile V2000 import/export
- MDL Molfile V3000 import/export, including `HILITE`
- SDF import/export with SD data fields and mixed V2000/V3000 records
- RXN and RDF support
- CML molecule and reaction support
- InChI integration
- CDK-derived 2D layout
- depiction support for:
  - aromatic display modes
  - Markush legends
  - Sgroup brackets and labels
  - query atoms and query bonds
  - highlight rendering modes

## Intentional Differences

- The API is Swift-first instead of mirroring CDK’s Java interfaces.
- The package emphasizes a small set of façade APIs rather than exposing every
  internal implementation detail as a first-class integration surface.
- Rendering output is package-owned as SVG or neutral scene data, not as a Java
  renderer object graph.
- Host concerns such as Quick Look, Spotlight, document state, and app UI stay
  outside the package.

## macOS-Oriented Differences

Compared with upstream CDK, this port adds host-facing concerns that are useful
for native Apple apps:

- coordinated file access via `CDKFileAccess`
- renderer-neutral scene generation for host Metal/Core Graphics pipelines
- package boundaries that make reuse from app extensions practical

Those are integration conveniences around the chemistry core, not attempts to
reproduce the original Java runtime model.

## Remaining Scope Limits

The following CDK areas are still outside the package’s current parity target:

- generalized SMARTS and broad substructure-search modules
- fingerprint/similarity toolkits
- force-field and 3D generation workflows
- the broader long tail of CDK QSAR descriptors
- first-class reaction-scheme object models
- every configurable reader/writer mode from upstream CDK

For reaction CML specifically, reaction schemes and lists can be read, but they
are flattened to `[CDKReaction]` because the Swift package does not currently
define a separate reaction-scheme type.

## Parity Tracking

Parity-oriented metadata lives under:

- `Tests/CDKSwiftNativePortTests/**/port_metadata.json`

Those metadata files link the current Swift test surface back to the upstream
CDK 2.12 tests that informed the port.
