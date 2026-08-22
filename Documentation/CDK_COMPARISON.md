# Comparison With Original CDK

This document compares `CDKSwiftNativePort` with upstream Java CDK.

Upstream reference: <https://github.com/cdk/cdk>

`CDKSwiftNativePort` includes code derived from and inspired by the Chemistry
Development Kit (CDK). The current parity target for the supported surface in
this package is CDK `2.13` (`cdk-2.13`, commit
`3bc6efe84119e5da0d00804bd4cf5f3bdd9c3d6d`).

## CDK Reference Citations

- Willighagen et al. (2017), *The Chemistry Development Kit (CDK) v2.0: atom
  typing, depiction, molecular formulas, and substructure searching*.
  DOI: [10.1186/s13321-017-0220-4](https://doi.org/10.1186/s13321-017-0220-4)
- May and Steinbeck (2014), *Efficient ring perception for the Chemistry
  Development Kit*.
  DOI: [10.1186/1758-2946-6-3](https://doi.org/10.1186/1758-2946-6-3)
- Steinbeck et al. (2006), *Recent Developments of the Chemistry Development
  Kit (CDK) - An Open-Source Java Library for Chemo- and Bioinformatics*.
  DOI: [10.2174/138161206777585274](https://doi.org/10.2174/138161206777585274)
- Steinbeck et al. (2003), *The Chemistry Development Kit (CDK): An Open-Source
  Java Library for Chemo- and Bioinformatics*.
  DOI: [10.1021/ci025584y](https://doi.org/10.1021/ci025584y)

For the full attribution and redistribution notice, see [`../NOTICE.md`](../NOTICE.md).

## Scope Summary

- **Original CDK** is a broad Java chemistry toolkit with a large module set.
- **CDKSwiftNativePort** is a Swift-native port that focuses on the workflows
  needed by native Apple hosts: parsing, layout, depiction, identifiers,
  descriptors, and major file formats.

That does not mean every CDK module is ported.

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
- CXSMILES layers used by CDK 2.13-backed host workflows:
  - atom labels and atom values
  - coordinates
  - radicals
  - enhanced stereo metadata
  - fragment grouping
  - `ha:` / `hb:` highlights
  - ligand ordering
  - coordinate bonds (`C:`)
  - Markush `RG:` and link-node `LN:`
  - polymer/data/generic/positional-variation Sgroups
- MDL Molfile V2000 import/export
- MDL Molfile V3000 import/export, including `HILITE`
- MDL V2000/V3000 reacting-center status round-trips
- R0/bare-R Markush roots and partial R-group definitions
- SDF import/export with SD data fields and mixed V2000/V3000 records
- RXN and RDF support
- CML molecule and reaction support
- InChI integration
- Universal-SMILES AuxInfo numbering with mixed fixed-H `/F:` components
- official-reference InChI parity harness for a curated upstream InChI CI corpus
- CDK-derived 2D layout
- small-ring-size queries and coordinate-bonded sandwich-ring layout
- circular ECFP/FCFP fingerprints, count/folded bit fingerprints, and
  Tanimoto/Dice/Cosine similarity helpers
- molecule-vs-molecule substructure matching over the Swift `Molecule` graph
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

- generalized SMARTS query parsing and the broader SMARTS module family
- non-circular fingerprint families such as MACCS, PubChem, shortest-path, and
  atom-pairs fingerprints
- force-field and 3D generation workflows
- the broader long tail of CDK QSAR descriptors
- every configurable reader/writer mode from upstream CDK

Reaction CML now has native Swift hierarchy types for sets, lists, schemes,
and list entries, and the reader/writer preserve branching and step grouping
through hierarchy round-trips instead of collapsing everything to flat reaction
arrays. The package now also exposes native reaction, reaction-set, and
reaction-scheme manipulator layers for counts, reversal, mapped-object lookup,
inline reaction-to-molecule conversion, atom-property propagation, basic atom
configuration helpers, and topology/path extraction, plus RXN V2000
strict/relaxed counts-line behavior, RXN V3000 parsing/writing, and
multi-record reaction SMILES set import. That maintained reaction
model/IO/manipulator surface is now tracked by a zero-gap strict parity gate in
[`REACTION_PARITY.md`](REACTION_PARITY.md).

What remains outside that reaction parity claim is a different layer of CDK:

- the `base/reaction` reaction-engine, mechanism, and reaction-type modules

`RInChI` is no longer outside the maintained Swift surface. The package now
includes a native `RInChI` implementation plus a vendored upstream-derived
reference gate documented in [`RINCHI_PARITY.md`](RINCHI_PARITY.md).

## Parity Tracking

Parity-oriented metadata lives under:

- `Tests/CDKSwiftNativePortTests/**/port_metadata.json`

Those metadata files link the current Swift test surface back to the upstream
CDK 2.13 tests that informed the port.

For InChI specifically, reference-grade work is tracked separately in
[`INCHI_REFERENCE_GRADE.md`](INCHI_REFERENCE_GRADE.md). The runtime remains
Swift native, but the package now vendors an official InChI reference corpus
and uses it as a parity oracle.

Reaction hierarchy and reaction IO parity are tracked separately in
[`REACTION_PARITY.md`](REACTION_PARITY.md). That harness vendors an
upstream-derived reaction corpus and keeps a committed known-gap inventory plus
an opt-in strict parity gate for reaction CML, RXN V2000/V3000, RDF, reaction
SMILES, and the maintained reaction manipulator surface.
