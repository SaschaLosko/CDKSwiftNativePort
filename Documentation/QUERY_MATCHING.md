# Molecular query matching

Release 1.4.41 adds native molecular SMARTS expressions, direct variable-attachment
matching and ordered coordination stereo. The reference is CDK 2.13, commit
`3bc6efe84119e5da0d00804bd4cf5f3bdd9c3d6d`. This is a documented subset of the
upstream query surface, not a claim that the entire SMARTS module is ported.

```swift
let query = try CDKSmartsParser.parse("[C;$(*=O)]-[O,N]")
let target = try CDKSmilesParser(generateCoordinates: false).parseSmiles("CC(=O)N")
let result = CDKQueryMatcher.match(query, in: target)
// result.verdict == .match; atomMapping identifies a compatible embedding.
```

## Expressions and results

`CDKAtomExpression` supports element, aromatic/aliphatic, isotope, formal charge,
explicit/heavy/total degree, valence, attached/implicit hydrogen, ring membership,
ring count, smallest ring and ring-bond count predicates. AND, OR, NOT and anchored
recursive query graphs can be composed programmatically or read from molecular
SMARTS. Text precedence is NOT, high-precedence AND/adjacency, OR, low-precedence
AND. Molecular atom-map numbers are annotations, not matching constraints.

`CDKBondExpression` supports bond order, any bond, ring/chain membership and
Boolean composition. An omitted SMARTS bond means single-or-aromatic. Branches,
ring closures, disconnected patterns and recursive `$()` environments are supported.
Ring-count predicates use the package's minimum cycle basis. Aromaticity and
hydrogen perception use the existing package models; full upstream aromaticity
model parity is not asserted. Valence queries use a localized Kekule assignment.

Results are `match`, `noMatch` or `unresolved`. A match contains an atom-ID mapping.
A negative is emitted only after complete traversal. Missing required stereo,
cancellation, invalid/unsupported input and resource limits remain unresolved.
Recursive predicates use three-valued Boolean logic: negating an unresolved
stereo query does not prove a positive. Unconstrained stereo permits either
configuration. `isExhaustive` reports traversal separately from the verdict.

The default mode is substructure matching. `.wholeGraph` additionally requires
all target vertices and bonds to be covered, after resolving attachment placeholders.
Explicit hydrogen vertices remain graph vertices in this API. This is distinct
from the hydrogen-normalized, closed finite-family semantics of `CDKMarkushMatcher`.
Neither a substructure hit nor a generic query hit is automatically finite
R-group family membership, containment, or a result about a source document.

Defaults: 2,000 atoms per graph, one million shared search steps, a cooperative
30-second deadline, at most 16 recursive queries, and bounded expression/parser
depth. Query text is limited to 32 KiB and 1,000 atoms. Recursive and attachment
searches share the budget. Graph preprocessing and individual perception routines
are not preemptible; the deadline is not a hard real-time guarantee.

## Imported query graphs and variable attachments

`try CDKQueryPattern(molecule: imported)` adapts supported MDL/CX atom lists,
atom/bond query types, ring topology, substitution counts and REST_H constraints.
Unsupported decorators are rejected instead of discarded. R-group definitions
must be expanded first. Bare text labels are not executable definitions.

An `ExtMulticenter`/`m` Sgroup identifies an explicit `*` placeholder followed by
its possible attachment atom IDs. The matcher skips these placeholders during
the initial embedding and resolves their bonds directly against the target:

- one or both bond ends can vary;
- query-bond predicates constrain the selected bond;
- fixed and variable bonds cannot consume the same target edge;
- overlapping choice sets are backtracked, so a failed early choice does not
  hide a later valid placement;
- the witness maps every placeholder to its selected target endpoint;
- metadata inconsistencies and chemically constrained placeholders are rejected.

This implements positional variation, not simultaneous multicentre coordination
bonding. Target coordinate bonds require a separate bonding model and remain unsupported.

## Stereo

Tetrahedral and geometric stereo constrain each complete mapping. Tetrahedral
ligand order, reordered traversal and ordinary explicit/implicit hydrogen are
handled. Relative tetrahedral groups must invert together. Imported AND/OR group
labels are used for relative query matching; this does not enumerate or identify
an individual compound from a racemic mixture.

`CDKQueryAtomStereo` additionally models square-planar (SP1–3), trigonal-bipyramidal
(TB1–20) and octahedral (OH1–30) orderings. Matching uses the upstream permutation
classes, not the literal configuration number. Explicit partial coordination
queries may repeat the focus ID for unspecified ligands. Square-planar/octahedral
normalization and the upstream limited TB/OH conversion are supported when only
the permitted subset of ligands is constrained.

```swift
let query = try CDKSmartsParser.parse("[NH3][Pt@SP1]([NH3])(Cl)Cl")
let target = try CDKQueryStructure.parseSmiles("[NH3][Pt@SP2](Cl)([NH3])Cl")
let result = CDKQueryMatcher.match(query, in: target) // match: two cis representations
```

`CDKQueryStructure` retains the molecule and ordered stereo sidecar together,
including through Codable. Use its SMILES reader for coordination-stereo targets.
The ordinary `Molecule` SMILES reader models only tetrahedral atom parity and now
rejects unsupported classes rather than silently converting them to tetrahedra.
Do not export only the molecule portion and expect the sidecar stereo to survive.
The existing depiction and general molecule writers do not consume this sidecar.

## Explicit remaining boundaries

Reaction SMARTS, component-level grouping, allene/extended stereo, arbitrary
Boolean combinations of stereo primitives, optional directional bond expressions,
directional ring closures and additional vendor primitives are rejected. The
SMARTS parser requires explicit ligand counts for coordination classes; partial
coordination uses the programmatic stereo API. This release does not add general
polymers, natural-language class expansion or arbitrary nested CX definitions.

## Reference evidence

The tests include independent examples from the
[Daylight SMARTS specification](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html)
and the CDK stereo class documentation, positive/negative/unknown results,
recursive anchoring, explicit/implicit hydrogen, joint group inversion,
cross-geometry partial constraints, two-ended variable bonds and overlapping
attachment choices requiring backtracking. Every SP/TB/OH carrier permutation
(24 + 120 + 720) is checked against the attributed upstream table data, including
rejection of all other configuration classes for distinct ligands.

Sources: [Expr](https://github.com/cdk/cdk/blob/cdk-2.13/base/isomorphism/src/main/java/org/openscience/cdk/isomorphism/matchers/Expr.java),
[VarAttachFilter](https://github.com/cdk/cdk/blob/cdk-2.13/base/isomorphism/src/main/java/org/openscience/cdk/isomorphism/VarAttachFilter.java),
[QueryStereoFilter](https://github.com/cdk/cdk/blob/cdk-2.13/base/isomorphism/src/main/java/org/openscience/cdk/isomorphism/QueryStereoFilter.java),
[SquarePlanar](https://github.com/cdk/cdk/blob/cdk-2.13/base/core/src/main/java/org/openscience/cdk/stereo/SquarePlanar.java),
[TrigonalBipyramidal](https://github.com/cdk/cdk/blob/cdk-2.13/base/core/src/main/java/org/openscience/cdk/stereo/TrigonalBipyramidal.java),
[Octahedral](https://github.com/cdk/cdk/blob/cdk-2.13/base/core/src/main/java/org/openscience/cdk/stereo/Octahedral.java).
The package consumes reference data and contains no upstream-language source or runtime dependency.

Release qualification on 2026-11-02: the complete macOS and Apple-container
Ubuntu suites each passed 667 tests with zero failures. Sixteen focused query
tests exercise the cases above; an additional ordinary SMILES regression checks
biaryl bond preservation and stereo-class rejection. Strict swift-format lint
passed for every changed Swift file. No application/customer references or
upstream-language source/build artifacts were added.
