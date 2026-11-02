# Finite compound-to-family matching

Status: implemented locally in the unreleased package line. No host application
or Java runtime is required.

```swift
let parser = CDKSmilesParser(flavor: .cdkDefault)
let family = try parser.parseSmiles("C* |$;R1$,RG:_R1={O},{N}|")
let compound = try parser.parseSmiles("CO")
let result = CDKMarkushMatcher.match(compound, in: family)

switch result.verdict {
case .match:
    // result.witness contains a permitted member, repeat counts and atom mapping.
    break
case .noMatch:
    // Every permitted finite assignment was checked, with no compatible member.
    assert(result.isExhaustive)
case .unresolved:
    // Present result.reason/detail; never treat an incomplete search as absence.
    break
}
```

## Membership contract

This is membership in the closed family described by explicit finite alternatives,
occurrences, dependencies and supported repeat ranges. The complete compound
graph must match a generated member. A common substructure does not suffice.
All disconnected components, formal charges, isotope assignments, bond orders
and hydrogen counts participate. Neutral terminal ordinary hydrogen vertices
are normalized to attached hydrogen counts; isotopic hydrogen remains explicit.
D/T normalize to hydrogen isotopes. Aromatic and perceived equivalent Kekule
representations are compared using the package's existing aromaticity model.
There is no tautomer, protonation, salt-stripping or general chemical-standardization
step, and no identifier-only shortcut.

REST_H and dependency/occurrence rules are applied by enumeration. Whole-graph
matching prevents a non-hydrogen substituent from occupying an omitted site.
REST_H=false does not invent unlisted substituents: this API is closed finite
membership, not open-ended generic substructure matching.

A positive result has a concrete witness and can stop early. A negative requires
complete enumeration and graph search. Limits, cancellation, unsupported input
and missing required compound stereochemistry produce `unresolved`. A logically
empty, well-formed finite family can produce an exhaustive negative. Missing
substituent definitions are not an empty family.

`isExhaustive` describes traversal, independently of the verdict: an unresolved
stereo comparison may still have examined every assignment. Positive witnesses
do not claim exhaustive traversal. Counts refer to assignments, including
symmetry-equivalent duplicates, rather than unique chemical identities.

## Variable labels

- Numeric R labels are nonnegative Swift integers, including R42 and R100. They
  retain their original numbers. Bare R denotes R0.
- Named pseudoatom labels such as X, Y, B, Rd and A13 can be bound to explicitly
  supplied finite CX `RG:` alternatives. Actual elemental B and Y remain elements.
- `CDKMarkushVariableLabels.resolve` returns a normalized copy and an explicit
  `numberBySourceLabel` mapping. Named variables receive unused positive numbers
  in sorted-label order, avoiding existing numeric sites and logic keys.
- The matcher reports that mapping as `variableNumberBySourceLabel`. Input
  structures, IDs and source labels are not modified. Ambiguous numeric aliases,
  negative/out-of-range numbers and missing definitions remain unresolved.
- The lower-level numeric query adapter can consume the explicitly resolved copy.
  Numeric occurrence/dependency logic belongs to that numeric representation;
  textual named-variable conditions are not interpreted automatically.
- RGfile export retains its R0–R32 format boundary and rejects out-of-domain
  queries rather than silently renumbering them. Numeric labels above 32 can
  round-trip with the Markush-preserving CX writer. Named-label CX export is
  still separate work; callers must retain source labels and the resolution map
  if exporting an explicitly normalized numeric representation.

The M2S backbones 35 and 47 now adapt without the former number-limit failure.
They still lack executable substituent definitions and therefore cannot supply
membership verdicts on their own.

## Stereochemical policy

Declared family stereochemistry constrains membership. Unspecified family
stereochemistry is unconstrained; the API does not enumerate every stereoisomer.

- Tetrahedral centers compare ligand permutation parity under the graph mapping,
  not the literal `@`/`@@` token. Implicit hydrogen and explicit ordinary hydrogen
  references are normalized consistently. Reordered atom IDs and traversal order
  do not change the result.
- Ordinary alkene stereo compares the referenced substituents on both ends,
  including reversed endpoints and swapped reference substituents. Explicit
  cis/trans references and complete SMILES directional pairs are supported.
- Opposite declared configurations yield no-match if every compatible graph
  mapping fails. A compound missing a required stereo assignment is unresolved
  unless another permitted member gives a definite match. A known incompatible
  constraint can reject a mapping even when another constraint is unspecified.
- Tetrahedral depiction wedges must not create alkene stereo constraints.
- Enhanced AND/OR/racemic groups, wedge-only stereo without ordered ligands,
  incomplete alkene directions and unsupported stereo types remain unresolved.
  Coordinates alone are not used to infer absolute stereochemistry.

Focused ordering examples follow the [OpenSMILES stereochemistry specification](https://opensmiles.org/opensmiles.html#section3.8).
This is targeted conformance coverage, not a claim of complete upstream stereo
or generic-query parity. CX fragment merging now remaps tetrahedral ligand IDs
and retains pseudoatom/stereo properties before enumeration.

## Bounds, witnesses and supported scope

Default limits are 1,000 configurations, 2,000 atoms in each input/member,
1,000,000 shared enumeration/search steps and a 30-second cooperative monotonic
deadline. Task cancellation is checked during enumeration and graph traversal.
Both enumerators expose a throwing `checkpoint` callback for caller-owned bounds.
Normalization and individual graph operations are not preemptible, so the
cooperative deadline is not a hard real-time latency guarantee.

The matcher streams members instead of retaining the entire expanded family.
It compares atom attributes and connectivity with bounded graph backtracking,
checking stereochemistry for each complete compatible mapping. A witness
contains the generated member and repeat counts, plus a member-to-compound atom
ID map. Suppressed ordinary hydrogen vertices are omitted from that map; it is
not a mapping back to all source-family alternative IDs.

Supported families are the finite scaffold repeats and flat R-group definitions
accepted by `CDKMarkushEnumerator`. Textual definitions, chemical classes,
exclusions, nested variables, positional/multicenter variation, general polymers,
radicals and coordinate bonds remain unresolved. Family overlap, containment,
equivalence and repository search integration are separate capabilities.

## Qualification

The 19 focused matching tests cover positive/negative membership, witnesses,
CO/CN alternatives, all nine repeat-example members, source preservation,
limits/cancellation, whole-graph and hydrogen/isotope/charge behavior, labels,
REST_H/dependencies, malformed input and the stereochemical cases above.

Full standalone macOS and Apple-container Ubuntu suites: 630 tests each, zero
failures and zero skips. The macOS suite runs from an isolated package snapshot.
The removed host-checkout tests are replaced by package-scoped boundary checks.
Strict swift-format lint passes for all changed Swift files. Validation logs are
retained outside version control in `/private/tmp/cdk-markush-matching/`.
