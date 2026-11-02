# Finite compound-to-family matching

Status: available in release 1.4.38. The package is standalone and requires no
JVM runtime.

The separate [molecular query API](QUERY_MATCHING.md), available in 1.4.41,
adds SMARTS expressions, direct positional-variation matching and coordination
stereo. Its generic/substructure semantics do not change this API's closed finite
family contract or silently reinterpret unsupported families as ordinary queries.

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
- Ordinary alkene and neutral imine stereo compare the referenced substituents on both ends,
  including reversed endpoints and swapped reference substituents. Explicit
  cis/trans references and complete SMILES directional pairs are supported.
- Opposite declared configurations yield no-match if every compatible graph
  mapping fails. A compound missing a required stereo assignment is unresolved
  unless another permitted member gives a definite match. A known incompatible
  constraint can reject a mapping even when another constraint is unspecified.
- Tetrahedral depiction wedges must not create alkene stereo constraints.
- ABS tetrahedra preserve their drawn configuration. OR groups on the scaffold
  expand by joint inversion; different groups vary independently. A two-centre
  OR group produces two configurations, not four. The compound must specify one
  stereoisomer. Relative groups inside alternatives require explicit scope and
  remain unsupported. See the [RDKit enhanced stereo semantics](https://www.rdkit.org/docs/RDKit_Book.html#support-for-enhanced-stereochemistry).
- AND/racemic mixtures, wedge-only stereo without ordered ligands,
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
Graph search checks stereo for every compatible mapping. A witness includes
repeat counts, selected positional targets, inverted OR groups, R-site choices
and attachment IDs. `sourceByMemberAtomID` maps each generated atom to its original
scaffold/alternative atom, with the original site and choice for substituent atoms.
Repeat copies can share an original atom ID. `compoundAtomIDByMemberAtomID` maps
member atoms to the compound; suppressed ordinary H vertices are omitted.
Alternative indices are zero-based in the compiled query, not an arbitrary source
text numbering. Preserve the reviewed sidecar to explain excluded/nested choices.

## Finite family relations

`CDKMarkushRelations.compare(left, right, relation: .overlaps)` also accepts
`.containedIn`, `.equivalent` and `.disjoint`. Results are `holds`, `doesNotHold`
or `unresolved`. Overlap supplies a shared member; failed containment supplies an
uncovered source member. `direction` identifies which input produced the witness.
Equivalence checks both directions. One deadline and search-step budget covers
all traversals; the configuration limit applies to each family enumeration.
Partial stereo may leave containment unresolved. Existential overlap also tries
the opposite direction when a specified member can resolve missing stereo.

## Reviewed finite definitions and positional variation

Release 1.4.44 adds alternative-local scopes, explicit finite alkyl/alkoxy classes
and source-choice evidence. See [scoped finite definitions](FINITE_DEFINITIONS.md)
for the current API, schema compatibility and boundaries. Source prose is still
retained rather than automatically interpreted.

`CDKFiniteMarkushDefinition` is a Codable sidecar retaining source identity,
verbatim definition text, explicit fragment SMILES, excluded alternative indices,
occurrence rules, same-scaffold dependencies and REST_H. Compilation requires
`reviewed == true` and no unresolved conditions. It does not translate natural
language or expand unspecified chemical classes. Fragment SMILES start at their
attachment atom. Nested variables expand recursively with cycle/depth, aggregate
alternative and retained-atom bounds; replacement of that first attachment atom
is rejected. Explicit finite choices are not a claim to cover every chemical
class or source document.

Single-bond positional variation requires an explicit unlabelled dummy endpoint,
non-stereogenic targets with available hydrogen, and consistent metadata. Endpoint
IDs and chosen targets are retained. Ambiguous carbon placeholders, coordination
bonds and general multicentre bonding remain unsupported. The distinction follows
[Chemaxon position variation](https://docs.chemaxon.com/latest/marvin_position-variation-in-marvinsketch.html).

The corpus retains 103 source annotations with provenance/hashes. Ten finite
interpretations (M2S 5, 6, 54, 57, 69, 72, 75, 98, 99, 101) are executable.
Six more (8, 14, 37, 66, 68, 80) have reviewed unresolved reasons; 87 remain
unreviewed. Shared repeat labels, ring-atom replacement and omitted definitions
are explicitly distinguished from executable independent alternatives. This covers the supplied annotations, not omitted document context.
The tests include independent positive/negative compounds and full Standard
InChI checks of positive witnesses, neutral imine stereo,
nested/excluded alternatives and unresolved cases. General polymers, broad text
classes, radicals, mixtures and arbitrary nested CX definitions remain unsupported.

## Qualification

The release is qualified by complete macOS and Apple-container Ubuntu package
suites (650 tests on each platform, zero failures), strict formatting and package-boundary checks. Focused regressions cover
original-ID provenance, omission/repetition, finite relations, cancellation,
malformed metadata, pre-allocation atom limits and reviewed-source gates. Tests
consume attributed reference data without an upstream runtime.
