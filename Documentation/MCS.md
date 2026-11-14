# Portable maximum common substructure

`CDKMCS` implements a bounded, portable exact reference search inspired by the
[CDK 2.14.1 MCS contract](https://github.com/cdk/cdk/blob/cdk-2.14.1/base/isomorphism/src/main/java/org/openscience/cdk/isomorphism/MCS.java).
It is implemented in Swift and needs no upstream runtime or native search dependency.
The package's broader baseline remains CDK 2.13 with selected later changes.

```swift
let parser = CDKSmilesParser(generateCoordinates: false)
let query = try parser.parseSmiles("CC(=O)NC")
let target = try parser.parseSmiles("CCC(=O)N")
let result = try CDKMCS(query: query).search(in: target)
// Complete: four atoms and three bonds, with one maximum mapping.
if result.isMaximumProven {
    let witness = result.mappings.first
    // witness maps query atom/bond IDs to target IDs, not array offsets.
}
```

## Search and result contract

- Connected search maximizes atoms, then common bonds. Every compatible bond
  between mapped endpoints is common; this is not induced matching.
- `minimumBondsPerFragment: n` selects disconnected common-edge search, maximizing
  bonds then atoms. Every connected fragment must contain at least `n` common
  bonds. Isolated atoms are excluded. The option is not a maximum fragment count.
- Immutable query/options snapshots may be shared; each call owns its search state.
- `search` returns best-so-far mappings with `completed`, `timedOut`, `searchLimit`
  or `cancelled`. Only completed traversal proves the maximum score. An empty
  interrupted result cannot establish that nothing is common.
- The mapping cap limits retained equal-best witnesses, not traversal. A completed
  result can prove its maximum while `mappingsTruncated` is true. Use
  `allMaximumMappingsReturned` when all mappings matter.
- `match` requires maximum proof and throws `incomplete` otherwise. Invalid or
  unsupported graphs/policies throw; they are not negative chemical conclusions.
- Defaults are 1,000 mappings, 256 atoms per input, 1,000,000 search steps and
  30 seconds. Limits include predicate evaluation and relevant-ring preparation;
  cancellation is cooperative. The atom cap is at most 256 and the retained
  mapping cap at most 100,000. The bounded defaults intentionally differ from
  upstream's unbounded default duration and much larger graph admission limit.

## Matching semantics

Default atoms match by symbol (pseudo atoms only other identical symbols); bonds
match by stored order, including aromatic order. Perceive aromaticity consistently
before searching. `atomMatching: .element` makes a query pseudo atom a wildcard;
`.any`, bond-order relaxation, charge/isotope constraints and computed ring/chain
constraints are explicit options. Query isotope constraints apply where supplied.
Explicit hydrogens count as vertices; implicit hydrogens do not.

A `CDKQueryPattern` uses its own atom/bond expressions, including Boolean and
anchored recursive predicates. Top-level SMARTS stereo is not checked; recursive
predicates retain the query engine's stereo rules. Unresolved recursive stereo
is rejected, not interpreted as false. Molecule matching overrides are rejected
for pattern queries. Family definitions, polymers, link nodes and positional
variation require expansion; MCS is not compound-to-family membership.

`completeRings` requires every common ring atom and bond to participate in a
whole relevant ring in both graphs. Relevant rings are the union of minimum cycle
bases, not one arbitrary SSSR. The implementation enumerates simple cycles within
the shared budget and tests independence against strictly shorter cycles over
GF(2). Benzene can match a naphthalene ring; the six-membered perimeter of norbornane
is not a relevant ring. Ring preparation never mutates source flags.

`matchStereochemistry` compares ordered tetrahedral and double-bond configurations
as far as common neighbors fix them. Three common neighbors fix a tetrahedral
center. Both alkene ends and a neighbor of each fix an alkene. Unspecified stereo
on one side remains compatible. Explicit/implicit H carriers are normalized;
relative and racemic groups permit joint inversion, not independent inversion of
each center. This does not promise allene, atropisomer or coordination-stereo MCS.
Unsupported ordered elements can be rejected rather than silently reinterpreted.

## Evidence and algorithm boundary

The reference enumerates injective partial mappings with safe atom-count pruning,
then checks the common-edge graph, ring and stereo constraints. It is exponential;
it is useful for bounded comparisons and independent qualification, not a claim of
native-engine speed on large or highly symmetric chemistry. No native accelerator
is selected implicitly, and no heuristic result is upgraded to a proof.

The release passes 725 tests on macOS and Ubuntu. Its 16 MCS tests include:

- 35 selected upstream score/mapping-count fixtures from
  [MCSTest](https://github.com/cdk/cdk/blob/cdk-2.14.1/tool/smarts/src/test/java/org/openscience/cdk/isomorphism/MCSTest.java),
  plus the public amide example; fixture attribution retains LGPL-2.1-or-later.
- 512 connected/disconnected comparisons with an independently implemented
  exhaustive subset/permutation oracle over all four-vertex source topologies.
- Fused, bridged and spiro rings; relevant-cycle multiplicity; absolute/relative
  stereo, explicit/implicit H, permutations, unspecified and negative controls.
- Predicate/recursive SMARTS, source nonmutation, mapping-cap semantics, bounded
  larger graphs and active cancellation with retained witnesses.

This establishes those cases and the stated API contract. It is not a complete
upstream MCS test-suite port, all-molecule parity or a performance replacement for
specialized native engines. Larger corpus qualification and stronger branch-and-bound
acceleration can be added without weakening the result contract.
