# Native R-group enumeration

Available in CDKSwiftNativePort 1.4.35. The upstream reference is Java CDK 2.13
(`cdk-2.13`, commit `3bc6efe84119e5da0d00804bd4cf5f3bdd9c3d6d`), principally
`RGroupQuery.getAllConfigurations`, `RGroupQueryReaderTest`, and
`RGroupQueryManipulatorTest`.

```swift
let query = try CDKRGFileReader.read(text: rgfile)
let members = try query.getAllConfigurations(maximumConfigurations: 10_000)

// Or consume one configuration at a time without storing the entire family.
let delivered = try query.enumerateConfigurations(maximumConfigurations: 10_000) { member in
    consume(member)
    return true // false stops deliberately, after delivering this member
}
```

For a flattened molecule containing R-group definitions, first use
`CDKRGroupQueryManipulator.toRGroupQuery(molecule)`. CXSMILES `RG:` definitions,
explicit `_AP1`/`_AP2` dummy atoms and `LO:` ordering are supported, as are RGfile
`M LOG`, `M AAL`, and `M APO` representations.

## Semantics

- Exact counts, comma-separated alternatives, inclusive ranges, `<n` and `>n`
  constrain the number of occupied sites for each R number.
- Every allowed positional distribution and substituent assignment is generated.
  Symmetry-equivalent molecules remain separate, as in Java CDK. Output is
  deterministic in ascending R number, occurrence, and scaffold atom order;
  alternatives retain their supplied order. Java's HashMap iteration order is
  not part of the compatibility contract.
- Required-R-group conditions apply only when the requiring group is present.
  Multiple dependencies and cycles are supported. Unsatisfiable occurrences
  produce an empty result; missing definitions are errors.
- Attachments follow explicit AAL/LO order, otherwise scaffold atom order.
  The first bond uses AP1, later bonds use AP2, including CDK's permissive
  three-or-more-bond behavior. APO=3 uses the same real atom for AP1 and AP2.
- Explicit AP dummy atoms are removed, real attachment atoms retained, bond
  orders and endpoint orientation preserved, IDs made unique, and stereo
  reference IDs remapped. Substituent 2D coordinates are translated to the site.
- Omitted sites set `Atom.rGroupRestH` on their neighbors, including an explicit
  false value, and increase a supplied hydrogen count by the removed bond order.
  Unspecified hydrogen counts continue to be inferred from the resulting graph.
  `rGroupRestH` is a query-matching annotation backed by the existing atom
  properties dictionary, not an extra hydrogen count. It survives Codable.
- Detached groups are supported. Matching CDK's lenient RGfile enumeration,
  a missing APO leaves that fragment disconnected and removes the unresolved
  scaffold bond. Such output must not be presented as a validated connected
  compound. Explicit references to nonexistent atoms are errors.
- Input molecules are values and are not mutated. Source-indexed CX state and
  Sgroups are cleared on products, as appropriate for a newly assembled graph.
  Root data fields are retained as source annotations, not recomputed descriptors.

## Scope and resource controls

This is R-group configuration enumeration, not general patent Markush expansion
or compound-to-family matching. Repeat/link-node expansion, polymers, multicenter
attachments, nested R-group definitions and uninterpreted CX extensions are
rejected explicitly; none is silently flattened into a claimed concrete member.
Other generic atom/bond queries may remain in a configuration. Clients must check
whether a product is fully specified before computing compound properties.
REST_H-aware substructure matching is a separate capability.

The default maximum is 100,000 delivered configurations, configurable per call.
The array API throws `configurationLimitExceeded` without returning a partial
array. Streaming callers may already have received members before an error and
must distinguish deliberate early stopping from complete enumeration. There is
a 256-site recursion guard and task-cancellation checks throughout traversal.
The result limit bounds output, not wall-clock search time; use cancellation for
time budgets. No Java runtime or new package dependency is needed in production.

## Verification and reproducibility

The package tests retain the seven upstream RGfiles and reference expectations:
valid fixtures 1, 2, 3, 4, 5, and 7 yield **4, 12, 8, 2, 17,820, and 20**
configurations; fixture 6 rejects missing definitions.

On 2026-11-01 the official Maven Central CDK 2.13 distribution was executed
independently with `Scripts/RGroupParity.java`. For every configuration, the
oracle records sorted element/charge/REST_H atom entries and sorted
endpoint-element/bond-order entries. The Swift test recomputes each complete
histogram and its SHA-256, comparing against the compact committed
`java-cdk-2.13.json` reference. All **17,866** configurations agree. These
signatures are not graph-isomorphism certificates; separate Standard InChI
checks verify aromatic substituents, both two-point attachment orders, scaffold
and fragment stereochemistry, and flattened round trips.

Recreate the oracle outside the checkout's build products:

```sh
mvn -B -f Scripts/rgroup-parity-pom.xml dependency:build-classpath \
  -Dmdep.outputFile=/tmp/cdk-rgroup-classpath.txt
java --class-path "$(cat /tmp/cdk-rgroup-classpath.txt)" Scripts/RGroupParity.java \
  Tests/CDKSwiftNativePortTests/MDL/RGroupEnumerationFixtures > /tmp/cdk-rgroup-histograms.json
```

For each fixture, serialize the sorted histogram as `signature<TAB>count<LF>`
and calculate SHA-256. The raw JSON evidence is retained at
`/private/tmp/cdk-rgroup-java/histograms.json` for this run, SHA-256
`9b6015fb14fcdfcaed28e8e012dd8dcde1e91aa9075b3b1b53a28003e68b5109`.
The versioned compact references, fixture bytes and oracle source permit reruns.

Release qualification: full macOS and Apple-container Ubuntu suites each ran
605 tests, with two existing skips and no failures. The ten dedicated enumeration
tests passed on both platforms. Strict swift-format lint passed for changed Swift
files. Evidence: `/private/tmp/cdk-rgroup-macos-full.log` and
`/private/tmp/cdk-rgroup-linux-final.log`.
