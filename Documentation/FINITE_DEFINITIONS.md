# Scoped finite definitions and chemical classes

Release 1.4.44 extends the reviewed finite-definition API with lexical scopes,
deterministic chemical classes and compilation evidence. These are original Swift
capabilities layered on the ported graph and R-group machinery. They are not
claims of upstream feature parity or automatic interpretation of source prose.

## Local definitions

`Variable.localDefinitions` contains an explicit alternative index and a list of
variables whose meanings are local to that alternative. Local labels shadow
enclosing labels; sibling alternatives do not share local definitions. A referenced
enclosing definition keeps its own declaration scope. Cycles are detected by
definition identity, so a local variable can safely reuse an outer label.

Occurrence, exclusions, dependencies and REST_H apply within the relevant scaffold.
Dependencies across separate scaffolds remain unsupported. Excluded alternatives
and their local definitions stay in the source sidecar. Embedded CX definitions
inside fragment strings remain unsupported by this compiler; use explicit local
definitions instead. This is not unrestricted nested CX evaluation.

## Finite classes

`Variable.chemicalClasses` supplements explicit `alternatives`. Each
`CDKFiniteChemicalClass` specifies `kind` (`alkyl` or `alkoxy`), `minimumCarbons`
and `maximumCarbons`, with an inclusive range contained in 1–6.

Both classes are unsubstituted, acyclic and saturated. All constitutional carbon
skeletons and distinct attachment positions are included. Alkoxy attaches through
oxygen; its carbon range counts only the alkyl part. Each site must have exactly
one single-bond attachment. Class fragments have unconstrained stereochemistry;
the expansion is not a list of individually specified stereoisomers. Isotope
substitution, rings, unsaturation, optional substitution and broader ranges are
not included. Unsupported ranges are rejected, never clipped to the supported range.

The generator enumerates unordered rooted carbon trees with at most three
children per carbon. The independent reference counts for C1–C6 are 1, 1, 2, 4,
8 and 17, respectively: [OEIS A000598](https://oeis.org/A000598). Tests additionally
compare all four butyl attachment identities with independently written SMILES,
check unique full Standard InChIs for each size and exercise both enantiomers of
an unconstrained stereogenic attachment.

Exclusion indices continue to refer only to explicit source alternatives, not to
generated class members. Class-specific exclusions or substitution conditions must
remain unresolved until represented by a supported explicit interpretation.

## Evidence

`compileWithEvidence()` returns the molecule and `CDKFiniteDefinitionEvidence`.
The evidence maps each compiled alternative index to its source alternative,
verbatim wording, local scope, nested choices and generated class fragment.
`evidence.choices(for: configuration.assignments)` explains selections, including
omitted sites. Use assignments from that exact compilation. Scope paths use
zero-based alternative indices to identify lexical scopes; they are not JSON pointers.
Atom IDs in a nested choice refer to its local scaffold.

The original sidecar is authoritative source evidence. A generated combination
does not assert that a complete individual compound was explicitly disclosed.
Export the sidecar alongside the selected choices and comparison results.
The flattened molecule alone does not retain the source wording or exclusions.

## Serialization and bounds

Legacy sidecars without a schema version remain readable. Extended sidecars use
`schemaVersion: 2` and `reviewedDefinition` for review state. Their legacy
`reviewed` field is false, keeping older readers' automatic compilation gate
closed when those readers ignore new keys. Older editors cannot preserve the new
semantics and must not be used to edit these drafts. Unknown schema versions are
rejected. Drafts remain serializable before review.

Compilation requires review, nonblank source identity and wording, and no unresolved
conditions. It permits at most 128 declared variables, 1,000 explicit alternatives
and class declarations combined, 16 nested scope/expansion levels and 32 KiB per
fragment. Default aggregate expansion limits remain 1,000 retained alternatives
and 2,000 retained atoms, with cancellation checks. Replacement of a fragment's
first attachment atom remains unsupported. Compilation failure never establishes
negative membership.

## Qualification and next work

Full macOS and Apple-container Ubuntu suites each pass 677 tests with zero failures.
Fourteen focused tests include the existing 103-annotation review gate, independent
identities, class boundaries, lexical shadowing, inherited-scope behavior,
dependencies, omission evidence and schema compatibility. The ten executable M2S
interpretations remain ten: a supported term does not resolve the rest of a source.

The M2S annotations motivate further scoped work: record 0 includes haloalkyl and
other conditions; 33 requests C1–C10; 36 refers to another formula; 81 also has
positional-variation ambiguity; 92 has a substitution clause absent from its
shortened annotation. These records must not become executable merely because
some alkyl terms are now supported.

Next extensions are bounded substitution rules, precisely delimited aromatic and
heterocyclic classes, broader source qualification and measured direct matching.
Complete upstream feature coverage remains a long-term objective, tracked
separately from original extensions. Query writing and stereo gaps should be
qualified against specific upstream parsing, matching and writing behavior.
