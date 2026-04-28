# Reaction Parity

`CDKSwiftNativePort` now tracks a curated reaction parity corpus derived from
upstream CDK `2.12` reaction model, reaction-IO, and reaction-manipulator
tests.

## Scope

The vendored corpus lives under:

- `Tests/CDKSwiftNativePortTests/Reaction/UpstreamReference/`

It currently covers:

- reaction CML hierarchy parsing and round-tripping
- reaction schemes and step lists
- reaction molecule-set references
- RXN V2000 parsing and round-tripping, including strict/relaxed counts-line
  agent behavior
- RXN V3000 parsing and round-tripping for vendored upstream resources
- RDF parsing and RDF round-tripping for multi-reaction sets
- reaction SMILES parsing and SMILES round-tripping

In addition to the curated vendored corpus, the direct Swift regression suite now
also exercises:

- RXN V3000 agent handling and zero-reactant handling
- RXN V3000 parsing and writing
- RXN writer round-trips for V2000 and V3000 reaction blocks
- multi-record reaction SMILES import into native reaction sets
- native reaction manipulator utilities, scheme-manipulator utilities, and
  inline reaction conversion helpers

The executable gate is:

- `Tests/CDKSwiftNativePortTests/Reaction/ReactionReferenceParityTests.swift`

## Current Status

- The committed known-gap inventory is `0` gaps across the current curated
  `22`-case reaction corpus.
- Strict reaction parity can therefore be required today for the maintained
  reaction model/IO/manipulator surface.

## Scope Boundary

The strict gate tracks the upstream-derived reaction surface that
`CDKSwiftNativePort` currently claims parity for:

- reaction hierarchy/model types
- reaction CML, RXN, RDF, and reaction-SMILES IO
- reaction manipulator, reaction-set manipulator, and reaction-scheme
  manipulator utilities

It does not currently claim parity for separate upstream subsystems such as:

- the `base/reaction` reaction-engine and reaction-type algorithms

`RInChI` is tracked separately in [`RINCHI_PARITY.md`](RINCHI_PARITY.md).

## Workflow

Refresh the committed inventory after reviewing intentional changes:

```bash
CDK_REACTION_WRITE_GAP_INVENTORY=1 swift test --package-path CDKSwiftNativePort --filter ReactionReferenceParityTests
```

Run the strict gate:

```bash
CDK_REACTION_REQUIRE_UPSTREAM_PARITY=1 swift test --package-path CDKSwiftNativePort --filter ReactionReferenceParityTests
```

Normal package tests compare the live gaps against the committed inventory so
inventory drift is visible even when strict mode is not enabled.
