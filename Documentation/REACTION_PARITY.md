# Reaction Parity

`CDKSwiftNativePort` now tracks a curated reaction parity corpus derived from
upstream CDK `2.12` reaction tests.

## Scope

The vendored corpus lives under:

- `Tests/CDKSwiftNativePortTests/Reaction/UpstreamReference/`

It currently covers:

- reaction CML hierarchy parsing and round-tripping
- reaction schemes and step lists
- reaction molecule-set references
- RXN parsing and RXN round-tripping for single reactions
- RDF parsing and RDF round-tripping for multi-reaction sets
- reaction SMILES parsing and SMILES round-tripping

The executable gate is:

- `Tests/CDKSwiftNativePortTests/Reaction/ReactionReferenceParityTests.swift`

## Current Status

- The committed known-gap inventory is `0` gaps across the current curated
  reaction corpus.
- Strict reaction parity can therefore be required today for that vendored
  surface.

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
