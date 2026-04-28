# RInChI Parity

`CDKSwiftNativePort` now includes a native Swift `RInChI` surface and an
executable vendored reference gate derived from upstream CDK `2.12`
`storage/rinchi` resources.

## Scope

The current Swift implementation covers:

- `CDKRInChIGenerator`
- `CDKRInChIToReaction`
- `CDKRInChIDecomposition`
- native Long/Short/Web `RInChIKey` derivation
- vendored upstream reaction-reference fixtures for RXN-based `RInChI` output

The vendored resource corpus lives under:

- `Tests/CDKSwiftNativePortTests/RInChI/UpstreamReference/`

The executable parity gate is:

- `Tests/CDKSwiftNativePortTests/RInChI/RInChIReferenceParityTests.swift`

## Current Status

- The vendored upstream-derived `RInChI` corpus is currently `0` gaps across
  `21` curated reaction cases.
- The direct Swift regression suite also covers nil-input handling,
  decomposition, `RInChI`-to-reaction reconstruction, equilibrium forcing, and
  `AuxInfo` exposure from the native InChI layer.

## Scope Boundary

The current gate tracks the maintained `RInChI` generation and reconstruction
surface that is actually exposed by the Swift package.

It does not claim byte-for-byte parity for every upstream `RInChI` internal
implementation detail outside the vendored executable reference cases.

## Workflow

Run the strict vendored gate:

```bash
swift test --package-path CDKSwiftNativePort --filter RInChIReferenceParityTests
```
