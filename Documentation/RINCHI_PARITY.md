# RInChI Parity

`CDKSwiftNativePort` includes a native Swift `RInChI` surface for reaction
identifiers. The implementation follows the official IUPAC RInChI reference
design from `https://github.com/IUPAC-InChI/RInChI` and keeps the InChI part
delegated to the vendored official IUPAC libinchi bridge.

## Scope

The current Swift implementation covers:

- `CDKRInChIGenerator`
- `CDKRInChIToReaction`
- `CDKRInChIDecomposition`
- Long/Short/Web `RInChIKey` derivation
- vendored upstream reaction-reference fixtures for RXN-based `RInChI` output

The key split is intentional:

- Long RInChIKey component InChIKeys are generated through libinchi.
- Short and Web RInChIKeys use the RInChI reference SHA-256/base-26 hashing
  algorithm implemented in Swift, because that hashing is part of RInChI
  itself rather than a substitute InChI implementation.

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
- The official InChI strict gate also passes with RInChI using libinchi for
  component InChIKey generation.

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
