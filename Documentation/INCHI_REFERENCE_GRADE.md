# Reference-Grade InChI Plan

This document tracks the work required to make the `CDKSwiftNativePort` InChI
surface reference-grade while keeping the runtime Swift native.

## Target

The package must eventually match the official InChI reference implementation
for:

- generated InChI strings
- generated InChIKeys
- InChI-to-structure reconstruction behavior
- round-trip stability on the supported structure surface

The runtime target remains native Swift. The official InChI implementation is
used as a test oracle, not as a runtime dependency.

## Current State

- `CDKInChINativeGenerator` is still the package-owned Swift implementation.
- `CDKInChIToStructure` is still a partial Swift parser.
- InChIKey derivation now follows the official major/minor block split and
  base-26 encoding in native Swift instead of the earlier pseudo-key path.
- The generator now handles disconnected official-style multicomponent output
  for repeated simple fragments and metal-disconnected salts on the supported
  subset, instead of collapsing those cases back into a single aggregate
  formula.
- Parsed InChI molecules now carry a signature-guarded source cache so
  unchanged parser output re-exports the exact original InChI while edited
  structures still fall back to native regeneration.
- The parser now accepts repeated hydrogen multiplicity tokens such as
  `2*1H3` and `72*1H` without warning.
- The package now vendors an official-reference parity corpus derived from the
  upstream InChI CI regression data.
- The current vendored snapshot is pinned to official InChI commit
  `94c871858c6f15ec686c0417265f4644366ef989`.
- The current curated corpus contains `48` cases:
  - `6` elemental official edge cases
  - `6` isotopic official edge cases
  - `12` neutral organic structure cases
  - `6` charged cases
  - `6` double-bond stereo cases
  - `6` tetrahedral stereo cases
  - `6` multicomponent cases
- The parity harness can run in two modes:
  - default mode compares the current implementation against a committed
    known-gap inventory
  - strict mode (`CDK_INCHI_REQUIRE_REFERENCE_GRADE=1`) requires zero known
    gaps
- The current measured known-gap inventory is `51` gaps across that corpus:
  - `25` exact InChI mismatches
  - `25` InChIKey mismatches
  - `1` parser warning case
- Measured completion against the original `128`-gap baseline is `60.2%`.

## Milestones

1. Official reference harness
   - Vendor a curated official InChI CI corpus with exact upstream InChI and
     InChIKey outputs.
   - Add executable Swift parity tests and exact gap reporting.
   - Keep a committed known-gap inventory so progress is measurable without
     pretending the current implementation is reference-grade.

2. Generator parity
   - Replace heuristic layer assembly with a Swift port of the official
     canonicalization and emission pipeline.
   - Close exact-string mismatches before treating InChIKey parity as done.

3. InChIKey parity
   - Keep the official-equivalent key derivation pipeline in Swift aligned
     with the remaining generator exact-string work.

4. Parser parity
   - Replace partial layer parsing and ignored-token behavior with a full Swift
     port of the official structure-reconstruction path.

5. Release gate
   - Keep strict reference parity as the release criterion for claiming
     reference-grade InChI support.

## Refresh Workflow

Regenerate the vendored official corpus from a local checkout of the official
InChI repository:

```bash
python3 Tools/import_official_inchi_reference.py --source /path/to/InChI
```

Refresh the current known-gap inventory after reviewing intentional changes:

```bash
CDK_INCHI_WRITE_GAP_INVENTORY=1 swift test --package-path CDKSwiftNativePort --filter OfficialReferenceParityTests
```

Run the strict gate:

```bash
CDK_INCHI_REQUIRE_REFERENCE_GRADE=1 swift test --package-path CDKSwiftNativePort --filter OfficialReferenceParityTests
```
