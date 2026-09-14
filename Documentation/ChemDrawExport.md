# Native ChemDraw export

`CDKFileExporter.writeData(molecule:as:)` and
`writeData(reaction:as:)` support `.cdx` and `.cdxml`. The hierarchy overload
writes one page per flattened reaction. The molecule-array overload writes one
page per molecule. Use the Data or URL API for binary CDX; the String API rejects
CDX rather than corrupting it through a text encoding. CDXML also supports the
existing String API. `CDKChemDrawWriter` exposes the underlying array writers.

The implementation is Swift and Foundation only. It follows the public
[Revvity CDX/CDXML specification](https://github.com/Glysade/chemdraw/tree/v1.0.14/docs).
No Indigo, RDKit, ChemDraw SDK, or C++ runtime dependency is added.

Supported: editable 2D molecular graphs, elements, isotopes, formal charge,
explicit hydrogen counts, aromatic bonds, wedge/hash bonds, explicit tetrahedral
ligand order, geometric cis/trans bonds, Unicode titles, multiple pages, reaction
participants/agents, common reaction arrows, and paired atom-map references.
CDX uses little-endian fields, a 22-byte prefix followed by the implicit six-byte
document record, extended property lengths, and UTF-8 font runs.

Exports contain chemical structures and titles, not an archival copy of arbitrary
SDF data fields or source-document annotations. Atom number labels are exported;
reaction map pairs are written into both automatic-neutral and manual mapping
properties. Reaction hierarchy nesting is flattened into pages.

The writer rejects unsupported chemistry explicitly: 3D coordinates, query atoms
and bonds, pseudo/R-group definitions, polymer/Sgroups, radicals, atom aliases,
valence overrides, enhanced/racemic stereo, coordinate bonds, reaction-center
annotations, non-unit stoichiometry, no-go arrows, and reaction-level CXSMILES annotations. Use MOL V3000, RGfile, or
CML as appropriate for those cases. Imported 2D wedge drawings are preserved;
explicit ligand parity is converted to outgoing wedges using the drawn geometry.
Degenerate geometry that cannot represent stereo fails rather than discarding it.
This release does not add CDX/CDXML readers.

## Validation

```sh
CDK_CHEMDRAW_EVIDENCE_DIR=/tmp/chemdraw-fixtures swift test --filter ChemDrawWriterTests
python3 Tools/validate_chemdraw_interop.py \
  --indigo-library /path/to/libindigo.dylib \
  --fixtures /tmp/chemdraw-fixtures
```

On 2026-09-14, Indigo 1.42.0 (`f856b9a9faf1edc0ff94be7aabde15faa450f187`)
independently read all 36 outputs from 17 molecule cases and one reaction case.
Canonical identities match, including both tetrahedral enantiomers, ring stereo,
multiple centers, unspecified double bonds, E/Z pairs, conjugated systems,
aromatic heterocycles, salts, and isotopes. Reaction participants and agents match.
Indigo does not read reaction mapping properties, so those references are checked
separately by the bounds-checked test-only CDX decoder. Tests also cover invalid
graphs, unsupported data, extended lengths, XML escaping, and atomic binary URL
writes. Full package suites pass on macOS and Linux: 596 tests, two platform skips,
zero failures on each platform. Actual desktop ChemDraw visual/open verification
is not part of this evidence.
