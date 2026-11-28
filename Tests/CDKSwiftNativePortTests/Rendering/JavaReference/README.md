# Java CDK depiction reference

Reference data for [Java CDK 2.14.1](https://github.com/cdk/cdk/tree/cdk-2.14.1)
(`4253aa225f474a2a5b26ad9ed3ecf73262355693`). It was first generated with Java CDK 2.13
(`3bc6efe84119e5da0d00804bd4cf5f3bdd9c3d6d`) using OpenJDK 17 on macOS. On 2026-11-28 the
official `cdk-2.14.1.jar` release asset (SHA-256
`9b7b9b6b64b7cb0aeb9475ac3584a7ce4226fd922a9ebe8ba3da7b7653af2985`) was run on the
same workflows below, next to `cdk-2.13.jar`. Every SVG, the regenerated MOL and
RXN files (apart from their timestamp lines) and every bond display form were
identical, so the recorded measurements stand unchanged for 2.14.1.
The inputs are small purpose-built fixtures. The reference records measurements
and display semantics from Java output.

- Create `participants.rxn` from `CCO.OCC1=CC=CC=C1>>CCOCC1=CC=CC=C1` using
  `SmilesParser`, `StructureDiagramGenerator.generateCoordinates`, and `MDLRXNWriter`.
  Read the file with `MDLRXNReader` and retain the implicit hydrogen counts from
  the parsed reaction for the geometry measurements.
- Depict the reaction at zoom 1, 2 and 3 using `DepictionGenerator.withZoom` and
  `toSvgStr(Depiction.UNITS_PX)`. Record the length and inherited stroke width
  of `mol1bnd1`, the unlabelled C–C single bond in the first reactant.
- Set each `IReaction.Direction` and depict at default zoom. Between the
  reactant and product SVG groups, record centerline endpoints, stroke widths,
  and the vertices of each filled path subpath. Record the open retrosynthetic
  path as two connected strokes. Infer the nominal drawing frame from the
  tips and shafts. Use the head's horizontal extent as the arrow height,
  or three times the lane separation for the headless form.
- Create `carbon-13.mol` from `[13CH3][NH2+]CC(=O)[O-]` using the same molecule
  parsing and coordinate generation, then `MDLV2000Writer`. The default Java SVG
  shows atom 1 with a superscript 13, including when ordinary carbon labels
  are implicit.

The fitted legibility samples use the same coordinates and Java
`withSize(width, height).withMargin(24).withFillToFit()` at 200 × 100, 560 ×
418, and 960 × 720 px. The ethanol (`CCO`) and stereo (`C[C@H](O)Cl`)
fixtures use the same MOL writer and reader workflow as the charged isotope.
The reference compares inherited SVG stroke widths, perpendicular distance
between parallel stroke centerlines, outline ink heights, and the minimum
distance from the terminal C–O bond endpoint to the OH outline bounds. The
widest hashed stereo stroke is measured perpendicular to the wedge axis.

The Swift suite checks proportional zoom, arrow meaning, and isotope visibility,
plus fitted single-bond weight, label ink height and clearance, parallel
carbonyl and aromatic strokes, and hashed stereo spread. Fitted comparisons
normalize by stroke width or glyph height because the application caps the
visual prominence of sparse molecules and uses compact-size stroke floors.
Glyph metrics run on Core Text platforms. Linux runs the geometry checks
using the package's portable text size estimates.

`stereochemistry.json` adds the two enantiomers of
`[13CH3][C@H](F)C(=O)[O-]`, the difluoroethene E/Z pair, a crowded charged
isotope fixture, and two reactions pairing those isomers. Generate coordinates
with Java `StructureDiagramGenerator` directly from each SMILES and depict with
`withSize`, `withMargin(12)`, `withFillToFit`, and `withAtomNumbers` at 224 × 176
and 560 × 418 px. For each SVG molecule group, count the atom outlines, bond
lines, and filled wedge paths. Exclude reaction arrows and atom-number
annotations from those counts. Read the solid or hashed bond display from
`IBond.getDisplay()`. For each configured double bond, compare the signed side
of the two `IDoubleBondStereochemistry` carrier atoms against the bond axis.
The recorded counts and display forms agree at both canvas sizes.

The Swift tests check stereobond display form, E/Z geometry and line topology,
reaction participant labels, and crowded label clearance with a 32-point
viewport margin. Layout and individual hatch spacing can differ from Java.
