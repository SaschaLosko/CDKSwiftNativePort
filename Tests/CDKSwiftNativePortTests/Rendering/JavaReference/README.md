# Java CDK depiction reference

Reference data generated using [Java CDK 2.13](https://github.com/cdk/cdk/tree/cdk-2.13)
at `3bc6efe84119e5da0d00804bd4cf5f3bdd9c3d6d` with OpenJDK 17 on macOS.
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
  reactant and product SVG groups, record the horizontal shafts, head locations
  and diagonal strokes crossing the shaft away from its ends.
- Create `carbon-13.mol` from `[13CH3][NH2+]CC(=O)[O-]` using the same molecule
  parsing and coordinate generation, then `MDLV2000Writer`. The default Java SVG
  shows atom 1 with a superscript 13, including when ordinary carbon labels
  are implicit.

The Swift suite tests proportional zoom, arrow meaning and isotope visibility.
Font metrics, absolute coordinates and arrowhead fill can vary by renderer.
