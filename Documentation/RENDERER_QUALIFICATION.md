# 2D reaction renderer qualification

The renderer resolves shared stroke and symbol sizes from the median fitted bond
length across reaction participants. Interactive zoom scales the reaction
coordinates, strokes, and symbol sizes together. Very small canvases use a
compact initial size. Projected atom coordinates retain a common scale across
reactants, agents, and products. Separators use the nearest available position
when a participant label occupies their usual location.

Explicit isotope masses use superscripts in the atom label. An isotope-labelled
carbon remains visible when the style hides ordinary carbon symbols. Automatic
reaction arrows distinguish forward, backward, bidirectional, no-go,
retrosynthetic, resonance, and undirected reactions.

The [Java CDK 2.13 reference data](../Tests/CDKSwiftNativePortTests/Rendering/JavaReference/README.md)
records measurements and arrow semantics from three Java zoom levels and seven
reaction directions, plus isotope visibility in one labelled molecule. Run the
offline reference suite with:

```sh
swift test --filter JavaDepictionReferenceTests
```

The Swift checks cover glyph visibility, arrow topology, and proportional zoom
at document and compact sizes. Typography, arrowhead fill, and raster edge
coverage have separate qualification requirements.

The complete package suite ran 744 tests with zero failures on ARM64 macOS and
ARM64 Ubuntu for this milestone.
