# 2D renderer qualification

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

The first milestone covers glyph visibility, arrow topology, and proportional
zoom at document and compact sizes.

The second milestone shares Core Text ink metrics across placement, clipping,
drawing, and vector outlines. Fitted stroke sizing keeps a readable floor for
compact structures, retains the sparse-molecule cap, and scales parallel-stroke
separation with bond width. The Java reference now covers measured ink height,
terminal bond clearance, and single, carbonyl, aromatic, and hashed stereo
strokes at 200 × 100, 560 × 418, and 960 × 720.

The complete package suite ran 746 tests on ARM64 macOS and 745 on ARM64
Ubuntu, with zero failures. The Core Text ink test runs on Apple platforms.

The third milestone matches the seven automatic arrow shapes measured from
Java CDK 2.13 SVG output. It checks polygon vertices, shaft and crossing
segments, stroke width, the absence of fill on the retrosynthetic and
undirected forms, and the fitted scene's zoom and rotation transforms.
Render hosts consume the same filled scene polygons alongside strokes.
