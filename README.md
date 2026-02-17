# CDKPortExperimental

Swift package containing the CDK-derived chemistry core used by ChemSketcher.

## Integration

Add as a Swift Package dependency and import `CDKPortExperimental`.

## Package Layout

- `Sources/CDKPortExperimental/Molecule.swift`
- `Sources/CDKPortExperimental/CDK/Layout/StructureDiagramGenerator.swift`
- `Sources/CDKPortExperimental/CDK/Rendering/StandardGenerator.swift`
- `Sources/CDKPortExperimental/CDK/InChI/InChIGeneratorFactory.swift`
- `Sources/CDKPortExperimental/CDK/InChI/InChIToStructure.swift`
- `Sources/CDKPortExperimental/CDK/Smiles/SmiFlavor.swift`
- `Sources/CDKPortExperimental/CDK/Smiles/SmilesParserFactory.swift`
- `Sources/CDKPortExperimental/CDK/Smiles/SmilesParser.swift`
- `Sources/CDKPortExperimental/CDK/Smiles/CxSmilesParser.swift`
- `Sources/CDKPortExperimental/CDK/Smiles/SmilesReaction.swift`
- `Sources/CDKPortExperimental/CDK/Smiles/SmilesReactionParser.swift`

## Scope

- CDK-style structure diagram generation and overlap refinement
- CDK-like depiction styling conventions
- CDK-style SMILES/CXSMILES/reaction parsing
- CDK-like InChI parsing layers (`/c`, `/h`, `/q`, `/i`, `/t`, `/m`)

## Tests

- `Tests/CDKPortExperimentalTests/InChI`
- `Tests/CDKPortExperimentalTests/Layout`
- `Tests/CDKPortExperimentalTests/Smiles`

Run with:

```bash
swift test
```

## License

This package is distributed under the GNU Lesser General Public License v2.1 (or later), see `LICENSE`.
