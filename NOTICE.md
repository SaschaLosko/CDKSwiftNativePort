# NOTICE

`CDKSwiftNativePort` includes code derived from and inspired by the
Chemistry Development Kit (CDK), developed by the CDK project contributors.

Upstream project:
- https://github.com/cdk/cdk

Reference target used for parity work in this port:
- CDK 2.12

`CDKSwiftNativePort` also vendors the official IUPAC InChI source distribution
as the `IUPACInChI` SwiftPM C target.

Upstream project:
- https://github.com/IUPAC-InChI/InChI

Vendored source:
- IUPAC InChI `v1.07.5`
- Git commit `11a87982bb518f57ac013f0b258c283655e1ea1d`
- License copy: `Sources/IUPACInChI/LICENSE-IUPAC-InChI.txt`
- External contributor notice: `Sources/IUPACInChI/External-contributors-IUPAC-InChI.txt`

Reference citations:

- Willighagen et al. (2017), The Chemistry Development Kit (CDK) v2.0:
  atom typing, depiction, molecular formulas, and substructure searching.
  DOI: [10.1186/s13321-017-0220-4](https://doi.org/10.1186/s13321-017-0220-4)
- May and Steinbeck (2014), Efficient ring perception for the Chemistry
  Development Kit. DOI: [10.1186/1758-2946-6-3](https://doi.org/10.1186/1758-2946-6-3)
- Steinbeck et al. (2006), Recent Developments of the Chemistry Development
  Kit (CDK) - An Open-Source Java Library for Chemo- and Bioinformatics.
  DOI: [10.2174/138161206777585274](https://doi.org/10.2174/138161206777585274)
- Steinbeck et al. (2003), The Chemistry Development Kit (CDK): An
  Open-Source Java Library for Chemo- and Bioinformatics.
  DOI: [10.1021/ci025584y](https://doi.org/10.1021/ci025584y)

Licensing:
- This repository is distributed under LGPL-2.1-or-later (`LICENSE`).
- The vendored IUPAC InChI source files are distributed under the MIT license
  terms preserved in `Sources/IUPACInChI/LICENSE-IUPAC-InChI.txt`.
- When redistributing binaries or modified source, ensure LGPL obligations
  (including notice and source availability requirements) are satisfied.
