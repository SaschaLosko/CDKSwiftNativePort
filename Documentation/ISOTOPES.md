# Isotope data and mass calculations

`CDKIsotopes` provides immutable isotope values, exact isotopic mass (Da), natural
abundance (percent), the most abundant natural isotope and natural average mass.
Lookups accept atomic number; `mass(element:massNumber:average:)` also accepts an
element symbol. Invalid/unknown isotopes return nil from lookup. A missing or
invalid bundled dataset throws `CDKIsotopeError.dataUnavailable`, distinct from
absent physical data. No natural average or major isotope is invented for an
element whose dataset supplies no natural abundance.

The checked-in `Resources/Isotopes.json` contains 3,171 records derived from
[CDK 2.14.1 isotopes.dat](https://github.com/cdk/cdk/blob/cdk-2.14.1/base/core/src/main/resources/org/openscience/cdk/config/data/isotopes.dat).
The binary source SHA-256 is
`878d90ac8a9b47406b08c768629d30ad408d2380bf1b301158a4e516a870a562`.
The 92 average entries and 109 major/common overrides retain the published
decimal precision of the [tagged static tables](https://github.com/cdk/cdk/blob/cdk-2.14.1/base/core/src/main/java/org/openscience/cdk/config/Isotopes.java)
(Git blob `b10ab6bc74cf77da3988cb780eb5f3c7b533af8e`). This is a versioned
upstream-parity dataset, not a claim to the latest physical reference values.

CDK attributes these values to the Blue Obelisk Data Repository, version 10.
[BODR 10 is CC0](https://github.com/egonw/bodr/blob/BODR-10/bodr/COPYING);
the license is retained alongside the packaged data. The CDK adaptation credits
Christoph Steinbeck (2001–2007), Egon Willighagen (2013, 2016) and John Mayfield
(2026), under LGPL-2.1-or-later, as retained in the Swift source and package license.

For regeneration, decode the big-endian binary bucket header, total/natural
counts, then the natural and additional isotope groups. Each group stores count,
atomic number and starting mass-number offset; mass numbers increment within a
group. Scale unsigned mass integers by 1e9 below mass number 4, 1e8 below 40,
otherwise 1e7. Abundance integers scale by 1e7. Overlay the tagged major/common
decimal entries, sort by atomic number/mass number and retain averages as listed.
The JSON is human-readable; no runtime conversion or network access is needed.

The dataset is decoded and indexed once using Swift's thread-safe static
initialization. The checked-in `Scripts/benchmark_isotopes.swift` probes the
public API in a fresh process. On the qualification Apple Silicon Mac with
Swift 6.4.2, an optimized probe linked to the Debug library measured 19.8 ms for
first lookup including loading, then 16.7 ms for 200,000 mixed common/rare
lookups. These are one-run observations, not a release throughput guarantee.
The measured cost does not justify a second fast-path cache or more lazy loading.

`CDKExactMassDescriptor.calculateValidated` and
`CDKMolecularWeightDescriptor.calculateValidated` throw when a required mass is
unavailable. Existing nonthrowing `calculate` methods return NaN for that case;
callers must display an unavailable result or omit nonfinite JSON values.
Explicit isotope labels use exact mass in both descriptors; unlabeled atoms use
major-isotope or average mass respectively. Implicit H is included by default.
Rule-of-Five evaluation remains partial when mass is unavailable.

CTfile mass-difference reference numbers remain in their dedicated MDL table.
They describe a file encoding convention and must not be derived from these
physical masses. Molecular formula output remains an elemental composition;
isotope labels are retained on the source graph and in supporting file formats.

Validation includes all-record uniqueness/order/lookups; common/rare expected
values; carbon-13 and heavy-water isotopologues; implicit-H inclusion; unavailable
natural and explicit mass; partial Rule-of-Five results; and rare-element CML
classification. The full suite passes 709 tests on macOS and Linux.
