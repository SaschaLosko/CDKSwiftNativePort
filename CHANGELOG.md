# Changelog

All notable changes to this project are documented in this file.

The format is based on Keep a Changelog and this project follows Semantic Versioning.

## [Unreleased]

- No entries yet.

## [1.0.0] - 2026-02-19

### Added

- Unified import/export dispatch APIs (`CDKFileImporter`, `CDKFileExporter`) across MDL, SMILES, InChI, MOL2, PDB, XYZ, CML, RXN, RDF, SVG.
- CDK-derived layout/depiction pipeline, including scene generation and style controls.
- Native Swift InChI parse/generate path and InChIKey generation service.
- Molecular property service and descriptor set including XLogP and Rule-of-Five helpers.
- Expanded CDK parity-oriented test suite and metadata tests.
- Package boundary guard tests to prevent app-level coupling from entering package sources.
- Release documentation set (`Documentation/API.md`, `CONTRIBUTING.md`, `SECURITY.md`).

### Changed

- Removed package source coupling to app-specific identifiers and app-branded defaults.
- Replaced app-scoped content-type identifiers in package IO metadata with package-neutral/common identifiers.
- Cleaned renderer dependency surface by removing `AppKit` import from package rendering implementation.

### Notes

- This package contains CDK-derived functionality and remains licensed under LGPL-2.1-or-later.
