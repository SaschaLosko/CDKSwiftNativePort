# Architecture Overview

`CDKSwiftNativePort` keeps CDK-derived chemistry functionality inside a
standalone Swift package with a narrow, host-consumable API surface.

## Goals

Complete upstream CDK feature coverage is the long-term objective. Current parity
claims remain limited to verified surfaces. Original finite-definition and family
reasoning extensions are documented separately from upstream-derived algorithms;
they retain the same standalone, portable package boundary.

- own chemistry implementations independently of any host application
- expose a Swift-first chemistry core for parsing, layout, depiction, and IO
- make the package reusable from apps, extensions, command-line tools, and tests
- keep host-specific UI and bundle policy outside the chemistry core

## Package Responsibilities

The package owns:

- molecule and reaction model types
- format readers and writers
- SMILES, CXSMILES, InChI, MDL, CML, RXN, and RDF parsing/generation
- 2D coordinate generation
- depiction preprocessing
- SVG export
- renderer-neutral scene generation for host drawing backends
- identifier and descriptor services

## Host Responsibilities

The host app or extension owns:

- window, document, and session lifecycle
- SwiftUI/AppKit/UIKit composition
- `MTKView`, Core Graphics, or other concrete drawing loops
- Quick Look and Spotlight extension entry points
- entitlements, bundle identifiers, signing, and App Store metadata
- persistence, bookmarks, and higher-level app policy

## Top-Level Source Layout

- `Sources/CDKSwiftNativePort/{Molecule,Atom,Bond,ChemistryTypes,ChemError}.swift`
  - molecular graph helpers, individual value types, stereochemistry/Sgroup types and errors
- `Sources/CDKSwiftNativePort/CDK/Common`
  - identifier and molecular property services
- `Sources/CDKSwiftNativePort/CDK/Smiles`
  - SMILES, CXSMILES, reaction SMILES parsing and generation
  - CX state is a standalone value model. The internal layer decoder reads text;
    metadata application and R-group graph assembly are separate responsibilities.
    `CDKCxSmilesParser` retains the public entry points and malformed-layer title
    fallback. Decoding never publishes a partially accepted metadata state.
- `Sources/CDKSwiftNativePort/CDK/InChI`
  - InChI parse/generate facades
- `Sources/CDKSwiftNativePort/CDK/Layout`
  - `Depiction2DGenerator` entry point and CDK-derived structure diagram generation
  - graph/ring perception and separate coordination, Markush, ring-system, chain,
    ring-placement, refinement and geometry stages; these helpers remain module-internal
- `Sources/CDKSwiftNativePort/CDK/Rendering`
  - depiction preprocessing, SVG generation, scene generation, highlight/query/Sgroup rendering
  - Metal scene values, preparation cache, fitting metrics and bond primitives
    have separate owners. Bond primitives use immutable per-render metrics;
    viewport transforms stay outside the shared preparation cache. The existing
    lock protects both cache entries and LRU order, while preparation happens
    outside that lock. The public scene builder and rendered chemistry are unchanged.
- `Sources/CDKSwiftNativePort/CDK/IO`
  - reader/writer implementations and unified import/export dispatch
- `Sources/CDKSwiftNativePort/CDK/QSAR`
  - descriptor implementations and property helpers

## Public API Shape

The package deliberately favors a small set of host-facing entry points:

- `CDKFileImporter`
- `CDKFileExporter`
- `Depiction2DGenerator`
- `CDKDepictionGenerator`
- `CDKMetalDepictionSceneBuilder`
- `CDKMetalReactionDepictionSceneBuilder`
- `CDKMoleculeIdentifierService`
- `CDKMoleculePropertyService`

Lower-level readers, writers, and helpers stay public where they are useful, but
host code usually does not need to bind to every internal port layer.

## Rendering Model

The rendering stack is split into three layers:

1. model preparation
   - layout generation and depiction preprocessing normalize coordinates,
     labels, Sgroups, highlights, and query visuals
2. package-owned outputs
   - `CDKDepictionGenerator` produces SVG
   - Metal scene builders produce geometry and labels for host renderers
3. host-owned drawing
   - a host turns the generated scene into on-screen pixels

This keeps the chemistry logic package-owned while letting the app keep control
over view hierarchy, gestures, animation, and GPU pipeline details.

## Boundary Contract

### Forbidden inside package sources

- host application imports or symbols
- Quick Look or Spotlight implementation code
- bundle IDs, entitlements, or app restoration logic
- `AppKit`, `UIKit`, or `SwiftUI` dependencies

### Allowed and expected in host code

- importing `CDKSwiftNativePort`
- turning scene data into actual drawing
- extension targets that call package APIs
- app-specific wrappers such as document/view-model layers

## Standalone Package Contract

All maintained sources, documentation, tooling and tests are package-scoped.
They must not name or depend on a particular consuming application or customer.
There is no Java source, JVM runtime or Maven build requirement. The official
native C InChI implementation is vendored in its own package target.

Consumers link the package product rather than compiling its sources directly.
App configuration and integration tests belong in the consuming repository.
Package tests must run from an isolated checkout without sibling projects.

## Boundary Enforcement

- Source tests reject application UI, Quick Look and Spotlight imports.
- Manifest tests reject external local-package paths.
- File checks reject Java source and compiled Java artifacts in maintained files.
- Package validation uses only this checkout and its bundled reference fixtures.

## Testing Strategy

The package test layout is organized by capability rather than by app feature:

- `Smiles`
- `MDL`
- `CML`
- `Rendering`
- `Layout`
- `Common`
- `QSAR`
- boundary tests

Coverage focuses on:

- reader/writer round trips
- parity-oriented cases ported from upstream CDK 2.13
- depiction behavior in both SVG and scene-generation paths
- boundary enforcement between the package and the host app

## R-group enumeration (1.4.35)

Native `CDKRGroupQuery.getAllConfigurations(maximumConfigurations:)` and
`enumerateConfigurations(maximumConfigurations:_:)` implement CDK 2.13 R-group
configuration semantics without a Java dependency. See the [API, limits and
Java parity evidence](RGROUP_ENUMERATION.md). Repeat-unit expansion and family-aware search
remain separate capabilities.
