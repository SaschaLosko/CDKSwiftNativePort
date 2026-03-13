# macOS-Specific Notes

`CDKSwiftNativePort` currently declares `macOS 14+` in its package manifest and
is designed to integrate cleanly into sandboxed native desktop hosts.

## 1. Sandboxed File Access Support

`CDKFileImporter.readMolecules(from:)` and `CDKFileImporter.readReaction(from:)`
call `startAccessingSecurityScopedResource()` for file URLs and stop access
afterward.

Implication:

- the package works cleanly with user-selected file URLs and security-scoped
  bookmarks in App Sandbox hosts
- host apps do not need to wrap every read path with their own temporary access
  handling when they already have a valid security-scoped URL

## 2. Native Rendering Surfaces for macOS Hosts

The package exposes rendering outputs that fit common macOS app stacks:

- `CDKStandardGenerator.draw(...)`
  - for SwiftUI `GraphicsContext` drawing
- `CDKMetalDepictionSceneBuilder`
  - for host-owned Metal molecule views
- `CDKMetalReactionDepictionSceneBuilder`
  - for host-owned Metal reaction views and hit-testing
- `CDKDepictionGenerator.toSVG(...)`
  - for export, previews, snapshots, and web-style embedding

Important: the Metal scene builders do not own a `MTKView` or a Metal pipeline.
They provide scene data only, which keeps the package reusable and keeps host UI
policy outside the chemistry core.

## 3. App Store-Friendly Boundary

The package intentionally excludes host extension wiring and app-level policies:

- no Quick Look extension implementation
- no Spotlight importer or indexer implementation
- no bundle identifiers or entitlements
- no window restoration or document-session logic

This separation is useful for App Store submission because host-review concerns
stay in the host app and its extensions, while chemistry logic remains packaged
as a standalone library.

## 4. No AppKit Dependency in the Chemistry Core

The package uses Swift / Foundation / CoreGraphics and selected Apple frameworks
that support the current rendering and text-measurement paths, but it does not
depend on `AppKit` or AtomLens code.

That keeps the chemistry core easier to test and reuse across multiple macOS
targets without pulling in host-app implementation.

## 5. Performance-Oriented Building Blocks

- depiction scene builders cache prepared depiction data for redraw-heavy flows
- scene output is renderer-neutral, which lets hosts optimize their own drawing
  backend
- Markush legend handling is computed once at scene-build time rather than being
  hard-coded into host drawing code

## 6. Recommended Host Architecture

- Use `CDKSwiftNativePort` for chemistry parsing, layout, identifiers,
  descriptors, and depiction scene generation.
- Keep all extension targets (Quick Look, Spotlight, metadata importers) in host
  app targets.
- Keep sandbox policy, bookmarks, persistence, and app-document coordination in
  host code.

## 7. What Is Not macOS-Specific

Although the package targets macOS today, most of the chemistry model and IO
logic is not AtomLens-specific and not tied to host-app state. The most
macOS-oriented parts are:

- security-scoped URL access helpers
- SwiftUI `GraphicsContext` renderer
- scene-generation APIs intended for native host rendering pipelines

If broader platform support is needed later, that should be treated as an
explicit compatibility project instead of an assumed guarantee.
