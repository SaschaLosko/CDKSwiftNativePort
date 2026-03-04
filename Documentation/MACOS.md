# macOS-Specific Notes

`CDKSwiftNativePort` is packaged for `macOS 14+` and is designed to integrate cleanly into sandboxed desktop hosts.

## Highlights

## 1. Sandboxed File Access Support

`CDKFileImporter.readMolecules(from:)` and `CDKFileImporter.readReaction(from:)` call `startAccessingSecurityScopedResource()` for file URLs and stop access afterward.

Implication:

- works cleanly with security-scoped bookmarks and user-selected file URLs in App Sandbox hosts.

## 2. Native Rendering Paths for macOS Hosts

The package provides multiple rendering targets that map well to macOS app stacks:

- `CDKStandardGenerator.draw(...)` for SwiftUI `GraphicsContext` usage.
- `CDKMetalDepictionSceneBuilder` / `CDKMetalReactionDepictionSceneBuilder` for custom Metal rendering pipelines.
- `CDKDepictionGenerator.toSVG(...)` for vector export and snapshots.

## 3. Host-Integration Boundary (App Store Friendly)

The package intentionally excludes host extension wiring and app-level policies:

- no Quick Look extension implementation
- no Spotlight importer/indexer implementation
- no app bundle IDs/entitlements/window-state logic

This separation keeps the chemistry core reusable and makes App Store host-review concerns local to the app project.

## 4. Recommended Host Architecture on macOS

- Use `CDKSwiftNativePort` for chemistry parsing, layout, identifiers, properties, and depiction scene generation.
- Keep all extension targets (Quick Look, Spotlight, metadata importers) in host app targets.
- Keep sandbox/entitlement decisions in host code, not inside this package.

## 5. Performance-Oriented Building Blocks

- Scene builder paths precompute and cache depiction data for repeated redraw scenarios.
- Renderer-agnostic scene data lets hosts optimize their own drawing backend (Metal, SwiftUI, offscreen export).

## 6. Platform Scope

Current package manifest declares:

- platform: `macOS(.v14)`

If broader platform support is required later, expansion should be treated as an explicit compatibility project.
