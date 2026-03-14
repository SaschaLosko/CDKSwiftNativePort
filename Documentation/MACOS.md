# macOS-Specific Notes

`CDKSwiftNativePort` currently declares `macOS 14+` and is designed to work
cleanly inside sandboxed native macOS hosts.

## 1. Sandboxed File Access

`CDKFileImporter` and related helpers route file reads through `CDKFileAccess`,
which provides:

- security-scoped resource access for file URLs
- coordinated reads via `NSFileCoordinator`
- shared text/data loading helpers for apps and extension targets

That is useful for:

- App Sandbox document-based apps
- Quick Look targets
- metadata importers
- helper tools that receive security-scoped URLs

## 2. Renderer-Neutral Outputs

The package does not own a window, view, or GPU pipeline. Instead it provides:

- `CDKDepictionGenerator`
  - SVG export for snapshots, previews, and web-style embedding
- `CDKMetalDepictionSceneBuilder`
  - molecule scene data for host Metal/Core Graphics renderers
- `CDKMetalReactionDepictionSceneBuilder`
  - reaction scene data and participant hit-testing

This is a good fit for macOS because the host can use whichever presentation
stack it wants without moving chemistry logic into app code.

## 3. App Store-Friendly Boundary

For App Store preparation, the key property is separation:

- the package contains chemistry logic
- the app and its extensions contain platform UI and bundle behavior

The package intentionally does not include:

- Quick Look extension entry points
- Spotlight importer entry points
- bundle identifiers
- entitlements
- document/session restoration logic

That keeps review-facing concerns where Apple expects them: in the app and its
extensions rather than in the chemistry library.

## 4. No App Framework Coupling in the Chemistry Core

The package uses `Foundation`, `CoreGraphics`, and related low-level Apple
frameworks required for the chemistry and rendering pipeline. It does not
depend on:

- `AtomLens`
- `AppKit`
- `UIKit`
- `SwiftUI`

That makes the package easier to test, publish, and reuse across multiple macOS
targets.

## 5. Good Extension Reuse Characteristics

Because the chemistry core stays package-owned and headless, the same package
can be shared by:

- the main app
- Quick Look preview and thumbnail targets
- Spotlight/metadata importers
- command-line utilities
- future helper or conversion tools

The extension target should call package APIs directly instead of copying
chemistry code into the extension bundle.

## 6. Recommended Host Architecture

- keep parsing, layout, depiction, export, identifiers, and descriptors inside
  `CDKSwiftNativePort`
- keep `MTKView`, Core Graphics, or view-model logic in the host
- keep extension lifecycle and metadata schema wiring in the host
- keep persistence, bookmarks, and sandbox policy in the host

## 7. Practical macOS Highlights

The most useful package characteristics for a macOS app today are:

- security-scoped file handling via `CDKFileAccess`
- a reusable chemistry core that can be shared by the app and extensions
- scene generation that lets the host control rendering and interaction policy
- no dependency on AtomLens or on higher-level app frameworks
