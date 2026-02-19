# Contributing

Thanks for contributing to `CDKSwiftNativePort`.

## Development Setup

1. Install Xcode with Swift 5.9+ toolchain.
2. Clone the repository.
3. Run:

```bash
swift test
```

## Contribution Rules

1. Keep package boundaries clean.
2. Do not add app-level dependencies into package source code.
3. Keep CDK-derived logic in CDK-like structure under `Sources/CDKSwiftNativePort/CDK/...`.
4. Add or update tests for every functional change.
5. Keep public API changes intentional and documented in `CHANGELOG.md` and `Documentation/API.md`.

## Package Boundary Policy

Package source code must not depend on:
- Spotlight APIs (`CoreSpotlight`)
- Quick Look APIs (`QuickLook`, `QuickLookThumbnailing`)
- app-specific bundle identifiers, app names, or app settings

If integration with an app is needed, that integration belongs in the app repository, not in this package.

## Testing

Run the complete suite before opening a PR:

```bash
swift test
```

The suite includes:
- parser/writer round-trip tests
- CDK parity metadata tests
- depiction/layout behavior tests
- package boundary guard tests

## Pull Request Checklist

- [ ] Tests pass locally (`swift test`)
- [ ] New behavior has test coverage
- [ ] No app-level coupling introduced
- [ ] Documentation updated (`README.md`, `Documentation/API.md`, `CHANGELOG.md`)
- [ ] License/attribution implications checked for CDK-derived code
