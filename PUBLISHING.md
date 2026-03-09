# Publishing Guide

This guide describes a typical professional GitHub release flow for `CDKSwiftNativePort`.

## 1) Pre-release Quality Gate

From repository root:

```bash
swift test
```

Required before release:
- tests pass
- docs updated (`README.md`, `Documentation/API.md`, `Documentation/INTEGRATION.md`, `Documentation/CDK_COMPARISON.md`, `Documentation/MACOS.md`, `CHANGELOG.md`)
- package boundary policy remains satisfied
- license and notice files present (`LICENSE`, `NOTICE.md`)

## 2) Prepare Release Commit

1. Ensure `CHANGELOG.md` has a finalized section for the target version.
2. Set a shell variable for the release version and reuse it in the commands below.
3. Commit release prep changes:

```bash
VERSION=1.1.1
git add .
git commit -m "Prepare v${VERSION} release"
```

## 3) Tag the Release

```bash
git tag -a "${VERSION}" -m "CDKSwiftNativePort ${VERSION}"
```

## 4) Configure Remote (first time)

Using GitHub CLI:

```bash
gh repo create SaschaLosko/CDKSwiftNativePort --public --source . --remote origin --push
```

or manually:

```bash
git remote add origin git@github.com:SaschaLosko/CDKSwiftNativePort.git
git push -u origin main
```

## 5) Push Tags

```bash
git push origin --tags
```

## 6) Create GitHub Release

Use GitHub UI or CLI:

```bash
gh release create "${VERSION}" --title "CDKSwiftNativePort ${VERSION}" --notes-file CHANGELOG.md
```

Recommended release body:
- summary of major additions
- breaking changes (if any)
- migration notes
- known limitations

## 7) Consumer Integration

Consumers add:

```text
https://github.com/SaschaLosko/CDKSwiftNativePort.git
```

and pin to:
- exact tag `1.1.1`, or
- semver range `from: "1.1.1"`

## 8) Post-release

- Create `Unreleased` section entries in `CHANGELOG.md` for next cycle.
- Track parity/bugfix tasks as issues.
