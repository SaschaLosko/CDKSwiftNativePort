# Publishing Guide

This guide describes a typical professional GitHub release flow for `CDKSwiftNativePort`.

## 1) Pre-release Quality Gate

From repository root:

```bash
swift test
```

Required before release:
- tests pass
- docs updated (`README.md`, `Documentation/API.md`, `CHANGELOG.md`)
- package boundary policy remains satisfied
- license and notice files present (`LICENSE`, `NOTICE.md`)

## 2) Prepare Release Commit

1. Ensure `CHANGELOG.md` has a finalized section for the target version.
2. Commit release prep changes:

```bash
git add .
git commit -m "Prepare v1.0.0 release"
```

## 3) Tag the Release

```bash
git tag -a 1.0.0 -m "CDKSwiftNativePort 1.0.0"
```

## 4) Configure Remote (first time)

Using GitHub CLI:

```bash
gh repo create <your-org-or-user>/CDKSwiftNativePort --public --source . --remote origin --push
```

or manually:

```bash
git remote add origin git@github.com:<your-org-or-user>/CDKSwiftNativePort.git
git push -u origin main
```

## 5) Push Tags

```bash
git push origin --tags
```

## 6) Create GitHub Release

Use GitHub UI or CLI:

```bash
gh release create 1.0.0 --title "CDKSwiftNativePort 1.0.0" --notes-file CHANGELOG.md
```

Recommended release body:
- summary of major additions
- breaking changes (if any)
- migration notes
- known limitations

## 7) Consumer Integration

Consumers add:

```text
https://github.com/<your-org-or-user>/CDKSwiftNativePort.git
```

and pin to:
- exact tag `1.0.0`, or
- semver range `from: "1.0.0"`

## 8) Post-release

- Create `Unreleased` section entries in `CHANGELOG.md` for next cycle.
- Track parity/bugfix tasks as issues.
