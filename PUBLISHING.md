# Publishing Guide

This guide describes a practical GitHub release flow for
`CDKSwiftNativePort`.

## 1. Pre-release Quality Gate

From the package root:

```bash
swift test
```

Required before publishing:

- the package test suite passes on macOS and on at least one Linux environment
- `README.md` and the documentation set match the current public API
- Linux installation and usage examples have been reviewed against the current
  release candidate
- `CHANGELOG.md` is updated
- package boundary checks remain green
- `LICENSE` and `NOTICE.md` are present

Recommended Linux verification target:

- Ubuntu 26.04 LTS with Swift 6.3.2

If you are working in the first-party monorepo, also verify that AtomLens still
builds against the package product:

```bash
xcodebuild -project /path/to/AtomLens.xcodeproj -scheme AtomLens -configuration Debug -sdk macosx CODE_SIGNING_ALLOWED=NO build
```

## 2. Review the Docs Set

Before release, make sure these documents are current:

- `README.md`
- `Documentation/API.md`
- `Documentation/INTEGRATION.md`
- `Documentation/LINUX.md`
- `Documentation/ARCHITECTURE.md`
- `Documentation/CDK_COMPARISON.md`
- `Documentation/MACOS.md`
- `CHANGELOG.md`

## 3. Prepare the Release Commit

Choose the release version and commit the release-prep changes:

```bash
VERSION=1.4.2
git add .
git commit -m "Prepare v${VERSION} release"
```

## 4. Tag the Release

```bash
git tag -a "${VERSION}" -m "CDKSwiftNativePort ${VERSION}"
```

## 5. Push the Branch and Tag

If the repository already exists:

```bash
git push origin main
git push origin --tags
```

If this is the first push, configure the remote first:

```bash
gh repo create SaschaLosko/CDKSwiftNativePort --public --source . --remote origin --push
```

## 6. Create the GitHub Release

Using GitHub CLI:

```bash
gh release create "${VERSION}" --title "CDKSwiftNativePort ${VERSION}" --notes-file CHANGELOG.md
```

Recommended release notes content:

- major additions
- compatibility notes
- breaking changes, if any
- migration notes
- known limitations

## 7. Consumer Guidance

Consumers should add:

```text
https://github.com/SaschaLosko/CDKSwiftNativePort.git
```

and pin either:

- an exact tag, or
- a semver range from the latest published release

If an important feature exists only on `main`, keep that explicit in the docs
instead of silently recommending an unreleased revision.

## 8. After the Release

- add a fresh `Unreleased` section in `CHANGELOG.md`
- track follow-up parity or bugfix work as issues
- refresh the documentation when new public APIs or supported formats change
