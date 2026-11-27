# CLAUDE.md

@AGENTS.md

## Notes for Claude

`AGENTS.md` (imported above) is the package's standing contract.

- Run the full suite with `swift test` on macOS and on Linux (ARM64 Ubuntu)
  before publishing a release. Keep Linux build artifacts outside the checkout,
  for example with `--scratch-path /tmp/cdk-swift-native-port-build`.
- Lint changed Swift files strictly with this package's `.swift-format`
  configuration before committing.
- A release is an immutable tag pushed to GitHub plus a GitHub release with a
  `CHANGELOG.md` entry. See `PUBLISHING.md`.
- Keep this repository free of consumer-specific names and paths, including in
  this file.
