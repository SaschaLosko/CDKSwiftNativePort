# Standalone library contract

- Maintain CDKSwiftNativePort as a self-contained Swift package. Its only native
  implementation dependency is the bundled C InChI target.
- Keep source, tests, documentation, scripts and metadata independent of named
  consuming applications, customers and private presentation context.
- Do not add Java source, compiled Java artifacts, Maven configuration or JVM
  requirements. External upstream comparisons may contribute attributed reference
  data; the package's tests must consume that data without an upstream runtime.
- Build and test from this checkout. Do not inspect or depend on sibling app
  repositories in package tests or use hard-coded developer checkout paths.
- Validate package changes with the full Swift test suite on macOS and Linux;
  use an isolated checkout when changing package boundaries. Apply this package's
  `.swift-format` configuration and lint changed Swift files strictly.
- Document supported chemistry and unresolved cases explicitly. Do not equate
  partial enumeration, missing definitions or resource limits with a negative
  membership result. Preserve upstream licenses and dataset attribution.
