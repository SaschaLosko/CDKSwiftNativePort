#!/usr/bin/env python3
"""Optional independent read-back gate. Indigo is a test tool, never a dependency."""
import argparse
import ctypes
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indigo-library", required=True, type=Path)
    parser.add_argument("--fixtures", required=True, type=Path)
    args = parser.parse_args()
    library = ctypes.CDLL(str(args.indigo_library.resolve()))
    library.indigoAllocSessionId.restype = ctypes.c_uint64
    library.indigoSetSessionId.argtypes = [ctypes.c_uint64]
    library.indigoReleaseSessionId.argtypes = [ctypes.c_uint64]
    for name in ("indigoGetLastError", "indigoCanonicalSmiles"):
        getattr(library, name).restype = ctypes.c_char_p
    for name in ("indigoLoadMoleculeFromBuffer", "indigoLoadReactionFromBuffer"):
        getattr(library, name).argtypes = [ctypes.c_char_p, ctypes.c_int]
    session = library.indigoAllocSessionId()
    library.indigoSetSessionId(session)
    failures = 0
    checks = 0
    try:
        fixtures = sorted(args.fixtures.glob("*.smi"))
        if not fixtures:
            raise SystemExit("No fixtures; run ChemDrawWriterTests with CDK_CHEMDRAW_EVIDENCE_DIR first.")
        for path in fixtures:
            reaction = path.stem == "reaction"
            loader = library.indigoLoadReactionFromBuffer if reaction else library.indigoLoadMoleculeFromBuffer
            source = path.read_bytes()
            original = loader(source, len(source))
            if original < 0:
                raise RuntimeError(library.indigoGetLastError())
            if reaction:
                # Indigo 1.42 does not read ReactionStepAtomMap. The Swift binary
                # decoder tests separately verify map pairs and object references.
                library.indigoClearAAM(original)
            expected = library.indigoCanonicalSmiles(original)
            for extension in (".cdx", ".cdxml"):
                data = path.with_suffix(extension).read_bytes()
                loaded = loader(data, len(data))
                actual = library.indigoCanonicalSmiles(loaded) if loaded >= 0 else None
                matches = actual is not None and expected == actual
                checks += 1
                failures += not matches
                print(f"{path.stem}{extension}: {'PASS' if matches else 'FAIL'}")
                if not matches:
                    print(f"  expected={expected!r}, actual={actual!r}, error={library.indigoGetLastError()!r}")
                if loaded >= 0:
                    library.indigoFree(loaded)
            library.indigoFree(original)
    finally:
        library.indigoReleaseSessionId(session)
    print(f"{checks} checks, {failures} failures")
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
