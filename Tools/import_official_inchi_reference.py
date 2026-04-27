#!/usr/bin/env python3

import argparse
import gzip
import json
import re
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Callable


SELECTION_LIMITS = {
    "inchi_elemental": 6,
    "inchi_isotopic": 6,
    "mcule_neutral_organic": 12,
    "mcule_charged": 6,
    "mcule_double_stereo": 6,
    "mcule_tetra_stereo": 6,
    "mcule_multicomponent": 6,
}


@dataclass
class ReferenceRecord:
    dataset: str
    record_id: str
    molfile: str
    properties: dict[str, str]
    expected_inchi: str
    expected_inchi_key: str
    expected_aux_info: str
    expected_exit: int
    expected_message: str
    formula: str
    inchi_length: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Vendor a curated official InChI reference corpus for CDKSwiftNativePort."
    )
    parser.add_argument(
        "--source",
        required=True,
        help="Path to a checked-out official InChI repository."
    )
    parser.add_argument(
        "--output",
        default=(
            "/Users/sascha/Development/AtomLens/"
            "CDKSwiftNativePort/Tests/CDKSwiftNativePortTests/InChI/OfficialReference/"
            "official_reference_cases.json"
        ),
        help="Output JSON path."
    )
    return parser.parse_args()


def git_commit(repo: Path) -> str:
    head = repo.joinpath(".git", "HEAD").read_text(encoding="utf-8").strip()
    if head.startswith("ref: "):
        ref = head.removeprefix("ref: ").strip()
        return repo.joinpath(".git", ref).read_text(encoding="utf-8").strip()
    return head


def parse_sdf_records(path: Path, id_property: str) -> list[tuple[str, str, dict[str, str]]]:
    text = gzip.open(path, "rt", encoding="utf-8", errors="ignore").read()
    blocks = [block.strip("\n") for block in text.split("$$$$") if block.strip()]
    records: list[tuple[str, str, dict[str, str]]] = []

    for block in blocks:
        if "M  END" not in block:
            continue

        molfile_part, property_part = block.split("M  END", 1)
        molfile = normalize_molfile_for_cdk(f"{molfile_part}M  END\n")
        properties: dict[str, str] = {}

        current_key: str | None = None
        current_lines: list[str] = []
        for line in property_part.splitlines():
            match = re.match(r">\s*<([^>]+)>", line)
            if match:
                if current_key is not None:
                    properties[current_key] = "\n".join(current_lines).strip()
                current_key = match.group(1).strip()
                current_lines = []
            elif current_key is not None:
                current_lines.append(line)
        if current_key is not None:
            properties[current_key] = "\n".join(current_lines).strip()

        record_id = (
            properties.get(id_property)
            or properties.get("ID")
            or properties.get("Mcule_ID")
            or properties.get("Index")
        )
        if not record_id:
            raise ValueError(f"Could not locate record identifier in {path}.")

        records.append((record_id.strip(), molfile, properties))

    return records


def normalize_molfile_for_cdk(molfile: str) -> str:
    lines = molfile.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    counts_index = None
    for index, line in enumerate(lines):
        upper = line.upper()
        if "V2000" in upper or "V3000" in upper:
            counts_index = index
            break
    if counts_index is None:
        return molfile

    while counts_index < 3:
        lines.insert(counts_index, "")
        counts_index += 1

    return "\n".join(lines).rstrip("\n") + "\n"


def load_reference_results(path: Path) -> dict[str, dict[str, object]]:
    connection = sqlite3.connect(path)
    cursor = connection.cursor()
    try:
        rows = cursor.execute("select molfile_id, result from results").fetchall()
        return {molfile_id.strip(): json.loads(result) for molfile_id, result in rows}
    finally:
        connection.close()


def inchi_formula(inchi: str) -> str:
    if inchi.startswith("InChI=1S/"):
        payload = inchi[len("InChI=1S/"):]
    elif inchi.startswith("InChI=1/"):
        payload = inchi[len("InChI=1/"):]
    else:
        payload = inchi
    return payload.split("/", 1)[0]


def build_records(dataset: str,
                  records: list[tuple[str, str, dict[str, str]]],
                  reference_results: dict[str, dict[str, object]]) -> list[ReferenceRecord]:
    out: list[ReferenceRecord] = []
    for record_id, molfile, properties in records:
        result = reference_results.get(record_id)
        if result is None:
            raise ValueError(f"Missing reference result for {dataset}:{record_id}.")

        expected_inchi = str(result["inchi"])
        out.append(
            ReferenceRecord(
                dataset=dataset,
                record_id=record_id,
                molfile=molfile,
                properties=properties,
                expected_inchi=expected_inchi,
                expected_inchi_key=str(result["key"]),
                expected_aux_info=str(result.get("aux", "")),
                expected_exit=int(result.get("exit", 0)),
                expected_message=str(result.get("message", "")),
                formula=inchi_formula(expected_inchi),
                inchi_length=len(expected_inchi),
            )
        )
    return out


def choose(records: list[ReferenceRecord],
           predicate: Callable[[ReferenceRecord], bool],
           limit: int,
           used_ids: set[tuple[str, str]]) -> list[ReferenceRecord]:
    selected: list[ReferenceRecord] = []
    for record in sorted(records, key=lambda item: (item.inchi_length, item.record_id)):
        key = (record.dataset, record.record_id)
        if key in used_ids or not predicate(record):
            continue
        used_ids.add(key)
        selected.append(record)
        if len(selected) == limit:
            break
    if len(selected) != limit:
        raise ValueError(f"Could not fill curated selection of size {limit}.")
    return selected


def main() -> None:
    args = parse_args()

    source = Path(args.source).expanduser().resolve()
    output = Path(args.output).expanduser().resolve()
    ci_root = source / "INCHI-1-TEST" / "tests" / "test_library" / "data" / "ci"

    mcule_records = build_records(
        dataset="mcule",
        records=parse_sdf_records(ci_root / "mcule.sdf.gz", "Mcule_ID"),
        reference_results=load_reference_results(ci_root / "mcule.sdf.regression_reference.sqlite"),
    )
    inchi_records = build_records(
        dataset="inchi",
        records=parse_sdf_records(ci_root / "inchi.sdf.gz", "ID"),
        reference_results=load_reference_results(ci_root / "inchi.sdf.regression_reference.sqlite"),
    )

    used_ids: set[tuple[str, str]] = set()
    selected: list[tuple[str, ReferenceRecord]] = []

    selected.extend((
        "inchi_elemental",
        record
    ) for record in choose(
        inchi_records,
        predicate=lambda record: record.record_id.startswith("_Elements.#"),
        limit=SELECTION_LIMITS["inchi_elemental"],
        used_ids=used_ids,
    ))
    selected.extend((
        "inchi_isotopic",
        record
    ) for record in choose(
        inchi_records,
        predicate=lambda record: "/i" in record.expected_inchi,
        limit=SELECTION_LIMITS["inchi_isotopic"],
        used_ids=used_ids,
    ))
    selected.extend((
        "mcule_neutral_organic",
        record
    ) for record in choose(
        mcule_records,
        predicate=lambda record: (
            record.formula.startswith("C")
            and "H" in record.formula
            and all(tag not in record.expected_inchi for tag in ("/q", "/p", "/i", "/t", "/b"))
            and "." not in record.formula
        ),
        limit=SELECTION_LIMITS["mcule_neutral_organic"],
        used_ids=used_ids,
    ))
    selected.extend((
        "mcule_charged",
        record
    ) for record in choose(
        mcule_records,
        predicate=lambda record: "/q" in record.expected_inchi or "/p" in record.expected_inchi,
        limit=SELECTION_LIMITS["mcule_charged"],
        used_ids=used_ids,
    ))
    selected.extend((
        "mcule_double_stereo",
        record
    ) for record in choose(
        mcule_records,
        predicate=lambda record: "/b" in record.expected_inchi,
        limit=SELECTION_LIMITS["mcule_double_stereo"],
        used_ids=used_ids,
    ))
    selected.extend((
        "mcule_tetra_stereo",
        record
    ) for record in choose(
        mcule_records,
        predicate=lambda record: "/t" in record.expected_inchi,
        limit=SELECTION_LIMITS["mcule_tetra_stereo"],
        used_ids=used_ids,
    ))
    selected.extend((
        "mcule_multicomponent",
        record
    ) for record in choose(
        mcule_records,
        predicate=lambda record: "." in record.formula,
        limit=SELECTION_LIMITS["mcule_multicomponent"],
        used_ids=used_ids,
    ))

    payload = {
        "schemaVersion": 1,
        "suite": "Official InChI Reference CI Corpus",
        "referenceRepository": "https://github.com/IUPAC-InChI/InChI",
        "referenceCommit": git_commit(source),
        "referenceSources": [
            "INCHI-1-TEST/tests/test_library/data/ci/inchi.sdf.gz",
            "INCHI-1-TEST/tests/test_library/data/ci/inchi.sdf.regression_reference.sqlite",
            "INCHI-1-TEST/tests/test_library/data/ci/mcule.sdf.gz",
            "INCHI-1-TEST/tests/test_library/data/ci/mcule.sdf.regression_reference.sqlite",
        ],
        "selectionLimits": SELECTION_LIMITS,
        "cases": [
            {
                "id": f"{record.dataset}:{record.record_id}",
                "dataset": record.dataset,
                "recordID": record.record_id,
                "selectionReason": selection_reason,
                "expectedInChI": record.expected_inchi,
                "expectedInChIKey": record.expected_inchi_key,
                "expectedAuxInfo": record.expected_aux_info,
                "expectedExitCode": record.expected_exit,
                "expectedMessage": record.expected_message,
                "molfile": record.molfile,
            }
            for selection_reason, record in selected
        ],
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
