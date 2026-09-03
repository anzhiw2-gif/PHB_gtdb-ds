#!/usr/bin/env python3
"""Create provenance-bound ePhaZ control panels.

The six annotation-only accessions that are strongly iPhaZ-like remain in a
separate challenge FASTA. They are never counted as ePhaZ positives.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


SUPPORTED_CLASSES = {
    "ePhaZ_curated_core",
    "ePhaZ_architecture_remote",
    "iPhaZ_like_challenge",
}
EPHAZ_PREFIX = "e-phaz"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        missing = required - fields
        if missing:
            raise ValueError(f"{path} missing columns: {sorted(missing)}")
        return [{(key or "").strip(): (value or "").strip() for key, value in row.items()} for row in reader]


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header = None
    sequence: list[str] = []

    def finish() -> None:
        nonlocal header, sequence
        if header is None:
            return
        accession = header.split("|", 1)[0].strip()
        value = "".join(sequence).strip().upper()
        if not accession or not value or accession in records:
            raise ValueError(f"invalid or duplicate FASTA record: {accession!r}")
        records[accession] = (header, value)

    for raw in path.read_text(encoding="utf-8-sig").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            header, sequence = line[1:].strip(), []
        else:
            if header is None:
                raise ValueError(f"sequence precedes FASTA header: {path}")
            sequence.append(line)
    finish()
    if not records:
        raise ValueError(f"empty FASTA: {path}")
    return records


def _write_fasta(path: Path, accessions: list[str], records: dict[str, tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for accession in accessions:
            header, sequence = records[accession]
            handle.write(f">{header}\n{sequence}\n")


def _file_record(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": sha256_file(path)}


def reclassify_controls(controls_tsv: Path | str, positive_faa: Path | str, class_manifest: Path | str, outdir: Path | str) -> dict[str, object]:
    controls_tsv, positive_faa, class_manifest, outdir = map(Path, (controls_tsv, positive_faa, class_manifest, outdir))
    controls = read_tsv(controls_tsv, {"accession", "label", "query_group"})
    records = read_fasta(positive_faa)
    class_rows = read_tsv(class_manifest, {"accession", "control_class"})
    classes: dict[str, str] = {}
    rationales: dict[str, str] = {}
    for row in class_rows:
        accession = row["accession"]
        if not accession or accession in classes:
            raise ValueError(f"duplicate/empty class accession: {accession!r}")
        control_class = row["control_class"]
        if control_class not in SUPPORTED_CLASSES:
            raise ValueError(f"unsupported control class for {accession}: {control_class}")
        classes[accession] = control_class
        rationales[accession] = row.get("rationale", "")

    ephaz_rows = [row for row in controls if row["label"].lower() == "positive" and row["query_group"].lower().startswith(EPHAZ_PREFIX)]
    ephaz_accessions = {row["accession"] for row in ephaz_rows}
    if ephaz_accessions != set(classes):
        missing = sorted(ephaz_accessions - set(classes))
        extra = sorted(set(classes) - ephaz_accessions)
        raise ValueError(f"class manifest does not exactly cover ePhaZ positives; missing={missing}, extra={extra}")
    if not ephaz_accessions.issubset(records):
        raise ValueError("positive FASTA is missing classified ePhaZ accession(s)")

    panels = {
        "ePhaZ_curated_core": sorted(acc for acc, cls in classes.items() if cls == "ePhaZ_curated_core"),
        "ePhaZ_architecture_remote": sorted(acc for acc, cls in classes.items() if cls == "ePhaZ_architecture_remote"),
        "iPhaZ_like_challenge": sorted(acc for acc, cls in classes.items() if cls == "iPhaZ_like_challenge"),
    }
    if not panels["ePhaZ_curated_core"]:
        raise ValueError("at least one ePhaZ_curated_core control is required")

    outdir.mkdir(parents=True, exist_ok=False)
    _write_fasta(outdir / "ephaz_positive_controls.faa", panels["ePhaZ_curated_core"] + panels["ePhaZ_architecture_remote"], records)
    _write_fasta(outdir / "iPhaZ_like_challenge.faa", panels["iPhaZ_like_challenge"], records)

    with (outdir / "ephaz_control_classes.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["accession", "control_class", "rationale"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for accession in sorted(classes):
            writer.writerow({"accession": accession, "control_class": classes[accession], "rationale": rationales.get(accession, "")})

    outputs = ["ephaz_positive_controls.faa", "iPhaZ_like_challenge.faa", "ephaz_control_classes.tsv"]
    manifest = {
        "schema_version": 1,
        "classification_policy": "six annotation-only iPhaZ-like accessions are challenge controls, never ePhaZ positives",
        "inputs": {name: _file_record(path) for name, path in (("controls_tsv", controls_tsv), ("positive_faa", positive_faa), ("class_manifest", class_manifest))},
        "counts": {key: len(value) for key, value in panels.items()},
        "outputs": {name: _file_record(outdir / name) for name in outputs},
    }
    (outdir / "control_governance.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"outdir": str(outdir), "counts": manifest["counts"], "outputs": outputs + ["control_governance.json"]}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--controls", required=True)
    parser.add_argument("--positive-faa", required=True)
    parser.add_argument("--class-manifest", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args(argv)
    try:
        reclassify_controls(args.controls, args.positive_faa, args.class_manifest, args.outdir)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
