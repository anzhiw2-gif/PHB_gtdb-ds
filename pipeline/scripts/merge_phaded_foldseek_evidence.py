#!/usr/bin/env python3
"""Merge accession-bound Foldseek structure comparisons into a PhaDED matrix.

Foldseek evidence is candidate-only structural homology evidence.  It must not
be interpreted as PHB/PHA degradation phenotype validation and does not resolve
profile, Pfam architecture, or family ambiguity by itself.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Iterable, Mapping


PHENOTYPE_BOUNDARY = (
    "Foldseek structure comparison is candidate-only structural evidence and "
    "does not validate PHB/PHA degradation phenotype."
)

FOLDSEEK_FIELDS = [
    "foldseek_status",
    "foldseek_reference_model_id",
    "foldseek_reference_class",
    "foldseek_fident",
    "foldseek_alnlen",
    "foldseek_alntmscore",
    "foldseek_qtmscore",
    "foldseek_ttmscore",
    "foldseek_rmsd",
    "foldseek_structural_evidence_band",
    "foldseek_interpretation",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[Mapping[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _status_from_band(band: str) -> str:
    normalized = (band or "").strip().lower()
    if normalized in {"moderate", "moderate_high", "high"}:
        return "foldseek_support_candidate_only"
    if normalized in {"conflict", "discordant"}:
        return "foldseek_conflict_candidate_only"
    if normalized:
        return "foldseek_review_candidate_only"
    return "foldseek_review_candidate_only"


def read_foldseek_summary(path: Path) -> dict[str, dict[str, str]]:
    """Read structural_evidence_summary.tsv keyed by candidate/accession id."""
    rows = read_tsv(path)
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        accession = (row.get("candidate_id") or row.get("accession") or row.get("source_accession") or "").strip()
        if not accession:
            raise ValueError("Foldseek summary contains an empty candidate_id/accession")
        if accession in indexed:
            raise ValueError(f"duplicate Foldseek candidate_id/accession: {accession}")
        band = row.get("structural_evidence_band", "").strip()
        indexed[accession] = {
            "foldseek_status": _status_from_band(band),
            "foldseek_reference_model_id": row.get("reference_model_id", "").strip() or "pending",
            "foldseek_reference_class": row.get("class", "").strip() or "pending",
            "foldseek_fident": row.get("fident", "").strip(),
            "foldseek_alnlen": row.get("alnlen", "").strip(),
            "foldseek_alntmscore": row.get("alntmscore", "").strip(),
            "foldseek_qtmscore": row.get("qtmscore", "").strip(),
            "foldseek_ttmscore": row.get("ttmscore", "").strip(),
            "foldseek_rmsd": row.get("rmsd", "").strip(),
            "foldseek_structural_evidence_band": band or "pending",
            "foldseek_interpretation": row.get("interpretation", "").strip() or "pending",
        }
    return indexed


def merge_foldseek_into_matrix(
    matrix_rows: list[Mapping[str, str]],
    foldseek_rows: Mapping[str, Mapping[str, str]],
) -> list[dict[str, str]]:
    """Add Foldseek columns while preserving matrix order and row count."""
    accessions = [str(row.get("accession", "")).strip() for row in matrix_rows]
    if any(not accession for accession in accessions) or len(accessions) != len(set(accessions)):
        raise ValueError("matrix contains empty or duplicate accessions")
    matrix_set = set(accessions)
    outside = sorted(set(foldseek_rows) - matrix_set)
    if outside:
        raise ValueError(f"Foldseek summary contains accession outside matrix: {outside[0]}")

    merged: list[dict[str, str]] = []
    pending = {
        "foldseek_status": "foldseek_not_tested_full_library",
        "foldseek_reference_model_id": "pending",
        "foldseek_reference_class": "pending",
        "foldseek_fident": "",
        "foldseek_alnlen": "",
        "foldseek_alntmscore": "",
        "foldseek_qtmscore": "",
        "foldseek_ttmscore": "",
        "foldseek_rmsd": "",
        "foldseek_structural_evidence_band": "not_tested_full_library",
        "foldseek_interpretation": "not tested in bounded accession-bound Foldseek comparison",
    }
    for row in matrix_rows:
        accession = str(row["accession"]).strip()
        item = dict(row)
        item.update(foldseek_rows.get(accession, pending))
        merged.append(item)
    return merged


def _fieldnames(rows: list[Mapping[str, str]]) -> list[str]:
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    for field in FOLDSEEK_FIELDS:
        if field not in fields:
            fields.append(field)
    return fields


def merge_files(matrix_path: Path, foldseek_path: Path, output_path: Path) -> dict[str, object]:
    matrix = read_tsv(matrix_path)
    foldseek = read_foldseek_summary(foldseek_path)
    merged = merge_foldseek_into_matrix(matrix, foldseek)
    fields = _fieldnames(merged)
    write_tsv(output_path, merged, fields)
    counts = Counter(row["foldseek_status"] for row in merged)
    manifest = {
        "schema_version": "1.0",
        "status": "completed_candidate_only",
        "candidate_count": len(merged),
        "accession_bound_foldseek_count": len(foldseek),
        "foldseek_status_counts": dict(sorted(counts.items())),
        "inputs": {
            "matrix": {
                "path": str(matrix_path.resolve()),
                "size": matrix_path.stat().st_size,
                "sha256": sha256(matrix_path),
            },
            "foldseek_summary": {
                "path": str(foldseek_path.resolve()),
                "size": foldseek_path.stat().st_size,
                "sha256": sha256(foldseek_path),
            },
        },
        "output": {
            "path": str(output_path.resolve()),
            "size": output_path.stat().st_size,
            "sha256": sha256(output_path),
        },
        "boundary": {
            "phenotype": PHENOTYPE_BOUNDARY,
            "classification": "Foldseek columns are an additional evidence layer and do not change subtype_call by themselves.",
        },
    }
    manifest_path = output_path.parent / "foldseek_merge_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    run_dir = output_path.parent.parent
    run_dir.mkdir(parents=True, exist_ok=True)
    contract = {
        "schema_version": "1.0",
        "run_id": run_dir.name,
        "status": "completed_candidate_only",
        "inputs": manifest["inputs"],
        "outputs": {
            "merged_matrix": manifest["output"],
            "merge_manifest": {
                "path": str(manifest_path.resolve()),
                "size": manifest_path.stat().st_size,
                "sha256": sha256(manifest_path),
            },
        },
        "phenotype_boundary": PHENOTYPE_BOUNDARY,
        "classification_boundary": manifest["boundary"]["classification"],
    }
    (run_dir / "input_contract.json").write_text(
        json.dumps(contract, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    (run_dir / "logs").mkdir(parents=True, exist_ok=True)
    (run_dir / "logs" / "foldseek_merge.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--foldseek-summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    print(json.dumps(merge_files(args.matrix, args.foldseek_summary, args.output), sort_keys=True, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
