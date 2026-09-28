#!/usr/bin/env python3
"""Freeze the 42-sequence PhaDED competition-review panel.

This is an input/provenance builder only. It does not infer function or run a
phylogenetic program. Four disjoint buckets are selected from the authoritative
priority tables and each sequence is paired with all PhaDED training controls
from the relevant superfamily. PhaC-like competition records are kept as
challenge taxa and are never silently promoted to PhaZ positives.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


BUCKETS = (
    "explicit_extracellular",
    "phac_vs_phaz",
    "structure_assignment_exception",
    "ephaz_competition_review",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    accession = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if line.startswith(">"):
            if accession is not None:
                if accession in records or not chunks:
                    raise ValueError(f"duplicate or empty FASTA record: {accession}")
                records[accession] = "".join(chunks)
            accession = line[1:].split(None, 1)[0]
            chunks = []
        elif line:
            if accession is None:
                raise ValueError("sequence before FASTA header")
            chunks.append(line)
    if accession is not None:
        if accession in records or not chunks:
            raise ValueError(f"duplicate or empty FASTA record: {accession}")
        records[accession] = "".join(chunks)
    return records


def bucket(row: dict[str, str], dossier: dict[str, dict[str, str]]) -> str | None:
    accession = row["accession"]
    label = dossier.get(accession, {}).get("provisional_label", "")
    if label == "A_retain_high_priority_architecture":
        return "explicit_extracellular"
    if label == "E_hold_synthase_like_structure":
        return "phac_vs_phaz"
    if label == "E_hold_structure_assignment_review":
        return "structure_assignment_exception"
    if row.get("review_stage") == "ePhaZ_competition_review" or row.get("layer") == "ePhaZ_competition_review":
        return "ephaz_competition_review"
    return None


def prepare(priority: Path, dossier_path: Path, fasta: Path, reference: Path, output_dir: Path) -> dict[str, object]:
    priority_rows = read_tsv(priority)
    dossier_rows = read_tsv(dossier_path)
    dossier = {row["accession"]: row for row in dossier_rows}
    sequence = read_fasta(fasta)
    refs = read_fasta(reference)
    selected: list[dict[str, str]] = []
    seen: set[str] = set()
    for row in priority_rows:
        accession = row.get("accession", "")
        current = bucket(row, dossier)
        if current is None:
            continue
        if accession in seen:
            raise ValueError(f"candidate appears in multiple buckets: {accession}")
        if accession not in sequence:
            raise ValueError(f"candidate missing from FASTA: {accession}")
        seen.add(accession)
        record = dict(row)
        record["review_bucket"] = current
        record["sequence_sha256"] = hashlib.sha256(sequence[accession].encode("ascii")).hexdigest()
        record["challenge_role"] = "phaC_like_challenge" if current == "phac_vs_phaz" else "phaZ_competition_candidate"
        record["phenotype_boundary"] = "candidate_only_not_phenotype_validation"
        selected.append(record)
    counts = Counter(row["review_bucket"] for row in selected)
    expected = {"explicit_extracellular": 4, "phac_vs_phaz": 12, "structure_assignment_exception": 8, "ephaz_competition_review": 18}
    if dict(counts) != expected:
        raise ValueError(f"unexpected bucket counts: {dict(counts)}")
    selected.sort(key=lambda row: (BUCKETS.index(row["review_bucket"]), int(row.get("review_order", "0")), row["accession"]))
    output_dir.mkdir(parents=True, exist_ok=True)
    panel_path = output_dir / "competition_candidates.tsv"
    fields = sorted({key for row in selected for key in row})
    with panel_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(selected)
    accessions_path = output_dir / "competition_accessions.txt"
    accessions_path.write_text("\n".join(row["accession"] for row in selected) + "\n", encoding="ascii")
    candidate_fasta = output_dir / "competition_candidates.faa"
    with candidate_fasta.open("w", encoding="ascii", newline="\n") as handle:
        for row in selected:
            accession = row["accession"]
            handle.write(f">{accession}\n{sequence[accession]}\n")
    manifest = {
        "status": "prepared_candidate_only",
        "counts": dict(counts),
        "total": len(selected),
        "inputs": {str(path.name): {"path": str(path.resolve()), "size": path.stat().st_size, "sha256": sha256(path)} for path in (priority, dossier_path, fasta, reference)},
        "outputs": {str(path.name): {"path": str(path.resolve()), "size": path.stat().st_size, "sha256": sha256(path)} for path in (panel_path, accessions_path, candidate_fasta)},
        "reference_control_count": len(refs),
        "reference_control_note": "Reference sequences are controls for placement. Challenge sequences are retained as challenge taxa and are not PhaZ-positive evidence.",
        "phenotype_boundary": "PhaDED, domain, motif, SignalP and phylogenetic evidence are candidate-only and do not prove PHB/PHA degradation.",
    }
    (output_dir / "competition_input_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--priority", type=Path, required=True)
    parser.add_argument("--dossier", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.priority, args.dossier, args.fasta, args.reference, args.output_dir)
    print(json.dumps({"status": result["status"], "counts": result["counts"], "total": result["total"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
