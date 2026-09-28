#!/usr/bin/env python3
"""Prepare the frozen candidate FASTA and profile inputs for a mapping run."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "size": path.stat().st_size, "sha256": sha256(path), "status": "verified"}


def read_fasta(path: Path, records: dict[str, str]) -> None:
    current: str | None = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current is not None:
                sequence = "".join(chunks)
                if current in records and records[current] != sequence:
                    raise ValueError(f"conflicting FASTA sequence for accession: {current}")
                records[current] = sequence
            current = line[1:].split()[0]
            chunks = []
        else:
            if current is None:
                raise ValueError(f"sequence before FASTA header: {path}")
            chunks.append(line)
    if current is not None:
        sequence = "".join(chunks)
        if current in records and records[current] != sequence:
            raise ValueError(f"conflicting FASTA sequence for accession: {current}")
        records[current] = sequence


def prepare(run_dir: Path, source_run: Path, profile_run: Path) -> dict[str, object]:
    inputs = run_dir / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    layer_path = source_run / "results" / "layers_retry_02" / "protein_layers.tsv"
    rows = list(csv.DictReader(layer_path.open(encoding="utf-8-sig", newline=""), delimiter="\t"))
    seen: set[str] = set()
    for row in rows:
        accession = row.get("accession", "").strip()
        if not accession or accession in seen:
            raise ValueError(f"duplicate or missing candidate accession: {accession!r}")
        seen.add(accession)

    sequences: dict[str, str] = {}
    for name in ("ePhaZ_tier1.faa", "ePhaZ_tier2.faa", "iPhaZ_tier1.faa", "iPhaZ_tier2.faa"):
        read_fasta(source_run / "inputs" / name, sequences)
    missing = sorted(seen - set(sequences))
    if missing:
        raise ValueError(f"candidate sequences missing from tier FASTA: {len(missing)}")

    fasta = inputs / "candidate_union.faa"
    with fasta.open("w", encoding="ascii", newline="\n") as handle:
        for row in rows:
            accession = row["accession"]
            handle.write(f">{accession}\n{sequences[accession]}\n")
    source_map = inputs / "candidate_source_map.tsv"
    with source_map.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    profile_manifest = profile_run / "results" / "profile_build_01" / "profile_manifest.tsv"
    shutil.copy2(profile_manifest, inputs / "profile_manifest.tsv")
    trained: list[str] = []
    with profile_manifest.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["model_status"] == "trained":
                trained.append(row["profile_id"])
    profile_dir = profile_run / "results" / "profile_build_01" / "profiles"
    for profile_id in trained:
        shutil.copy2(profile_dir / f"{profile_id}.hmm", inputs / f"{profile_id}.hmm")

    contract_path = run_dir / "input_contract.json"
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    contract["status"] = "prepared_candidate_only"
    contract["mapping_policy"]["score_gap_rule"] = (
        "min_score_gap=1.0 bits; declared before HMMER results; "
        "ties and gaps below threshold remain ambiguous"
    )
    contract["mapping_policy"]["trained_profile_count"] = len(trained)
    contract["mapping_policy"]["trained_profile_ids"] = trained
    contract["inputs"]["candidate_source_map"] = file_record(source_map)
    contract["inputs"]["candidate_union_fasta"] = file_record(fasta)
    contract["inputs"]["profile_manifest"] = file_record(inputs / "profile_manifest.tsv")
    contract["inputs"]["trained_profiles"] = [file_record(inputs / f"{profile_id}.hmm") for profile_id in trained]
    contract["candidate_counts"] = {
        "layer_rows": len(rows),
        "unique_accessions": len(seen),
        "sequence_records": len(sequences),
        "candidate_fasta_records": len(rows),
        "missing_sequences": len(missing),
    }
    contract_path.write_text(json.dumps(contract, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return {"candidates": len(rows), "sequences": len(sequences), "trained_profiles": len(trained), "missing": len(missing)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--source-run", type=Path, required=True)
    parser.add_argument("--profile-run", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.run_dir, args.source_run, args.profile_run), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
