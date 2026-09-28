#!/usr/bin/env python3
"""Prepare an isolated, candidate-only archaeal PhaZh1 run."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import urllib.request
from pathlib import Path
from typing import Iterable, Mapping


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from run_context import create_run_layout, sha256_file, write_input_contract  # noqa: E402


REQUIRED_ACCESSIONS = {"I3RBH0", "M1XPT2"}


def _write_fasta(path: Path, records: Iterable[Mapping[str, object]], family: str) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            accession = str(record["accession"])
            sequence = "".join(str(record["sequence"]).split()).upper()
            organism = str(record.get("organism", "unknown")).replace("|", "/")
            handle.write(f">{accession}|{family}|{organism}\n{sequence}\n")


def _validate_records(records: list[Mapping[str, object]], *, require_curated: bool = False) -> None:
    seen = set()
    for record in records:
        accession = str(record.get("accession", "")).strip()
        sequence = "".join(str(record.get("sequence", "")).split())
        if not accession or accession in seen:
            raise ValueError(f"invalid or duplicate accession: {accession}")
        if not sequence or any(char not in "ACDEFGHIKLMNPQRSTVWYX" for char in sequence.upper()):
            raise ValueError(f"invalid sequence for {accession}")
        seen.add(accession)
    missing = REQUIRED_ACCESSIONS - seen
    if require_curated and missing:
        raise ValueError(f"missing required accessions: {', '.join(sorted(missing))}")


def prepare_run(repo_root: Path | str, run_id: str, records: list[Mapping[str, object]], negative: list[Mapping[str, object]]) -> Path:
    """Create the new run and return its path; no production files are touched."""
    _validate_records(records, require_curated=True)
    if negative:
        _validate_records([{"accession": item["accession"], "sequence": item["sequence"]} for item in negative])
    repo_root = Path(repo_root)
    run = create_run_layout(repo_root, run_id)
    seed = run / "inputs" / "seed_panel.faa"
    neg = run / "inputs" / "negative_panel.faa"
    manifest = run / "inputs" / "reference_manifest.tsv"
    _write_fasta(seed, records, "PhaZh1_like_candidate")
    _write_fasta(neg, negative, "negative_control")
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["accession", "organism", "evidence", "sequence_length"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow({"accession": record["accession"], "organism": record.get("organism", ""), "evidence": record.get("evidence", ""), "sequence_length": len("".join(str(record["sequence"]).split()))})
    metadata = {
        "run_id": run_id,
        "status": "candidate-only",
        "required_accessions": sorted(REQUIRED_ACCESSIONS),
        "formal_scan_authorized": False,
        "formal_registry_modified": False,
        "inputs": {"seed_panel": sha256_file(seed), "negative_panel": sha256_file(neg), "reference_manifest": sha256_file(manifest)},
    }
    (run / "run_metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    write_input_contract(run, run_id=run_id, inputs={"seed_panel": seed, "negative_panel": neg, "reference_manifest": manifest})
    return run


def fetch_uniprot_records(accessions: Iterable[str], raw_dir: Path) -> tuple[list[dict[str, object]], dict[str, Path]]:
    """Fetch accession JSON records and retain exact raw responses for provenance."""
    raw_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, object]] = []
    raw_paths: dict[str, Path] = {}
    for accession in accessions:
        accession = str(accession).strip()
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
        with urllib.request.urlopen(url, timeout=60) as response:
            payload = response.read()
        raw_path = raw_dir / f"{accession}.json"
        raw_path.write_bytes(payload)
        raw_paths[accession] = raw_path
        obj = json.loads(payload.decode("utf-8"))
        sequence = (obj.get("sequence") or {}).get("value", "")
        organism = (obj.get("organism") or {}).get("scientificName", "")
        description = obj.get("proteinDescription") or {}
        recommended = (description.get("recommendedName") or {}).get("fullName", {}).get("value", "")
        entry_type = obj.get("entryType", "")
        records.append({
            "accession": accession,
            "sequence": sequence,
            "organism": organism,
            "evidence": f"UniProt REST {url}; entry_type={entry_type}; protein_name={recommended}",
        })
    return records, raw_paths


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--seed-json", type=Path, required=True)
    parser.add_argument("--negative-json", type=Path, required=True)
    args = parser.parse_args()
    records = json.loads(args.seed_json.read_text(encoding="utf-8"))
    negative = json.loads(args.negative_json.read_text(encoding="utf-8"))
    print(prepare_run(args.repo_root, args.run_id, records, negative))


if __name__ == "__main__":
    main()
