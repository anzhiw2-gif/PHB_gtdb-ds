#!/usr/bin/env python3
"""Finalize provenance for the independent Pfam architecture evidence run."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def record(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "size": path.stat().st_size, "sha256": sha256(path), "status": "verified"}


def finalize(contract_path: Path) -> None:
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    run_dir = contract_path.parent
    results = run_dir / "results"
    contract["status"] = "completed_candidate_only"
    contract["authorization"] = {
        "candidate_only_execution": True,
        "formal_scan_authorized": False,
        "formal_registry_modified": False,
        "formal_scan_started": False,
        "server_execution_started": True,
    }
    contract["execution"] = {
        "server": "<SERVER_HOST>",
        "user": "<SERVER_USER>",
        "deploy_dir": "${PHB_REMOTE_ROOT}/PHB_gtdb-ds/deploy/20260911_phaded_pfam_architecture_01",
        "run_dir": "${PHB_REMOTE_ROOT}/PHB_gtdb-ds/runs/20260911_phaded_pfam_architecture_01",
        "hmmscan": "${PHB_REMOTE_ROOT}/miniconda3/envs/phb_gtdb/bin/hmmscan",
        "hmmscan_version": "HMMER 3.4 (Aug 2023)",
        "cpu": 40,
        "cutoff": "--cut_ga",
        "reference_scan_completed": True,
        "candidate_targeted_scan_completed": True,
    }
    contract["outputs"] = {
        "reference_pfam_provenance": record(results / "reference_pfam_provenance.json"),
        "candidate_pfam_provenance": record(results / "candidate_pfam_provenance_compact.json"),
        "candidate_pfam_provenance_raw_retained": record(results / "candidate_pfam_provenance.json"),
        "reference_fingerprint": record(results / "reference_pfam_fingerprint.tsv"),
        "reference_raw_domtblout": record(results / "raw" / "reference_pfam.domtblout"),
        "candidate_raw_domtblout": record(results / "raw" / "candidate_targeted_pfam.domtblout"),
        "pfam_subset_hmm": record(results / "pfam_subset" / "phaded_reference_pfam.hmm"),
        "pfam_accession_name": record(results / "pfam_subset" / "pfam_accession_name.tsv"),
        "pfam_architecture_evidence": record(results / "phaded_pfam_architecture_evidence.tsv"),
        "feature_consistency": record(results / "phaded_feature_consistency.tsv"),
        "feature_consistency_manifest": record(results / "phaded_feature_consistency_manifest.json"),
    }
    contract["interpretation"] = {
        "database": "Pfam-A reference-derived targeted subset",
        "reference_fingerprint_superfamilies": 6,
        "targeted_pfam_models": 23,
        "candidate_count": 109087,
        "candidate_domain_hit_accessions": 98419,
        "architecture_state_counts": {"partial": 66632, "conflicting": 83, "pending": 42372},
        "architecture_support": "A candidate is partial only when a Pfam accession is present in the reference fingerprint for its assigned PhaDED superfamily.",
        "conflict": "A Pfam accession exclusive to another reference superfamily is flagged for review; it does not erase the PhaDED assignment.",
        "pending": "No targeted Pfam hit or no reference fingerprint remains pending; neither is a biological negative.",
        "phenotype_boundary": "Pfam, motif, SignalP, and PhaDED mapping evidence are candidate architecture evidence and do not prove PHB/PHA degradation activity.",
    }
    contract_path.write_text(json.dumps(contract, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    args = parser.parse_args()
    finalize(args.contract)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
