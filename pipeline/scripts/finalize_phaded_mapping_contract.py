#!/usr/bin/env python3
"""Bind completed candidate-only mapping outputs to the run contract."""

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


def finalize(run_dir: Path, deploy_manifest: Path | None = None) -> None:
    contract_path = run_dir / "input_contract.json"
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    results = run_dir / "results"
    inputs = run_dir / "inputs"
    contract["status"] = "completed_candidate_only"
    contract["authorization"]["server_execution_started"] = True
    contract["authorization"]["formal_scan_authorized"] = False
    contract["authorization"]["formal_registry_modified"] = False
    contract["authorization"]["formal_scan_started"] = False
    contract["execution"] = {
        "server": "<SERVER_HOST>",
        "user": "<SERVER_USER>",
        "deploy_dir": "${PHB_REMOTE_ROOT}/PHB_gtdb-ds/deploy/20260909_phaded_candidate_mapping_01",
        "run_dir": "${PHB_REMOTE_ROOT}/PHB_gtdb-ds/runs/20260909_phaded_candidate_mapping_01",
        "hmmsearch_cpu": 40,
        "hmmsearch_profiles_sequential": True,
        "hmmsearch_max_threads": 40,
        "hmmsearch_version": "HMMER 3.4 (Aug 2023)",
    }
    contract["outputs"] = {
        "profile_scores": record(results / "phaded_profile_scores.tsv"),
        "assignment_table": record(results / "phaded_assignment.tsv"),
        "mapping_summary": record(results / "mapping_summary.json"),
        "mapping_reconciliation": record(results / "mapping_reconciliation.json"),
        "hmmsearch_provenance": record(results / "hmmsearch_provenance.json"),
        "feature_annotation": record(results / "phaded_feature_annotation.tsv"),
        "assignment_evidence": record(results / "phaded_assignment_evidence.tsv"),
        "feature_input_manifest": record(inputs / "feature_input_manifest.json"),
        "feature_annotation_manifest": record(results / "feature_annotation_manifest.json"),
        "assignment_evidence_manifest": record(results / "assignment_evidence_manifest.json"),
    }
    contract["interpretation"] = {
        "candidate_only": True,
        "phenotype_assignment_forbidden": True,
        "trained_profile_scope": "Only 9 trained profiles were used; reference_only families were not used as calibrated HMMs.",
        "unassigned_scope": "unassigned_PhaDED_like means no score was reported by the trained profile set at the declared reporting setting; it is not a biological negative.",
    }
    contract_path.write_text(json.dumps(contract, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if deploy_manifest is not None:
        manifest = json.loads(deploy_manifest.read_text(encoding="utf-8"))
        manifest["status"] = "executed_candidate_only"
        manifest["server_execution_started"] = True
        manifest["formal_scan_authorized"] = False
        manifest["formal_registry_modified"] = False
        manifest["formal_scan_started"] = False
        manifest["execution"] = contract["execution"]
        manifest["bound_inputs"] = {
            "input_contract_sha256": sha256(contract_path),
            "candidate_union_fasta_sha256": sha256(inputs / "candidate_union.faa"),
            "profile_manifest_sha256": sha256(inputs / "profile_manifest.tsv"),
        }
        deploy_manifest.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--deploy-manifest", type=Path)
    args = parser.parse_args()
    finalize(args.run_dir, args.deploy_manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
