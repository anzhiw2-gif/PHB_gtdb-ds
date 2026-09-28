#!/usr/bin/env python3
"""Apply conservative, candidate-only architecture review labels."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def classify_candidate(row: dict[str, str]) -> dict[str, str]:
    reasons = []
    if row.get("architecture_evidence") == "conflicting":
        reasons.append("PFAM_OTHER_SUPERFAMILY")
    if row.get("sequence_integrity") == "possible_N_truncation":
        reasons.append("POSSIBLE_N_TRUNCATION")
    if row.get("structure_quality") == "low" or float(row.get("pdb_low_fraction") or 0) >= 0.45:
        reasons.append("LOW_STRUCTURE_CONFIDENCE")
    if reasons and "PFAM_OTHER_SUPERFAMILY" in reasons:
        decision = "hold_architecture_conflict"
    elif reasons:
        decision = "hold_gene_model_or_structure"
    else:
        decision = "retain_candidate_comparison"
    return {
        "source_accession": row.get("source_accession", ""),
        "review_decision": decision,
        "decision_reason": ";".join(reasons) or "concordant_candidate_architecture",
        "phenotype_boundary": "candidate_only",
    }


def review(input_path: Path, output_dir: Path) -> dict[str, object]:
    with input_path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    decisions = []
    for row in rows:
        decision = classify_candidate(row)
        decision.update({"merged_priority": row.get("merged_priority", ""), "superfamily": row.get("superfamily", ""), "architecture_evidence": row.get("architecture_evidence", ""), "pfam_accessions": row.get("pfam_accessions", ""), "interpro_secondary_state": row.get("interpro_secondary_state", ""), "structure_quality": row.get("structure_quality", ""), "pdb_low_fraction": row.get("pdb_low_fraction", ""), "sequence_integrity": row.get("sequence_integrity", "")})
        decisions.append(decision)
    output_dir.mkdir(parents=True, exist_ok=True)
    table = output_dir / "candidate_architecture_review.tsv"
    fields = list(decisions[0]) if decisions else []
    with table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(decisions)
    counts = Counter(row["review_decision"] for row in decisions)
    summary = output_dir / "candidate_architecture_review_summary.tsv"
    with summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["review_decision", "candidate_count"]); writer.writerows(sorted(counts.items()))
    manifest = {"schema_version": "1.0", "status": "completed_candidate_only", "candidate_count": len(decisions), "input": {"path": str(input_path.resolve()), "size": input_path.stat().st_size, "sha256": sha256(input_path)}, "outputs": {"review": {"path": str(table.resolve()), "size": table.stat().st_size, "sha256": sha256(table)}, "summary": {"path": str(summary.resolve()), "size": summary.stat().st_size, "sha256": sha256(summary)}}, "review_decision_counts": dict(sorted(counts.items())), "phenotype_boundary": "Architecture review labels are candidate-only and do not establish PHB/PHA degradation or curate training positives."}
    (output_dir / "candidate_architecture_review_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    run_dir = output_dir.parent
    (run_dir / "logs").mkdir(parents=True, exist_ok=True)
    (run_dir / "inputs").mkdir(parents=True, exist_ok=True)
    (run_dir / "input_contract.json").write_text(json.dumps({"schema_version": "1.0", "run_id": run_dir.name, "status": "completed_candidate_only", "inputs": {"review_table": manifest["input"]}, "outputs": manifest["outputs"], "phenotype_boundary": manifest["phenotype_boundary"]}, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(review(args.input, args.output_dir), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
