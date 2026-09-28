#!/usr/bin/env python3
"""Prepare a fail-closed external PhaDED panel acquisition table."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from collections import defaultdict


FIELDS = (
    "accession", "role", "evidence_status", "family_binding_status", "bound_family", "training_accession",
    "pmid", "doi", "source_database", "sequence_sha256", "independence_status",
    "eligibility", "new_family_call_created",
)


def classify_candidate(role: str, binding_status: str, evidence_status: str, training_accession: bool = False) -> dict[str, str]:
    role = role.strip()
    binding_status = binding_status.strip()
    evidence_status = evidence_status.strip()
    resolved = binding_status.startswith("resolved_existing_family")
    if role == "independent_positive" and training_accession:
        eligibility = "not_eligible_training_accession"
    elif role == "independent_positive":
        eligibility = "eligible_positive_pending_independence" if resolved and evidence_status == "experimental_positive" else "not_eligible_unresolved_family"
    elif role == "family_resolved_negative":
        eligibility = "eligible_negative_pending_validation" if resolved and evidence_status == "experimental_negative" else "not_eligible_unresolved_family"
    else:
        eligibility = "eligible_challenge" if resolved else "challenge_unresolved_family"
    return {
        "role": role,
        "eligibility": eligibility,
        "independence_status": "pending_family_resolution" if not resolved else "pending_leakage_check",
        "new_family_call_created": "false",
    }


def bind_existing_family(accession: str, ledger_path: str | Path) -> tuple[str, str]:
    exact: dict[str, str] = {}
    by_base: dict[str, set[str]] = defaultdict(set)
    with Path(ledger_path).open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            acc = row.get("accession", "").strip()
            family = row.get("phaded_family_id", "").strip()
            if acc and family:
                exact.setdefault(acc, family)
                by_base[acc.split(".", 1)[0]].add(family)
    accession = accession.strip()
    if accession in exact:
        return "resolved_existing_family", exact[accession]
    families = by_base.get(accession.split(".", 1)[0], set())
    if len(families) == 1:
        return "resolved_existing_family_versionless", next(iter(families))
    if len(families) > 1:
        return "ambiguous_existing_family", ""
    return "unresolved", ""


def write_panel(path: str | Path, rows: list[dict[str, str]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for source in rows:
            training_accession = source.get("training_accession", "").strip().lower() == "true"
            row = {field: source.get(field, "") for field in FIELDS}
            row.update(classify_candidate(row["role"], row["family_binding_status"], row["evidence_status"], training_accession))
            writer.writerow(row)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--ledger", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.input.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("external panel input is empty")
    if args.ledger:
        training = {}
        manifest = Path("runs/20260915_phaded_reference_panel_acquisition_02/inputs/profile_manifest.tsv")
        if manifest.is_file():
            with manifest.open(encoding="utf-8-sig", newline="") as handle:
                for profile in csv.DictReader(handle, delimiter="\t"):
                    family = profile.get("phaded_family_id", "").strip()
                    training[family] = {x.strip().split(".", 1)[0] for x in profile.get("training_accessions", "").split(";") if x.strip()}
        for row in rows:
            if not row.get("bound_family", "").strip():
                binding, family = bind_existing_family(row.get("accession", ""), args.ledger)
                row["family_binding_status"] = binding
                row["bound_family"] = family
            row["training_accession"] = "true" if row.get("bound_family", "").strip() and row.get("accession", "").split(".", 1)[0] in training.get(row.get("bound_family", "").strip(), set()) else "false"
    write_panel(args.output, rows)
    summary = {
        "status": "completed_candidate_only",
        "input": {"path": str(args.input.resolve()), "sha256": _sha256(args.input)},
        "output": {"path": str(args.output.resolve()), "sha256": _sha256(args.output)},
        "records": len(rows),
        "new_family_call_created": False,
    }
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
