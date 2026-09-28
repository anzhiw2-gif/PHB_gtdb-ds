#!/usr/bin/env python3
"""Reconcile external PhaDED controls without creating new family calls.

The input panel is evidence, not a classifier.  An accession is bound only to
an existing ledger family (exactly or by an unambiguous versionless accession).
Records that cannot be bound to one of the 37 reference-only families remain
unresolved and cannot make a profile calibration-eligible.

Rule #3 (axis-anchored credit): the family-bound counters above only see records
that bind to the profile's **own** family, so an operator may additionally pass the
discrimination-axis panel (``--axis-panel`` + ``--axis-definitions``).  Only then are
the extra ``axis_anchored_*`` columns appended, and they are computed by
``credit_phaded_axis_anchored_controls`` with a fail-closed audit requirement.  The
legacy columns and their values are never touched by that enrichment; without the two
extra inputs the readiness output is byte-for-byte what it was before this change.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
from collections import Counter, defaultdict
from pathlib import Path


PANEL_FIELDS = {"accession", "panel", "pmid", "doi", "sequence_length", "normalized_sequence_sha256", "source_sha256", "independence_status", "decision"}
AXIS_PANEL_REQUIRED = {
    "axis_id", "axis_end", "role", "accession", "evidence_status", "anchor_evidence_kind",
    "family_binding_status", "bound_family", "bound_reference_only_profile_id",
    "anchor_justification", "source_doi_or_pmid",
}
AXIS_READINESS_FIELDS = (
    "axis_anchored_negative_count", "axis_anchored_challenge_count",
    "axis_anchor_basis_status", "axis_anchored_negative_credit_ids",
    "axis_anchored_challenge_credit_ids",
)
DEFAULT_AXIS_CREDIT_SCRIPT = "credit_phaded_axis_anchored_controls.py"


def _load_credit_module(script_name: str = DEFAULT_AXIS_CREDIT_SCRIPT):
    """Load the axis-credit module from this script's own directory."""
    path = Path(__file__).resolve().parent / script_name
    spec = importlib.util.spec_from_file_location(Path(script_name).stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
PROFILE_FIELDS = {"profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id", "model_status"}
LEDGER_FIELDS = {"accession", "phaded_family_id", "phaded_superfamily", "evidence_status"}
OUTPUT_FIELDS = (
    "accession", "panel", "panel_role", "bound_family", "bound_superfamily",
    "family_binding_status", "reference_only_profile_id", "evidence_status",
    "pmid", "doi", "sequence_length", "normalized_sequence_sha256",
    "source_sha256", "independence_status", "decision", "calibration_eligibility",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read(path: Path, required: set[str]) -> list[dict[str, str]]:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"input is not a regular file: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f"{path} missing columns: {','.join(missing)}")
        return list(reader)


def _base_accession(accession: str) -> str:
    return accession.strip().split(".", 1)[0]


def _panel_role(panel: str) -> str:
    value = panel.strip().lower()
    if value == "independent_experimental_positive":
        return "independent_positive"
    if value in {"intracellular_non_ephaz_negative", "mcl_pha_non_phb_negative"}:
        return "family_resolved_negative_candidate"
    return "challenge_control"


def _binding(accession: str, exact: dict[str, dict[str, str]], by_base: dict[str, list[dict[str, str]]]) -> tuple[str, dict[str, str] | None]:
    if accession in exact:
        return "resolved_existing_family", exact[accession]
    matches = by_base.get(_base_accession(accession), [])
    if len(matches) == 1:
        return "resolved_existing_family_versionless", matches[0]
    if len(matches) > 1:
        return "ambiguous_existing_family", None
    return "unresolved", None


def reconcile(profile_path: str | Path, ledger_path: str | Path, panel_path: str | Path, output_dir: str | Path, *, expected_reference_only_count: int | None = 37, axis_panel: str | Path | None = None, axis_definitions: str | Path | None = None) -> dict[str, object]:
    profile_path, ledger_path, panel_path, output_dir = map(Path, (profile_path, ledger_path, panel_path, output_dir))
    if (axis_panel is None) != (axis_definitions is None):
        raise ValueError("axis_anchored credit requires both axis_panel and axis_definitions")
    profiles = _read(profile_path, PROFILE_FIELDS)
    ledger = _read(ledger_path, LEDGER_FIELDS)
    panel = _read(panel_path, PANEL_FIELDS)
    reference_only = [row for row in profiles if row["model_status"].strip() == "reference_only"]
    if expected_reference_only_count is not None and len(reference_only) != expected_reference_only_count:
        raise ValueError(f"expected {expected_reference_only_count} reference-only profiles, observed {len(reference_only)}")
    profile_by_family = {row["phaded_family_id"].strip(): row for row in reference_only if row["profile_kind"].strip() == "family" and row["phaded_family_id"].strip()}
    exact: dict[str, dict[str, str]] = {}
    by_base: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in ledger:
        accession = row["accession"].strip()
        if not accession:
            continue
        exact.setdefault(accession, row)
        by_base[_base_accession(accession)].append(row)
    seen: set[str] = set()
    records: list[dict[str, str]] = []
    for source in sorted(panel, key=lambda row: row["accession"].strip()):
        accession = source["accession"].strip()
        if not accession:
            raise ValueError("panel contains empty accession")
        if accession in seen:
            raise ValueError(f"duplicate panel accession: {accession}")
        seen.add(accession)
        binding_status, bound = _binding(accession, exact, by_base)
        family = (bound or {}).get("phaded_family_id", "").strip()
        superfamily = (bound or {}).get("phaded_superfamily", "").strip()
        target = profile_by_family.get(family)
        role = _panel_role(source["panel"])
        target_id = target["profile_id"].strip() if target else ""
        is_bound_reference_only = bool(target_id and binding_status.startswith("resolved"))
        if role == "independent_positive":
            evidence = "experimental_positive"
            eligible = "eligible_positive" if is_bound_reference_only else "not_eligible_unresolved_family"
        elif role == "family_resolved_negative_candidate":
            evidence = "experimental_negative" if source["panel"].strip() == "mcl_pha_non_phb_negative" else "challenge_control"
            eligible = "eligible_negative" if is_bound_reference_only and evidence == "experimental_negative" else "not_eligible_unresolved_or_nonformal"
        else:
            evidence = "challenge_control"
            eligible = "eligible_challenge" if is_bound_reference_only else "challenge_unresolved_family"
        records.append({
            "accession": accession,
            "panel": source["panel"].strip(),
            "panel_role": role,
            "bound_family": family,
            "bound_superfamily": superfamily,
            "family_binding_status": binding_status,
            "reference_only_profile_id": target_id,
            "evidence_status": evidence,
            "pmid": source.get("pmid", "").strip(),
            "doi": source.get("doi", "").strip(),
            "sequence_length": source.get("sequence_length", "").strip(),
            "normalized_sequence_sha256": source.get("normalized_sequence_sha256", "").strip(),
            "source_sha256": source.get("source_sha256", "").strip(),
            "independence_status": source.get("independence_status", "").strip(),
            "decision": source.get("decision", "").strip(),
            "calibration_eligibility": eligible,
        })

    by_profile: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in records:
        if row["reference_only_profile_id"]:
            by_profile[row["reference_only_profile_id"]].append(row)
    readiness: list[dict[str, object]] = []
    for profile in sorted(reference_only, key=lambda row: row["profile_id"].strip()):
        profile_id = profile["profile_id"].strip()
        family = profile.get("phaded_family_id", "").strip()
        existing_positive = sum(row.get("phaded_family_id", "").strip() == family and row.get("evidence_status", "").strip() == "experimental_positive" for row in ledger)
        bound_rows = by_profile.get(profile_id, [])
        positive = sum(row["calibration_eligibility"] == "eligible_positive" for row in bound_rows)
        negative = sum(row["calibration_eligibility"] == "eligible_negative" for row in bound_rows)
        challenge = sum(row["calibration_eligibility"] == "eligible_challenge" for row in bound_rows)
        decision = "reference_only_insufficient_panel" if profile["profile_kind"].strip() == "family" else "planned_not_run_superfamily_requires_family_resolution"
        readiness.append({
            "profile_id": profile_id,
            "profile_kind": profile["profile_kind"].strip(),
            "phaded_superfamily": profile["phaded_superfamily"].strip(),
            "phaded_family_id": family,
            "existing_experimental_positive": existing_positive,
            "external_bound_positive": positive,
            "external_bound_negative": negative,
            "external_bound_challenge": challenge,
            "heldout_positive": 0,
            "family_resolved_negative_status": "sufficient" if negative else "missing",
            "challenge_status": "sufficient" if challenge else "missing",
            "calibration_decision": decision,
            "new_family_call_created": "false",
        })

    axis_credit: dict[str, object] | None = None
    if axis_panel is not None:
        credit_module = _load_credit_module()
        axis_rows = _read(Path(axis_panel), AXIS_PANEL_REQUIRED)
        definitions = credit_module.load_axis_definitions(axis_definitions)
        family_superfamily: dict[str, str] = {}
        for row in profiles:
            family = row.get("phaded_family_id", "").strip()
            superfamily = row.get("phaded_superfamily", "").strip()
            if not family or not superfamily:
                continue
            previous = family_superfamily.get(family)
            if previous is not None and previous != superfamily:
                raise ValueError(f"{family}: conflicting superfamily {previous!r} vs {superfamily!r}")
            family_superfamily[family] = superfamily
        axis_credit = credit_module.attribute_axis_credit(axis_rows, definitions, reference_only, family_superfamily)
        for row in readiness:
            payload = axis_credit["profiles"].get(row["profile_id"], {})
            row["axis_anchored_negative_count"] = payload.get("axis_anchored_negative_count", 0)
            row["axis_anchored_challenge_count"] = payload.get("axis_anchored_challenge_count", 0)
            row["axis_anchor_basis_status"] = payload.get("axis_anchor_basis_status", "not_audited")
            row["axis_anchored_negative_credit_ids"] = ";".join(payload.get("axis_anchored_negative_credit_ids", []))
            row["axis_anchored_challenge_credit_ids"] = ";".join(payload.get("axis_anchored_challenge_credit_ids", []))

    output_dir.mkdir(parents=True, exist_ok=True)
    panel_out = output_dir / "panel_binding.tsv"
    readiness_out = output_dir / "profile_calibration_readiness.tsv"
    report_out = output_dir / "panel_reconciliation.json"
    for path in (panel_out, readiness_out, report_out):
        if path.exists():
            raise FileExistsError(path)
    with panel_out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(records)
    readiness_fields = list(readiness[0]) if readiness else ["profile_id"]
    with readiness_out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=readiness_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(readiness)
    binding_counts = Counter(row["family_binding_status"] for row in records)
    role_counts = Counter(row["panel_role"] for row in records)
    decision_counts = Counter(row["calibration_decision"] for row in readiness)
    report = {
        "schema_version": 1,
        "status": "completed_candidate_only",
        "inputs": {name: {"path": str(path.resolve()), "size": path.stat().st_size, "sha256": _sha256(path)} for name, path in (("profiles", profile_path), ("ledger", ledger_path), ("panel", panel_path))},
        "outputs": {name: {"path": str(path.resolve()), "size": path.stat().st_size, "sha256": _sha256(path)} for name, path in (("panel_binding", panel_out), ("profile_calibration_readiness", readiness_out))},
        "profile_count": len(reference_only),
        "panel_record_count": len(records),
        "binding_status_counts": dict(sorted(binding_counts.items())),
        "panel_role_counts": dict(sorted(role_counts.items())),
        "calibration_decision_counts": dict(sorted(decision_counts.items())),
        "summary": {
            "new_family_call_created": False,
            "reference_only_profiles_calibrated": 0,
            "external_records_bound_to_reference_only": sum(bool(row["reference_only_profile_id"]) for row in records),
            "external_records_unresolved_or_existing_trained_family": sum(not bool(row["reference_only_profile_id"]) for row in records),
            "calibrated_candidate_model": decision_counts.get("calibrated_candidate_model", 0),
            "reference_only_insufficient_panel": decision_counts.get("reference_only_insufficient_panel", 0),
        },
        "boundaries": [
            "Unresolved external accessions remain evidence records and are not assigned to a reference-only family.",
            "Annotation-only and taxonomy-only records are not formal negatives.",
            "No reference-only profile is calibration-eligible without independent positives, family-resolved negatives, challenge controls, leakage checks, and held-out recovery.",
            "All outputs are candidate-only sequence evidence and do not prove PHB/PHA degradation phenotype.",
        ],
    }
    if axis_credit is not None:
        report["axis_anchored_credit"] = {
            "convention": "axis_anchored_credit_v1",
            "legacy_counters_unchanged": True,
            "legacy_fields": [
                "existing_experimental_positive", "external_bound_positive",
                "external_bound_negative", "external_bound_challenge", "heldout_positive",
                "family_resolved_negative_status", "challenge_status",
            ],
            "source_files": {
                "axis_panel": {"path": str(Path(axis_panel).resolve()), "size": Path(axis_panel).stat().st_size, "sha256": _sha256(Path(axis_panel))},
                "axis_definitions": {"path": str(Path(axis_definitions).resolve()), "size": Path(axis_definitions).stat().st_size, "sha256": _sha256(Path(axis_definitions))},
            },
            "counts": axis_credit["counts"],
            "profiles_with_axis_credit": sum(
                bool(payload["axis_anchored_negative_count"] or payload["axis_anchored_challenge_count"])
                for payload in axis_credit["profiles"].values()
            ),
            "boundaries": [
                "Axis-anchored credit is an additional, separately named counting channel; it never feeds the family-bound counters.",
                "A record only counts when axis, anchor_justification and source DOI/PMID are all present and its evidence kind is pre-registered for that axis end.",
                "A record only credits a profile that holds a family-resolved anchor on the opposite end of the same axis.",
                "Axis-anchored credit is candidate-only counting evidence and does not prove PHB/PHA degradation phenotype.",
            ],
        }
    report_out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return {"records": records, "readiness": readiness, "summary": report["summary"], "report": report, "axis_credit": axis_credit}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", type=Path, required=True)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-reference-only-count", type=int, default=37)
    parser.add_argument("--axis-panel", type=Path, default=None, help="optional discrimination-axis panel for axis_anchored credit")
    parser.add_argument("--axis-definitions", type=Path, default=None, help="optional axis_reachability.json carrying the frozen axis definitions")
    args = parser.parse_args()
    result = reconcile(args.profiles, args.ledger, args.panel, args.output_dir, expected_reference_only_count=args.expected_reference_only_count, axis_panel=args.axis_panel, axis_definitions=args.axis_definitions)
    print(json.dumps(result["report"]["calibration_decision_counts"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
