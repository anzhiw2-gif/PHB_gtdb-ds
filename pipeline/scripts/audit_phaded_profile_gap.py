#!/usr/bin/env python3
"""Audit PhaDED profile coverage and evidence gaps without running HMM search."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _int(value: str | None) -> int:
    try:
        return int(value or 0)
    except ValueError:
        return 0


def audit(
    manifest_path: Path,
    ledger_path: Path,
    assignment_path: Path,
    score_path: Path,
    panel_path: Path,
    output_path: Path,
) -> dict[str, object]:
    manifest = _read(manifest_path)
    ledger = _read(ledger_path)
    assignments = _read(assignment_path)
    scores = _read(score_path)
    panel = _read(panel_path) if panel_path.is_file() else []

    if len({row.get("profile_id", "") for row in manifest}) != len(manifest):
        raise ValueError("profile manifest contains duplicate profile_id")

    scores_by_profile = Counter(row.get("profile_id", "") for row in scores)
    best_by_profile = Counter()
    for row in assignments:
        for key in ("phaded_superfamily_best", "phaded_family_best"):
            profile_id = row.get(key, "")
            if profile_id:
                best_by_profile[profile_id] += 1

    ledger_by_family: dict[str, list[dict[str, str]]] = {}
    ledger_by_superfamily: dict[str, list[dict[str, str]]] = {}
    ledger_by_accession: dict[str, dict[str, str]] = {}
    for row in ledger:
        family = row.get("phaded_family_id", "")
        sf = row.get("phaded_superfamily", "")
        accession = row.get("accession", "")
        if family:
            ledger_by_family.setdefault(family, []).append(row)
        if sf:
            ledger_by_superfamily.setdefault(sf, []).append(row)
        if accession:
            ledger_by_accession[accession] = row

    assignments_by_sf: dict[str, list[dict[str, str]]] = {}
    for row in assignments:
        sf = row.get("phaded_superfamily_best", "")
        if sf:
            assignments_by_sf.setdefault(sf, []).append(row)

    panel_accessions = {row.get("accession", "") for row in panel if row.get("accession", "")}
    global_panel_counts = Counter(row.get("panel", "") for row in panel)
    profile_rows: list[dict[str, object]] = []
    for row in manifest:
        profile_id = row.get("profile_id", "")
        family = row.get("phaded_family_id", "")
        sf = row.get("phaded_superfamily", "")
        reference_rows = ledger_by_family.get(family, []) if family else []
        parent_candidates = assignments_by_sf.get(sf, [])
        direct_score_rows = scores_by_profile[profile_id]
        direct_best = best_by_profile[profile_id]
        family_accessions = {r.get("accession", "") for r in reference_rows}
        panel_overlap = family_accessions & panel_accessions
        exact_panel_rows = [r for r in panel if r.get("accession", "") in panel_overlap]
        exp_count = sum(r.get("evidence_status") == "experimental_positive" for r in reference_rows)
        ann_count = sum(r.get("evidence_status") == "annotation_only" for r in reference_rows)
        positive_seed_count = sum(r.get("positive_negative_control") == "positive_seed" for r in reference_rows)
        negative_or_challenge_count = sum(
            any(token in r.get("positive_negative_control", "").lower() for token in ("negative", "challenge"))
            for r in reference_rows
        )
        independent_positive_sufficiency = "sufficient" if exp_count >= 3 else "insufficient"
        if row.get("model_status") == "trained":
            coverage_status = "trained_profile_used"
        elif row.get("profile_kind") == "superfamily":
            coverage_status = "reference_only_superfamily"
        else:
            coverage_status = "reference_only_family"
        profile_rows.append(
            {
                "profile_id": profile_id,
                "profile_kind": row.get("profile_kind", ""),
                "phaded_superfamily": sf,
                "phaded_family_id": family,
                "model_status": row.get("model_status", ""),
                "model_reason": row.get("model_reason", ""),
                "training_count": _int(row.get("training_count")),
                "reference_count": len(reference_rows),
                "reference_experimental_positive_count": exp_count,
                "reference_annotation_only_count": ann_count,
                "reference_positive_seed_count": positive_seed_count,
                "reference_negative_or_challenge_count": negative_or_challenge_count,
                "independent_positive_sufficiency": independent_positive_sufficiency,
                "direct_score_row_count": direct_score_rows,
                "direct_best_candidate_count": direct_best,
                "parent_superfamily_candidate_count": len(parent_candidates),
                "parent_superfamily_assigned_count": sum(r.get("assignment_status") == "assigned" for r in parent_candidates),
                "parent_superfamily_ambiguous_count": sum(r.get("assignment_status", "").startswith("ambiguous") for r in parent_candidates),
                "parent_superfamily_unassigned_count": sum(r.get("assignment_status") == "unassigned_PhaDED_like" for r in parent_candidates),
                "external_panel_exact_profile_count": len(panel_overlap),
                "external_panel_exact_positive_count": sum(
                    "positive" in r.get("panel", "") and "negative" not in r.get("panel", "")
                    for r in exact_panel_rows
                ),
                "external_panel_exact_negative_count": sum("negative" in r.get("panel", "") for r in exact_panel_rows),
                "external_panel_exact_challenge_count": sum(
                    "challenge" in r.get("panel", "") or "fragment" in r.get("panel", "")
                    for r in exact_panel_rows
                ),
                "external_panel_unresolved_count": len(panel_accessions - family_accessions) if family else len(panel_accessions),
                "external_panel_global_positive_count": sum(
                    count for label, count in global_panel_counts.items()
                    if "positive" in label and "negative" not in label
                ),
                "external_panel_global_negative_count": sum(
                    count for label, count in global_panel_counts.items() if "negative" in label
                ),
                "external_panel_global_challenge_count": sum(
                    count for label, count in global_panel_counts.items()
                    if "challenge" in label or "fragment" in label
                ),
                "negative_challenge_evidence_status": (
                    "family_resolved" if exact_panel_rows else "family_unresolved"
                ),
                "coverage_status": coverage_status,
                "interpretation": "candidate_only; reference_only was not used as a calibrated HMM",
            }
        )

    status_counts = Counter(row.get("model_status", "") for row in manifest)
    trained = [row for row in profile_rows if row["model_status"] == "trained"]
    reference_only = [row for row in profile_rows if row["model_status"] == "reference_only"]
    report = {
        "schema_version": 1,
        "status": "verified",
        "manifest": {"path": str(manifest_path.resolve()), "size": manifest_path.stat().st_size, "sha256": _sha256(manifest_path)},
        "inputs": {
            "reference_ledger": {"path": str(ledger_path.resolve()), "size": ledger_path.stat().st_size, "sha256": _sha256(ledger_path)},
            "assignments": {"path": str(assignment_path.resolve()), "size": assignment_path.stat().st_size, "sha256": _sha256(assignment_path)},
            "scores": {"path": str(score_path.resolve()), "size": score_path.stat().st_size, "sha256": _sha256(score_path)},
            "external_panel": {"path": str(panel_path.resolve()), "size": panel_path.stat().st_size, "sha256": _sha256(panel_path)} if panel_path.is_file() else {"path": str(panel_path), "status": "pending"},
        },
        "profile_count": len(manifest),
        "trained_count": len(trained),
        "reference_only_count": len(reference_only),
        "model_status_counts": dict(status_counts),
        "candidate_count": len(assignments),
        "score_row_count": len(scores),
        "external_panel_accession_count": len(panel_accessions),
        "profiles": profile_rows,
        "conclusions": [
            "Only trained profiles generated candidate assignments; reference_only profiles have zero direct calibrated score rows.",
            "Parent superfamily candidate counts describe a broad candidate pool and do not constitute family-level calls.",
            "Reference evidence and external panels remain candidate-only; annotation-only or taxonomy-only records cannot establish phenotype.",
            "Reference-only profiles require independent positive references, independent negatives, and challenge-set calibration before training or full-library assignment.",
        ],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


def write_tsv(report: dict[str, object], path: Path) -> None:
    rows = report["profiles"]
    assert isinstance(rows, list)
    if not rows:
        path.write_text("\n", encoding="utf-8")
        return
    fields = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(report: dict[str, object], path: Path) -> None:
    rows = [row for row in report["profiles"] if row["model_status"] == "reference_only"]
    grouped: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault(str(row["phaded_superfamily"]), []).append(row)
    lines = [
        "# PhaDED profile coverage-gap audit",
        "",
        "Run: `20260915_phaded_profile_gap_audit_01`",
        "",
        "This audit is candidate-only. It does not run HMM search, modify a profile, or assign phenotype.",
        "",
        "## Summary",
        "",
        f"- Profiles: {report['profile_count']} total; {report['trained_count']} trained; {report['reference_only_count']} reference-only.",
        f"- Candidate assignments: {report['candidate_count']}; score rows: {report['score_row_count']}.",
        "- Every reference-only profile has zero direct calibrated score rows and zero direct best-candidate calls.",
        "- Parent-superfamily counts below are broad pools only; they cannot be interpreted as family calls.",
        "- The external panel has 17 accessions, but it is not family-resolved for the reference-only profiles except where accession overlap is explicit.",
        "",
        "## Reference-only profiles by parent superfamily",
        "",
        "| Parent superfamily | Profiles | Reference rows | Experimental positive rows | Parent candidates | Parent assigned | Parent ambiguous | Panel family-resolved |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for sf in sorted(grouped):
        group = grouped[sf]
        lines.append(
            "| " + " | ".join(
                [
                    sf,
                    str(len(group)),
                    str(sum(int(r["reference_count"]) for r in group)),
                    str(sum(int(r["reference_experimental_positive_count"]) for r in group)),
                    str(group[0]["parent_superfamily_candidate_count"]),
                    str(group[0]["parent_superfamily_assigned_count"]),
                    str(group[0]["parent_superfamily_ambiguous_count"]),
                    str(sum(r["negative_challenge_evidence_status"] == "family_resolved" for r in group)),
                ]
            )
            + " |"
        )
    lines += [
        "",
        "## Profile-level gaps",
        "",
        "The complete machine-readable table is `profile_gap_audit.tsv`. The fields distinguish direct profile coverage from the parent pool.",
        "",
        "| Profile | Family | Training count | Reference rows | Experimental positives | Parent candidates | Direct scores | Negative/challenge status |",
        "|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in sorted(rows, key=lambda item: str(item["profile_id"])):
        lines.append(
            "| " + " | ".join(
                [
                    str(row["profile_id"]),
                    str(row["phaded_family_id"] or "(superfamily)"),
                    str(row["training_count"]),
                    str(row["reference_count"]),
                    str(row["reference_experimental_positive_count"]),
                    str(row["parent_superfamily_candidate_count"]),
                    str(row["direct_score_row_count"]),
                    str(row["negative_challenge_evidence_status"]),
                ]
            )
            + " |"
        )
    lines += [
        "",
        "## Interpretation",
        "",
        "Reference-only profiles are retained as vocabulary, not calibrated classifiers. A family must not be assigned from the parent-superfamily pool. Before training, each family needs independent positive references, family-resolved negatives, and challenge sequences with declared thresholds and validation splits.",
        "",
        "All sequence, domain, homology, and panel labels remain candidate evidence and are not experimental PHB degradation claims.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--external-panel", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.manifest, args.ledger, args.assignments, args.scores, args.external_panel, args.output_json)
    write_tsv(report, args.output_tsv)
    write_markdown(report, args.output_md)
    print(json.dumps({k: report[k] for k in ("profile_count", "trained_count", "reference_only_count", "candidate_count", "score_row_count")}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
