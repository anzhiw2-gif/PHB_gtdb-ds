#!/usr/bin/env python3
"""Prepare an accession-bound acquisition manifest for reference-only PhaDED profiles.

This stage only declares what must be acquired.  It does not retrieve sequences,
interpret annotations, train profiles, calibrate thresholds, or assign any new
family.  Every reference-only profile receives one row for each of the three
independent panel roles so that missing evidence remains explicit and auditable.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
from typing import Iterable


REQUIRED_PROFILE_COLUMNS = {
    "profile_id",
    "profile_kind",
    "phaded_superfamily",
    "phaded_family_id",
    "model_status",
}
ROLE_ORDER = (
    "independent_positive",
    "family_resolved_negative",
    "challenge_control",
)
OUTPUT_FIELDS = (
    "profile_id",
    "profile_kind",
    "family",
    "superfamily",
    "role",
    "accession",
    "source",
    "evidence_status",
    "sequence_sha256",
    "retrieval_status",
    "independence_status",
    "notes",
)


class PanelPreparationError(ValueError):
    """Raised when a profile manifest cannot be converted safely."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_profile_manifest(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.is_symlink():
        raise PanelPreparationError(f"profile manifest is not a regular file: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        missing = sorted(REQUIRED_PROFILE_COLUMNS - fields)
        if missing:
            raise PanelPreparationError("profile manifest missing columns: " + ",".join(missing))
        rows = list(reader)
    if not rows:
        raise PanelPreparationError("profile manifest is empty")

    seen: set[str] = set()
    for line_number, row in enumerate(rows, start=2):
        profile_id = row.get("profile_id", "").strip()
        if not profile_id:
            raise PanelPreparationError(f"line {line_number}: profile_id is empty")
        if profile_id in seen:
            raise PanelPreparationError(f"duplicate profile_id: {profile_id}")
        if row.get("profile_kind", "").strip() not in {"family", "superfamily"}:
            raise PanelPreparationError(f"{profile_id}: profile_kind must be family or superfamily")
        if not row.get("phaded_superfamily", "").strip():
            raise PanelPreparationError(f"{profile_id}: phaded_superfamily is empty")
        if row.get("profile_kind", "").strip() == "family" and not row.get("phaded_family_id", "").strip():
            raise PanelPreparationError(f"{profile_id}: family profile lacks phaded_family_id")
        seen.add(profile_id)
    return rows


def _role_source(role: str) -> str:
    return {
        "independent_positive": "UniProt reviewed; NCBI Protein; DED/PhaDED literature",
        "family_resolved_negative": "UniProt/NCBI adjacent-family controls with accession-bound family evidence",
        "challenge_control": "UniProt/NCBI structural-domain competitors and incomplete-sequence controls",
    }[role]


def _role_notes(role: str, profile_kind: str) -> str:
    if role == "independent_positive":
        return (
            "Acquire >=3 independent accessions with explicit experimental PHA/PHB result; "
            "deduplicate sequence and publication lineage before calibration."
        )
    if role == "family_resolved_negative":
        family_clause = (
            "Family identity must be accession-bound; do not infer it from parent-superfamily scores."
            if profile_kind == "family"
            else "Profile is superfamily-only; retain this as a control target and do not create a family call."
        )
        return (
            "Acquire same-superfamily, family-resolved non-target controls with explicit negative or "
            "non-matching function evidence. "
            + family_clause
        )
    return (
        "Acquire challenge controls spanning adjacent families, PhaC/PhaY/MaoC/generic esterases, "
        "domain-similar non-targets, and fragments; report challenge failures separately from biological negatives."
    )


def _build_rows(profiles: Iterable[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for profile in sorted(profiles, key=lambda item: item["profile_id"].strip()):
        profile_id = profile["profile_id"].strip()
        profile_kind = profile["profile_kind"].strip()
        family = profile.get("phaded_family_id", "").strip()
        superfamily = profile["phaded_superfamily"].strip()
        for role in ROLE_ORDER:
            rows.append(
                {
                    "profile_id": profile_id,
                    "profile_kind": profile_kind,
                    "family": family,
                    "superfamily": superfamily,
                    "role": role,
                    # Empty accession is intentional: acquisition has not run.
                    "accession": "",
                    "source": _role_source(role),
                    "evidence_status": "planned_not_run",
                    "sequence_sha256": "pending",
                    "retrieval_status": "pending",
                    "independence_status": "pending",
                    "notes": _role_notes(role, profile_kind),
                }
            )
    return rows


def _write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def prepare(
    profile_manifest: str | Path,
    output_dir: str | Path,
    *,
    expected_reference_only_count: int | None = 37,
) -> dict[str, object]:
    """Create a three-role acquisition manifest for reference-only profiles.

    ``expected_reference_only_count`` is configurable for fixture tests, but the
    CLI defaults to 37 so an accidental partial manifest fails closed in a real
    run.  No output file is overwritten.
    """
    profile_manifest = Path(profile_manifest)
    output_dir = Path(output_dir)
    all_profiles = _read_profile_manifest(profile_manifest)
    reference_only = [
        row for row in all_profiles if row.get("model_status", "").strip() == "reference_only"
    ]
    if expected_reference_only_count is not None and len(reference_only) != expected_reference_only_count:
        raise PanelPreparationError(
            f"expected {expected_reference_only_count} reference-only profiles, observed {len(reference_only)}"
        )
    rows = _build_rows(reference_only)
    output_dir.mkdir(parents=True, exist_ok=True)
    panel_path = output_dir / "reference_panel_acquisition.tsv"
    report_path = output_dir / "reference_panel_acquisition.json"
    for path in (panel_path, report_path):
        if path.exists():
            raise FileExistsError(f"refusing to overwrite existing output: {path}")
    _write_tsv(panel_path, rows)

    role_counts = {role: sum(row["role"] == role for row in rows) for role in ROLE_ORDER}
    report: dict[str, object] = {
        "schema_version": 1,
        "status": "planned_not_run",
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "profile_manifest": {
            "path": str(profile_manifest.resolve()),
            "size": profile_manifest.stat().st_size,
            "sha256": sha256_file(profile_manifest),
        },
        "reference_only_profile_count": len(reference_only),
        "record_count": len(rows),
        "role_counts": role_counts,
        "outputs": {
            "reference_panel_acquisition_tsv": {
                "path": str(panel_path.resolve()),
                "size": panel_path.stat().st_size,
                "sha256": sha256_file(panel_path),
            }
        },
        "calibration_boundary": (
            "Targets are acquisition placeholders only. Until independent positives, family-resolved "
            "negatives, challenge controls, leakage checks, and held-out calibration are complete, "
            "reference-only profiles remain untrained and cannot produce new family assignments."
        ),
    }
    # The report itself is intentionally not self-hashed: a self-referential
    # digest cannot be stable.  The outer run manifest records this JSON file.
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-reference-only-count", type=int, default=37)
    args = parser.parse_args()
    report = prepare(
        args.profile_manifest,
        args.output_dir,
        expected_reference_only_count=args.expected_reference_only_count,
    )
    print(
        json.dumps(
            {
                "status": report["status"],
                "reference_only_profile_count": report["reference_only_profile_count"],
                "record_count": report["record_count"],
                "role_counts": report["role_counts"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
