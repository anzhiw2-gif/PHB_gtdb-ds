#!/usr/bin/env python3
"""Validate and score a separate four-seed ePhaZ subtype calibration panel."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path


DECISIONS = {"calibrated_candidate_model", "reference_only_insufficient_panel", "rejected_mixed_or_nonseparable"}
ALLOWED_EVIDENCE = {"experimental_positive", "experimental_negative", "challenge_control"}
REQUIRED = {
    "subtype", "accession", "role", "substrate_class", "localization", "phaded_superfamily",
    "genus", "sequence_sha256", "evidence_status", "threshold", "coverage_threshold",
    "heldout_hit", "negative_hit",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"panel manifest is not a regular file: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(REQUIRED - set(reader.fieldnames or []))
        if missing:
            raise ValueError("panel manifest missing columns: " + ",".join(missing))
        rows = list(reader)
    if not rows:
        raise ValueError("panel manifest is empty")
    return rows


def _expected_substrate(subtype: str) -> str | None:
    label = subtype.lower().replace("-", "_")
    if "dphamcl" in label or "mcl" in label:
        return "MCL"
    if "dphascl" in label or "phaz7" in label or "scl" in label:
        return "SCL"
    return None


def _truth(value: str) -> bool:
    return value.strip().lower() in {"true", "1", "yes", "hit"}


def calibrate_subtype_models(panel_manifest: str | Path, output_dir: str | Path) -> dict[str, object]:
    panel_manifest, output_dir = Path(panel_manifest), Path(output_dir)
    rows = _read(panel_manifest)
    seen_accessions: set[str] = set()
    train_ids: set[str] = set()
    heldout_ids: set[str] = set()
    by_subtype: dict[str, list[dict[str, str]]] = defaultdict(list)
    for index, row in enumerate(rows, start=2):
        subtype = row["subtype"].strip()
        accession = row["accession"].strip()
        if not subtype or not accession:
            raise ValueError(f"row {index}: subtype and accession are required")
        if accession in seen_accessions:
            roles = {existing["role"] for existing in rows if existing["accession"].strip() == accession}
            if "train_positive" in roles and "heldout_positive" in roles:
                raise ValueError(f"train-heldout leakage: {accession}")
            raise ValueError(f"duplicate panel accession: {accession}")
        seen_accessions.add(accession)
        role = row["role"].strip()
        evidence_status = row["evidence_status"].strip()
        if evidence_status not in ALLOWED_EVIDENCE:
            raise ValueError(f"invalid evidence_status: {accession}")
        if role in {"train_positive", "heldout_positive"} and evidence_status != "experimental_positive":
            raise ValueError(f"train/heldout evidence_status must be experimental_positive: {accession}")
        digest = row["sequence_sha256"].strip().lower()
        if len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest):
            raise ValueError(f"invalid sequence_sha256: {accession}")
        if role == "train_positive":
            train_ids.add(accession)
        elif role == "heldout_positive":
            heldout_ids.add(accession)
        by_subtype[subtype].append(row)
    overlap = train_ids & heldout_ids
    if overlap:
        raise ValueError("train-heldout leakage: " + ",".join(sorted(overlap)))
    decisions: dict[str, str] = {}
    metrics: list[dict[str, object]] = []
    for subtype, group in sorted(by_subtype.items()):
        substrates = {row["substrate_class"].strip().upper() for row in group}
        localizations = {row["localization"].strip().lower() for row in group}
        superfamilies = {row["phaded_superfamily"].strip() for row in group}
        expected = _expected_substrate(subtype)
        if expected and substrates != {expected}:
            raise ValueError(f"substrate-class mismatch: {subtype} expects {expected}, observed {sorted(substrates)}")
        if len(substrates) > 1 or len(localizations) > 1 or len(superfamilies) > 1:
            raise ValueError(f"mixed substrate/localization/PhaDED boundary: {subtype}")
        thresholds = {(row["threshold"].strip(), row["coverage_threshold"].strip()) for row in group}
        if len(thresholds) != 1:
            raise ValueError(f"threshold change within subtype: {subtype}")
        train = [row for row in group if row["role"] == "train_positive"]
        heldout = [row for row in group if row["role"] == "heldout_positive"]
        negatives = [row for row in group if row["role"] in {"formal_negative", "challenge_control"}]
        genera = {row["genus"].strip() for row in train + heldout if row["genus"].strip()}
        heldout_recovered = all(_truth(row["heldout_hit"]) for row in heldout) if heldout else False
        negative_hits = sum(_truth(row["negative_hit"]) for row in negatives)
        eligible = len(train) >= 3 and bool(heldout) and len(genera) >= 3 and heldout_recovered and negative_hits == 0
        decision = "calibrated_candidate_model" if eligible else "reference_only_insufficient_panel"
        decisions[subtype] = decision
        metrics.append({
            "subtype": subtype, "train_positive": len(train), "heldout_positive": len(heldout),
            "genera": len(genera), "heldout_recovered": str(heldout_recovered).lower(),
            "negative_hits": negative_hits, "threshold": next(iter(thresholds))[0],
            "coverage_threshold": next(iter(thresholds))[1], "decision": decision,
        })
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    metrics_path = output_dir / "leave_one_out_metrics.tsv"
    with metrics_path.open("w", encoding="utf-8", newline="") as handle:
        fields = list(metrics[0]) if metrics else ["subtype", "decision"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(metrics)
    decision_path = output_dir / "subtype_model_decision.json"
    report = {
        "status": "candidate-only", "decisions": decisions,
        "input_manifest": {"path": str(panel_manifest.resolve()), "sha256": _sha256(panel_manifest)},
        "metrics": str(metrics_path.resolve()), "formal_registry_modified": False,
        "historical_ephaz_model_overwritten": False,
    }
    decision_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    print(json.dumps(calibrate_subtype_models(args.panel_manifest, args.output_dir), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
