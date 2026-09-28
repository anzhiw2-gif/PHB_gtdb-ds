#!/usr/bin/env python3
"""Merge full-library SignalP 6 predictions into the PhaDED feature table."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Iterable, Mapping


SIGNAL_TYPES = {"SP", "LIPO", "TAT", "TATLIPO", "PILIN", "OTHER"}


def parse_signalp_results(
    path: Path,
    expected_accessions: Iterable[str],
    excluded_accessions: set[str],
) -> dict[str, dict[str, str]]:
    expected = [value.strip() for value in expected_accessions if value.strip()]
    if len(expected) != len(set(expected)):
        raise ValueError("expected accession list contains duplicates or empty values")
    expected_set = set(expected)
    predictions: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8-sig") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            accession, signal_class = fields[0].strip(), fields[1].strip()
            if accession not in expected_set:
                raise ValueError(f"SignalP accession outside expected set: {accession}")
            if accession in predictions:
                raise ValueError(f"duplicate SignalP accession: {accession}")
            if signal_class not in SIGNAL_TYPES:
                raise ValueError(f"unknown SignalP prediction class: {signal_class}")
            predictions[accession] = {
                "signalp_class": signal_class,
                "signalp_status": "signalp_supported",
            }
    if not excluded_accessions.issubset(expected_set):
        raise ValueError("SignalP exclusion contains accession outside expected set")
    for accession in expected:
        if accession in excluded_accessions:
            predictions[accession] = {
                "signalp_class": "pending",
                "signalp_status": "signalp_tool_input_excluded",
            }
        elif accession not in predictions:
            predictions[accession] = {
                "signalp_class": "pending",
                "signalp_status": "signalp_not_reported",
            }
    return predictions


def merge_signalp_into_features(
    feature_rows: list[Mapping[str, str]],
    signalp_rows: Mapping[str, Mapping[str, str]],
) -> list[dict[str, str]]:
    accessions = [str(row.get("accession", "")).strip() for row in feature_rows]
    if any(not accession for accession in accessions) or len(accessions) != len(set(accessions)):
        raise ValueError("feature table contains empty or duplicate accessions")
    if set(signalp_rows) != set(accessions):
        raise ValueError("SignalP and feature accession sets differ")
    merged: list[dict[str, str]] = []
    for row in feature_rows:
        accession = str(row["accession"]).strip()
        item = dict(row)
        item["signalp_class"] = signalp_rows[accession]["signalp_class"]
        item["signalp_evidence"] = signalp_rows[accession]["signalp_status"]
        merged.append(item)
    return merged


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--signalp", type=Path, required=True)
    parser.add_argument("--excluded", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    with args.features.open(encoding="utf-8-sig", newline="") as handle:
        features = list(csv.DictReader(handle, delimiter="\t"))
    with args.excluded.open(encoding="utf-8-sig", newline="") as handle:
        excluded = {
            row.get("accession", "").strip()
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("accession", "").strip()
        }
    signalp = parse_signalp_results(args.signalp, [row.get("accession", "") for row in features], excluded)
    merged = merge_signalp_into_features(features, signalp)
    fields: list[str] = []
    for row in merged:
        for field in row:
            if field not in fields:
                fields.append(field)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(merged)
    summary = {
        "status": "completed_candidate_only",
        "candidate_count": len(merged),
        "signalp_status_counts": {},
        "inputs": {
            "features": {"path": str(args.features.resolve()), "sha256": _sha256(args.features)},
            "signalp": {"path": str(args.signalp.resolve()), "sha256": _sha256(args.signalp)},
            "excluded": {"path": str(args.excluded.resolve()), "sha256": _sha256(args.excluded)},
        },
        "output": {"path": str(args.output.resolve()), "sha256": _sha256(args.output)},
        "phenotype_boundary": "SignalP is localization evidence only and does not prove PHB/PHA degradation phenotype.",
    }
    for row in merged:
        status = row["signalp_evidence"]
        summary["signalp_status_counts"][status] = summary["signalp_status_counts"].get(status, 0) + 1
    (args.output.parent / "signalp_merge_manifest.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
