#!/usr/bin/env python3
"""Candidate-only joint classifier for archaeal PhaZh1-like patatin hits."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def _bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def _markers(value: object) -> set[str]:
    return {item.strip() for item in str(value or "").replace(",", ";").split(";") if item.strip()}


def parse_domtblout(lines: list[str]) -> list[dict[str, object]]:
    """Parse standard HMMER domtblout lines and calculate target coverage."""
    parsed = []
    for raw in lines:
        if not raw or raw.startswith("#"):
            continue
        fields = raw.split()
        if len(fields) < 22:
            continue
        try:
            target_length = int(fields[2])
            ali_from, ali_to = int(fields[17]), int(fields[18])
            coverage = round((ali_to - ali_from + 1) / target_length, 6)
            parsed.append({"target": fields[0], "target_length": target_length, "evalue": float(fields[6]), "bitscore": float(fields[7]), "coverage": coverage})
        except (ValueError, ZeroDivisionError):
            continue
    return parsed


def audit_fasta(domtblout: Path, fasta: Path, context_tsv: Path, output_tsv: Path, summary_tsv: Path) -> Counter:
    """Create four-layer archaeal audit rows from probe hits and neighborhood context."""
    with context_tsv.open(newline="", encoding="utf-8") as handle:
        context = {row.get("hit_locus", ""): row for row in csv.DictReader(handle, delimiter="\t")}
    hits = parse_domtblout(domtblout.read_text(encoding="utf-8", errors="replace").splitlines())
    rows = []
    for hit in hits:
        genome, _, locus = hit["target"].partition("|")
        seq_len = hit["target_length"]
        c = context.get(locus, {})
        markers = c.get("nearby_markers", "")
        rows.append({"genome": genome, "locus": locus, "domain": "Archaea", "coverage": hit["coverage"], "evalue": hit["evalue"], "complete": str(seq_len >= 200).lower(), "nearby_markers": markers, "patatin_conflict": str(bool(_markers(markers) & {"PhaC", "PhaE", "phasin"})).lower()})
    input_tsv = output_tsv.with_suffix(".input.tsv")
    with input_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["genome", "locus", "domain", "coverage", "evalue", "complete", "nearby_markers", "patatin_conflict"], delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    return classify_table(input_tsv, output_tsv, summary_tsv)


def classify_record(row: dict[str, object]) -> str:
    """Return one conservative candidate-only layer for one hit."""
    if str(row.get("domain", "")).strip().lower() not in {"archaea", "d__archaea"}:
        return "non_PhaZh1_patatin"
    if _bool(row.get("patatin_conflict", False)):
        return "non_PhaZh1_patatin"
    try:
        coverage = float(row.get("coverage", 0.0))
        evalue = float(row.get("evalue", "inf"))
    except (TypeError, ValueError):
        return "archaeal_patatin_exploratory"
    if evalue > 1e-5 or coverage < 0.6:
        return "archaeal_patatin_exploratory"
    markers = _markers(row.get("nearby_markers", ""))
    mobilization = bool(markers & {"BdhA", "PhaJ"})
    contextual = bool(markers & {"PhaC", "PhaE", "phasin"})
    if coverage >= 0.8 and _bool(row.get("complete", False)) and mobilization:
        return "PhaZh1_like_high"
    if _bool(row.get("complete", False)) and (mobilization or contextual):
        return "PhaZh1_like_review"
    return "archaeal_patatin_exploratory"


def classify_table(input_tsv: Path, output_tsv: Path, summary_tsv: Path) -> Counter:
    with input_tsv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"domain", "coverage", "evalue", "complete", "nearby_markers", "patatin_conflict"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"missing classifier columns: {', '.join(sorted(missing))}")
        rows = list(reader)
    counts: Counter = Counter()
    for row in rows:
        row["layer"] = classify_record(row)
        counts[row["layer"]] += 1
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("w", newline="", encoding="utf-8") as handle:
        fields = list(rows[0].keys()) if rows else [*sorted(required), "layer"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with summary_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["layer", "records"])
        writer.writerows(sorted(counts.items()))
    return counts


def build_input_from_tier(tier_tsv: Path, context_tsv: Path, output_tsv: Path) -> None:
    """Join archaeal patatin tier records to neighborhood evidence conservatively."""
    with tier_tsv.open(newline="", encoding="utf-8") as handle:
        tier = list(csv.DictReader(handle, delimiter="\t"))
    with context_tsv.open(newline="", encoding="utf-8") as handle:
        context = list(csv.DictReader(handle, delimiter="\t"))
    by_locus = {row.get("hit_locus", ""): row for row in context}
    rows = []
    for row in tier:
        locus = row.get("locus", "")
        evidence = by_locus.get(locus, {})
        taxonomy = row.get("gtdb_taxonomy", "")
        rows.append({
            "genome": row.get("genome", ""),
            "locus": locus,
            "domain": "Archaea" if taxonomy.startswith("d__Archaea;") else "Bacteria",
            "coverage": row.get("coverage", "0"),
            "evalue": row.get("evalue", "inf"),
            "complete": row.get("complete", "false"),
            "nearby_markers": evidence.get("nearby_markers", ""),
            "patatin_conflict": evidence.get("patatin_conflict", "false"),
        })
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["genome", "locus", "domain", "coverage", "evalue", "complete", "nearby_markers", "patatin_conflict"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    classify_table(args.input, args.output, args.summary)


if __name__ == "__main__":
    main()
