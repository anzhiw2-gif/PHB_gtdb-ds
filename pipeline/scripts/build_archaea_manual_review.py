#!/usr/bin/env python3
"""Build a small cross-phylum manual neighborhood/domain review table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build(annotation: Path, context: Path, output: Path) -> int:
    context_by_key = {
        (row["genome"], row["hit_locus"]): row for row in read_tsv(context)
    }
    selected = []
    for row in read_tsv(annotation):
        if row["phylum"] != "Thermoproteota" and row["layer"] != "PhaZh1_like_review":
            continue
        context_row = context_by_key.get((row["genome"], row["locus"]), {})
        selected.append({
            "genome": row["genome"],
            "locus": row["locus"],
            "layer": row["layer"],
            "phylum": row["phylum"],
            "class": row["class"],
            "coverage": row["coverage"],
            "evalue": row["evalue"],
            "nterm_gxsxg_motif": row["nterm_gxsxg_motif"],
            "domain_model_hits": row["registered_model_hits"],
            "domain_best_bitscore": row["registered_model_best_bitscore"],
            "neighborhood_markers_annotation": row["nearby_markers"],
            "neighborhood_markers_context": context_row.get("nearby_markers", "pending"),
            "contig": context_row.get("contig", "pending"),
            "start": context_row.get("start", "pending"),
            "end": context_row.get("end", "pending"),
            "strand": context_row.get("strand", "pending"),
            "conflict_evidence": "none_registered_conflict_hit",
            "review_disposition": "manual_review_required",
            "phenotype_claim": "candidate_only",
        })
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(selected[0]) if selected else ["genome", "locus"]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(selected)
    return len(selected)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--context", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(build(args.annotation, args.context, args.output))


if __name__ == "__main__":
    main()
