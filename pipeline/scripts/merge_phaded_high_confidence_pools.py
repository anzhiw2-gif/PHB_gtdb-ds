#!/usr/bin/env python3
"""Merge pool-internal and pool-external high-confidence candidate sets.

Adopts the combined high-confidence candidate set:

    pool_internal (36,559) + pool_external strong-tier (4,131) = 40,690.

The pool-internal source is the nucleophile-annotated high-confidence set
(``runs/20260919_phaded_nucleophile_unify_01/results/high_confidence_candidates.tsv``,
36,611 rows) MINUS the 52 SignalP-export Cys accessions demoted to hold
(``runs/20260920_phaded_three_gaps_01/results/cys_demotion/hold_candidates_addendum_52.tsv``).
Using this source (rather than the 15-column revised demotion TSV) preserves the
``nucleophile_type`` / ``nucleophile_family`` columns for the pool-internal rows.

Fail-closed invariants:
  - pool-internal and pool-external share the identical column set;
  - the two sets have zero overlapping accessions;
  - final counts are exactly 36,559 / 4,131 / 40,690.

Boundary: candidate-only. This merge introduces no new scoring, no registry
change, no family call and no deletion; it only unions two frozen outputs and
adds a ``pool_origin`` provenance column.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

POOL_IN = "pool_in"
POOL_EXTERNAL = "pool_external"
EXPECTED_POOL_IN = 36559
EXPECTED_POOL_EXTERNAL = 4131
EXPECTED_TOTAL = 40690


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or ()), list(reader)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def merge(pool_in_final: list[dict], pool_external: list[dict], columns: list[str]):
    """Union two high-confidence pools with a ``pool_origin`` provenance column.

    Returns ``(out_columns, merged_rows)``. Raises ``ValueError`` on any
    overlapping accession between the two pools or on a malformed column list.
    """
    if not columns or columns[0] != "accession":
        raise ValueError("columns must be non-empty and start with 'accession'")

    in_acc = {r["accession"] for r in pool_in_final}
    ext_acc = {r["accession"] for r in pool_external}
    overlap = in_acc & ext_acc
    if overlap:
        sample = sorted(overlap)[:5]
        raise ValueError(
            f"overlapping accessions between pools: {sample!r} ... ({len(overlap)} total)"
        )

    out_columns = ["accession", "pool_origin", *columns[1:]]
    origin_rank = {POOL_IN: 0, POOL_EXTERNAL: 1}
    merged = [dict(r, pool_origin=POOL_IN) for r in pool_in_final]
    merged += [dict(r, pool_origin=POOL_EXTERNAL) for r in pool_external]
    merged.sort(key=lambda r: (origin_rank[r["pool_origin"]], r["accession"]))
    return out_columns, merged


def merge_files(pool_in_hc: Path, pool_in_demote: Path, pool_external_hc: Path):
    """Read the three frozen inputs and return merged rows with invariants."""
    in_cols, in_rows = read_rows(pool_in_hc)
    dem_cols, dem_rows = read_rows(pool_in_demote)
    ext_cols, ext_rows = read_rows(pool_external_hc)

    if in_cols != ext_cols:
        raise ValueError(
            f"column mismatch:\n  pool_in : {in_cols}\n  pool_ext: {ext_cols}"
        )
    if "accession" not in dem_cols:
        raise ValueError(f"demote file missing accession column: {dem_cols}")

    demote = {r["accession"] for r in dem_rows}
    pool_in_final = [r for r in in_rows if r["accession"] not in demote]

    if len(pool_in_final) != EXPECTED_POOL_IN:
        raise ValueError(
            f"pool-in final rows {len(pool_in_final)} != expected {EXPECTED_POOL_IN}"
        )
    if len(ext_rows) != EXPECTED_POOL_EXTERNAL:
        raise ValueError(
            f"pool-external rows {len(ext_rows)} != expected {EXPECTED_POOL_EXTERNAL}"
        )

    out_cols, merged = merge(pool_in_final, ext_rows, in_cols)
    if len(merged) != EXPECTED_TOTAL:
        raise ValueError(f"merged total {len(merged)} != expected {EXPECTED_TOTAL}")

    stats = {
        "pool_in_hc_input": len(in_rows),
        "demoted": len(dem_rows),
        "pool_in_final": len(pool_in_final),
        "pool_external": len(ext_rows),
        "merged_total": len(merged),
        "overlap": 0,
    }
    return out_cols, merged, stats


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pool-in-hc", type=Path, required=True)
    parser.add_argument("--pool-in-demote", type=Path, required=True)
    parser.add_argument("--pool-external-hc", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()

    out_cols, merged, stats = merge_files(
        args.pool_in_hc, args.pool_in_demote, args.pool_external_hc
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=out_cols, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(merged)

    manifest = {
        "schema_version": "1.0",
        "run_id": args.run_id,
        "status": "completed_candidate_only",
        "purpose": "merge pool-internal and pool-external high-confidence candidate sets (40,690)",
        "counts": stats,
        "columns": out_cols,
        "inputs": {
            "pool_in_hc": {"path": str(args.pool_in_hc), "sha256": sha256(args.pool_in_hc)},
            "pool_in_demote": {"path": str(args.pool_in_demote), "sha256": sha256(args.pool_in_demote)},
            "pool_external_hc": {"path": str(args.pool_external_hc), "sha256": sha256(args.pool_external_hc)},
        },
        "output": {"path": str(args.output), "sha256": sha256(args.output)},
        "boundary": (
            "candidate-only; union of two frozen high-confidence outputs; no new "
            "scoring, no registry change, no family call, no deletion."
        ),
    }
    args.manifest.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
