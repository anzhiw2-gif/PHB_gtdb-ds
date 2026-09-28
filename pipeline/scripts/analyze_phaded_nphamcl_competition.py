#!/usr/bin/env python3
"""Score the nPHAMCL competition panel against the 987 candidates.

Reads a BLAST tabular result (``-outfmt 6 qseqid sseqid pident length evalue
bitscore``) of the 987 candidates against the competition panel built by
``build_phaded_nphamcl_competition_panel.py`` and answers the P4 question:

    for each candidate, is its best hit the nPHAMCL positive anchor (8YNV_A) or
    one of the measured hard competitors?

A candidate whose best hit is a competitor cannot be called nPHAMCL-like on the
evidence available: the competitor fold explains the sequence at least as well.
Such a candidate belongs to ``function_unresolved`` -- it is NOT deleted,
demoted or excluded.

Candidate-only boundary: this scores a review.  It changes no disposition by
itself, and a panel score is not a phenotype claim.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

DISPOSITION_SUPPORTED = "nphamcl_like_supported"
DISPOSITION_COMPETED = "competed_by_measured_confounder"
DISPOSITION_UNRESOLVED = "function_unresolved_no_anchor_hit"
DISPOSITION_MISSING = "no_hit_to_panel"

ANCHOR = "8YNV_A"


def read_blast(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"BLAST output is not a regular file: {path}")
    fields = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"]
    rows: list[dict[str, str]] = []
    with path.open(encoding="utf-8-sig", errors="strict") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < len(fields):
                raise ValueError(f"malformed BLAST row: {line!r}")
            row = dict(zip(fields, parts[: len(fields)]))
            try:
                row["evalue_f"] = float(row["evalue"])
                row["bitscore_f"] = float(row["bitscore"])
            except ValueError as exc:
                raise ValueError(f"non-numeric BLAST score in {line!r}") from exc
            rows.append(row)
    return rows


def classify(
    blast_rows: list[dict[str, str]],
    *,
    anchor: str = ANCHOR,
    competitor_ids: list[str] | None = None,
    control_ids: list[str] | None = None,
    supported_label: str = DISPOSITION_SUPPORTED,
    competed_label: str = DISPOSITION_COMPETED,
    unresolved_label: str = DISPOSITION_UNRESOLVED,
) -> tuple[list[dict[str, str]], dict]:
    """Per-candidate best hit against the anchor and against the competitors.

    The three labels are arguments so the same panel logic can score a different
    candidate pool (the with-lipase deferred pool, for instance) without pretending
    its members are nPHAMCL candidates.
    """
    competitor_ids = competitor_ids or []
    control_ids = control_ids or []
    by_query: dict[str, list[dict[str, str]]] = {}
    for row in blast_rows:
        by_query.setdefault(row["qseqid"], []).append(row)

    records: list[dict[str, str]] = []
    for query, rows in sorted(by_query.items()):
        anchor_rows = [r for r in rows if r["sseqid"] == anchor]
        competitor_rows = [r for r in rows if r["sseqid"] in competitor_ids]
        control_rows = [r for r in rows if r["sseqid"] in control_ids]
        best_anchor = min(anchor_rows, key=lambda r: r["evalue_f"]) if anchor_rows else None
        best_competitor = (
            min(competitor_rows, key=lambda r: r["evalue_f"]) if competitor_rows else None
        )
        best_control = (
            min(control_rows, key=lambda r: r["evalue_f"]) if control_rows else None
        )
        best_overall = min(rows, key=lambda r: r["evalue_f"])
        if best_anchor is None and best_competitor is None:
            disposition = unresolved_label
        elif best_competitor is None:
            disposition = supported_label
        elif best_anchor is None:
            disposition = competed_label
        elif best_competitor["evalue_f"] < best_anchor["evalue_f"]:
            disposition = competed_label
        else:
            disposition = supported_label
        margin = ""
        if best_anchor is not None and best_competitor is not None:
            anchor_e = best_anchor["evalue_f"]
            competitor_e = best_competitor["evalue_f"]
            # A very strong hit underflows to E = 0.0, so the ratio cannot simply be
            # divided. The three cases are reported rather than special-cased into a
            # number that would look like data:
            #   both 0        -> indistinguishable at this precision
            #   anchor 0 only -> the anchor is infinitely better
            #   competitor 0  -> the competitor is infinitely better
            if anchor_e == 0.0 and competitor_e == 0.0:
                margin = "both_below_precision"
            elif anchor_e == 0.0:
                margin = "inf"
            elif competitor_e == 0.0:
                margin = "0"
            else:
                margin = f"{competitor_e / anchor_e:.6g}"
        records.append({
            "accession": query,
            "best_overall_subject": best_overall["sseqid"],
            "anchor_subject": anchor if best_anchor else "",
            "anchor_evalue": best_anchor["evalue"] if best_anchor else "pending",
            "anchor_bitscore": best_anchor["bitscore"] if best_anchor else "pending",
            "competitor_subject": best_competitor["sseqid"] if best_competitor else "",
            "competitor_evalue": best_competitor["evalue"] if best_competitor else "pending",
            "competitor_bitscore": best_competitor["bitscore"] if best_competitor else "pending",
            "control_best_evalue": best_control["evalue"] if best_control else "pending",
            "competitor_over_anchor_evalue_ratio": margin,
            "disposition": disposition,
        })
    summary = summarize(records, competitor_ids, control_ids,
                        competed_label=competed_label, supported_label=supported_label)
    summary["queries_seen"] = len(by_query)
    return records, summary


def summarize(
    records: list[dict[str, str]],
    competitor_ids: list[str],
    control_ids: list[str],
    *,
    competed_label: str = DISPOSITION_COMPETED,
    supported_label: str = DISPOSITION_SUPPORTED,
) -> dict:
    counts: dict[str, int] = {}
    for record in records:
        counts[record["disposition"]] = counts.get(record["disposition"], 0) + 1
    total = len(records)
    competed = counts.get(competed_label, 0)
    supported = counts.get(supported_label, 0)
    return {
        "candidates_scored": total,
        "disposition_counts": counts,
        "competed_share": (competed / total) if total else None,
        "supported_share": (supported / total) if total else None,
        "hard_competitors_used": list(competitor_ids),
        "specificity_controls_used": list(control_ids),
        "interpretation": (
            "A candidate whose best panel hit is a measured confounder cannot be called "
            "nPHAMCL-like on this evidence; it belongs to function_unresolved and is NOT "
            "deleted, demoted or excluded. A candidate nearer the anchor is nPHAMCL-like at "
            "this resolution only, which is sequence-level and does not settle the fold."
        ),
        "boundary": (
            "Sequence-level competition only. The structural stage (Foldseek/TM-score "
            "against the same panel) requires predicted or experimental structures for the "
            "candidates and is a separate, heavier run."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blast", required=True, type=Path)
    parser.add_argument("--competitors", required=True, nargs="+")
    parser.add_argument("--controls", required=True, nargs="+")
    parser.add_argument("--anchor", default=ANCHOR)
    parser.add_argument("--supported-label", default=DISPOSITION_SUPPORTED,
                        help="disposition name for a candidate whose best panel hit is the anchor; "
                             "override it when scoring a different candidate pool")
    parser.add_argument("--competed-label", default=DISPOSITION_COMPETED,
                        help="disposition name for a candidate whose best panel hit is a competitor")
    parser.add_argument("--unresolved-label", default=DISPOSITION_UNRESOLVED)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)

    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise ValueError(f"refusing to write into a non-empty output dir: {args.out_dir}")
    records, summary = classify(
        read_blast(args.blast),
        anchor=args.anchor,
        competitor_ids=args.competitors,
        control_ids=args.controls,
        supported_label=args.supported_label,
        competed_label=args.competed_label,
        unresolved_label=args.unresolved_label,
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    fields = list(records[0].keys()) if records else ["accession"]
    with (args.out_dir / "nphamcl_competition_scores.tsv").open(
        "w", encoding="utf-8", newline="",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    (args.out_dir / "nphamcl_competition_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n",
    )
    print(json.dumps({
        "candidates_scored": summary["candidates_scored"],
        "disposition_counts": summary["disposition_counts"],
        "competed_share": summary["competed_share"],
        "out_dir": str(args.out_dir),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
