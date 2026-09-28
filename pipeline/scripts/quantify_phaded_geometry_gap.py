#!/usr/bin/env python3
"""Quantify the geometry gap F14 flagged, without re-running its chain.

F14's adapter joins the v1 merge and the demotion addendum for the catalytic
domain verdict, and reported 668 accessions whose geometry stayed ``unresolved``
because neither table covers them.  The frozen catalytic-domain-type table
(109,087 rows) DOES cover them; it simply was not in the declared input set.

This module answers the question that decides whether a re-run is worth it:
**what would those rows become?**  It changes nothing itself.

Candidate-only boundary: reading two frozen tables and counting. No disposition
changes, and the catalog is not rewritten.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

UNRESOLVED = "unresolved"
ALLOWED = ("type1_verified", "type2_verified", "other", UNRESOLVED)


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"input is not a regular file: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"empty TSV: {path}")
    return rows


def normalize_verdict(value: str) -> str:
    """Map a frozen geometry verdict onto the schema's type vocabulary.

    Unrecognised verdicts (for example ``undetermined_no_oxyanion``) stay
    ``unresolved``: the frozen layer uses them to mean "geometry could not
    decide", which is not a verified type.
    """
    text = (value or "").strip()
    if text in ("type1_verified", "type2_verified"):
        return text
    return UNRESOLVED


def quantify(
    adapted_rows: list[dict[str, str]],
    geometry_rows: list[dict[str, str]],
) -> tuple[list[dict[str, str]], dict]:
    geometry = {
        (row.get("accession") or "").strip(): (row.get("catalytic_domain_type") or "").strip()
        for row in geometry_rows
        if (row.get("accession") or "").strip()
    }
    if len(geometry) != len(geometry_rows):
        raise ValueError("geometry table holds duplicate accessions")

    gap_rows: list[dict[str, str]] = []
    closed = Counter()
    still_unresolved = Counter()
    for row in adapted_rows:
        accession = (row.get("accession") or "").strip()
        if not accession:
            raise ValueError("adapted row without an accession")
        current = (row.get("catalytic_domain_type") or "").strip()
        if current != UNRESOLVED:
            continue
        verdict = geometry.get(accession)
        if verdict is None:
            gap_rows.append({
                "accession": accession,
                "current": current,
                "geometry_verdict": "absent_from_the_geometry_table",
                "would_become": UNRESOLVED,
                "closes_the_gap": "no",
            })
            still_unresolved["absent_from_the_geometry_table"] += 1
            continue
        would_be = normalize_verdict(verdict)
        closes = would_be != UNRESOLVED
        gap_rows.append({
            "accession": accession,
            "current": current,
            "geometry_verdict": verdict,
            "would_become": would_be,
            "closes_the_gap": "yes" if closes else "no",
        })
        if closes:
            closed[would_be] += 1
        else:
            still_unresolved[verdict or "<empty>"] += 1

    summary = {
        "adapted_rows": len(adapted_rows),
        "unresolved_before": len(gap_rows),
        "closed_by_the_geometry_table": sum(closed.values()),
        "closed_breakdown": dict(closed),
        "still_unresolved": sum(still_unresolved.values()),
        "still_unresolved_breakdown": dict(still_unresolved),
        "closure_rate": (sum(closed.values()) / len(gap_rows)) if gap_rows else None,
        "interpretation": (
            "A re-run of the adapter with the geometry table added as a third source would move "
            "these rows out of 'unresolved'. Unrecognised frozen verdicts stay unresolved on "
            "purpose: they mean the geometry could not decide, which is not a verified type."
        ),
        "boundary": (
            "Quantification only: no catalog is rewritten, no disposition changes, and the "
            "deferred layer is untouched."
        ),
    }
    return gap_rows, summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--adapted-input", required=True, type=Path,
                        help="F14's results/phaded_v2_catalog_input.tsv")
    parser.add_argument("--geometry", required=True, type=Path,
                        help="frozen runs/20260917_phaded_catalytic_domain_type_01/results/"
                             "catalytic_domain_type.tsv")
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)

    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise ValueError(f"refusing to write into a non-empty output dir: {args.out_dir}")
    adapted = read_tsv(args.adapted_input)
    geometry = read_tsv(args.geometry)
    rows, summary = quantify(adapted, geometry)
    summary["geometry_rows"] = len(geometry)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    fields = ["accession", "current", "geometry_verdict", "would_become", "closes_the_gap"]
    with (args.out_dir / "geometry_gap_rows.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (args.out_dir / "geometry_gap_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n",
    )
    print(json.dumps({
        "unresolved_before": summary["unresolved_before"],
        "closed_by_the_geometry_table": summary["closed_by_the_geometry_table"],
        "closed_breakdown": summary["closed_breakdown"],
        "still_unresolved": summary["still_unresolved"],
        "out_dir": str(args.out_dir),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
