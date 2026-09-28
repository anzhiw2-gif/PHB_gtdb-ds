#!/usr/bin/env python3
"""Score the P4 structural survey tranche and test whether the sequence margin predicts it.

Reads the Foldseek comparison of the survey's predicted folds against the same
panel the pilot used (the experimental 8YNV structure plus the three measured
competitor predictions) and answers the question the pilot raised:

    does the sequence-level panel margin predict the STRUCTURAL margin?

The pilot found that both extremes of the sequence margin point at the confounder,
which - if general - means the sequence-level ``nphamcl_like_supported`` label
carries no structural warranty anywhere.  This module measures that directly,
with Spearman's rank correlation over a margin-stratified sample and a plain count
of which side wins.

The plan's pre-registered structural criterion is TM >= 0.5 with a supporting
E-value; it is applied here per comparison, and the *margin* between the anchor and
the competitor is reported alongside it, because the pilot showed both can clear
0.5 while the competitor still wins.

Candidate-only boundary: scores a review. Deletes, demotes and promotes nothing,
and a predicted fold is not a phenotype claim.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

ANCHOR_PREFIX = "anchor_8YNV"
COMPETITOR_PREFIX = "competitor_"
TMSOURCE = "qtmscore"
PRE_REGISTERED_TM = 0.5


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"input is not a regular file: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"empty TSV: {path}")
    return rows


def read_foldseek(path: Path) -> list[dict[str, str]]:
    """Foldseek ``easy-search`` output: query, target, alntm, qtm, ttm, prob, evalue, lddt."""
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"Foldseek output is not a regular file: {path}")
    fields = ["query", "target", "alntmscore", "qtmscore", "ttmscore", "prob", "evalue", "lddt"]
    rows: list[dict[str, str]] = []
    with path.open(encoding="utf-8-sig", errors="strict") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < len(fields):
                raise ValueError(f"malformed Foldseek row: {line!r}")
            row = dict(zip(fields, parts[: len(fields)]))
            try:
                row["tm"] = float(row[TMSOURCE])
                row["evalue_f"] = float(row["evalue"])
            except ValueError as exc:
                raise ValueError(f"non-numeric Foldseek score in {line!r}") from exc
            rows.append(row)
    if not rows:
        raise ValueError(f"Foldseek output is empty: {path}")
    return rows


def rank(values: list[float]) -> list[float]:
    """Average ranks, so ties do not distort Spearman."""
    order = sorted(range(len(values)), key=lambda index: values[index])
    ranks = [0.0] * len(values)
    position = 0
    while position < len(order):
        end = position
        while end + 1 < len(order) and values[order[end + 1]] == values[order[position]]:
            end += 1
        average = (position + end) / 2.0 + 1.0
        for index in range(position, end + 1):
            ranks[order[index]] = average
        position = end + 1
    return ranks


def spearman(first: list[float], second: list[float]) -> float | None:
    if len(first) != len(second):
        raise ValueError("Spearman needs equal-length samples")
    if len(first) < 3:
        return None
    left, right = rank(first), rank(second)
    n = len(left)
    mean_left, mean_right = sum(left) / n, sum(right) / n
    covariance = sum((a - mean_left) * (b - mean_right) for a, b in zip(left, right))
    spread_left = sum((a - mean_left) ** 2 for a in left) ** 0.5
    spread_right = sum((b - mean_right) ** 2 for b in right) ** 0.5
    if spread_left == 0 or spread_right == 0:
        return None
    return covariance / (spread_left * spread_right)


def candidate_of(query: str) -> str:
    """Foldseek query name back to its accession.

    ColabFold names outputs ``<accession>_unrelaxed_rank_00N_...pdb``; the
    accession itself contains ``|`` and ``.`` so it must not be re-split on those.
    """
    marker = "_unrelaxed_rank_"
    index = query.find(marker)
    if index < 0:
        raise ValueError(f"unrecognised prediction filename: {query!r}")
    return query[:index]


def score(
    foldseek_rows: list[dict[str, str]],
    composition_rows: list[dict[str, str]],
) -> tuple[list[dict[str, str]], dict]:
    margins: dict[str, str] = {}
    for row in composition_rows:
        accession = (row.get("accession") or "").strip()
        if not accession:
            raise ValueError("composition row without an accession")
        margins[accession] = (row.get("panel_margin_competitor_minus_anchor_bits") or "").strip()

    by_candidate: dict[str, list[dict[str, str]]] = {}
    for row in foldseek_rows:
        by_candidate.setdefault(candidate_of(row["query"]), []).append(row)

    records: list[dict[str, str]] = []
    for accession in sorted(by_candidate):
        rows = by_candidate[accession]
        anchor_rows = [row for row in rows if row["target"].startswith(ANCHOR_PREFIX)]
        competitor_rows = [row for row in rows if row["target"].startswith(COMPETITOR_PREFIX)]
        if not anchor_rows:
            raise ValueError(f"{accession} has no comparison against the anchor")
        if not competitor_rows:
            raise ValueError(f"{accession} has no comparison against any competitor")
        best_anchor = max(anchor_rows, key=lambda row: row["tm"])
        best_competitor = max(competitor_rows, key=lambda row: row["tm"])
        structural_margin = best_anchor["tm"] - best_competitor["tm"]
        sequence_margin = margins.get(accession)
        if sequence_margin is None:
            raise ValueError(f"{accession} is not in the survey composition table")
        records.append({
            "accession": accession,
            "anchor_target": best_anchor["target"],
            "anchor_tm": f"{best_anchor['tm']:.4f}",
            "anchor_evalue": best_anchor["evalue"],
            "competitor_target": best_competitor["target"],
            "competitor_tm": f"{best_competitor['tm']:.4f}",
            "competitor_evalue": best_competitor["evalue"],
            "structural_margin_anchor_minus_competitor": f"{structural_margin:.4f}",
            "sequence_margin_competitor_minus_anchor_bits": sequence_margin,
            "anchor_clears_preregistered_tm": "yes" if best_anchor["tm"] >= PRE_REGISTERED_TM else "no",
            "structural_winner": "anchor" if structural_margin > 0 else "competitor",
        })

    sequence_values: list[float] = []
    structural_values: list[float] = []
    for record in records:
        try:
            sequence_values.append(float(record["sequence_margin_competitor_minus_anchor_bits"]))
            structural_values.append(
                -float(record["structural_margin_anchor_minus_competitor"])
            )
        except ValueError:
            continue
    anchor_wins = sum(1 for record in records if record["structural_winner"] == "anchor")
    competitor_wins = len(records) - anchor_wins
    correlation = spearman(sequence_values, structural_values)
    summary = {
        "candidates_scored": len(records),
        "structural_winner_counts": {"anchor": anchor_wins, "competitor": competitor_wins},
        "share_won_by_the_competitor": competitor_wins / len(records) if records else None,
        "anchor_clears_preregistered_tm": sum(
            1 for record in records if record["anchor_clears_preregistered_tm"] == "yes"
        ),
        "preregistered_tm_threshold": PRE_REGISTERED_TM,
        "spearman_sequence_margin_vs_structural_margin": correlation,
        "pairs_used_for_correlation": len(sequence_values),
        "interpretation": (
            "The sequence margin is 'competitor minus anchor bits', the structural margin is "
            "'anchor minus competitor TM', so a POSITIVE Spearman means the sequence-level ranking "
            "predicts the structural one. A correlation near zero means the sequence label carries "
            "no structural warranty, which is what the two-candidate pilot suggested."
        ),
        "boundary": (
            "Predicted folds, not experimental structures; a tranche of the stratified sample, not "
            "a census. Scores a review and changes no disposition."
        ),
    }
    return records, summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--foldseek", required=True, type=Path)
    parser.add_argument("--composition", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)

    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise ValueError(f"refusing to write into a non-empty output dir: {args.out_dir}")
    records, summary = score(read_foldseek(args.foldseek), read_tsv(args.composition))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    fields = list(records[0].keys()) if records else ["accession"]
    with (args.out_dir / "survey_structural_scores.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    (args.out_dir / "survey_structural_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n",
    )
    print(json.dumps({
        "candidates_scored": summary["candidates_scored"],
        "structural_winner_counts": summary["structural_winner_counts"],
        "share_won_by_the_competitor": summary["share_won_by_the_competitor"],
        "spearman": summary["spearman_sequence_margin_vs_structural_margin"],
        "out_dir": str(args.out_dir),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
