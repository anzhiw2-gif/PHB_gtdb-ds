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

    ColabFold names outputs ``<accession>_unrelaxed_rank_00N_...pdb``, and it
    **replaces ``|`` with ``_``** when it writes those filenames. The accession
    itself contains ``|`` and ``.``, so it must not be re-split on those, and the
    recovered name is the *normalised* one - matching against the composition table
    goes through :func:`normalize_name`.
    """
    marker = "_unrelaxed_rank_"
    index = query.find(marker)
    if index < 0:
        raise ValueError(f"unrecognised prediction filename: {query!r}")
    return query[:index]


def normalize_name(accession: str) -> str:
    """The form ColabFold uses in filenames: ``|`` becomes ``_``."""
    return accession.replace("|", "_")


#: Two-tailed t at alpha = 0.05, by degrees of freedom. A lookup rather than a
#: single value, because the requirement grows sharply as n falls: at df=3 the
#: critical coefficient is 0.878, at df=14 it is 0.497, and quoting one number for
#: all n would understate the bar exactly where it matters most.
T_CRITICAL_ALPHA_05 = {
    1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365,
    8: 2.306, 9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179, 13: 2.160, 14: 2.145,
    15: 2.131, 16: 2.120, 17: 2.110, 18: 2.101, 19: 2.093, 20: 2.086,
    22: 2.074, 24: 2.064, 26: 2.056, 28: 2.048, 30: 2.042, 40: 2.021,
    60: 2.000, 120: 1.980,
}

#: Below this structural-margin spread (in TM), a rank correlation is ranking
#: differences smaller than a predicted-fold comparison can resolve. Measured on
#: the real survey the spread was 0.024 TM, which is why the coefficient there must
#: not be read as a relationship.
STRUCTURAL_SPAN_RESOLUTION_TM = 0.05


def spearman_critical_value(n: int, alpha_table: dict = T_CRITICAL_ALPHA_05) -> float | None:
    """The |rho| a two-tailed test at alpha = 0.05 requires for this n.

    Uses the t-approximation rho = t / sqrt(n - 2 + t^2) with the tabulated t for
    the nearest available degrees of freedom, erring toward the larger t when the
    exact df is not tabulated so the bar is never understated.
    """
    if n < 4:
        return None
    df = n - 2
    available = sorted(alpha_table)
    candidate = next((key for key in available if key >= df), available[-1])
    t = alpha_table[candidate]
    return t / (df + t * t) ** 0.5


def model_rank(query: str) -> int:
    """The ``rank_00N`` index of a prediction filename; 0 when it cannot be read."""
    marker = "_unrelaxed_rank_"
    index = query.find(marker)
    if index < 0:
        return 0
    digits = query[index + len(marker):index + len(marker) + 3]
    try:
        return int(digits)
    except ValueError:
        return 0


def select_models(
    by_candidate: dict[str, list[dict[str, str]]],
    selection: str,
) -> dict[str, list[dict[str, str]]]:
    """Reduce each candidate's rows to one representative model.

    ``rank_001`` (the default) keeps the comparison self-consistent: the anchor and
    the competitor are scored against the SAME model. Taking each side's best
    independently - the ``max`` option - can mix two different models of one
    candidate, which is not a single structural comparison and inflates both sides.
    A candidate with no rank_001 row falls back to all its rows rather than being
    dropped.
    """
    if selection not in ("rank_001", "max"):
        raise ValueError(f"unknown model selection {selection!r}")
    reduced: dict[str, list[dict[str, str]]] = {}
    for accession, rows in by_candidate.items():
        if selection == "max":
            reduced[accession] = rows
            continue
        ranked = [row for row in rows if model_rank(row["query"]) == 1]
        reduced[accession] = ranked or rows
    return reduced


def index_composition(composition_rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    """Map normalised accession -> composition row, failing closed on collisions.

    Two distinct accessions can normalise to the same filename (``A|B`` and
    ``A_B``), and silently picking one would attribute a structure to the wrong
    candidate. That is refused rather than guessed.
    """
    by_normalized: dict[str, dict[str, str]] = {}
    for row in composition_rows:
        accession = (row.get("accession") or "").strip()
        if not accession:
            raise ValueError("composition row without an accession")
        key = normalize_name(accession)
        if key in by_normalized:
            raise ValueError(
                f"accessions {by_normalized[key]['accession']!r} and {accession!r} both "
                f"normalise to {key!r}; the structure could not be attributed safely"
            )
        by_normalized[key] = row
    if not by_normalized:
        raise ValueError("the composition table is empty")
    return by_normalized


def score(
    foldseek_rows: list[dict[str, str]],
    composition_rows: list[dict[str, str]],
    model_selection: str = "rank_001",
) -> tuple[list[dict[str, str]], dict]:
    members = index_composition(composition_rows)

    raw_by_candidate: dict[str, list[dict[str, str]]] = {}
    for row in foldseek_rows:
        raw_by_candidate.setdefault(candidate_of(row["query"]), []).append(row)
    by_candidate = select_models(raw_by_candidate, model_selection)

    # The composition table is the authority on who is in the survey. A Foldseek
    # query outside it is not a survey member (a panel structure compared against
    # itself, for instance) and is counted rather than treated as an error; a
    # member with no rows at all IS an error, so a failed prediction cannot hide.
    ignored = sorted(name for name in by_candidate if name not in members)

    records: list[dict[str, str]] = []
    for key in sorted(members):
        rows = by_candidate.get(key)
        composition = members[key]
        if not rows:
            raise ValueError(
                f"{composition['accession']} has no Foldseek comparison at all; the prediction "
                "may have failed"
            )
        anchor_rows = [row for row in rows if row["target"].startswith(ANCHOR_PREFIX)]
        competitor_rows = [row for row in rows if row["target"].startswith(COMPETITOR_PREFIX)]
        if not anchor_rows:
            raise ValueError(f"{composition['accession']} has no comparison against the anchor")
        if not competitor_rows:
            raise ValueError(f"{composition['accession']} has no comparison against any competitor")
        best_anchor = max(anchor_rows, key=lambda row: row["tm"])
        best_competitor = max(competitor_rows, key=lambda row: row["tm"])
        structural_margin = best_anchor["tm"] - best_competitor["tm"]
        sequence_margin = (composition.get("panel_margin_competitor_minus_anchor_bits") or "").strip()
        records.append({
            "accession": composition["accession"].strip(),
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

    sequence_span = (max(sequence_values) - min(sequence_values)) if sequence_values else None
    structural_span = (max(structural_values) - min(structural_values)) if structural_values else None
    critical = spearman_critical_value(len(sequence_values))
    significant = (
        correlation is not None and critical is not None and abs(correlation) > critical
    )
    # Interpretability is about RESOLUTION, not significance: a coefficient computed
    # by ranking a band narrower than a predicted-fold comparison can resolve is
    # measuring the sample, not the biology, however large it happens to be.
    interpretable = structural_span is not None and structural_span >= STRUCTURAL_SPAN_RESOLUTION_TM
    if not interpretable:
        interpretation = (
            "NOT interpretable as a relationship. The structural margin spans only "
            f"{structural_span if structural_span is None else round(structural_span, 4)} TM "
            f"across a sequence margin spanning "
            f"{sequence_span if sequence_span is None else round(sequence_span, 2)} bits, which is "
            f"below the {STRUCTURAL_SPAN_RESOLUTION_TM} TM resolution floor: the ranking the "
            "coefficient consumes is dominated by differences smaller than a predicted-fold "
            "comparison can resolve. Read the spans and the winner counts instead - and note this "
            "guard is about resolution, independent of whether the coefficient is significant."
        )
    else:
        interpretation = (
            "The sequence margin is 'competitor minus anchor bits', the structural margin is "
            "'anchor minus competitor TM', so a POSITIVE Spearman means the sequence-level ranking "
            "predicts the structural one. Read it together with spearman_is_significant and the "
            "reported spans: a coefficient below the alpha=0.05 threshold for this n is inside the "
            "null range, and one computed over a narrow structural band is not evidence either way."
        )
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
        "sequence_margin_span_bits": sequence_span,
        "structural_margin_span_tm": structural_span,
        "structural_span_resolution_floor_tm": STRUCTURAL_SPAN_RESOLUTION_TM,
        "spearman_significance_threshold_alpha_05": critical,
        "spearman_is_significant": significant,
        "spearman_is_interpretable": interpretable,
        "model_selection": model_selection,
        "non_survey_queries_ignored": ignored,
        "non_survey_query_count": len(ignored),
        "interpretation": interpretation,
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
    parser.add_argument(
        "--model-selection", default="rank_001", choices=("rank_001", "max"),
        help="rank_001 (default) scores the anchor and the competitor against the SAME model; "
             "max lets each side pick its best model independently, which mixes models",
    )
    args = parser.parse_args(argv)

    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise ValueError(f"refusing to write into a non-empty output dir: {args.out_dir}")
    records, summary = score(read_foldseek(args.foldseek), read_tsv(args.composition),
                             args.model_selection)
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
