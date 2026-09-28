#!/usr/bin/env python3
"""P4 structural survey: a margin-stratified tranche of the 987 nPHAMCL candidates.

The two-candidate pilot showed something that needs testing rather than
believing: **both** ends of the sequence-level margin pointed structurally at the
measured confounder, including the candidate with the largest anchor margin in the
whole set.  If that generalises, the sequence-level ``nphamcl_like_supported``
label carries no structural warranty anywhere in the 987.

Testing it needs candidates spread ACROSS the margin range, not just its two
extremes.  A proportional sample would not do: 90% of the 987 sit on the negative
side, so a proportional draw would barely touch the competitor-favoured tail that
the pilot found decisive.  This module therefore stratifies by **margin decile**
and draws a fixed number per decile, so every part of the range is represented by
construction.

Scope: this is a TRANCHE of the plan's structural stage, not the whole survey. The
pilot measured ~6-8 min per sequence including MSA generation, so all 987 would be
roughly 100-130 GPU-hours on one card; the tranche is sized to be answerable
within one background job and is stated as such.

Candidate-only boundary: selects sequences for structural review. Deletes,
demotes and promotes nothing, and a fold comparison is not a phenotype claim.
"""

from __future__ import annotations

import argparse
import csv
import json
import random
from pathlib import Path

NPHAMCL_SUPERFAMILY = "intracellular nPHAMCL"
DEFAULT_SEED = 20260928
DECILES = 10


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"input is not a regular file: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"empty TSV: {path}")
    return rows


def read_fasta(path: Path) -> dict[str, str]:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"FASTA is not a regular file: {path}")
    records: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8", errors="strict") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(chunks)
                current = line[1:].split()[0]
                if current in records:
                    raise ValueError(f"duplicate FASTA header: {current}")
                chunks = []
            else:
                if current is None:
                    raise ValueError("sequence before the first header")
                chunks.append(line)
    if current is not None:
        records[current] = "".join(chunks)
    if not records:
        raise ValueError(f"empty FASTA: {path}")
    return records


STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
STRIPPED = "* \t\r\n"
AMBIGUOUS_AA = frozenset("XBZJUO")


def sanitize_sequence(name: str, sequence: str) -> str:
    """Strip ``*``/whitespace, fail closed on anything else non-standard.

    Learned the hard way: the candidate FASTA carries a trailing stop codon on
    some records and ColabFold rejects it outright, so the survey must not repeat
    the pilot's wasted run.
    """
    cleaned = "".join(character for character in sequence if character not in STRIPPED)
    bad = sorted({character for character in cleaned if character not in STANDARD_AA})
    unexpected = [character for character in bad if character not in AMBIGUOUS_AA]
    if unexpected:
        raise ValueError(
            f"sequence {name} carries unsupported residue(s) {unexpected}; "
            "sanitise the source FASTA rather than guessing"
        )
    if not cleaned:
        raise ValueError(f"sequence {name} is empty after sanitisation")
    return cleaned


def scored_candidates(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Candidates with a usable panel margin, as floats, sorted by margin."""
    scored: list[dict[str, str]] = []
    for row in rows:
        accession = (row.get("accession") or "").strip()
        if not accession:
            raise ValueError("score row without an accession")
        try:
            anchor = float(row.get("anchor_bitscore") or "")
            competitor = float(row.get("competitor_bitscore") or "")
        except ValueError:
            continue
        scored.append({**row, "accession": accession,
                       "margin": competitor - anchor})
    if not scored:
        raise ValueError("no candidate carries a usable panel margin")
    scored.sort(key=lambda row: (row["margin"], row["accession"]))
    return scored


def decile_strata(scored: list[dict[str, str]]) -> list[list[dict[str, str]]]:
    """Split the sorted candidates into deciles, largest margin last.

    The split is on POSITION, not on quantile edges, so ties cannot collapse a
    stratum to nothing; each decile holds ceil/floor of n/10 members.
    """
    total = len(scored)
    if total < DECILES:
        raise ValueError(f"need at least {DECILES} candidates to form deciles, got {total}")
    base, remainder = divmod(total, DECILES)
    strata: list[list[dict[str, str]]] = []
    start = 0
    for index in range(DECILES):
        size = base + (1 if index < remainder else 0)
        strata.append(scored[start:start + size])
        start += size
    if sum(len(stratum) for stratum in strata) != total:
        raise ValueError("decile split lost or duplicated a candidate")
    return strata


def draw(
    strata: list[list[dict[str, str]]],
    per_decile: int,
    seed: int,
) -> list[dict[str, str]]:
    if per_decile < 1:
        raise ValueError("per-decile draw must be at least 1")
    rng = random.Random(seed)
    drawn: list[dict[str, str]] = []
    for index, stratum in enumerate(strata):
        if per_decile > len(stratum):
            raise ValueError(
                f"decile {index} holds {len(stratum)} candidates, fewer than the "
                f"{per_decile} requested"
            )
        picks = rng.sample(stratum, per_decile)
        for pick in picks:
            drawn.append({**pick, "margin_decile": index})
    if len({row["accession"] for row in drawn}) != len(drawn):
        raise ValueError("the stratified draw repeated a candidate")
    return drawn


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--competition-scores", required=True, type=Path)
    parser.add_argument("--candidate-union", required=True, type=Path)
    parser.add_argument("--merged", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--per-decile", type=int, default=3)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--expect-scored", type=int, default=987)
    args = parser.parse_args(argv)

    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise ValueError(f"refusing to write into a non-empty output dir: {args.out_dir}")

    scored = scored_candidates(read_tsv(args.competition_scores))
    if args.expect_scored is not None and len(scored) != args.expect_scored:
        raise ValueError(
            f"{len(scored)} scored candidates does not match the declared {args.expect_scored}; "
            "the universe changed, so the stratification would not be the declared one"
        )
    superfamilies = {
        (row.get("accession") or "").strip(): (row.get("superfamily") or "").strip()
        for row in read_tsv(args.merged)
    }
    union = read_fasta(args.candidate_union)
    strata = decile_strata(scored)
    drawn = draw(strata, args.per_decile, args.seed)

    members: dict[str, str] = {}
    provenance: list[dict[str, str]] = []
    for row in drawn:
        accession = row["accession"]
        superfamily = superfamilies.get(accession, "")
        if superfamily != NPHAMCL_SUPERFAMILY:
            raise ValueError(
                f"survey candidate {accession} is {superfamily!r}, not {NPHAMCL_SUPERFAMILY!r}"
            )
        if accession not in union:
            raise ValueError(f"survey candidate {accession} has no sequence in the candidate union")
        members[accession] = sanitize_sequence(accession, union[accession])
        provenance.append({
            "accession": accession,
            "margin_decile": str(row["margin_decile"]),
            "panel_margin_competitor_minus_anchor_bits": f"{row['margin']:.4f}",
            "anchor_bitscore": row.get("anchor_bitscore", ""),
            "anchor_evalue": row.get("anchor_evalue", ""),
            "competitor_subject": row.get("competitor_subject", ""),
            "competitor_bitscore": row.get("competitor_bitscore", ""),
            "sequence_stage_disposition": row.get("disposition", ""),
        })

    args.out_dir.mkdir(parents=True, exist_ok=True)
    with (args.out_dir / "survey_targets.faa").open("w", encoding="ascii", newline="\n") as handle:
        for name, sequence in members.items():
            handle.write(f">{name}\n{sequence}\n")
    fields = list(provenance[0].keys())
    with (args.out_dir / "survey_composition.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(provenance)
    margins = [row["margin"] for row in scored]
    summary = {
        "scored_candidates": len(scored),
        "deciles": DECILES,
        "per_decile": args.per_decile,
        "sample_size": len(members),
        "seed": args.seed,
        "margin_range": {"min": min(margins), "max": max(margins)},
        "margin_median": margins[len(margins) // 2],
        "share_negative_margin": sum(1 for value in margins if value < 0) / len(margins),
        "stratum_sizes": [len(stratum) for stratum in strata],
        "draw_per_decile": {
            str(index): sorted(row["accession"] for row in drawn if row["margin_decile"] == index)
            for index in range(DECILES)
        },
        "why_deciles_not_proportional": (
            f"{sum(1 for value in margins if value < 0) / len(margins):.1%} of the candidates sit on "
            "the anchor-favoured side, so a proportional draw would barely touch the "
            "competitor-favoured tail that the pilot found decisive. A fixed draw per decile "
            "represents every part of the margin range by construction."
        ),
        "scope": (
            "A TRANCHE of the plan's structural stage, not the whole survey: the pilot measured "
            "~6-8 min per sequence including MSA generation, so all 987 would be roughly "
            "100-130 GPU-hours on one card."
        ),
        "boundary": (
            "Selects sequences for structural review. Deletes, demotes and promotes nothing, and "
            "a fold comparison is not a phenotype claim."
        ),
    }
    (args.out_dir / "survey_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n",
    )
    print(json.dumps({
        "scored_candidates": len(scored),
        "sample_size": len(members),
        "seed": args.seed,
        "per_decile": args.per_decile,
        "out_dir": str(args.out_dir),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
