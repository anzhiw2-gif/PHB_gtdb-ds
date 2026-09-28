#!/usr/bin/env python3
"""Assemble the smallest defensible structural pilot for packet P4.

The frozen evidence measured single-sequence prediction as unusable for
comparison (AF2 pLDDT 36.1, same-sequence TM 0.20), so the structural stage needs
MSA-based prediction.  This module builds the smallest set that can demonstrate
the pipeline end to end, with every member chosen by frozen measurement rather
than by taste:

* two candidates taken from the two EXTREMES of the frozen pilot's blastp-vs-8YNV
  table -- the strongest hit (E=1.16e-122, bitscore 341) and the weakest
  (E~8.3e-4, bitscore ~26).  If the structural stage cannot separate those two,
  it cannot separate anything;
* the three HARD COMPETITORS of the panel (P24640 / Q02104 / Q88N36), which are
  exactly the sequences ``hfam_4`` failed to reject.

8YNV chain A is the positive reference and comes from the PDB rather than from a
prediction, so the comparison has one experimental anchor.

Candidate-only boundary: this assembles a pilot.  It asserts nothing about any
candidate and changes no disposition.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

ANCHOR_ACCESSION = "8YNV_A"

STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
#: Characters that are stripped rather than rejected.  ``*`` is a stop codon and
#: the frozen candidate FASTA carries it on some records; whitespace is not part
#: of a sequence at all.  The project already learned this: the pool-external
#: pool ships a ``.noasterisk`` variant for exactly this reason.
STRIPPED = "* \t\r\n"
AMBIGUOUS_AA = frozenset("XBZJUO")


def sanitize_sequence(name: str, sequence: str) -> tuple[str, list[str]]:
    """Strip ``*``/whitespace and fail closed on anything else non-standard.

    Predictors reject a stop codon outright (ColabFold: "Invalid character in the
    sequence: *"), so a pilot MUST clean its input rather than discover this after
    an hour of MSA generation.  Ambiguous residues are kept and reported, because
    dropping them would silently shorten the sequence.
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
    return cleaned, bad

PILOT_HITS = "pilot_vs_8ynv.tsv"
NPHAMCL_SUPERFAMILY = "intracellular nPHAMCL"


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


def write_fasta(path: Path, records: dict[str, str]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n{sequence}\n")


def pilot_extremes_from_scores(score_rows: list[dict[str, str]]) -> tuple[dict[str, str], dict[str, str]]:
    """The two extreme candidates of the P4 sequence-level competition.

    The P4 score table holds, per candidate, the best panel hit against the
    nPHAMCL anchor (8YNV_A) and against the hard competitors.  The extremes of
    the **margin** are the natural pilot ends: the most anchor-favoured candidate
    and the most competitor-favoured one.  If a structural stage cannot separate
    those two, it cannot separate anything.
    """
    parsed: list[tuple[float, dict[str, str]]] = []
    for row in score_rows:
        accession = (row.get("accession") or "").strip()
        if not accession:
            raise ValueError("score row without an accession")
        try:
            anchor = float(row.get("anchor_bitscore") or "")
            competitor = float(row.get("competitor_bitscore") or "")
        except ValueError:
            # An unscored row (e.g. "pending" on either side) carries no margin and
            # therefore cannot be a pilot extreme. It is skipped, never guessed at,
            # and the "fewer than two" check below still fails closed if nothing
            # usable remains.
            continue
        if anchor != anchor or competitor != competitor:  # NaN check
            continue
        parsed.append((competitor - anchor, {**row, "accession": accession,
                                            "margin": f"{competitor - anchor:.4f}"}))
    if len(parsed) < 2:
        raise ValueError("the P4 score table holds fewer than two scored candidates")
    parsed.sort(key=lambda item: item[0])
    most_anchor_favoured = parsed[0][1]
    most_competitor_favoured = parsed[-1][1]
    if most_anchor_favoured["accession"] == most_competitor_favoured["accession"]:
        raise ValueError("the P4 score table holds a single candidate; no extremes to take")
    return most_anchor_favoured, most_competitor_favoured


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--competition-scores", required=True, type=Path,
                        help="P4 nphamcl_competition_scores.tsv (per-candidate panel margins)")
    parser.add_argument("--candidate-union", required=True, type=Path)
    parser.add_argument("--merged", required=True, type=Path)
    parser.add_argument("--challenge-fasta", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--competitors", nargs="+", default=["P24640", "Q02104", "Q88N36"])
    args = parser.parse_args(argv)

    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise ValueError(f"refusing to write into a non-empty output dir: {args.out_dir}")

    anchor_favoured, competitor_favoured = pilot_extremes_from_scores(
        read_tsv(args.competition_scores),
    )
    merged = {
        (row.get("accession") or "").strip(): (row.get("superfamily") or "").strip()
        for row in read_tsv(args.merged)
    }
    union = read_fasta(args.candidate_union)
    challenge = read_fasta(args.challenge_fasta)

    members: dict[str, str] = {}
    provenance: dict[str, dict[str, str]] = {}
    for label, hit in (("most_anchor_favoured", anchor_favoured),
                       ("most_competitor_favoured", competitor_favoured)):
        accession = hit["accession"]
        if accession not in union:
            raise ValueError(f"pilot candidate {accession} has no sequence in the candidate union")
        superfamily = merged.get(accession, "")
        if superfamily != NPHAMCL_SUPERFAMILY:
            raise ValueError(
                f"pilot candidate {accession} is {superfamily!r}, not {NPHAMCL_SUPERFAMILY!r}"
            )
        cleaned, ambiguous = sanitize_sequence(accession, union[accession])
        members[accession] = cleaned
        provenance[accession] = {
            "role": label,
            "superfamily": superfamily,
            "panel_margin_competitor_minus_anchor_bits": hit["margin"],
            "anchor_bitscore": hit.get("anchor_bitscore", ""),
            "anchor_evalue": hit.get("anchor_evalue", ""),
            "competitor_subject": hit.get("competitor_subject", ""),
            "competitor_bitscore": hit.get("competitor_bitscore", ""),
            "competitor_evalue": hit.get("competitor_evalue", ""),
            "sequence_stage_disposition": hit.get("disposition", ""),
        }
        if ambiguous:
            provenance[accession]["ambiguous_residues_kept"] = "".join(ambiguous)
    for name in args.competitors:
        if name not in challenge:
            raise ValueError(f"competitor {name} has no sequence in the challenge FASTA")
        cleaned, ambiguous = sanitize_sequence(name, challenge[name])
        members[name] = cleaned
        provenance[name] = {"role": "hard_competitor", "superfamily": "",
                            "panel_margin_competitor_minus_anchor_bits": "",
                            "anchor_bitscore": "", "anchor_evalue": "",
                            "competitor_subject": "", "competitor_bitscore": "",
                            "competitor_evalue": "", "sequence_stage_disposition": ""}

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_fasta(args.out_dir / "pilot_prediction_targets.faa", members)
    (args.out_dir / "pilot_composition.json").write_text(
        json.dumps({
            "members": provenance,
            "member_count": len(members),
            "anchor": ANCHOR_ACCESSION,
            "anchor_source": "experimental PDB structure (downloaded), not a prediction",
            "why_these_candidates": (
                "the two extremes of the P4 sequence-level panel margins: the most "
                "anchor-favoured candidate and the most competitor-favoured one. If the "
                "structural stage cannot separate those two, it cannot separate anything."
            ),
            "boundary": (
                "Pilot assembly only. Asserts nothing about any candidate and changes no "
                "disposition."
            ),
        }, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n",
    )
    print(json.dumps({
        "members": sorted(members),
        "member_count": len(members),
        "most_anchor_favoured": anchor_favoured["accession"],
        "most_competitor_favoured": competitor_favoured["accession"],
        "out_dir": str(args.out_dir),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
