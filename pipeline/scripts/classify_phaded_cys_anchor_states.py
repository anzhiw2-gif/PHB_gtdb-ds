#!/usr/bin/env python3
"""P3: recast the Cys anchor evidence into substitution / truncation / uncertain.

Packet P3 requires the 29,974 Cys candidates to be mapped against the literature
anchor while **explicitly distinguishing substitution, truncation and uncertain**
states, and requires that uncertainty is never written as a negative.

The frozen layer already holds the *pattern strength* per candidate
(``cys_nucleophile_motif`` in {anchor_exact_vcq,
anchor_context_hydrophobic_c_q, minus1_hydrophobic_only, not_detected,
no_cysteine}) but not the plan's state vocabulary.  This module adds that
vocabulary with rules derived from the single documented anchor rather than from
taste:

* anchor: PhaZ1 **Cys183** (*R. eutropha* H16), context ``VCQ``, PMID 16233560,
  site-directed mutagenesis — C183A and D355A/H388Q abolish activity while
  **C183S retains it**.
* so the anchor POSITION is 183, and a protein shorter than that cannot host it.

The four states, each with its rule stated in the code and in the output:

``pattern_present``
    the anchor context is reproduced (exact ``VCQ`` or the documented loose
    ``hydrophobic-C-Q``).
``substitution``
    a Cys is present with the documented hydrophobic -1 but the +1 is not Q --
    the residue is there and its context differs from the anchor.
``truncation``
    the sequence is shorter than the anchor position, so the anchor cannot be
    present at all.
``uncertain``
    everything else, including a full-length protein with **no** cysteine.
    Uncertainty is NOT a negative: C183S retains activity, so an absent or
    unmatched pattern cannot be read as "not a Cys-type enzyme".

Candidate-only boundary: this is sequence-context evidence against ONE literature
anchor.  It is not catalytic-activity verification, not a family call, and not a
phenotype claim, and it changes no disposition.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

ANCHOR_ACCESSION = "CAJ92291.1"
ANCHOR_POSITION = 183
ANCHOR_CONTEXT = "VCQ"
ANCHOR_PMID = "16233560"

TIER_PATTERN = ("anchor_exact_vcq", "anchor_context_hydrophobic_c_q")
TIER_SUBSTITUTION = ("minus1_hydrophobic_only",)

STATE_PATTERN = "pattern_present"
STATE_SUBSTITUTION = "substitution"
STATE_TRUNCATION = "truncation"
STATE_UNCERTAIN = "uncertain"

CYS_SUPERFAMILY = "intracellular nPHASCL without lipase box"

CAVEAT = (
    "Sequence context against a single literature anchor (PhaZ1 C183, PMID "
    f"{ANCHOR_PMID}). NOT catalytic-activity verification, NOT a family call. "
    "'uncertain' is not a negative: C183S retains activity, so an absent or "
    "unmatched pattern must never be read as evidence against a Cys-type enzyme."
)


def read_fasta_lengths(path: Path) -> dict[str, int]:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"FASTA is not a regular file: {path}")
    lengths: dict[str, int] = {}
    current: str | None = None
    count = 0
    with path.open(encoding="utf-8", errors="strict") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    lengths[current] = count
                current = line[1:].split()[0]
                if current in lengths:
                    raise ValueError(f"duplicate FASTA header: {current}")
                count = 0
            else:
                count += len(line)
    if current is not None:
        lengths[current] = count
    if not lengths:
        raise ValueError(f"empty FASTA: {path}")
    return lengths


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"TSV is not a regular file: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"empty TSV: {path}")
    return rows


def anchor_state(tier: str, n_cysteine: int, length: int) -> tuple[str, str]:
    """Return (state, reason). Rules are ordered and each names its evidence."""
    if length < ANCHOR_POSITION:
        return STATE_TRUNCATION, (
            f"length {length} < anchor position {ANCHOR_POSITION}: the anchor cannot be present"
        )
    if tier in TIER_PATTERN:
        return STATE_PATTERN, f"tier {tier}: documented anchor context reproduced"
    if tier in TIER_SUBSTITUTION:
        return STATE_SUBSTITUTION, (
            f"tier {tier}: cysteine present with the hydrophobic -1 but +1 is not Q"
        )
    if tier == "no_cysteine":
        return STATE_UNCERTAIN, (
            "no cysteine in a full-length protein: the nucleophile residue is absent, which is "
            "neither a clean substitution call nor decidable from context alone"
        )
    if tier == "not_detected":
        return STATE_UNCERTAIN, (
            f"tier not_detected with n_cysteine={n_cysteine}: no hydrophobic -1 cysteine; "
            "context alone cannot decide"
        )
    return STATE_UNCERTAIN, f"unrecognised tier {tier!r}: reported rather than forced into a call"


def classify(
    evidence_rows: list[dict[str, str]],
    lengths: dict[str, int],
    target_accessions: set[str] | None = None,
) -> tuple[list[dict[str, str]], dict]:
    out: list[dict[str, str]] = []
    missing_length: list[str] = []
    for row in evidence_rows:
        accession = (row.get("accession") or "").strip()
        if not accession:
            raise ValueError("evidence row without an accession")
        if target_accessions is not None and accession not in target_accessions:
            continue
        tier = (row.get("cys_nucleophile_motif") or "").strip()
        try:
            n_cysteine = int((row.get("n_cysteine") or "0").strip() or 0)
        except ValueError as exc:
            raise ValueError(f"non-integer n_cysteine for {accession}") from exc
        length = lengths.get(accession)
        if length is None:
            missing_length.append(accession)
            state, reason = STATE_UNCERTAIN, "sequence length not found; cannot test the anchor position"
        else:
            state, reason = anchor_state(tier, n_cysteine, length)
        out.append({
            "accession": accession,
            "pool": (row.get("pool") or "").strip(),
            "superfamily": (row.get("pool") or "").strip(),
            "length": str(length) if length is not None else "pending",
            "n_cysteine": str(n_cysteine),
            "frozen_tier": tier,
            "anchor_state": state,
            "state_rule": reason,
            "evidence_basis": (row.get("evidence_basis") or "").strip(),
            "catalytic_activity_verified": (row.get("catalytic_activity_verified") or "").strip(),
        })
    summary = summarize(out)
    summary["accessions_missing_sequence_length"] = missing_length
    return out, summary


def summarize(rows: list[dict[str, str]]) -> dict:
    counts: dict[str, int] = {}
    by_tier: dict[str, dict[str, int]] = {}
    for row in rows:
        counts[row["anchor_state"]] = counts.get(row["anchor_state"], 0) + 1
        by_tier.setdefault(row["frozen_tier"], {})
        by_tier[row["frozen_tier"]][row["anchor_state"]] = (
            by_tier[row["frozen_tier"]].get(row["anchor_state"], 0) + 1
        )
    total = len(rows)
    return {
        "candidates": total,
        "state_counts": counts,
        "state_shares": {key: value / total for key, value in counts.items()} if total else {},
        "tier_to_state": by_tier,
        "anchor": {
            "accession": ANCHOR_ACCESSION,
            "position": ANCHOR_POSITION,
            "context": ANCHOR_CONTEXT,
            "pmid": ANCHOR_PMID,
            "evidence_type": "site_directed_mutagenesis",
        },
        "rule_order": [
            f"truncation: length < {ANCHOR_POSITION}",
            f"pattern_present: tier in {list(TIER_PATTERN)}",
            f"substitution: tier in {list(TIER_SUBSTITUTION)}",
            "uncertain: every other tier, including no_cysteine at full length",
        ],
        "caveat": CAVEAT,
        "boundary": (
            "Adds a state vocabulary to frozen evidence; deletes, demotes and promotes nothing, "
            "and makes no family call."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--motif-evidence", required=True, type=Path)
    parser.add_argument("--candidate-union", required=True, type=Path)
    parser.add_argument(
        "--extra-fasta", type=Path, action="append", default=[],
        help="additional sequence source (e.g. the pool-external pool); repeatable. "
             "The pool-external Cys candidates are not in candidate_union.faa, so without "
             "this their sequence length is pending and their state is uncertain.",
    )
    parser.add_argument("--merged", required=True, type=Path,
                        help="merged table, to restrict to the Cys candidates")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--expect-candidates", type=int, default=29974)
    args = parser.parse_args(argv)

    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise ValueError(f"refusing to write into a non-empty output dir: {args.out_dir}")

    merged = read_rows(args.merged)
    targets = {
        (row.get("accession") or "").strip()
        for row in merged
        if (row.get("nucleophile_type") or "").strip() == "cys"
    }
    lengths = read_fasta_lengths(args.candidate_union)
    extra_sources: dict[str, str] = {}
    for extra in args.extra_fasta:
        extra_lengths = read_fasta_lengths(extra)
        for name, length in extra_lengths.items():
            if name in lengths:
                raise ValueError(f"sequence {name} appears in more than one FASTA source")
            lengths[name] = length
            extra_sources[name] = str(extra)
    rows, summary = classify(read_rows(args.motif_evidence), lengths, targets)
    summary["extra_fasta_sources"] = [str(path) for path in args.extra_fasta]
    summary["sequences_from_extra_sources"] = len(extra_sources)
    summary["cys_candidates_from_merged"] = len(targets)
    summary["declared_expectation"] = args.expect_candidates
    summary["expectation_matches"] = len(targets) == args.expect_candidates
    summary["candidates_scored"] = len(rows)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys()) if rows else ["accession"]
    with (args.out_dir / "cys_anchor_states.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (args.out_dir / "cys_anchor_states_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n",
    )
    print(json.dumps({
        "cys_candidates": len(targets),
        "candidates_scored": len(rows),
        "state_counts": summary["state_counts"],
        "expectation_matches": summary["expectation_matches"],
        "out_dir": str(args.out_dir),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
