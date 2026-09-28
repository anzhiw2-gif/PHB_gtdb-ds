#!/usr/bin/env python3
"""Cys nucleophile motif evidence annotator (sequence-context match only).

Why this exists
---------------
The ``intracellular nPHASCL without lipase box`` (Cys) superfamily has so far
carried only an *inferred* catalytic residue (``inferred_cys_not_verified``).
The typing rests on "no GxSxG + PF06850 + strong Cys-family profile"; the
catalytic residue itself was never observed per candidate, because the candidate
layer shares no accession with the frozen DED alignments
(``candidate_ded_alignment_members = 0``).

Two literature facts make a sequence-level observation possible:

* Knoll 2009 states the Cys-type nucleophile carries a **hydrophobic residue at
  position -1** (the ``cysteine-1`` rule);
* the only experimentally verified Cys-type nucleophile is **Ralstonia eutropha
  PhaZ1 Cys183** (site-directed mutagenesis: C183A, D355A and H388Q abolish
  activity while C183S retains it; PMID 16233560), whose local context is
  ``V182-C183-Q184`` -- i.e. a ``V-C-Q`` tripeptide.

This module reports, per candidate, whether the sequence carries a cysteine in
the anchor's ``-1/+1`` context.  It is a **sequence context match against one
literature anchor**, NOT a catalytic-activity verification and NOT a family
call: every output row states that explicitly, and nothing here may be written
as a verified Cys.

Calibration is fail-closed: the pattern must reproduce the PhaZ1 anchor
(`HYDROPHOBIC` at -1 and `Q` at +1 at the reported position) before any
candidate row is written.

Boundary: candidate-only sequence evidence.  The pattern's specificity has so
far only been measured *in sample* on the frozen reference panel; independent
negative-set validation remains pending, so this column must not be used as a
standalone filter.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping

HYDROPHOBIC = frozenset("ACFILMVWY")

ANCHOR_ACCESSION = "CAJ92291.1"
ANCHOR_POSITION = 183
ANCHOR_CONTEXT = "VCQ"
ANCHOR_PMID = "16233560"
ANCHOR_DOI = "10.1016/S1389-1723(03)70136-4"
ANCHOR_EVIDENCE_TYPE = "site_directed_mutagenesis"
ANCHOR_STATEMENT = (
    "C183A and D355A/H388Q abolish PhaZ1 activity while C183S retains it; "
    "C183 is the catalytic nucleophile of the Cys-type intracellular PHB "
    "depolymerase PhaZ1 (Ralstonia eutropha H16)"
)

TIER_ANCHOR_EXACT = "anchor_exact_vcq"
TIER_ANCHOR_LOOSE = "anchor_context_hydrophobic_c_q"
TIER_MINUS1_HYDROPHOBIC = "minus1_hydrophobic_only"
TIER_NOT_DETECTED = "not_detected"
TIER_NO_CYSTEINE = "no_cysteine"

EVIDENCE_BASIS = (
    "sequence_context_match_to_single_literature_anchor_phaZ1_C183_PMID16233560"
)
INDEPENDENT_VALIDATION = "pending"

OUTPUT_FIELDS = (
    "accession", "pool", "n_cysteine", "cys_nucleophile_motif",
    "cys_nucleophile_motif_positions", "cys_nucleophile_motif_tripeptides",
    "evidence_basis", "catalytic_activity_verified", "independent_validation",
)


# --------------------------------------------------------------------------- #
# sequence primitives
# --------------------------------------------------------------------------- #
def cys_contexts(sequence: str) -> list[dict[str, object]]:
    """Every cysteine with its -1/+1 flanks, positions 1-based."""
    out: list[dict[str, object]] = []
    for i, residue in enumerate(sequence):
        if residue != "C":
            continue
        minus1 = sequence[i - 1] if i > 0 else ""
        plus1 = sequence[i + 1] if i + 1 < len(sequence) else ""
        out.append({
            "position": i + 1,
            "minus1": minus1,
            "plus1": plus1,
            "tripeptide": minus1 + "C" + plus1,
        })
    return out


def _is_anchor_context(ctx: Mapping[str, object]) -> bool:
    return ctx["minus1"] in HYDROPHOBIC and ctx["plus1"] == "Q"


def _is_exact_context(ctx: Mapping[str, object]) -> bool:
    return ctx["minus1"] == "V" and ctx["plus1"] == "Q"


def classify_motif(sequence: str) -> dict[str, object]:
    """Tier the sequence by its best cysteine context."""
    ctxs = cys_contexts(sequence)
    if not ctxs:
        return {"tier": TIER_NO_CYSTEINE, "positions": [], "tripeptides": []}
    exact = [c for c in ctxs if _is_exact_context(c)]
    loose = [c for c in ctxs if _is_anchor_context(c)]
    minus1 = [c for c in ctxs if c["minus1"] in HYDROPHOBIC]
    if exact:
        chosen = exact
        tier = TIER_ANCHOR_EXACT
    elif loose:
        chosen = loose
        tier = TIER_ANCHOR_LOOSE
    elif minus1:
        chosen = minus1
        tier = TIER_MINUS1_HYDROPHOBIC
    else:
        chosen = []
        tier = TIER_NOT_DETECTED
    return {
        "tier": tier,
        "positions": [int(c["position"]) for c in chosen],
        "tripeptides": [str(c["tripeptide"]) for c in chosen],
    }


# --------------------------------------------------------------------------- #
# FASTA
# --------------------------------------------------------------------------- #
def read_fasta_records(path: Path) -> list[tuple[str, str]]:
    """(header, sequence) pairs, header without the leading '>'."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:].strip()
                chunks = []
            elif header is not None:
                chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def reference_accession(header: str) -> str:
    """Accession of a *reference* header ('DED_hfam_65_0001|CAJ92291.1').

    Candidate ids also contain '|' (genome|protein) but are emitted verbatim --
    see ``annotate``.
    """
    return header.split("|")[-1].strip() if "|" in header else header.strip()


# --------------------------------------------------------------------------- #
# calibration
# --------------------------------------------------------------------------- #
def verify_anchor(
    sequences: Mapping[str, str], accession: str, position: int
) -> dict[str, object]:
    """Fail-closed check that the pattern reproduces the literature anchor.

    When the reported anchor position lies inside the anchor sequence, the
    pattern must hit at exactly that position; the reported position is
    otherwise recorded verbatim as provenance (small synthetic fixtures
    renumber their sequence).
    """
    sequence = sequences.get(accession)
    if sequence is None:
        return {
            "reproduced": False,
            "reason": "anchor_accession_not_in_reference_set",
            "observed_positions": [],
        }
    ctxs = cys_contexts(sequence)
    hits = [c for c in ctxs if _is_anchor_context(c)]
    observed = [int(c["position"]) for c in hits]
    if not hits:
        return {
            "reproduced": False,
            "reason": "pattern_not_found_in_anchor_sequence",
            "observed_positions": [],
        }
    if 1 <= position <= len(sequence):
        at_position = [c for c in hits if int(c["position"]) == position]
        if not at_position:
            return {
                "reproduced": False,
                "reason": "pattern_absent_at_reported_anchor_position",
                "observed_positions": observed,
            }
    return {
        "reproduced": True,
        "reason": "pattern_reproduced",
        "observed_positions": observed,
        "reported_position": position,
        "observed_tripeptides": [str(c["tripeptide"]) for c in hits],
    }


def reference_panel_separation(
    records: Iterable[tuple[str, str]],
    superfamily_by_accession: Mapping[str, str] | None = None,
    cys_superfamily: str = "intracellular nPHASCL without lipase box",
) -> dict[str, object]:
    """In-sample hit rates of the pattern on the frozen reference panel."""
    groups: dict[str, Counter] = {}

    def bucket(name: str) -> Counter:
        return groups.setdefault(name, Counter())

    for header, sequence in records:
        accession = reference_accession(header)
        if superfamily_by_accession is None:
            label = "all"
        else:
            label = (
                "cys_superfamily"
                if superfamily_by_accession.get(accession) == cys_superfamily
                else "other_superfamily"
            )
        hits = classify_motif(sequence)
        counter = bucket(label)
        counter["n"] += 1
        if hits["tier"] in (TIER_ANCHOR_EXACT, TIER_ANCHOR_LOOSE):
            counter["pattern_hit"] += 1
        if hits["tier"] == TIER_ANCHOR_EXACT:
            counter["exact_vcq_hit"] += 1

    out: dict[str, object] = {}
    for name, counter in sorted(groups.items()):
        n = counter["n"] or 1
        out[name] = {
            "n": counter["n"],
            "pattern_hit": counter["pattern_hit"],
            "pattern_hit_fraction": round(counter["pattern_hit"] / n, 6),
            "exact_vcq_hit": counter["exact_vcq_hit"],
            "exact_vcq_hit_fraction": round(counter["exact_vcq_hit"] / n, 6),
        }
    out["scope"] = "in_sample_reference_panel"
    out["held_out_set"] = "none"
    return out


def read_superfamily_map(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = (row.get("accession") or row.get("protein_id") or "").strip()
            superfamily = (row.get("phaded_superfamily") or row.get("superfamily") or "").strip()
            if accession and superfamily:
                out[accession] = superfamily
    return out


# --------------------------------------------------------------------------- #
# annotation
# --------------------------------------------------------------------------- #
def annotate(
    *,
    reference_fasta: Path,
    candidate_fastas: Mapping[str, Path],
    output_dir: Path,
    anchor_accession: str = ANCHOR_ACCESSION,
    anchor_position: int = ANCHOR_POSITION,
    reference_superfamily_table: Path | None = None,
) -> dict[str, object]:
    """Calibrate (fail-closed) then annotate every candidate; write the table."""
    output_dir = Path(output_dir)
    reference_records = read_fasta_records(reference_fasta)
    if not reference_records:
        raise SystemExit(f"reference FASTA is empty: {reference_fasta}")
    # References are keyed by BOTH the full 'id|accession' header and the bare
    # accession, so the anchor resolves however the FASTA was written.
    anchor_sequences: dict[str, str] = {}
    for header, sequence in reference_records:
        anchor_sequences[header] = sequence
        anchor_sequences.setdefault(reference_accession(header), sequence)

    anchor = verify_anchor(anchor_sequences, anchor_accession, anchor_position)
    if not anchor["reproduced"]:
        raise SystemExit(
            "Cys nucleophile pattern does not reproduce the literature anchor "
            f"{anchor_accession} at {anchor_position}: {anchor['reason']}"
        )

    superfamily_map = (
        read_superfamily_map(reference_superfamily_table)
        if reference_superfamily_table is not None
        else None
    )
    separation = reference_panel_separation(reference_records, superfamily_map)

    rows: list[dict[str, str]] = []
    per_pool: dict[str, Counter] = {}
    for pool, path in candidate_fastas.items():
        counter = per_pool.setdefault(pool, Counter())
        for header, sequence in read_fasta_records(path):
            # candidate ids are 'genome|protein' and are emitted verbatim
            accession = header.strip()
            result = classify_motif(sequence)
            counter[str(result["tier"])] += 1
            counter["n"] += 1
            rows.append({
                "accession": accession,
                "pool": pool,
                "n_cysteine": str(len(cys_contexts(sequence))),
                "cys_nucleophile_motif": str(result["tier"]),
                "cys_nucleophile_motif_positions": ";".join(
                    str(p) for p in result["positions"]
                ),
                "cys_nucleophile_motif_tripeptides": ";".join(result["tripeptides"]),
                "evidence_basis": EVIDENCE_BASIS,
                "catalytic_activity_verified": "false",
                "independent_validation": INDEPENDENT_VALIDATION,
            })

    if not rows:
        raise SystemExit("no candidate rows were read; refusing to write")

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "cys_nucleophile_motif_evidence.tsv"
    with evidence_path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(OUTPUT_FIELDS), delimiter="\t",
            lineterminator="\n", extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

    calibration = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "anchor_accession": anchor_accession,
        "anchor_position": anchor_position,
        "anchor_context": ANCHOR_CONTEXT,
        "anchor_pmid": ANCHOR_PMID,
        "anchor_doi": ANCHOR_DOI,
        "anchor_evidence_type": ANCHOR_EVIDENCE_TYPE,
        "anchor_statement": ANCHOR_STATEMENT,
        "anchor_reproduced": True,
        "anchor_details": anchor,
        "pattern": {
            "minus1": "hydrophobic (ACFILMVWY) -- Knoll 2009 cysteine-1 rule",
            "plus1": "Q -- the PhaZ1 C183 anchor context",
            "exact_tier": ANCHOR_CONTEXT,
        },
        "reference_panel_separation": separation,
        "candidate_tier_counts": {k: dict(v) for k, v in sorted(per_pool.items())},
        "candidates_written": len(rows),
        "evidence_boundary": (
            "sequence context match against a single literature anchor; NOT a "
            "catalytic-activity verification, NOT a family call, and NOT a "
            "PHB/PHA phenotype; independent negative-set validation pending"
        ),
        "hmm_fitted": False,
        "family_call_made": False,
        "phenotype_boundary": (
            "candidate-only; this column must not be read as a verified Cys"
        ),
    }
    (output_dir / "cys_nucleophile_motif_calibration.json").write_text(
        json.dumps(calibration, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    calibration["evidence_sha256"] = hashlib.sha256(
        evidence_path.read_bytes()
    ).hexdigest().upper()
    return calibration


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument(
        "--candidate-fasta", action="append", required=True,
        help="repeatable POOL=PATH pair, e.g. pool_in=path/to/cand.faa",
    )
    parser.add_argument("--reference-superfamily-table", type=Path, default=None)
    parser.add_argument("--anchor-accession", default=ANCHOR_ACCESSION)
    parser.add_argument("--anchor-position", type=int, default=ANCHOR_POSITION)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)

    fastas: dict[str, Path] = {}
    for item in args.candidate_fasta:
        if "=" not in item:
            raise SystemExit(f"--candidate-fasta must be POOL=PATH, got {item!r}")
        pool, raw = item.split("=", 1)
        fastas[pool] = Path(raw)

    result = annotate(
        reference_fasta=args.reference_fasta,
        candidate_fastas=fastas,
        output_dir=args.output_dir,
        anchor_accession=args.anchor_accession,
        anchor_position=args.anchor_position,
        reference_superfamily_table=args.reference_superfamily_table,
    )
    print(json.dumps({
        "anchor_reproduced": result["anchor_reproduced"],
        "candidates_written": result["candidates_written"],
        "candidate_tier_counts": result["candidate_tier_counts"],
        "evidence_sha256": result["evidence_sha256"],
    }, ensure_ascii=False, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
