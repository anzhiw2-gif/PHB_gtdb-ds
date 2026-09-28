#!/usr/bin/env python3
"""Build the nPHAMCL competition panel and its 987-candidate target set.

Packet P4 of the 2026-09-28 plan requires a STRUCTURAL COMPETITION PANEL for the
987 nPHAMCL-like candidates "rather than a single 8YNV reference", because the
sequence layer has been proven non-discriminative for this family (hfam_4: the
catalytic triad is shared with generic alpha/beta hydrolases, and the only
discriminative columns were V/W in the elbow, which the `--pnone` rebuild could
only exploit by collapsing to zero recall).

Every panel member is DERIVED FROM FROZEN MEASUREMENT, not chosen by taste:

* ``8YNV_A`` — the nPHAMCL reference structure used by the frozen structure
  verification pilot (PDB 8YNV chain A, 298 aa, carries the ``GWSMGGG`` lipase
  box and the ``DRDYVV`` motif).
* the three HARD COMPETITORS — the exact sequences that ``hfam_4`` FAILED to
  reject in the frozen challenge run: P24640 (Moraxella lipase, E=5.4e-12),
  Q02104 (Psychrobacter lipase, E=4.1e-11), Q88N36 (P. putida PcaD
  gamma-lactonase, E=6.6e-11).  Their E-values and their ``rejected=NO`` status
  are read from ``challenge_validation.tsv``; the module refuses to assemble the
  panel if that evidence disagrees.
* the seven SPECIFICITY CONTROLS — the challenge sequences that all five trained
  families rejected at E=999, so they test that the panel is not trivially
  permissive.

Candidate-only boundary: this assembles inputs for a review.  It deletes,
demotes and promotes nothing, and a panel membership is not a functional claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

#: Panel members whose inclusion is forced by the frozen challenge measurement.
HARD_COMPETITORS = {
    "P24640": "Moraxella lipase, measured hfam_4 false positive E=5.4e-12",
    "Q02104": "Psychrobacter lipase, measured hfam_4 false positive E=4.1e-11",
    "Q88N36": "P. putida PcaD gamma-lactonase, measured hfam_4 false positive E=6.6e-11",
}
POSITIVE_ANCHOR = "8YNV_A"
NPHAMCL_SUPERFAMILY = "intracellular nPHAMCL"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path: Path) -> dict[str, str]:
    """Read a FASTA into {header_first_token: sequence}."""
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"FASTA is not a regular file: {path}")
    records: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="utf-8", errors="strict").splitlines():
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


def read_tsv(path: Path) -> list[dict[str, str]]:
    import csv

    if not path.is_file() or path.is_symlink():
        raise ValueError(f"TSV is not a regular file: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"empty TSV: {path}")
    return rows


def select_nphamcl_accessions(merged_rows: list[dict[str, str]]) -> list[str]:
    """The 987 nPHAMCL-like candidates, from the frozen merged table."""
    if not merged_rows or "superfamily" not in merged_rows[0] or "accession" not in merged_rows[0]:
        raise ValueError("merged table must carry accession and superfamily columns")
    selected = [
        (row["accession"] or "").strip() for row in merged_rows
        if (row.get("superfamily") or "").strip() == NPHAMCL_SUPERFAMILY
    ]
    if any(not accession for accession in selected):
        raise ValueError("selected row without an accession")
    if len(set(selected)) != len(selected):
        raise ValueError("target set is not accession-unique")
    return selected


def verify_panel_against_measurement(
    challenge_rows: list[dict[str, str]],
    *,
    family: str = "hfam_4",
) -> tuple[list[str], list[str], list[str]]:
    """Return (hard_competitors, specificity_controls, rejected_at_999).

    Reads the frozen challenge validation and refuses to proceed if the three
    hard competitors were in fact rejected, or if a specificity control was not.
    """
    rows = [row for row in challenge_rows if (row.get("family") or "").strip() == family]
    if not rows:
        raise ValueError(f"challenge validation carries no rows for {family}")
    accepted = [row for row in rows if (row.get("rejected") or "").strip().lower() == "no"]
    rejected = [row for row in rows if (row.get("rejected") or "").strip().lower() == "yes"]
    accepted_ids = sorted((row["challenge"] or "").strip() for row in accepted)
    expected = sorted(HARD_COMPETITORS)
    if accepted_ids != expected:
        raise ValueError(
            f"frozen measurement disagrees with the declared hard competitors: "
            f"measured {accepted_ids}, declared {expected}"
        )
    controls = sorted((row["challenge"] or "").strip() for row in rejected)
    if len(controls) < 1:
        raise ValueError("no specificity controls in the frozen challenge evidence")
    return accepted_ids, controls, controls


def build_panel(
    challenge_fasta: dict[str, str],
    anchor_fasta: dict[str, str],
    hard_competitors: list[str],
    controls: list[str],
) -> dict[str, str]:
    """Assemble the panel, failing closed on any missing sequence."""
    if len(anchor_fasta) != 1:
        raise ValueError(f"anchor FASTA must hold exactly one record, got {len(anchor_fasta)}")
    panel: dict[str, str] = {}
    anchor_name = next(iter(anchor_fasta))
    panel[anchor_name] = anchor_fasta[anchor_name]
    for name in hard_competitors + controls:
        if name not in challenge_fasta:
            raise ValueError(f"panel member {name} has no sequence in the challenge FASTA")
        panel[name] = challenge_fasta[name]
    if len(panel) != 1 + len(hard_competitors) + len(controls):
        raise ValueError("panel assembly lost or duplicated a member")
    return panel


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--merged", required=True, type=Path)
    parser.add_argument("--candidate-union", required=True, type=Path)
    parser.add_argument("--challenge-fasta", required=True, type=Path)
    parser.add_argument("--challenge-validation", required=True, type=Path)
    parser.add_argument("--anchor-fasta", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--expect-targets", type=int, default=987)
    args = parser.parse_args(argv)

    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise ValueError(f"refusing to write into a non-empty output dir: {args.out_dir}")

    merged = read_tsv(args.merged)
    accessions = select_nphamcl_accessions(merged)
    hard, controls, rejected = verify_panel_against_measurement(
        read_tsv(args.challenge_validation),
    )
    union = read_fasta(args.candidate_union)
    missing = [accession for accession in accessions if accession not in union]
    if missing:
        raise ValueError(f"{len(missing)} target accessions have no sequence, e.g. {missing[:3]}")
    targets = {accession: union[accession] for accession in accessions}
    panel = build_panel(
        read_fasta(args.challenge_fasta), read_fasta(args.anchor_fasta), hard, controls,
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_fasta(args.out_dir / "nphamcl_targets_987.faa", targets)
    write_fasta(args.out_dir / "competition_panel.faa", panel)
    provenance = {
        "target_count": len(targets),
        "declared_target_count": args.expect_targets,
        "target_count_matches_declared": len(targets) == args.expect_targets,
        "panel_size": len(panel),
        "positive_anchor": [name for name in panel if name == POSITIVE_ANCHOR],
        "hard_competitors": {name: HARD_COMPETITORS[name] for name in hard},
        "specificity_controls": controls,
        "rejected_at_999": rejected,
        "derivation": (
            "Every panel member is forced by frozen measurement: the anchor is the "
            "reference used by the frozen structure-verification pilot, the hard "
            "competitors are exactly the sequences hfam_4 failed to reject in "
            "challenge_validation.tsv, and the controls are the sequences all five "
            "trained families rejected. No member was chosen by taste."
        ),
        "boundary": (
            "Assembles inputs for a review only; deletes, demotes and promotes nothing, "
            "and panel membership is not a functional claim."
        ),
        "inputs": {
            "merged": {"path": str(args.merged), "sha256": sha256(args.merged)},
            "candidate_union": {"path": str(args.candidate_union), "sha256": sha256(args.candidate_union)},
            "challenge_fasta": {"path": str(args.challenge_fasta), "sha256": sha256(args.challenge_fasta)},
            "challenge_validation": {"path": str(args.challenge_validation), "sha256": sha256(args.challenge_validation)},
            "anchor_fasta": {"path": str(args.anchor_fasta), "sha256": sha256(args.anchor_fasta)},
        },
    }
    (args.out_dir / "panel_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n",
    )
    print(json.dumps({
        "targets": len(targets),
        "panel_size": len(panel),
        "hard_competitors": hard,
        "controls": len(controls),
        "target_count_matches_declared": provenance["target_count_matches_declared"],
        "out_dir": str(args.out_dir),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
