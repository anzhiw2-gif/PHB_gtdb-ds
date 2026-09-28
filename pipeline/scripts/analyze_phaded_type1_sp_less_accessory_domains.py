#!/usr/bin/env python3
"""P2: accessory-domain check for the SP-less type-1 candidates.

Decision question (packet P2 of the 2026-09-28 redesign plan, sharpened by the P1
result): the candidate layer holds 48,038 type-1 candidates on the
secretion-signal rule alone; 11,610 of them also carry `type1_verified`
catalytic-domain geometry and a trained-profile hit.  P1 measured that only
46.8% of the *inherited*-label extracellular references and 40.8% of the
inherited type-1 references are predicted to carry a secretion signal, while
100% of the experimentally documented extracellular references are.  So the
signal-peptide gate cannot decide these candidates, and the reference layer
cannot be used to justify holding them.

This module asks the INDEPENDENT question: do those candidates carry the type-1
accessory/connector domains (FN3 / Ig-like / TSP3 / CHB) -- an extracellular
line of evidence that does NOT run through SignalP?

This is an OPTIONAL-module annotation, never a pass/hold criterion (the frozen
layer's own header states this).  Nothing here deletes, demotes or promotes a
candidate; it adds evidence flags.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

#: The project's authoritative type-1 accessory marker set, copied verbatim from
#: runs/20260920_phaded_three_gaps_01/plan/build_type1_accessory_layer.py so the
#: two layers cannot drift apart.
MARKERS = {
    "fn3": ("PF00041", "SM00060", "cd00063", "SSF49265", "PS50853", "G3DSA:2.60.40.10"),
    "ig_like": ("PF16403", "PF17957"),
    "tsp3": ("PF02412",),
    "chb": ("PF13290",),
}
MARKER_ORDER = ("fn3", "ig_like", "tsp3", "chb")

TYPE1_SUPERFAMILY = "extracellular dPHASCL type 1"
LOCALIZATION_HOLD_REASON = "localization_conflict"
TYPE1_GEOMETRY = "type1_verified"

STATUS_YES = "yes"
STATUS_NO = "no"
STATUS_PENDING = "pending"


def tokens(*fields: str) -> set[str]:
    """Split comma/semicolon/pipe/space separated signature fields into tokens."""
    out: set[str] = set()
    for field in fields:
        if not field:
            continue
        for part in field.replace(",", " ").replace(";", " ").replace("|", " ").split():
            out.add(part.strip())
    return out


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"input is not a regular file: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"empty input: {path}")
    return rows


def select_sp_less_type1(hold_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """The held SP-less type-1 target set, by the documented three keys.

    superfamily == extracellular dPHASCL type 1
      AND hold_reason == localization_conflict   (the held-on-signal subset)
      AND catalytic_domain_type == type1_verified (type-1 geometry confirmed)
    """
    missing = [
        column for column in ("accession", "superfamily", "hold_reason", "catalytic_domain_type")
        if column not in (hold_rows[0] if hold_rows else {})
    ]
    if missing:
        raise ValueError(f"hold table missing required columns: {','.join(missing)}")
    selected = [
        row for row in hold_rows
        if (row.get("superfamily") or "").strip() == TYPE1_SUPERFAMILY
        and (row.get("hold_reason") or "").strip() == LOCALIZATION_HOLD_REASON
        and (row.get("catalytic_domain_type") or "").strip() == TYPE1_GEOMETRY
    ]
    accessions = [(row.get("accession") or "").strip() for row in selected]
    if any(not accession for accession in accessions):
        raise ValueError("selected row without an accession")
    if len(set(accessions)) != len(accessions):
        duplicates = [a for a, n in Counter(accessions).items() if n > 1]
        raise ValueError(f"target set is not accession-unique: {duplicates[:5]}")
    return selected


def accessory_hit(signatures: str, pfam: str, interpro: str) -> dict[str, str]:
    """Presence of each accessory domain from the frozen signature fields."""
    observed = tokens(signatures, pfam, interpro)
    result: dict[str, str] = {}
    for name in MARKER_ORDER:
        result[name] = STATUS_YES if observed.intersection(MARKERS[name]) else STATUS_NO
    return result


def classify_rows(
    target_rows: list[dict[str, str]],
    matrix_by_accession: dict[str, dict[str, str]],
) -> tuple[list[dict[str, str]], dict]:
    joined: list[dict[str, str]] = []
    absent_from_matrix: list[str] = []
    for row in target_rows:
        accession = (row.get("accession") or "").strip()
        matrix_row = matrix_by_accession.get(accession)
        if matrix_row is None:
            absent_from_matrix.append(accession)
            hits = {name: STATUS_PENDING for name in MARKER_ORDER}
            pfam = interpro = signatures = ""
            in_matrix = "no"
        else:
            in_matrix = "yes"
            pfam = (matrix_row.get("pfam_accessions") or "").strip()
            interpro = (matrix_row.get("interpro_signatures") or "").strip()
            signatures = (matrix_row.get("interpro_signatures") or "").strip()
            hits = accessory_hit(signatures, pfam, interpro)
        present = [name for name in MARKER_ORDER if hits[name] == STATUS_YES]
        joined.append({
            "accession": accession,
            "genome": (row.get("genome") or "").strip(),
            "superfamily": (row.get("superfamily") or "").strip(),
            "hold_reason": (row.get("hold_reason") or "").strip(),
            "catalytic_domain_type": (row.get("catalytic_domain_type") or "").strip(),
            "signalp_class": (row.get("signalp_class") or "").strip(),
            "in_matrix": in_matrix,
            "has_fn3": hits["fn3"],
            "has_ig_like": hits["ig_like"],
            "has_tsp3": hits["tsp3"],
            "has_chb": hits["chb"],
            "n_accessory_domains": str(len(present)),
            "accessory_domains": ";".join(present),
            "pfam_accessions": pfam,
            "interpro_signatures": interpro,
        })
    summary = summarize(joined)
    summary["target_rows_absent_from_matrix"] = absent_from_matrix
    return joined, summary


def summarize(joined: list[dict[str, str]]) -> dict:
    total = len(joined)
    any_domain = [row for row in joined if int(row["n_accessory_domains"]) > 0]
    resolved = [row for row in joined if row["in_matrix"] == "yes"]
    per_domain = {
        name: sum(1 for row in joined if row[f"has_{name}"] == STATUS_YES)
        for name in MARKER_ORDER
    }
    combination_counts = Counter(row["accessory_domains"] for row in joined)
    return {
        "target_count": total,
        "in_matrix_count": len(resolved),
        "with_any_accessory_domain": len(any_domain),
        "with_any_accessory_domain_rate": (len(any_domain) / total) if total else None,
        "per_domain_counts": per_domain,
        "combination_counts": dict(combination_counts.most_common()),
        "interpretation_note": (
            "Accessory-domain presence is an OPTIONAL-module annotation and is NEVER a "
            "pass/hold criterion; absence of these domains is NOT evidence of an "
            "intracellular protein, and presence is NOT proof of extracellular activity."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hold", required=True, type=Path,
                        help="frozen hold_candidates.tsv (71,860 rows)")
    parser.add_argument("--matrix", required=True, type=Path,
                        help="frozen per-candidate evidence matrix carrying pfam_accessions "
                             "and interpro_signatures (109,087 rows)")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--expect-targets", type=int, default=11610,
                        help="declared target-set size; a mismatch is reported, not hidden")
    args = parser.parse_args(argv)

    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise ValueError(f"refusing to write into a non-empty output dir: {args.out_dir}")

    hold_rows = read_tsv(args.hold)
    matrix_rows = read_tsv(args.matrix)
    matrix_by_accession = {
        (row.get("accession") or "").strip(): row for row in matrix_rows
    }
    if len(matrix_by_accession) != len(matrix_rows):
        raise ValueError("matrix has duplicate accessions")
    target_rows = select_sp_less_type1(hold_rows)
    joined, summary = classify_rows(target_rows, matrix_by_accession)

    mismatch = None
    if args.expect_targets is not None and len(target_rows) != args.expect_targets:
        mismatch = {
            "declared": args.expect_targets,
            "observed": len(target_rows),
            "note": "the documented 11,610 figure is not reproduced by this derivation",
        }
    summary["declared_target_count"] = args.expect_targets
    summary["target_count_mismatch"] = mismatch

    args.out_dir.mkdir(parents=True, exist_ok=True)
    fields = list(joined[0].keys())
    with (args.out_dir / "type1_sp_less_accessory_domains.tsv").open(
        "w", encoding="utf-8", newline="",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(joined)
    with (args.out_dir / "type1_sp_less_accessory_summary.json").open(
        "w", encoding="utf-8", newline="",
    ) as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({
        "target_count": summary["target_count"],
        "with_any_accessory_domain": summary["with_any_accessory_domain"],
        "rate": summary["with_any_accessory_domain_rate"],
        "per_domain_counts": summary["per_domain_counts"],
        "target_count_mismatch": mismatch,
        "out_dir": str(args.out_dir),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
