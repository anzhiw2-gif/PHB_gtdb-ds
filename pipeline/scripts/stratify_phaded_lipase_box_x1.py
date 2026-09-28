#!/usr/bin/env python3
"""Stratify PhaDED candidates by the lipase-box x1 residue (Knoll 2009 criterion).

Pure analysis: this script only *reads* an existing candidate evidence table and
writes new derived tables. It runs no HMM/profile/structural search, starts no
server job, and modifies no historical run, registry or formal scan result.

Scientific basis (Knoll M, et al. The PHA Depolymerase Engineering Database.
BMC Bioinformatics 2009;10:89, PMC2666664): the x1 position -- the residue
immediately after the catalytic serine/cysteine in the lipase box GxSxG -- is
typically a *polar* residue in lipases and esterases, whereas PHA depolymerases
almost exclusively carry a *hydrophobic* x1. A non-hydrophobic x1 is therefore a
lipase/esterase-confusion suspicion for a candidate, never a biological negative.

Evidence boundary: an x1 call, hydrophobic or not, denotes candidate homology or
functional potential only. It does not establish a validated PHB/PHA degradation
phenotype, and no candidate may be deleted on the basis of this axis alone.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path


PHENOTYPE_BOUNDARY = (
    "candidate_only; lipase-box x1 hydrophobicity is primary-sequence evidence of "
    "candidate homology or functional potential and is not a validated PHB/PHA "
    "degradation phenotype"
)
CITATION = "Knoll M, et al. BMC Bioinformatics 2009;10:89 (PhaDED; PMC2666664)"

# Knoll 2009 criterion: PHA depolymerase lipase boxes carry a hydrophobic x1.
HYDROPHOBIC_RESIDUES = frozenset("ACFILMVWY")

# Polar/charged residues -- the lipase/esterase-like side of the criterion.
NON_HYDROPHOBIC_RESIDUES = frozenset("DEHKNQRST")

# Gly and Pro have no polar side chain; they are non-hydrophobic (so they cannot
# satisfy the PHA depolymerase x1 criterion) but they are annotated separately
# instead of being pooled with the polar residues.
GLYCINE_RESIDUE = "G"
PROLINE_RESIDUE = "P"
SEPARATELY_ANNOTATED_RESIDUES = frozenset({GLYCINE_RESIDUE, PROLINE_RESIDUE})

STANDARD_AMINO_ACIDS = frozenset(
    HYDROPHOBIC_RESIDUES | NON_HYDROPHOBIC_RESIDUES | SEPARATELY_ANNOTATED_RESIDUES
)

X1_CLASSES = ("hydrophobic", "non_hydrophobic", "not_tested")

# Descriptive side-chain grouping. Only HYDROPHOBIC_RESIDUES drives the
# classification; the group column exists so that Gly and Pro stay individually
# identifiable inside the non_hydrophobic class.
_RESIDUE_GROUPS = {
    **{residue: "hydrophobic" for residue in HYDROPHOBIC_RESIDUES},
    "D": "acidic",
    "E": "acidic",
    "K": "basic",
    "R": "basic",
    "H": "basic",
    "N": "polar_uncharged",
    "Q": "polar_uncharged",
    "S": "polar_uncharged",
    "T": "polar_uncharged",
    GLYCINE_RESIDUE: "glycine",
    PROLINE_RESIDUE: "proline",
}

_RESIDUE_NAMES = {
    "A": "alanine",
    "C": "cysteine",
    "D": "aspartate",
    "E": "glutamate",
    "F": "phenylalanine",
    "G": "glycine",
    "H": "histidine",
    "I": "isoleucine",
    "K": "lysine",
    "L": "leucine",
    "M": "methionine",
    "N": "asparagine",
    "P": "proline",
    "Q": "glutamine",
    "R": "arginine",
    "S": "serine",
    "T": "threonine",
    "V": "valine",
    "W": "tryptophan",
    "Y": "tyrosine",
}

CAP_HYDROPHOBIC = "high_candidate_only"
CAP_NON_HYDROPHOBIC = "moderate_candidate_only"
CAP_NOT_TESTED = "not_assessed_by_x1_axis"

CONFOUNDER_FLAG_NON_HYDROPHOBIC = "suspected_lipase_esterase_confounder"
CONFOUNDER_FLAG_HYDROPHOBIC = "none"
CONFOUNDER_FLAG_NOT_TESTED = "not_assessed"

# No x1 class is a biological negative: "not_tested" is missing evidence, and a
# non-hydrophobic x1 is a suspicion that demands stronger discriminating
# evidence -- not an exclusion.
NEGATIVE_CLASSES: frozenset[str] = frozenset()

CAPS = {
    "hydrophobic": CAP_HYDROPHOBIC,
    "non_hydrophobic": CAP_NON_HYDROPHOBIC,
    "not_tested": CAP_NOT_TESTED,
}

CONFOUNDER_FLAGS = {
    "hydrophobic": CONFOUNDER_FLAG_HYDROPHOBIC,
    "non_hydrophobic": CONFOUNDER_FLAG_NON_HYDROPHOBIC,
    "not_tested": CONFOUNDER_FLAG_NOT_TESTED,
}

STRATIFICATION_COLUMNS = [
    "accession",
    "lipase_box_x1",
    "x1_class",
    "confounder_flag",
    "recommended_confidence_cap",
    "x1_group",
]
CONFOUNDER_COLUMNS = [
    "accession",
    "genome",
    "lipase_box_x1",
    "x1_group",
    "confounder_flag",
    "recommended_confidence_cap",
    "motif_panel_status",
]
BREAKDOWN_COLUMNS = [
    "residue",
    "residue_name",
    "x1_class",
    "x1_group",
    "count",
    "fraction_of_total",
    "fraction_of_tested",
]
SUMMARY_COLUMNS = [
    "metric",
    "count",
    "fraction_of_total",
    "fraction_of_tested",
    "interpretation",
]

REQUIRED_EVIDENCE_COLUMNS = ("accession", "lipase_box_x1")


def _residue_token(residue: object) -> str:
    return "" if residue is None else str(residue).strip()


def classify_x1(residue: object) -> str:
    """Return ``hydrophobic``/``non_hydrophobic``/``not_tested`` for one x1 value.

    Fail-closed: anything that is not a standard uppercase amino-acid letter is
    rejected rather than silently folded into a class.
    """
    token = _residue_token(residue)
    if token == "":
        return "not_tested"
    if token not in STANDARD_AMINO_ACIDS:
        raise ValueError(f"unsupported lipase_box_x1 residue: {residue!r}")
    return "hydrophobic" if token in HYDROPHOBIC_RESIDUES else "non_hydrophobic"


def x1_group(residue: object) -> str:
    """Return the descriptive side-chain group used to keep Gly/Pro distinguishable."""
    token = _residue_token(residue)
    if token == "":
        return "not_tested"
    if token not in STANDARD_AMINO_ACIDS:
        raise ValueError(f"unsupported lipase_box_x1 residue: {residue!r}")
    return _RESIDUE_GROUPS[token]


def _validated_class(x1_class: str, table: dict[str, str], what: str) -> str:
    if x1_class not in table:
        raise ValueError(f"unknown x1 class for {what}: {x1_class!r}")
    return table[x1_class]


def recommended_confidence_cap(x1_class: str) -> str:
    """Return the confidence cap that this single axis imposes."""
    return _validated_class(x1_class, CAPS, "recommended_confidence_cap")


def confounder_flag(x1_class: str) -> str:
    """Return the lipase/esterase-confusion flag for an x1 class."""
    return _validated_class(x1_class, CONFOUNDER_FLAGS, "confounder_flag")


def is_biological_negative(x1_class: str) -> bool:
    """Return whether an x1 class may count as a biological negative (never)."""
    if x1_class not in X1_CLASSES:
        raise ValueError(f"unknown x1 class: {x1_class!r}")
    return x1_class in NEGATIVE_CLASSES


def read_evidence(path: Path) -> list[dict[str, str]]:
    """Read and validate the candidate evidence table (read-only)."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        missing = [name for name in REQUIRED_EVIDENCE_COLUMNS if name not in fields]
        if missing:
            raise ValueError(f"evidence table missing required columns: {', '.join(missing)}")
        for line_number, row in enumerate(reader, start=2):
            accession = (row.get("accession") or "").strip()
            if not accession:
                raise ValueError(f"empty accession at line {line_number}")
            if accession in seen:
                raise ValueError(f"duplicate accession: {accession}")
            seen.add(accession)
            residue = _residue_token(row.get("lipase_box_x1"))
            try:
                x1_class = classify_x1(residue)
            except ValueError as error:
                raise ValueError(f"{error} (accession {accession}, line {line_number})") from error
            rows.append(
                {
                    "accession": accession,
                    "lipase_box_x1": residue,
                    "x1_class": x1_class,
                    "x1_group": x1_group(residue),
                    "confounder_flag": confounder_flag(x1_class),
                    "recommended_confidence_cap": recommended_confidence_cap(x1_class),
                    "genome": (row.get("genome") or "").strip(),
                    "motif_panel_status": (row.get("motif_panel_status") or "").strip(),
                }
            )
    if not rows:
        raise ValueError(f"evidence table has no data rows: {path}")
    return rows


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _fraction(numerator: int, denominator: int) -> str:
    if denominator <= 0:
        return ""
    return f"{numerator / denominator:.6f}"


def _write_tsv(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def stratify(evidence_path: Path, output_dir: Path) -> dict[str, object]:
    """Write the x1 stratification tables and return an auditable manifest."""
    evidence_path = Path(evidence_path)
    output_dir = Path(output_dir)
    rows = sorted(read_evidence(evidence_path), key=lambda row: row["accession"])

    total = len(rows)
    counts = Counter(row["x1_class"] for row in rows)
    hydrophobic = counts["hydrophobic"]
    non_hydrophobic = counts["non_hydrophobic"]
    not_tested = counts["not_tested"]
    tested = hydrophobic + non_hydrophobic
    negative_calls = sum(1 for row in rows if is_biological_negative(row["x1_class"]))
    confounders = [row for row in rows if row["confounder_flag"] == CONFOUNDER_FLAG_NON_HYDROPHOBIC]
    residue_counts = Counter(row["lipase_box_x1"] for row in rows)

    output_dir.mkdir(parents=True, exist_ok=True)
    _write_tsv(
        output_dir / "x1_stratification.tsv",
        STRATIFICATION_COLUMNS,
        [{column: row[column] for column in STRATIFICATION_COLUMNS} for row in rows],
    )
    _write_tsv(
        output_dir / "confounder_candidates.tsv",
        CONFOUNDER_COLUMNS,
        [{column: row[column] for column in CONFOUNDER_COLUMNS} for row in confounders],
    )

    breakdown_rows = []
    for residue in sorted(STANDARD_AMINO_ACIDS):
        count = residue_counts[residue]
        if count == 0:
            continue
        breakdown_rows.append(
            {
                "residue": residue,
                "residue_name": _RESIDUE_NAMES[residue],
                "x1_class": classify_x1(residue),
                "x1_group": x1_group(residue),
                "count": count,
                "fraction_of_total": _fraction(count, total),
                "fraction_of_tested": _fraction(count, tested),
            }
        )
    missing_count = residue_counts[""]
    if missing_count:
        breakdown_rows.append(
            {
                "residue": "not_tested",
                "residue_name": "missing_lipase_box_x1",
                "x1_class": "not_tested",
                "x1_group": "not_tested",
                "count": missing_count,
                "fraction_of_total": _fraction(missing_count, total),
                "fraction_of_tested": "",
            }
        )
    _write_tsv(output_dir / "x1_residue_breakdown.tsv", BREAKDOWN_COLUMNS, breakdown_rows)

    summary_rows = [
        {
            "metric": "total_candidates",
            "count": total,
            "fraction_of_total": _fraction(total, total),
            "fraction_of_tested": "",
            "interpretation": "candidate universe rows consumed from the evidence table",
        },
        {
            "metric": "tested_count",
            "count": tested,
            "fraction_of_total": _fraction(tested, total),
            "fraction_of_tested": _fraction(tested, tested),
            "interpretation": "candidates with a measured lipase-box x1 residue (denominator for x1 shares)",
        },
        {
            "metric": "hydrophobic",
            "count": hydrophobic,
            "fraction_of_total": _fraction(hydrophobic, total),
            "fraction_of_tested": _fraction(hydrophobic, tested),
            "interpretation": "x1 in A,C,F,I,L,M,V,W,Y; consistent with the PHA depolymerase x1 criterion",
        },
        {
            "metric": "non_hydrophobic",
            "count": non_hydrophobic,
            "fraction_of_total": _fraction(non_hydrophobic, total),
            "fraction_of_tested": _fraction(non_hydrophobic, tested),
            "interpretation": "x1 polar/charged or Gly/Pro; lipase/esterase confusion suspicion, not a negative",
        },
        {
            "metric": "not_tested",
            "count": not_tested,
            "fraction_of_total": _fraction(not_tested, total),
            "fraction_of_tested": "",
            "interpretation": "no lipase-box x1 measured; missing evidence, never counted as a negative",
        },
        {
            "metric": "hydrophobic_share_of_tested",
            "count": "",
            "fraction_of_total": "",
            "fraction_of_tested": _fraction(hydrophobic, tested),
            "interpretation": "share of tested candidates whose x1 is hydrophobic",
        },
        {
            "metric": "non_hydrophobic_share_of_tested",
            "count": "",
            "fraction_of_total": "",
            "fraction_of_tested": _fraction(non_hydrophobic, tested),
            "interpretation": "share of tested candidates whose x1 is non-hydrophobic",
        },
        {
            "metric": "confounder_flagged_count",
            "count": len(confounders),
            "fraction_of_total": _fraction(len(confounders), total),
            "fraction_of_tested": _fraction(len(confounders), tested),
            "interpretation": "rows written to confounder_candidates.tsv",
        },
        {
            "metric": "negative_call_count",
            "count": negative_calls,
            "fraction_of_total": _fraction(negative_calls, total),
            "fraction_of_tested": _fraction(negative_calls, tested),
            "interpretation": "x1 axis issues zero negative calls: not_tested is missing evidence and a "
            "non-hydrophobic x1 is a suspicion, so no candidate is excluded here",
        },
    ]
    _write_tsv(output_dir / "x1_summary.tsv", SUMMARY_COLUMNS, summary_rows)

    outputs = {}
    for name in (
        "x1_stratification.tsv",
        "x1_summary.tsv",
        "confounder_candidates.tsv",
        "x1_residue_breakdown.tsv",
    ):
        path = output_dir / name
        outputs[name] = {
            "path": str(path.resolve()),
            "size": path.stat().st_size,
            "sha256": _sha256(path),
        }

    manifest: dict[str, object] = {
        "schema_version": "1.0",
        "status": "completed_candidate_only",
        "analysis": "lipase_box_x1_hydrophobicity_stratification",
        "criterion": {
            "source": CITATION,
            "hydrophobic_x1": "".join(sorted(HYDROPHOBIC_RESIDUES)),
            "non_hydrophobic_x1": "".join(sorted(NON_HYDROPHOBIC_RESIDUES)),
            "annotated_separately": [GLYCINE_RESIDUE, PROLINE_RESIDUE],
            "glycine_group": x1_group(GLYCINE_RESIDUE),
            "proline_group": x1_group(PROLINE_RESIDUE),
        },
        "input": {
            "path": str(evidence_path.resolve()),
            "size": evidence_path.stat().st_size,
            "sha256": _sha256(evidence_path),
            "rows_read": total,
            "unique_accessions": len({row["accession"] for row in rows}),
        },
        "total_candidates": total,
        "tested_count": tested,
        "hydrophobic_count": hydrophobic,
        "non_hydrophobic_count": non_hydrophobic,
        "not_tested_count": not_tested,
        "confounder_count": len(confounders),
        "negative_call_count": negative_calls,
        "hydrophobic_share_of_tested": _fraction(hydrophobic, tested),
        "non_hydrophobic_share_of_tested": _fraction(non_hydrophobic, tested),
        "residue_counts": {key: value for key, value in sorted(residue_counts.items())},
        "confidence_caps": {
            "hydrophobic": CAP_HYDROPHOBIC,
            "non_hydrophobic": CAP_NON_HYDROPHOBIC,
            "not_tested": CAP_NOT_TESTED,
            "note": "axis-local cap; the effective cap for a candidate is the minimum across axes",
        },
        "totals_self_consistent": hydrophobic + non_hydrophobic + not_tested == total,
        "outputs": outputs,
        "phenotype_boundary": PHENOTYPE_BOUNDARY,
        "no_candidate_deletion": "No candidate is removed, folded or reclassified as negative by this axis.",
    }
    (output_dir / "x1_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def write_residue_policy(path: Path) -> None:
    """Write the versioned x1 residue policy table so the criterion is auditable."""
    rows: list[dict[str, object]] = []
    for residue in sorted(STANDARD_AMINO_ACIDS):
        note = ""
        if residue == GLYCINE_RESIDUE:
            note = "no side chain; non-hydrophobic, annotated separately from polar residues"
        elif residue == PROLINE_RESIDUE:
            note = "cyclic secondary amine; non-hydrophobic, annotated separately from polar residues"
        rows.append(
            {
                "residue": residue,
                "residue_name": _RESIDUE_NAMES[residue],
                "x1_class": classify_x1(residue),
                "x1_group": x1_group(residue),
                "hydrophobic_set_member": "yes" if residue in HYDROPHOBIC_RESIDUES else "no",
                "note": note,
                "citation": CITATION,
            }
        )
    rows.append(
        {
            "residue": "not_tested",
            "residue_name": "missing_lipase_box_x1",
            "x1_class": "not_tested",
            "x1_group": "not_tested",
            "hydrophobic_set_member": "no",
            "note": "no x1 measured; missing evidence, never a negative",
            "citation": CITATION,
        }
    )
    _write_tsv(
        path,
        [
            "residue",
            "residue_name",
            "x1_class",
            "x1_group",
            "hydrophobic_set_member",
            "note",
            "citation",
        ],
        rows,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument(
        "--create-run",
        action="store_true",
        help="create the run layout with pipeline/scripts/run_context.py before writing",
    )
    args = parser.parse_args(argv)

    run_dir = args.run_dir.resolve()
    if run_dir.parent.name != "runs":
        raise SystemExit(f"run directory must live under a runs/ directory: {run_dir}")

    scripts_dir = Path(__file__).resolve().parent
    if str(scripts_dir) not in sys.path:
        sys.path.insert(0, str(scripts_dir))
    from run_context import create_run_layout, sha256_file, write_input_contract  # noqa: E402

    if args.create_run:
        create_run_layout(run_dir.parent.parent, run_dir.name)
    for name in ("logs", "inputs", "results"):
        path = run_dir / name
        if not path.is_dir():
            raise SystemExit(f"run directory is incomplete, missing {name}/: {run_dir}")

    policy_path = run_dir / "inputs" / "x1_residue_policy.tsv"
    write_residue_policy(policy_path)

    manifest = stratify(args.evidence, run_dir / "results")

    log_lines = [
        f"run_id: {run_dir.name}",
        f"generated_at: {dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')}",
        f"evidence: {Path(args.evidence).resolve()}",
        f"evidence_sha256: {manifest['input']['sha256']}",
        f"total_candidates: {manifest['total_candidates']}",
        f"tested_count: {manifest['tested_count']}",
        f"hydrophobic_count: {manifest['hydrophobic_count']}",
        f"non_hydrophobic_count: {manifest['non_hydrophobic_count']}",
        f"not_tested_count: {manifest['not_tested_count']}",
        f"negative_call_count: {manifest['negative_call_count']}",
        f"totals_self_consistent: {manifest['totals_self_consistent']}",
        "server_execution_started: false",
    ]
    (run_dir / "logs" / "x1_stratification.log").write_text(
        "\n".join(log_lines) + "\n", encoding="utf-8"
    )

    run_manifest = {
        "run_id": run_dir.name,
        "status": "completed_candidate_only",
        "purpose": "x1 (lipase box) hydrophobicity stratification; pure analysis, zero new computation",
        "analysis": manifest["analysis"],
        "criteria_source": CITATION,
        "counts": {
            "total_candidates": manifest["total_candidates"],
            "tested_count": manifest["tested_count"],
            "hydrophobic_count": manifest["hydrophobic_count"],
            "non_hydrophobic_count": manifest["non_hydrophobic_count"],
            "not_tested_count": manifest["not_tested_count"],
            "confounder_count": manifest["confounder_count"],
            "negative_call_count": manifest["negative_call_count"],
        },
        "residue_counts": manifest["residue_counts"],
        "authorization": {
            "candidate_only_execution": True,
            "server_execution_started": False,
            "formal_scan_started": False,
            "formal_registry_modified": False,
            "historical_run_modified": False,
        },
        "source_snapshot": {
            "stratification_script": {
                "path": str(Path(__file__).resolve()),
                "size": Path(__file__).resolve().stat().st_size,
                "sha256": sha256_file(Path(__file__).resolve()),
            },
            "run_context": {
                "path": str((scripts_dir / "run_context.py").resolve()),
                "size": (scripts_dir / "run_context.py").stat().st_size,
                "sha256": sha256_file(scripts_dir / "run_context.py"),
            },
        },
        "outputs": manifest["outputs"],
        "phenotype_boundary": PHENOTYPE_BOUNDARY,
        "non_modification_statement": (
            "Read-only consumption of the 20260915 motif reconciliation evidence table; no input was "
            "rewritten and no historical run, registry or formal scan artifact was touched."
        ),
    }
    (run_dir / "run_manifest.json").write_text(
        json.dumps(run_manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    contract = write_input_contract(
        run_dir,
        run_id=run_dir.name,
        inputs={
            "motif_candidate_evidence": Path(args.evidence).resolve(),
            "x1_residue_policy": policy_path,
        },
    )

    print(
        json.dumps(
            {
                "run_id": run_dir.name,
                "total_candidates": manifest["total_candidates"],
                "tested_count": manifest["tested_count"],
                "hydrophobic_count": manifest["hydrophobic_count"],
                "non_hydrophobic_count": manifest["non_hydrophobic_count"],
                "not_tested_count": manifest["not_tested_count"],
                "negative_call_count": manifest["negative_call_count"],
                "totals_self_consistent": manifest["totals_self_consistent"],
                "input_contract_status": contract["status"],
                "status": "completed_candidate_only",
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
