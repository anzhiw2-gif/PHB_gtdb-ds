#!/usr/bin/env python3
"""Calibrate the PhaDED secretion-signal gate against the 723-reference panel.

Decision question (packet P1 of the 2026-09-28 redesign plan): the candidate
layer HOLDS a type-1 candidate unless SignalP predicts a secretion signal
(measured: 48,038 of 67,720 type-1 candidates = 71% are held on that rule
alone), while the reference layer was never tested by it -- the "extracellular"
label on a reference is DEFINITIONAL, inherited from the historical DED
superfamily name rather than measured.  This module measures how often SignalP
actually predicts a secretion signal for a reference whose extracellular label
is inherited, stratified by superfamily, by evidence grade and by whether an
experimental localization is documented.

Candidate-only boundaries: a SignalP prediction is transport-signal evidence and
is NOT experimental localization, and none of this proves a PHB/PHA degradation
phenotype.  Nothing here deletes, demotes or promotes a candidate.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

#: Historical DED superfamily -> the localization its NAME implies.  This is an
#: inherited label, not a measurement; the whole point of this module is to test
#: how well the inherited label survives an actual prediction.
SUPERFAMILY_LOCALIZATION = {
    "intracellular nPHASCL without lipase box": "intracellular",
    "intracellular nPHASCL with lipase box": "intracellular",
    "intracellular nPHAMCL": "intracellular",
    "extracellular dPHASCL type 1": "extracellular",
    "extracellular dPHASCL type 2": "extracellular",
    "extracellular dPHAMCL": "extracellular",
    "extracellular native-SCL/PhaZ7-like": "extracellular",
    "periplasmic PHA depolymerases": "periplasmic",
}

EXPORT_PREDICTIONS = frozenset({"SP", "LIPO", "TAT", "TATLIPO"})
NON_EXPORT_PREDICTIONS = frozenset({"OTHER", "PILIN"})

#: SignalP column layout of the ``--format txt`` prediction_results.txt.
PREDICTION_COLUMNS = (
    "OTHER", "SP(Sec/SPI)", "LIPO(Sec/SPII)", "TAT(Tat/SPI)",
    "TATLIPO(Tat/SPII)", "PILIN(Sec/SPIII)",
)


def parse_prediction_results(path: Path) -> dict[str, dict[str, str]]:
    """Read SignalP ``prediction_results.txt`` into {reference_id: record}.

    The project FASTA header convention is ``>reference_id|accession``, so the
    ID column is split on ``|`` and the record is keyed by ``reference_id``.
    """
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"SignalP prediction file is not a regular file: {path}")
    records: dict[str, dict[str, str]] = {}
    header: list[str] | None = None
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for line in handle:
            line = line.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            if line.startswith("#"):
                # The second comment line carries the column header.
                fields = [field.strip() for field in line.lstrip("#").split("\t")]
                if "ID" in fields and "Prediction" in fields:
                    header = fields
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(f"malformed SignalP row: {line!r}")
            raw_id = fields[0].strip()
            if not raw_id:
                raise ValueError(f"SignalP row without an ID: {line!r}")
            reference_id = raw_id.split("|")[0]
            if reference_id in records:
                raise ValueError(f"duplicate SignalP prediction for {reference_id}")
            record = {"reference_id": reference_id, "signalp_id": raw_id}
            if header is not None:
                for name, value in zip(header[1:], fields[1:]):
                    record[name] = value
            else:  # pragma: no cover - header is always present in SignalP 6 output
                record["Prediction"] = fields[1]
            records[reference_id] = record
    if not records:
        raise ValueError(f"no SignalP predictions parsed from {path}")
    return records


def read_reference_rows(path: Path) -> list[dict[str, str]]:
    """Read the curated reference ledger (one row per reference_id)."""
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"reference ledger is not a regular file: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"empty reference ledger: {path}")
    return rows


def historical_localization(superfamily: str) -> str:
    """Localization implied by the historical superfamily NAME (inherited)."""
    return SUPERFAMILY_LOCALIZATION.get(superfamily.strip(), "unmapped")


def read_experimental_localization(path: Path | None) -> dict[str, str]:
    """Documented (experimental) localization per reference_id, when available.

    The provenance amendment records ``reported_localization`` for the
    experimental positives.  Absent file -> empty map; the analysis then reports
    the inherited stratum only.
    """
    if path is None:
        return {}
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"provenance amendment is not a regular file: {path}")
    mapping: dict[str, str] = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            reference_id = (row.get("reference_id") or "").strip()
            reported = (row.get("reported_localization") or "").strip()
            if reference_id and reported:
                mapping[reference_id] = reported
    return mapping


def classify_rows(
    ledger_rows: list[dict[str, str]],
    predictions: dict[str, dict[str, str]],
    experimental_localization: dict[str, str] | None = None,
) -> tuple[list[dict[str, str]], dict]:
    """Join the ledger with the SignalP predictions and compute the strata."""
    experimental_localization = experimental_localization or {}
    joined: list[dict[str, str]] = []
    missing_prediction: list[str] = []
    unmapped_superfamily: list[str] = []
    for row in ledger_rows:
        reference_id = (row.get("reference_id") or "").strip()
        superfamily = (row.get("phaded_superfamily") or "").strip()
        inherited = historical_localization(superfamily)
        if inherited == "unmapped":
            unmapped_superfamily.append(reference_id)
        prediction_record = predictions.get(reference_id)
        if prediction_record is None:
            missing_prediction.append(reference_id)
            prediction = ""
            export = ""
        else:
            prediction = (prediction_record.get("Prediction") or "").strip()
            if prediction in EXPORT_PREDICTIONS:
                export = "true"
            elif prediction in NON_EXPORT_PREDICTIONS:
                export = "false"
            else:
                export = "unrecognised"
        joined.append({
            "reference_id": reference_id,
            "accession": (row.get("accession") or "").strip(),
            "phaded_superfamily": superfamily,
            "phaded_family_id": (row.get("phaded_family_id") or "").strip(),
            "inherited_localization": inherited,
            "documented_localization": experimental_localization.get(reference_id, ""),
            "experimental_evidence_grade": (row.get("experimental_evidence_grade") or "").strip(),
            "sequence_integrity": (row.get("sequence_integrity") or "").strip(),
            "signalp_prediction": prediction,
            "has_export_signal": export,
        })
    summary = summarize(joined)
    summary["rows_without_prediction"] = missing_prediction
    summary["rows_with_unmapped_superfamily"] = unmapped_superfamily
    return joined, summary


def _rate(numerator: int, denominator: int) -> float | None:
    if denominator <= 0:
        return None
    return numerator / denominator


def summarize(joined: list[dict[str, str]]) -> dict:
    """Cross-tabulate prediction class against the inherited label."""
    by_inherited: dict[str, dict[str, int]] = {}
    by_prediction: dict[str, int] = {}
    crosstab: dict[str, dict[str, int]] = {}
    for row in joined:
        inherited = row["inherited_localization"]
        prediction = row["signalp_prediction"] or "NO_PREDICTION"
        export = row["has_export_signal"]
        by_prediction[prediction] = by_prediction.get(prediction, 0) + 1
        bucket = by_inherited.setdefault(
            inherited, {"total": 0, "export_signal": 0, "no_export_signal": 0},
        )
        bucket["total"] += 1
        if export == "true":
            bucket["export_signal"] += 1
        elif export == "false":
            bucket["no_export_signal"] += 1
        crosstab.setdefault(inherited, {})
        crosstab[inherited][prediction] = crosstab[inherited].get(prediction, 0) + 1

    for inherited, bucket in by_inherited.items():
        bucket["export_signal_rate"] = _rate(bucket["export_signal"], bucket["total"])

    extracellular = by_inherited.get("extracellular", {"total": 0, "export_signal": 0})
    intracellular = by_inherited.get("intracellular", {"total": 0, "export_signal": 0})
    return {
        "total_references": len(joined),
        "prediction_class_counts": by_prediction,
        "by_inherited_localization": by_inherited,
        "crosstab_inherited_x_prediction": crosstab,
        "extracellular_export_signal_rate": _rate(
            extracellular["export_signal"], extracellular["total"],
        ),
        "intracellular_export_signal_rate": _rate(
            intracellular["export_signal"], intracellular["total"],
        ),
        "specificity_note": (
            "A secretion-signal prediction is transport evidence, not experimental "
            "localization; an absent prediction is NOT evidence of an intracellular "
            "protein (Knoll 2009 records an extracellular type-1 family member that "
            "lacks a signal peptide)."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", required=True, type=Path,
                        help="curated reference ledger (723 rows)")
    parser.add_argument("--signalp", required=True, type=Path,
                        help="SignalP prediction_results.txt")
    parser.add_argument("--provenance", type=Path, default=None,
                        help="optional positive provenance amendment (reported_localization)")
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)

    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise ValueError(f"refusing to write into a non-empty output dir: {args.out_dir}")
    ledger_rows = read_reference_rows(args.ledger)
    predictions = parse_prediction_results(args.signalp)
    experimental = read_experimental_localization(args.provenance)
    joined, summary = classify_rows(ledger_rows, predictions, experimental)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    fields = list(joined[0].keys())
    with (args.out_dir / "reference_signalp_stratification.tsv").open(
        "w", encoding="utf-8", newline="",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(joined)
    with (args.out_dir / "reference_signalp_calibration.json").open(
        "w", encoding="utf-8", newline="",
    ) as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({
        "rows": len(joined),
        "extracellular_export_signal_rate": summary["extracellular_export_signal_rate"],
        "intracellular_export_signal_rate": summary["intracellular_export_signal_rate"],
        "prediction_class_counts": summary["prediction_class_counts"],
        "out_dir": str(args.out_dir),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
