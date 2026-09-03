#!/usr/bin/env python3
"""Calibrate the two ePhaZ HMM layers against accessioned controls.

This intentionally does not fall back to the historical ``ePhaZ.hmm`` model.
Both layer HMMs and a layer manifest are required so the calibration record can
be traced back to the exact seed decision and evidence level.
"""

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


FAMILIES = {
    "ePhaZ_curated_core": "curated_hmm",
    "ePhaZ_broad_discovery": "broad_hmm",
}
DEFAULT_THRESHOLDS = (1e-2, 1e-3, 1e-5, 1e-8, 1e-10, 1e-15, 1e-20)
EPHAZ_PREFIX = "e-phaz"
UNKNOWN_EVIDENCE = {"", "pending", "unknown", "none", "na", "n/a"}


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path):
    candidate = Path(path).resolve()
    if not candidate.is_file() or candidate.stat().st_size == 0:
        raise ValueError(f"required non-empty file missing: {path}")
    return {"path": str(candidate), "size": candidate.stat().st_size, "sha256": sha256(candidate)}


def read_tsv(path, required):
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        missing = set(required) - fields
        if missing:
            raise ValueError(f"{path} missing columns: {sorted(missing)}")
        rows = list(reader)
    return rows


def read_fasta(path):
    records = {}
    header = None
    sequence = []
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    accession = header.split("|", 1)[0].strip()
                    if not accession or accession in records or not sequence:
                        raise ValueError(f"invalid or duplicate FASTA accession in {path}: {accession!r}")
                    records[accession] = (header, "".join(sequence))
                header = line[1:].strip()
                sequence = []
            else:
                if header is None:
                    raise ValueError(f"FASTA sequence precedes header in {path}")
                sequence.append(line)
    if header is not None:
        accession = header.split("|", 1)[0].strip()
        if not accession or accession in records or not sequence:
            raise ValueError(f"invalid or duplicate FASTA accession in {path}: {accession!r}")
        records[accession] = (header, "".join(sequence))
    if not records:
        raise ValueError(f"empty FASTA: {path}")
    return records


def _manifest_layer(row):
    return (row.get("layer") or row.get("assigned_layer") or row.get("target_layer") or "").strip()


def _manifest_evidence(row):
    return (row.get("evidence") or row.get("evidence_level") or row.get("evidence_status") or "").strip()


def read_layer_manifest(path):
    rows = read_tsv(path, {"accession"})
    mapping = {}
    for row in rows:
        accession = row["accession"].strip()
        if not accession or accession in mapping:
            raise ValueError(f"duplicate/empty accession in layer manifest: {accession!r}")
        mapping[accession] = {
            "layer": _manifest_layer(row),
            "evidence": _manifest_evidence(row),
            "architecture": (row.get("architecture") or row.get("architecture_status") or "").strip(),
            "length": (row.get("length") or "").strip(),
            "decision_reason": (row.get("decision_reason") or row.get("reason") or "").strip(),
        }
    return mapping


def _family_for_control(row):
    group = (row.get("query_group") or "").strip().lower()
    return group.startswith(EPHAZ_PREFIX)


def _has_explicit_evidence(value):
    normalized = value.strip().lower()
    return normalized not in UNKNOWN_EVIDENCE and not normalized.startswith("pending")


def parse_thresholds(value):
    try:
        thresholds = tuple(float(item) for item in value.split(",") if item.strip())
    except ValueError as exc:
        raise ValueError(f"invalid --thresholds: {value}") from exc
    if not thresholds or any(item <= 0 for item in thresholds):
        raise ValueError("--thresholds must contain positive numbers")
    return thresholds


def parse_hits(tbl_path, dom_path):
    hits = {}
    with open(tbl_path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            if len(fields) < 6:
                continue
            try:
                hits[fields[0].split("|", 1)[0]] = {"E": float(fields[4]), "cov": 0.0}
            except ValueError:
                continue
    with open(dom_path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            if len(fields) < 21:
                continue
            try:
                accession = fields[0].split("|", 1)[0]
                qlen = float(fields[5])
                # domtblout columns 20/21 (zero-based 19/20) are HMM
                # envelope coordinates; fields 15/16 are E-values.
                coverage = (int(fields[20]) - int(fields[19]) + 1) / qlen
            except (IndexError, ValueError, ZeroDivisionError):
                continue
            if accession in hits:
                hits[accession]["cov"] = max(hits[accession]["cov"], min(1.0, coverage))
    return hits


def run_hmmsearch(hmm, probe, workdir, hmmsearch_bin, cpu):
    os.makedirs(workdir, exist_ok=True)
    tbl = os.path.join(workdir, "hits.tbl")
    dom = os.path.join(workdir, "hits.dom")
    subprocess.run(
        [hmmsearch_bin, "--tblout", tbl, "--domtblout", dom, "-E", "1e-2", "--cpu", str(cpu), hmm, probe],
        check=True,
        capture_output=True,
        text=True,
    )
    if not os.path.isfile(tbl) or not os.path.isfile(dom):
        raise RuntimeError(f"hmmsearch did not produce both outputs for {hmm}")
    return parse_hits(tbl, dom)


def _metrics(pos_acc, neg_acc, hits, threshold, min_cov):
    def detected(accession):
        value = hits.get(accession)
        return value is not None and value["E"] <= threshold and value.get("cov", 0.0) >= min_cov

    tp = sum(detected(item) for item in pos_acc)
    fn = len(pos_acc) - tp
    fp = sum(detected(item) for item in neg_acc)
    tn = len(neg_acc) - fp
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    denominator = ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5
    mcc = (tp * tn - fp * fn) / denominator if denominator else 0.0
    return {
        "TP": tp, "FP": fp, "FN": fn, "TN": tn,
        "precision": round(precision, 6), "recall": round(recall, 6),
        "F1": round(f1, 6), "MCC": round(mcc, 6),
    }


def _write_probe(path, accessions, records):
    with open(path, "w", encoding="utf-8") as handle:
        for accession in accessions:
            header, sequence = records[accession]
            handle.write(f">{header}\n{sequence}\n")


def calibrate(curated_hmm, broad_hmm, controls, layer_manifest, outdir, thresholds=DEFAULT_THRESHOLDS, min_cov=0.0, hmmsearch_bin=None, cpu=4):
    curated_info = file_record(curated_hmm)
    broad_info = file_record(broad_hmm)
    if Path(curated_hmm).resolve() == Path(broad_hmm).resolve():
        raise ValueError("curated and broad HMMs must be distinct explicit files")
    controls_root = Path(controls).resolve()
    controls_tsv = controls_root / "controls.tsv"
    positive_faa = controls_root / "positive.faa"
    negative_faa = controls_root / "negative.faa"
    control_rows = read_tsv(controls_tsv, {"accession", "label", "query_group"})
    control_accessions = [row["accession"].strip() for row in control_rows]
    if any(not accession for accession in control_accessions):
        raise ValueError("controls.tsv contains an empty accession")
    if len(control_accessions) != len(set(control_accessions)):
        raise ValueError("controls.tsv contains duplicate accession")
    control_accessions = [row["accession"].strip() for row in control_rows]
    if any(not accession for accession in control_accessions):
        raise ValueError("controls.tsv contains an empty accession")
    if len(control_accessions) != len(set(control_accessions)):
        raise ValueError("controls.tsv contains duplicate accessions")
    manifest = read_layer_manifest(layer_manifest)
    pos_fasta = read_fasta(positive_faa)
    neg_fasta = read_fasta(negative_faa)
    positives = [row for row in control_rows if row["label"].strip().lower() == "positive" and _family_for_control(row)]
    negatives = [row for row in control_rows if row["label"].strip().lower() == "negative"]
    broad_rows = positives
    curated_rows = []
    for row in positives:
        info = manifest.get(row["accession"].strip())
        if info and info["layer"] == "ePhaZ_curated_core" and _has_explicit_evidence(info["evidence"]):
            curated_rows.append(row)
    if not broad_rows:
        raise ValueError("no ePhaZ positive controls found")
    if not curated_rows:
        raise ValueError("no evidence-backed ePhaZ_curated_core positive controls found")
    for row in broad_rows + negatives:
        accession = row["accession"].strip()
        records = pos_fasta if row["label"].strip().lower() == "positive" else neg_fasta
        if accession not in records:
            raise ValueError(f"control accession missing from FASTA: {accession}")

    hmmsearch_bin = hmmsearch_bin or shutil.which("hmmsearch")
    if not hmmsearch_bin:
        raise ValueError("hmmsearch executable not found; pass --hmmsearch-bin explicitly")
    os.makedirs(outdir, exist_ok=True)
    thresholds = tuple(thresholds)
    negative_accessions = [row["accession"].strip() for row in negatives]
    layer_rows = {
        "ePhaZ_curated_core": curated_rows,
        "ePhaZ_broad_discovery": broad_rows,
    }
    hmm_paths = {
        "ePhaZ_curated_core": curated_hmm,
        "ePhaZ_broad_discovery": broad_hmm,
    }
    fasta_by_label = {"positive": pos_fasta, "negative": neg_fasta}
    summary_rows = []
    hit_rows = []
    layer_metadata = {}
    with tempfile.TemporaryDirectory(prefix="ephaz_calibration_", dir=outdir) as work_root:
        for family, selected_rows in layer_rows.items():
            positive_accessions = [row["accession"].strip() for row in selected_rows]
            probe = os.path.join(work_root, f"{family}.faa")
            _write_probe(probe, positive_accessions, pos_fasta)
            _write_probe(probe + ".neg", negative_accessions, neg_fasta)
            with open(probe, "a", encoding="utf-8") as handle, open(probe + ".neg", encoding="utf-8") as negative_handle:
                handle.write(negative_handle.read())
            hits = run_hmmsearch(hmm_paths[family], probe, os.path.join(work_root, family), hmmsearch_bin, cpu)
            for accession, values in sorted(hits.items()):
                source = next((row for row in control_rows if row["accession"].strip() == accession), {})
                manifest_info = manifest.get(accession, {})
                hit_rows.append({
                    "family": family,
                    "accession": accession,
                    "is_positive": accession in positive_accessions,
                    "E_value": f"{values['E']:.6g}",
                    "cov": f"{values.get('cov', 0.0):.6f}",
                    "evidence": manifest_info.get("evidence", ""),
                    "control_layer": manifest_info.get("layer", ""),
                    "reviewed": source.get("reviewed", ""),
                })
            for threshold in thresholds:
                row = _metrics(positive_accessions, negative_accessions, hits, threshold, min_cov)
                summary_rows.append({"family": family, "threshold": f"{threshold:.0e}", "min_cov": f"{min_cov:.3f}", "positive_count": len(positive_accessions), "negative_count": len(negative_accessions), **row})
            layer_metadata[family] = {
                "hmm": file_record(hmm_paths[family]),
                "positive_controls": [
                    {
                        "accession": row["accession"].strip(),
                        "evidence": manifest.get(row["accession"].strip(), {}).get("evidence", ""),
                        "layer": manifest.get(row["accession"].strip(), {}).get("layer", ""),
                        "reviewed": row.get("reviewed", ""),
                    }
                    for row in selected_rows
                ],
                "negative_control_count": len(negative_accessions),
            }

    summary_path = os.path.join(outdir, "calibration_summary.tsv")
    summary_columns = ["family", "threshold", "min_cov", "positive_count", "negative_count", "TP", "FP", "FN", "TN", "precision", "recall", "F1", "MCC"]
    with open(summary_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)
    hits_path = os.path.join(outdir, "calibration_hits.tsv")
    hit_columns = ["family", "accession", "is_positive", "E_value", "cov", "evidence", "control_layer", "reviewed"]
    with open(hits_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=hit_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(hit_rows)
    metadata = {
        "schema_version": 1,
        "tool": "calibrate_ephaz_layers.py",
        "parameters": {"thresholds": list(thresholds), "min_cov": min_cov, "hmmsearch_bin": hmmsearch_bin, "cpu": cpu},
        "inputs": {
            "curated_hmm": curated_info,
            "broad_hmm": broad_info,
            "controls_tsv": file_record(controls_tsv),
            "positive_faa": file_record(positive_faa),
            "negative_faa": file_record(negative_faa),
            "layer_manifest": file_record(layer_manifest),
        },
        "layers": layer_metadata,
        "outputs": {"calibration_summary.tsv": file_record(summary_path), "calibration_hits.tsv": file_record(hits_path)},
    }
    metadata_path = os.path.join(outdir, "calibration_metadata.json")
    with open(metadata_path, "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    return metadata


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--curated-hmm", required=True)
    parser.add_argument("--broad-hmm", required=True)
    parser.add_argument("--controls", required=True, help="directory containing controls.tsv, positive.faa, negative.faa")
    parser.add_argument("--layer-manifest", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--thresholds", default=",".join(str(item) for item in DEFAULT_THRESHOLDS))
    parser.add_argument("--min-cov", type=float, default=0.0)
    parser.add_argument("--hmmsearch-bin", default=None)
    parser.add_argument("--cpu", type=int, default=4)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.min_cov < 0 or args.min_cov > 1:
        raise SystemExit("--min-cov must be between 0 and 1")
    if args.cpu < 1:
        raise SystemExit("--cpu must be positive")
    try:
        metadata = calibrate(
            args.curated_hmm,
            args.broad_hmm,
            args.controls,
            args.layer_manifest,
            args.outdir,
            thresholds=parse_thresholds(args.thresholds),
            min_cov=args.min_cov,
            hmmsearch_bin=args.hmmsearch_bin,
            cpu=args.cpu,
        )
    except (OSError, ValueError, RuntimeError) as exc:
        raise SystemExit(str(exc)) from exc
    print(f"layered calibration complete: {args.outdir} ({', '.join(metadata['layers'])})")


if __name__ == "__main__":
    main()
