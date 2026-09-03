#!/usr/bin/env python3
"""Run a bounded external ePhaZ panel calibration with separate denominators."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path


PANEL_FILES = {
    "PHB_experimental": "independent_experimental_positive.faa",
    "MCL-PHA_experimental": "mcl_pha_experimental_positive.faa",
    "intracellular_PHB": "intracellular_non_ephaz_negative.faa",
    "MCL-PHA_non_PHB": "mcl_pha_non_phb_negative.faa",
    "annotation_only": "annotation_only_near_neighbor_negative.faa",
    "fragment_challenge": "fragment_or_incomplete_negative.faa",
    "excluded_challenge": "ephaz_excluded_challenge.faa",
}
FORMAL_NEGATIVE = {"intracellular_PHB", "MCL-PHA_non_PHB", "annotation_only"}
POSITIVE = {"PHB_experimental", "MCL-PHA_experimental"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                accession = header.split("|", 1)[0].split()[0]
                if accession in records or not chunks:
                    raise ValueError(f"invalid or duplicate FASTA accession: {accession}")
                records[accession] = (header, "".join(chunks))
            header, chunks = line[1:], []
        elif header is None:
            raise ValueError(f"sequence precedes FASTA header in {path}")
        else:
            chunks.append(line.replace("-", "").replace(".", "").upper())
    if header is not None:
        accession = header.split("|", 1)[0].split()[0]
        if accession in records or not chunks:
            raise ValueError(f"invalid or duplicate FASTA accession: {accession}")
        records[accession] = (header, "".join(chunks))
    if not records:
        raise ValueError(f"empty FASTA: {path}")
    return records


def parse_tblout(lines: list[str]) -> dict[str, dict[str, float]]:
    hits: dict[str, dict[str, float]] = {}
    for raw in lines:
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        fields = raw.split()
        if len(fields) < 6:
            continue
        try:
            accession = fields[0].split("|", 1)[0]
            evalue, bitscore = float(fields[4]), float(fields[5])
        except (ValueError, IndexError):
            continue
        current = hits.get(accession)
        if current is None or bitscore > current["bitscore"]:
            hits[accession] = {"evalue": evalue, "bitscore": bitscore}
    return hits


def parse_domtblout(lines: list[str], lengths: dict[str, int]) -> dict[str, dict[str, float]]:
    coverage: dict[str, dict[str, float]] = {}
    for raw in lines:
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        fields = raw.split()
        if len(fields) < 21:
            continue
        try:
            accession = fields[0].split("|", 1)[0]
            target_length = lengths[accession]
            hmm_length = float(fields[5])
            hmm_from, hmm_to = int(fields[15]), int(fields[16])
            target_from, target_to = int(fields[17]), int(fields[18])
        except (KeyError, ValueError, IndexError, ZeroDivisionError):
            continue
        item = {
            "hmm_coverage": min(1.0, (hmm_to - hmm_from + 1) / max(hmm_length, 1.0)),
            "target_coverage": min(1.0, (target_to - target_from + 1) / max(target_length, 1)),
        }
        if accession not in coverage or item["hmm_coverage"] > coverage[accession]["hmm_coverage"]:
            coverage[accession] = item
    return coverage


def summarize_panel_hits(panels: dict[str, list[str]], hits: dict[str, dict[str, float]], threshold: float, min_hmm_coverage: float = 0.0) -> list[dict[str, object]]:
    rows = []
    for panel, accessions in panels.items():
        detected = sum(
            1 for acc in accessions
            if acc in hits and hits[acc]["evalue"] <= threshold and hits[acc].get("hmm_coverage", 0.0) >= min_hmm_coverage
        )
        row: dict[str, object] = {
            "panel": panel,
            "tested": len(accessions),
            "detected": detected,
            "threshold": threshold,
            "min_hmm_coverage": min_hmm_coverage,
            "formal_denominator": panel in FORMAL_NEGATIVE or panel in POSITIVE,
        }
        if panel in POSITIVE:
            row.update({"TP": detected, "FN": len(accessions) - detected, "FP": 0, "TN": 0, "challenge_detected": 0})
        elif panel in FORMAL_NEGATIVE:
            row.update({"TP": 0, "FN": 0, "FP": detected, "TN": len(accessions) - detected, "challenge_detected": 0})
        else:
            row.update({"TP": 0, "FN": 0, "FP": 0, "TN": 0, "challenge_detected": detected})
        rows.append(row)
    return rows


def summarize_sensitivity_grid(panels: dict[str, list[str]], hits: dict[str, dict[str, float]], thresholds: list[float], coverages: list[float]) -> list[dict[str, object]]:
    rows = []
    for threshold in thresholds:
        for coverage in coverages:
            rows.extend(summarize_panel_hits(panels, hits, threshold, coverage))
    return rows


def calibrate(models: dict[str, Path], panel_dir: Path, outdir: Path, threshold: float = 1e-5, min_hmm_coverage: float = 0.0, cpu: int = 4, hmmsearch_bin: str | None = None, sensitivity_thresholds: list[float] | None = None, sensitivity_coverages: list[float] | None = None) -> dict[str, object]:
    panel_records: dict[str, dict[str, tuple[str, str]]] = {}
    accession_panel: dict[str, str] = {}
    for panel, filename in PANEL_FILES.items():
        path = panel_dir / filename
        records = read_fasta(path)
        panel_records[panel] = records
        for accession in records:
            if accession in accession_panel:
                raise ValueError(f"accession appears in multiple panels: {accession}")
            accession_panel[accession] = panel
    all_records = {acc: record for records in panel_records.values() for acc, record in records.items()}
    panels = {panel: sorted(records) for panel, records in panel_records.items()}
    hmmsearch_bin = hmmsearch_bin or shutil.which("hmmsearch")
    if not hmmsearch_bin:
        raise ValueError("hmmsearch executable not found")
    outdir.mkdir(parents=True, exist_ok=False)
    probe = outdir / "panel_probe.faa"
    with probe.open("w", encoding="utf-8", newline="\n") as handle:
        for accession in sorted(all_records):
            header, sequence = all_records[accession]
            handle.write(f">{header}\n{sequence}\n")
    summary_rows, sensitivity_rows, hit_rows = [], [], []
    sensitivity_thresholds = sensitivity_thresholds or [1e-5, 1e-10, 1e-20]
    sensitivity_coverages = sensitivity_coverages or [0.0, 0.4, 0.6, 0.8]
    model_records = {}
    with tempfile.TemporaryDirectory(prefix="panel_calibration_", dir=outdir) as work:
        for model, hmm in models.items():
            hmm = Path(hmm)
            if not hmm.is_file() or hmm.stat().st_size == 0:
                raise ValueError(f"missing HMM: {hmm}")
            tbl, dom = Path(work) / f"{model}.tbl", Path(work) / f"{model}.dom"
            subprocess.run([hmmsearch_bin, "--tblout", str(tbl), "--domtblout", str(dom), "-E", "1e-2", "--cpu", str(cpu), str(hmm), str(probe)], check=True, capture_output=True, text=True)
            hits = parse_tblout(tbl.read_text(encoding="utf-8").splitlines())
            lengths = {acc: len(seq) for acc, (_, seq) in all_records.items()}
            for acc, cov in parse_domtblout(dom.read_text(encoding="utf-8").splitlines(), lengths).items():
                if acc in hits:
                    hits[acc].update(cov)
            for acc, value in sorted(hits.items()):
                hit_rows.append({"model": model, "accession": acc, "panel": accession_panel.get(acc, "unknown"), **value})
            for row in summarize_panel_hits(panels, hits, threshold, min_hmm_coverage):
                summary_rows.append({"model": model, **row})
            for row in summarize_sensitivity_grid(panels, hits, sensitivity_thresholds, sensitivity_coverages):
                sensitivity_rows.append({"model": model, **row})
            model_records[model] = {"path": str(hmm.resolve()), "bytes": hmm.stat().st_size, "sha256": sha256_file(hmm)}
    probe.unlink()
    with (outdir / "panel_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["model", "panel", "tested", "detected", "threshold", "min_hmm_coverage", "formal_denominator", "TP", "FN", "FP", "TN", "challenge_detected"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n"); writer.writeheader(); writer.writerows(summary_rows)
    with (outdir / "panel_hits.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["model", "accession", "panel", "evalue", "bitscore", "hmm_coverage", "target_coverage"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore"); writer.writeheader(); writer.writerows(hit_rows)
    with (outdir / "panel_sensitivity.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["model", "panel", "tested", "detected", "threshold", "min_hmm_coverage", "formal_denominator", "TP", "FN", "FP", "TN", "challenge_detected"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n"); writer.writeheader(); writer.writerows(sensitivity_rows)
    metadata = {
        "schema_version": 1,
        "run_type": "small_scale_external_panel_calibration",
        "parameters": {"threshold": threshold, "min_hmm_coverage": min_hmm_coverage, "sensitivity_thresholds": sensitivity_thresholds, "sensitivity_coverages": sensitivity_coverages, "cpu": cpu, "hmmsearch_bin": hmmsearch_bin},
        "models": model_records,
        "panels": {panel: {"file": str((panel_dir / filename).resolve()), "bytes": (panel_dir / filename).stat().st_size, "sha256": sha256_file(panel_dir / filename), "tested": len(accessions)} for panel, filename in PANEL_FILES.items() for accessions in [panels[panel]]},
        "outputs": {name: {"bytes": (outdir / name).stat().st_size, "sha256": sha256_file(outdir / name)} for name in ("panel_summary.tsv", "panel_sensitivity.tsv", "panel_hits.tsv")},
    }
    (outdir / "panel_calibration_metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return metadata


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--model", action="append", required=True, help="NAME=HMM_PATH; repeat per model")
    parser.add_argument("--threshold", type=float, default=1e-5)
    parser.add_argument("--min-hmm-coverage", type=float, default=0.0)
    parser.add_argument("--cpu", type=int, default=4)
    parser.add_argument("--hmmsearch-bin", default=None)
    args = parser.parse_args(argv)
    models = {}
    for item in args.model:
        if "=" not in item:
            parser.error("--model must be NAME=HMM_PATH")
        name, path = item.split("=", 1)
        if not name or name in models:
            parser.error(f"invalid or duplicate model name: {name}")
        models[name] = Path(path)
    calibrate(models, args.panel_dir, args.outdir, args.threshold, args.min_hmm_coverage, args.cpu, args.hmmsearch_bin)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
