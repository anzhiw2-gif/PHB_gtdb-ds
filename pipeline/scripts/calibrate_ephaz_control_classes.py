#!/usr/bin/env python3
"""Calibrate ePhaZ/iPhaZ models with evidence-class-separated controls."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path


SUPPORTED_CLASSES = {"ePhaZ_curated_core", "ePhaZ_architecture_remote", "iPhaZ_like_challenge", "negative"}


def validate_control_classes(classes: dict[str, str]) -> None:
    unknown = {value for value in classes.values() if value not in SUPPORTED_CLASSES}
    if unknown:
        raise ValueError(f"unsupported control class(es): {sorted(unknown)}")


def parse_tblout(lines) -> dict[str, dict[str, float]]:
    hits: dict[str, dict[str, float]] = {}
    for raw in lines:
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        fields = raw.split()
        if len(fields) < 6:
            continue
        try:
            accession = fields[0].split("|", 1)[0]
            evalue = float(fields[4])
            bitscore = float(fields[5])
        except (IndexError, ValueError):
            continue
        current = hits.get(accession)
        if current is None or bitscore > current["bitscore"]:
            hits[accession] = {"evalue": evalue, "bitscore": bitscore, "coverage": 0.0}
    return hits


def parse_domtblout(lines, lengths: dict[str, int]) -> dict[str, float]:
    coverage: dict[str, float] = {}
    for raw in lines:
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        fields = raw.split()
        if len(fields) < 21:
            continue
        try:
            accession = fields[0].split("|", 1)[0]
            start, end = int(fields[19]), int(fields[20])
            length = lengths[accession]
        except (KeyError, ValueError, IndexError):
            continue
        coverage[accession] = max(coverage.get(accession, 0.0), min(1.0, (end - start + 1) / max(length, 1)))
    return coverage


def summarize_by_class(classes: dict[str, str], hits: dict[str, dict[str, float]], threshold: float, min_cov: float = 0.0) -> list[dict[str, object]]:
    rows = []
    for control_class in sorted(SUPPORTED_CLASSES | set(classes.values())):
        accessions = [acc for acc, value in classes.items() if value == control_class]
        if not accessions:
            continue
        detected = sum(1 for acc in accessions if acc in hits and hits[acc]["evalue"] <= threshold and hits[acc].get("coverage", 0.0) >= min_cov)
        row = {"control_class": control_class, "tested": len(accessions), "detected": detected, "threshold": threshold, "min_cov": min_cov}
        if control_class == "negative":
            row.update({"TP": 0, "FN": 0, "FP": detected, "TN": len(accessions) - detected})
        elif control_class == "iPhaZ_like_challenge":
            row.update({"TP": 0, "FN": 0, "FP": 0, "TN": 0, "challenge_detected": detected})
        else:
            row.update({"TP": detected, "FN": len(accessions) - detected, "FP": 0, "TN": 0})
        rows.append(row)
    return rows


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _file_record(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": _sha256(path)}


def _read_fasta(path: Path) -> tuple[dict[str, tuple[str, str]], dict[str, int]]:
    records = {}
    header = None
    sequence = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith(">"):
            if header is not None:
                accession = header.split("|", 1)[0]
                if accession in records:
                    raise ValueError(f"duplicate FASTA accession: {accession}")
                records[accession] = (header, "".join(sequence))
            header, sequence = line[1:], []
        elif line:
            sequence.append(line)
    if header is not None:
        accession = header.split("|", 1)[0]
        if accession in records:
            raise ValueError(f"duplicate FASTA accession: {accession}")
        records[accession] = (header, "".join(sequence))
    if not records:
        raise ValueError(f"invalid FASTA: {path}")
    return records, {acc: len(seq) for acc, (_, seq) in records.items()}


def _run_model(hmm: Path, probe: Path, workdir: Path, hmmsearch_bin: str, cpu: int, lengths: dict[str, int]) -> dict[str, dict[str, float]]:
    workdir.mkdir(parents=True, exist_ok=True)
    tbl, dom = workdir / "hits.tbl", workdir / "hits.dom"
    subprocess.run([hmmsearch_bin, "--tblout", str(tbl), "--domtblout", str(dom), "-E", "1e-2", "--cpu", str(cpu), str(hmm), str(probe)], check=True, capture_output=True, text=True)
    hits = parse_tblout(tbl.read_text(encoding="utf-8").splitlines())
    coverage = parse_domtblout(dom.read_text(encoding="utf-8").splitlines(), lengths)
    for accession, value in hits.items():
        value["coverage"] = coverage.get(accession, 0.0)
    return hits


def calibrate_control_classes(ephaz_hmm: Path | str, iphaz_hmm: Path | str, positive_faa: Path | str, negative_faa: Path | str, class_manifest: Path | str, outdir: Path | str, threshold: float = 1e-5, min_cov: float = 0.0, hmmsearch_bin: str | None = None, cpu: int = 4, challenge_faa: Path | str | None = None) -> dict[str, object]:
    ephaz_hmm, iphaz_hmm, positive_faa, negative_faa, class_manifest, outdir = map(Path, (ephaz_hmm, iphaz_hmm, positive_faa, negative_faa, class_manifest, outdir))
    rows = []
    with class_manifest.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append(row)
    classes = {row["accession"]: row["control_class"] for row in rows}
    validate_control_classes(classes)
    pos, pos_lengths = _read_fasta(positive_faa)
    challenge, challenge_lengths = ({}, {}) if challenge_faa is None else _read_fasta(Path(challenge_faa))
    neg, neg_lengths = _read_fasta(negative_faa)
    if set(pos) & set(challenge) or set(pos) & set(neg) or set(challenge) & set(neg):
        raise ValueError("positive, challenge, and negative FASTA accession sets must be disjoint")
    classified_records = set(pos) | set(challenge)
    if set(classes) != classified_records:
        raise ValueError("class manifest must exactly cover positive and challenge FASTA accession sets")
    all_records = {**pos, **challenge, **neg}
    lengths = {**pos_lengths, **challenge_lengths, **neg_lengths}
    all_classes = {**classes, **{acc: "negative" for acc in neg}}
    probe = outdir / "control_probe.faa"
    outdir.mkdir(parents=True, exist_ok=False)
    with probe.open("w", encoding="utf-8") as handle:
        for accession in sorted(all_records):
            header, sequence = all_records[accession]
            handle.write(f">{header}\n{sequence}\n")
    hmmsearch_bin = hmmsearch_bin or shutil.which("hmmsearch")
    if not hmmsearch_bin:
        raise ValueError("hmmsearch executable not found")
    hit_rows = []
    summary_rows = []
    with tempfile.TemporaryDirectory(prefix="control_calibration_", dir=outdir) as work:
        for model, hmm in (("ePhaZ", ephaz_hmm), ("iPhaZ", iphaz_hmm)):
            hits = _run_model(hmm, probe, Path(work) / model, hmmsearch_bin, cpu, lengths)
            for accession, value in sorted(hits.items()):
                hit_rows.append({"model": model, "accession": accession, "control_class": all_classes[accession], "evalue": value["evalue"], "bitscore": value["bitscore"], "coverage": value["coverage"]})
            for summary in summarize_by_class(all_classes, hits, threshold, min_cov):
                summary_rows.append({"model": model, **summary})
    with (outdir / "calibration_by_class.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["model", "control_class", "tested", "detected", "threshold", "min_cov", "TP", "FN", "FP", "TN", "challenge_detected"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader(); writer.writerows(summary_rows)
    with (outdir / "calibration_hits.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["model", "accession", "control_class", "evalue", "bitscore", "coverage"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(hit_rows)
    metadata = {"schema_version": 1, "parameters": {"threshold": threshold, "min_cov": min_cov, "hmmsearch_bin": hmmsearch_bin, "cpu": cpu}, "inputs": {name: _file_record(path) for name, path in (("ephaz_hmm", ephaz_hmm), ("iphaz_hmm", iphaz_hmm), ("positive_faa", positive_faa), ("challenge_faa", Path(challenge_faa) if challenge_faa is not None else None), ("negative_faa", negative_faa), ("class_manifest", class_manifest)) if path is not None}, "outputs": {name: _file_record(outdir / name) for name in ("calibration_by_class.tsv", "calibration_hits.tsv")}}
    (outdir / "calibration_metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    probe.unlink()
    return metadata


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ephaz-hmm", required=True); parser.add_argument("--iphaz-hmm", required=True)
    parser.add_argument("--positive-faa", required=True); parser.add_argument("--challenge-faa", required=True); parser.add_argument("--negative-faa", required=True)
    parser.add_argument("--class-manifest", required=True); parser.add_argument("--outdir", required=True)
    parser.add_argument("--threshold", type=float, default=1e-5); parser.add_argument("--min-cov", type=float, default=0.0)
    parser.add_argument("--hmmsearch-bin", default=None); parser.add_argument("--cpu", type=int, default=4)
    args = parser.parse_args(argv)
    try:
        calibrate_control_classes(**vars(args))
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
