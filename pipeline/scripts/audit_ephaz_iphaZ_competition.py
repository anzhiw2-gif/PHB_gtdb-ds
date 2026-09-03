#!/usr/bin/env python3
"""Audit ePhaZ-vs-iPhaZ model competition among SignalP OTHER candidates."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path


def select_other_predictions(predictions: dict[str, str]) -> set[str]:
    return {accession for accession, prediction in predictions.items() if prediction.strip().upper() == "OTHER"}


def parse_tblout(lines) -> dict[str, dict[str, float]]:
    hits: dict[str, dict[str, float]] = {}
    for raw in lines:
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        fields = raw.split()
        if len(fields) < 6:
            continue
        try:
            accession = fields[0]
            evalue, bitscore = float(fields[4]), float(fields[5])
        except (ValueError, IndexError):
            continue
        if accession not in hits or bitscore > hits[accession]["bitscore"]:
            hits[accession] = {"evalue": evalue, "bitscore": bitscore}
    return hits


def parse_domtblout(lines, lengths: dict[str, int]) -> dict[str, dict[str, float]]:
    domains = {}
    for raw in lines:
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        fields = raw.split()
        if len(fields) < 21:
            continue
        try:
            accession = fields[0]
            hmm_length = float(fields[5])
            hmm_from, hmm_to = int(fields[15]), int(fields[16])
            target_from, target_to = int(fields[17]), int(fields[18])
            target_length = lengths[accession]
        except (KeyError, ValueError, IndexError, ZeroDivisionError):
            continue
        item = {
            "hmm_from": hmm_from,
            "hmm_to": hmm_to,
            "target_from": target_from,
            "target_to": target_to,
            "hmm_coverage": min(1.0, (hmm_to - hmm_from + 1) / max(hmm_length, 1.0)),
            "target_coverage": min(1.0, (target_to - target_from + 1) / max(target_length, 1)),
        }
        if accession not in domains or item["hmm_coverage"] > domains[accession]["hmm_coverage"]:
            domains[accession] = item
    return domains


def classify_competition(scores: dict[str, dict[str, dict[str, float]]], margin_bits: float = 20.0) -> dict[str, str]:
    result = {}
    for accession, values in scores.items():
        ephaz, iphaz = values.get("ephaz", {}), values.get("iphaz", {})
        if not ephaz and not iphaz:
            result[accession] = "no_reportable_hit"
        elif not iphaz:
            result[accession] = "ePhaZ_like"
        elif not ephaz:
            result[accession] = "iPhaZ_like"
        elif iphaz["bitscore"] - ephaz["bitscore"] >= margin_bits:
            result[accession] = "iPhaZ_like"
        elif ephaz["bitscore"] - iphaz["bitscore"] >= margin_bits:
            result[accession] = "ePhaZ_like"
        else:
            result[accession] = "ambiguous"
    return result


def _read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    records = {}
    header = None; sequence = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith(">"):
            if header is not None:
                records[header.split(None, 1)[0]] = (header, "".join(sequence))
            header, sequence = line[1:], []
        elif line:
            sequence.append(line)
    if header is not None:
        records[header.split(None, 1)[0]] = (header, "".join(sequence))
    if not records:
        raise ValueError(f"empty FASTA: {path}")
    return records


def _read_predictions(path: Path) -> dict[str, str]:
    predictions = {}
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) >= 2:
                predictions[fields[0]] = fields[1]
    return predictions


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _run(hmm: Path, probe: Path, workdir: Path, binary: str, cpu: int, lengths: dict[str, int]) -> dict[str, dict[str, float]]:
    workdir.mkdir(parents=True, exist_ok=True)
    tbl, dom = workdir / "hits.tbl", workdir / "hits.dom"
    subprocess.run([binary, "--tblout", str(tbl), "--domtblout", str(dom), "-E", "1e-2", "--cpu", str(cpu), str(hmm), str(probe)], check=True, capture_output=True, text=True)
    hits = parse_tblout(tbl.read_text(encoding="utf-8").splitlines())
    for accession, values in parse_domtblout(dom.read_text(encoding="utf-8").splitlines(), lengths).items():
        if accession in hits:
            hits[accession].update(values)
    return hits


def audit_competition(tier_faa: Path | str, signalp: Path | str, ephaz_hmm: Path | str, iphaz_hmm: Path | str, outdir: Path | str, hmmsearch_bin: str | None = None, cpu: int = 4, margin_bits: float = 20.0) -> dict[str, object]:
    tier_faa, signalp, ephaz_hmm, iphaz_hmm, outdir = map(Path, (tier_faa, signalp, ephaz_hmm, iphaz_hmm, outdir))
    records = _read_fasta(tier_faa)
    other = select_other_predictions(_read_predictions(signalp))
    selected = sorted(other & set(records))
    if not selected:
        raise ValueError("SignalP OTHER has no matching tier FASTA records")
    outdir.mkdir(parents=True, exist_ok=False)
    probe = outdir / "no_signal_peptide.faa"
    with probe.open("w", encoding="utf-8") as handle:
        for accession in selected:
            header, sequence = records[accession]
            handle.write(f">{header}\n{sequence}\n")
    hmmsearch_bin = hmmsearch_bin or shutil.which("hmmsearch")
    if not hmmsearch_bin:
        raise ValueError("hmmsearch executable not found")
    lengths = {accession: len(sequence) for accession, (_, sequence) in records.items()}
    with tempfile.TemporaryDirectory(prefix="competition_", dir=outdir) as work:
        ephaz = _run(ephaz_hmm, probe, Path(work) / "ephaz", hmmsearch_bin, cpu, lengths)
        iphaz = _run(iphaz_hmm, probe, Path(work) / "iphaz", hmmsearch_bin, cpu, lengths)
    scores = {acc: {"ephaz": ephaz.get(acc, {}), "iphaz": iphaz.get(acc, {})} for acc in selected}
    assignments = classify_competition(scores, margin_bits)
    with (outdir / "no_signal_competition.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["accession", "length", "ephaz_evalue", "ephaz_bitscore", "ephaz_hmm_from", "ephaz_hmm_to", "ephaz_target_from", "ephaz_target_to", "ephaz_hmm_coverage", "ephaz_target_coverage", "iphaz_evalue", "iphaz_bitscore", "iphaz_hmm_from", "iphaz_hmm_to", "iphaz_target_from", "iphaz_target_to", "iphaz_hmm_coverage", "iphaz_target_coverage", "assignment"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n"); writer.writeheader()
        for accession in selected:
            e, i = scores[accession]["ephaz"], scores[accession]["iphaz"]
            writer.writerow({"accession": accession, "length": len(records[accession][1]), "ephaz_evalue": e.get("evalue", ""), "ephaz_bitscore": e.get("bitscore", ""), "ephaz_hmm_from": e.get("hmm_from", ""), "ephaz_hmm_to": e.get("hmm_to", ""), "ephaz_target_from": e.get("target_from", ""), "ephaz_target_to": e.get("target_to", ""), "ephaz_hmm_coverage": e.get("hmm_coverage", ""), "ephaz_target_coverage": e.get("target_coverage", ""), "iphaz_evalue": i.get("evalue", ""), "iphaz_bitscore": i.get("bitscore", ""), "iphaz_hmm_from": i.get("hmm_from", ""), "iphaz_hmm_to": i.get("hmm_to", ""), "iphaz_target_from": i.get("target_from", ""), "iphaz_target_to": i.get("target_to", ""), "iphaz_hmm_coverage": i.get("hmm_coverage", ""), "iphaz_target_coverage": i.get("target_coverage", ""), "assignment": assignments[accession]})
    metadata = {"schema_version": 1, "parameters": {"hmmsearch_bin": hmmsearch_bin, "cpu": cpu, "margin_bits": margin_bits}, "counts": {label: sum(value == label for value in assignments.values()) for label in sorted(set(assignments.values()))}, "selected_other": len(selected), "inputs": {name: {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": _sha256(path)} for name, path in (("tier_faa", tier_faa), ("signalp", signalp), ("ephaz_hmm", ephaz_hmm), ("iphaz_hmm", iphaz_hmm))}, "outputs": {"no_signal_peptide.faa": {"path": str(probe.resolve()), "bytes": probe.stat().st_size, "sha256": _sha256(probe)}, "no_signal_competition.tsv": {"path": str((outdir / "no_signal_competition.tsv").resolve()), "bytes": (outdir / "no_signal_competition.tsv").stat().st_size, "sha256": _sha256(outdir / "no_signal_competition.tsv")}}}
    (outdir / "competition_metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return metadata


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tier-faa", required=True); parser.add_argument("--signalp", required=True)
    parser.add_argument("--ephaz-hmm", required=True); parser.add_argument("--iphaz-hmm", required=True)
    parser.add_argument("--outdir", required=True); parser.add_argument("--hmmsearch-bin", default=None)
    parser.add_argument("--cpu", type=int, default=4); parser.add_argument("--margin-bits", type=float, default=20.0)
    args = parser.parse_args(argv)
    try:
        audit_competition(**vars(args))
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
