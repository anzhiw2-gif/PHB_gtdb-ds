#!/usr/bin/env python3
"""Compare baseline and bridge-augmented ePhaZ HMMs with leave-one-out tests."""

from __future__ import annotations

import csv
import argparse
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Iterable


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def describe_file(path: Path) -> dict[str, object]:
    """Describe a non-empty preserved calibration artifact."""
    if not path.is_file() or not path.stat().st_size:
        raise RuntimeError(f"missing or empty calibration artifact: {path}")
    return {"path": str(path.resolve()), "size": path.stat().st_size, "sha256": sha256(path)}


def parse_fasta_accession(header: str) -> str:
    """Normalize UniProt/refseq-style headers to the accession used in ledgers."""
    token = header.lstrip("> ").split()[0]
    fields = token.split("|")
    if len(fields) >= 2 and fields[0].lower() in {"sp", "tr", "ref", "gi"}:
        return fields[1]
    return fields[0]


def read_fasta(path: str | Path) -> dict[str, str]:
    records: dict[str, str] = {}
    accession = None
    chunks: list[str] = []
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        if raw.startswith(">"):
            if accession is not None:
                if accession in records or not chunks:
                    raise ValueError(f"invalid FASTA record: {accession}")
                records[accession] = "".join(chunks)
            accession, chunks = parse_fasta_accession(raw[1:]), []
        elif raw.strip():
            chunks.append(raw.strip())
    if accession is not None:
        if accession in records or not chunks:
            raise ValueError(f"invalid FASTA record: {accession}")
        records[accession] = "".join(chunks)
    if not records:
        raise ValueError(f"empty FASTA: {path}")
    return records


def write_fasta(path: Path, records: dict[str, str], accessions: Iterable[str]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for accession in accessions:
            handle.write(f">{accession}\n{records[accession]}\n")


def parse_hits(tbl: Path) -> dict[str, float]:
    hits: dict[str, float] = {}
    for line in tbl.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split()
        if len(fields) >= 5:
            try:
                hits[parse_fasta_accession(fields[0])] = float(fields[4])
            except ValueError:
                pass
    return hits


def _run(command: list[str], **kwargs):
    return subprocess.run(command, check=True, capture_output=True, text=True, **kwargs)


def tool_version(executable: str) -> str:
    """Return the first version/help line without masking a failed executable."""
    flag = "--version" if Path(executable).name.lower().startswith("mafft") else "-h"
    result = subprocess.run([executable, flag], check=False, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"cannot determine tool version for {executable}: exit {result.returncode}")
    lines = (result.stdout + "\n" + result.stderr).splitlines()
    return next((line.strip() for line in lines if line.strip()), "pending")


def tool_record(executable: str) -> dict[str, object]:
    """Bind an executable file to the version used by this calibration."""
    path = Path(executable).resolve()
    if not path.is_file() or not path.stat().st_size:
        raise FileNotFoundError(f"tool is missing or empty: {path}")
    return {
        "path": str(path),
        "version": tool_version(str(path)),
        "size": path.stat().st_size,
        "sha256": sha256(path),
    }


def build_hmm(training_fasta: Path, output_hmm: Path, aligned: Path, mafft_bin: str, hmmbuild_bin: str) -> None:
    result = _run([mafft_bin, "--auto", str(training_fasta)])
    aligned.write_text(result.stdout, encoding="utf-8", newline="\n")
    _run([hmmbuild_bin, str(output_hmm), str(aligned)])
    if not output_hmm.is_file() or not output_hmm.stat().st_size:
        raise RuntimeError(f"hmmbuild did not create {output_hmm}")


def run_loo_calibration(
    core_fasta: str | Path, bridge_fasta: str | Path, challenge_fasta: str | Path,
    negative_fasta: str | Path, outdir: str | Path, *, thresholds: tuple[float, ...] = (1e-2, 1e-5, 1e-10, 1e-20),
    mafft_bin: str = "mafft", hmmbuild_bin: str = "hmmbuild", hmmsearch_bin: str = "hmmsearch",
    software_versions: dict[str, str] | None = None,
    software_records: dict[str, dict[str, object]] | None = None,
) -> dict[str, object]:
    core, bridge = read_fasta(core_fasta), read_fasta(bridge_fasta)
    challenge, negative = read_fasta(challenge_fasta), read_fasta(negative_fasta)
    if len(core) < 4 or not bridge:
        raise ValueError("need at least four core and one bridge record")
    outdir = Path(outdir)
    if outdir.exists():
        raise FileExistsError(outdir)
    outdir.mkdir(parents=True)
    folds_dir = outdir / "folds"
    folds_dir.mkdir()
    summary: list[dict[str, object]] = []
    model_artifacts: list[dict[str, object]] = []
    hit_rows: list[dict[str, object]] = []
    for holdout in core:
        train = [item for item in core if item != holdout]
        for name, training_ids in ((f"baseline_without_{holdout}", train), (f"bridge_augmented_without_{holdout}", train + list(bridge))):
                fold = folds_dir / name
                fold.mkdir()
                training = fold / "training.faa"
                combined = dict(core)
                combined.update(bridge)
                write_fasta(training, combined, training_ids)
                alignment = fold / "alignment.faa"
                hmm = fold / "model.hmm"
                build_hmm(training, hmm, alignment, mafft_bin, hmmbuild_bin)
                probe_records = dict(core)
                probe_records.update(bridge)
                probe_records.update(challenge)
                probe_records.update(negative)
                probe = fold / "probe.faa"
                write_fasta(probe, probe_records, [holdout, *challenge, *negative])
                tbl, dom = fold / "hits.tblout", fold / "hits.domtblout"
                search = _run([hmmsearch_bin, "--tblout", str(tbl), "--domtblout", str(dom), "-E", "1e-2", str(hmm), str(probe)])
                (fold / "hmmsearch.stdout.log").write_text(search.stdout or "", encoding="utf-8", newline="\n")
                (fold / "hmmsearch.stderr.log").write_text(search.stderr or "", encoding="utf-8", newline="\n")
                command = {"mafft": [mafft_bin, "--auto", str(training)], "hmmbuild": [hmmbuild_bin, str(hmm), str(alignment)], "hmmsearch": [hmmsearch_bin, "--tblout", str(tbl), "--domtblout", str(dom), "-E", "1e-2", str(hmm), str(probe)]}
                (fold / "commands.json").write_text(json.dumps(command, indent=2) + "\n", encoding="utf-8")
                hits = parse_hits(tbl)
                for accession in [holdout, *challenge, *negative]:
                    hit_rows.append({"model": name, "holdout": holdout, "accession": accession, "control_class": "held_out_core" if accession == holdout else ("iPhaZ_like_challenge" if accession in challenge else "negative"), "full_sequence_evalue": hits.get(accession, "no_hit")})
                for threshold in thresholds:
                    summary.append({
                        "model": name, "holdout": holdout, "threshold": f"{threshold:.0e}",
                        "holdout_hit": str(hits.get(holdout, float("inf")) <= threshold).lower(),
                        "challenge_hits": sum(value <= threshold for accession, value in hits.items() if accession in challenge),
                        "negative_hits": sum(value <= threshold for accession, value in hits.items() if accession in negative),
                    })
                model_artifacts.append({"model": name, "holdout": holdout, "training_accessions": training_ids, "artifacts": {"training_fasta": describe_file(training), "alignment_fasta": describe_file(alignment), "model_hmm": describe_file(hmm), "probe_fasta": describe_file(probe), "tblout": describe_file(tbl), "domtblout": describe_file(dom), "commands": describe_file(fold / "commands.json")}})
    fields = ["model", "holdout", "threshold", "holdout_hit", "challenge_hits", "negative_hits"]
    with (outdir / "loo_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary)
    with (outdir / "loo_hit_ledger.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["model", "holdout", "accession", "control_class", "full_sequence_evalue"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(hit_rows)
    software_versions = software_versions or {}
    software_records = software_records or {}
    metadata = {
        "schema_version": 1, "method": "leave_one_core_out", "fold_count": len(core),
        "thresholds": list(thresholds), "role_constraints": {"iPhaZ_like_challenge": "must_remain_zero_hits", "negative": "must_remain_zero_hits"},
        "inputs": {str(Path(path).name): describe_file(Path(path)) for path in (core_fasta, bridge_fasta, challenge_fasta, negative_fasta)},
        "models": model_artifacts,
        "software": {
            "mafft": software_records.get("mafft", {"path": mafft_bin, "version": software_versions.get("mafft", "pending"), "size": "pending", "sha256": "pending"}),
            "hmmbuild": software_records.get("hmmbuild", {"path": hmmbuild_bin, "version": software_versions.get("hmmbuild", "pending"), "size": "pending", "sha256": "pending"}),
            "hmmsearch": software_records.get("hmmsearch", {"path": hmmsearch_bin, "version": software_versions.get("hmmsearch", "pending"), "size": "pending", "sha256": "pending"}),
        },
    }
    (outdir / "loo_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return {"fold_count": len(core), "summary": str(outdir / "loo_summary.tsv")}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core", required=True)
    parser.add_argument("--bridge", required=True)
    parser.add_argument("--challenge", required=True)
    parser.add_argument("--negative", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--mafft", default="mafft")
    parser.add_argument("--hmmbuild", default="hmmbuild")
    parser.add_argument("--hmmsearch", default="hmmsearch")
    args = parser.parse_args(argv)
    result = run_loo_calibration(
        args.core, args.bridge, args.challenge, args.negative, args.outdir,
        mafft_bin=args.mafft, hmmbuild_bin=args.hmmbuild, hmmsearch_bin=args.hmmsearch,
        software_records={
            "mafft": tool_record(args.mafft),
            "hmmbuild": tool_record(args.hmmbuild),
            "hmmsearch": tool_record(args.hmmsearch),
        },
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
