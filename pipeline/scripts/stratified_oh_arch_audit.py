#!/usr/bin/env python3
"""Deterministic cross-phylum audit for OH and ArchPhaZ hydrolase candidates."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

FAMILIES = {"OH", "ArchPhaZ_hydrolase"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _read_fasta(path: Path) -> dict[str, str]:
    records, header, seq = {}, None, []
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(seq)
                header, seq = line[1:].split()[0], []
            elif header is not None:
                seq.append(line.strip())
        if header is not None:
            records[header] = "".join(seq)
    return records


def _coverage_tier(family: str, cov: float) -> str:
    if family == "OH":
        return "oh_high" if cov >= 0.8 else "oh_review"
    if cov >= 0.8:
        return "arch_high"
    if cov >= 0.6:
        return "arch_review"
    return "arch_exploratory"


def audit(tier_table: Path, hits: Path, fasta: Path | list[Path], outdir: Path, *, per_stratum: int = 5, seed: int = 20260903) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    tier_rows: dict[tuple[str, str], dict[str, str]] = {}
    with tier_table.open(newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("family") in FAMILIES:
                tier_rows[(row["genome"], row["family"])] = row

    hit_by_key: defaultdict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    with hits.open(newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("family") not in FAMILIES:
                continue
            genome = row.get("protein", "").split("|", 1)[0]
            if (genome, row["family"]) in tier_rows:
                hit_by_key[(genome, row["family"])].append(row)

    fasta_paths = fasta if isinstance(fasta, list) else [fasta]
    records = {}
    for path in fasta_paths:
        records.update(_read_fasta(path))

    candidates = []
    for key, row in tier_rows.items():
        hits_for_key = hit_by_key.get(key, [])
        if not hits_for_key:
            candidates.append({**row, "coverage": "", "E-value": "", "score": "", "protein": "", "coverage_tier": "missing_raw_hit"})
            continue
        raw_best = max(hits_for_key, key=lambda h: (float(h.get("cov", 0)), -float(h.get("E-value", "inf"))))
        available = [h for h in hits_for_key if h.get("protein", "") in records]
        best = max(available or hits_for_key, key=lambda h: (float(h.get("cov", 0)), -float(h.get("E-value", "inf"))))
        cov = float(best.get("cov", 0))
        candidates.append({**row, "coverage": f"{cov:.6g}", "E-value": best.get("E-value", ""), "score": best.get("score", ""),
                           "protein": best.get("protein", ""), "raw_best_protein": raw_best.get("protein", ""),
                           "raw_best_coverage": raw_best.get("cov", ""), "coverage_tier": _coverage_tier(row["family"], cov)})

    with (outdir / "phylum_counts.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["family", "coverage_tier", "domain", "phylum", "genomes"])
        counts = Counter()
        for row in candidates:
            taxonomy = row.get("gtdb_taxonomy", "")
            domain = taxonomy.split(";", 1)[0].removeprefix("d__") or "unknown"
            phylum = row.get("phylum") or (taxonomy.split(";")[1].removeprefix("p__") if ";" in taxonomy else "unknown")
            counts[(row["family"], row["coverage_tier"], domain, phylum)] += 1
        for key, value in sorted(counts.items()):
            writer.writerow([*key, value])

    strata: defaultdict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in candidates:
        domain = row.get("gtdb_taxonomy", "").split(";", 1)[0].removeprefix("d__") or "unknown"
        phylum = row.get("phylum") or "unknown"
        strata[(row["family"], row["coverage_tier"], phylum)].append(row)
    rng = random.Random(seed)
    selected = []
    for key, rows in sorted(strata.items()):
        rows = sorted(rows, key=lambda r: (r["genome"], r.get("protein", "")))
        if len(rows) > per_stratum:
            rows = rng.sample(rows, per_stratum)
        selected.extend(rows)

    sample_fields = ["genome", "family", "coverage_tier", "coverage", "E-value", "score", "protein", "gtdb_taxonomy", "phylum", "class", "sequence_status", "sequence_length"]
    missing = []
    with (outdir / "sample.tsv").open("w", newline="", encoding="utf-8") as fh, (outdir / "sample.faa").open("w", encoding="utf-8") as fa:
        writer = csv.DictWriter(fh, fieldnames=sample_fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in sorted(selected, key=lambda r: (r["family"], r["coverage_tier"], r.get("phylum", ""), r["genome"])):
            protein = row.get("protein", "")
            sequence = records.get(protein, "")
            row["sequence_status"] = "verified" if sequence else "missing"
            row["sequence_length"] = len(sequence) if sequence else ""
            if sequence:
                fa.write(f">{protein}\n{sequence}\n")
            else:
                missing.append(row)
            writer.writerow(row)
    with (outdir / "missing_sequences.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=sample_fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader(); writer.writerows(missing)

    summary = {
        "status": "completed", "scope": "read-only cross-phylum sample; candidate-only",
        "seed": seed, "per_stratum": per_stratum, "families": sorted(FAMILIES),
        "candidate_genomes": len(candidates), "sample_records": len(selected), "missing_sequences": len(missing),
        "strata": {"|".join(key): len(rows) for key, rows in sorted(strata.items())},
        "inputs": {
            "tier_table": {"path": str(tier_table), "size": tier_table.stat().st_size if tier_table.is_file() else None,
                           "sha256": sha256(tier_table) if tier_table.is_file() else None, "status": "verified" if tier_table.is_file() else "pending"},
            "raw_hits": {"path": str(hits), "size": hits.stat().st_size if hits.is_file() else None,
                         "sha256": sha256(hits) if hits.is_file() else None, "status": "verified" if hits.is_file() else "pending"},
            "validated_fasta": [{"path": str(path), "size": path.stat().st_size if path.is_file() else None,
                                 "sha256": sha256(path) if path.is_file() else None, "status": "verified" if path.is_file() else "pending"}
                                for path in fasta_paths],
        },
        "outputs": {name: {"path": str(outdir / name), "sha256": sha256(outdir / name)}
                    for name in ["phylum_counts.tsv", "sample.tsv", "sample.faa", "missing_sequences.tsv"]},
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    report = ["# OH / ArchPhaZ cross-phylum audit", "", "Status: `completed`", "", "Read-only, deterministic stratified sample; evidence remains `candidate-only`.", "",
              f"- Candidate genomes represented in tier table: `{len(candidates)}`", f"- Sample records: `{len(selected)}`", f"- Missing FASTA records: `{len(missing)}`",
              f"- Sampling seed: `{seed}`; maximum per family × coverage tier × phylum stratum: `{per_stratum}`", "",
              "The full phylum counts are in `phylum_counts.tsv`; sampled evidence is in `sample.tsv` and `sample.faa`. Domain composition must be read from the stratified table, not inferred from the sample alone.", "",
              "HMM/domain/sequence evidence is not experimental PHB degradation confirmation."]
    (outdir / "AUDIT_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tier-table", type=Path, required=True)
    parser.add_argument("--hits", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True, nargs="+")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--per-stratum", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260903)
    args = parser.parse_args()
    audit(args.tier_table, args.hits, args.fasta, args.outdir, per_stratum=args.per_stratum, seed=args.seed)


if __name__ == "__main__":
    main()
