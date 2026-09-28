#!/usr/bin/env python3
"""Read-only QA audit for completed formal scan 13 and tier processing."""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path


def _threshold(value: str) -> float:
    return float(value.replace("e-", "1e-"))


def _registry(path: Path) -> dict[str, tuple[float, float]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return {
            row["model"]: (_threshold(row["threshold"]), float(row["min_cov"]))
            for row in csv.DictReader(fh, delimiter="\t")
        }


def _sha256(path: Path) -> str:
    if not path.is_file():
        return "pending"
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _file_record(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {"path": str(path), "status": "pending", "size": None, "sha256": None}
    return {"path": str(path), "status": "verified", "size": path.stat().st_size, "sha256": _sha256(path)}


def _tier_stats(path: Path) -> tuple[dict[str, int], Counter[str], set[str]]:
    family_genomes: defaultdict[str, set[str]] = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            genome, family = row.get("genome", ""), row.get("family", "")
            if genome and family:
                family_genomes[family].add(genome)
    combos = Counter("+".join(sorted(fams)) for fams in (
        {family for family, genomes in family_genomes.items() if genome in genomes}
        for genome in set().union(*family_genomes.values())
    ))
    counts = {family: len(genomes) for family, genomes in family_genomes.items()}
    union = set().union(*family_genomes.values()) if family_genomes else set()
    return counts, combos, union


def _taxonomy_stats(path: Path | None, genomes: set[str]) -> dict[str, object]:
    if path is None or not path.is_file():
        return {"status": "pending", "genomes": len(genomes), "classified": None, "unclassified": None}
    found: set[str] = set()
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t", 1)
            if len(parts) < 2:
                continue
            accession = parts[0]
            if accession in genomes:
                found.add(accession)
            elif accession.startswith("RS_") and accession[3:] in genomes:
                found.add(accession[3:])
            elif accession.startswith("GB_") and accession[3:] in genomes:
                found.add(accession[3:])
    return {"status": "verified", "genomes": len(genomes), "classified": len(found), "unclassified": len(genomes - found)}


def _job_stats(parent: Path) -> dict[str, int]:
    joblog = parent / "logs" / "parallel.joblog"
    total = failed = 0
    if joblog.is_file():
        with joblog.open(encoding="utf-8", errors="replace") as fh:
            next(fh, None)
            for line in fh:
                fields = line.rstrip("\n").split(None, 8)
                if len(fields) >= 7:
                    total += 1
                    failed += int(fields[6]) != 0
    ok_files = len(list((parent / "results" / "task_status").glob("*.ok")))
    manifest = {}
    manifest_path = parent / "results" / "scan_manifest.json"
    if manifest_path.is_file():
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            manifest = {}
    return {
        "joblog_tasks": total,
        "joblog_failed": failed,
        "ok_markers": ok_files,
        "manifest_task_total": int(manifest.get("task_total", 0)),
        "manifest_task_completed": int(manifest.get("task_completed", 0)),
        "manifest_reused_tasks": int(manifest.get("reused_tasks", 0)),
        "manifest_status": manifest.get("status", "pending"),
    }


def audit(parent: Path, tier: Path, registry: Path, outdir: Path, taxonomy: Path | None = None,
          historical: Path | None = None) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    rules = _registry(registry)
    counts = Counter()
    accepted_genomes: defaultdict[str, set[str]] = defaultdict(set)
    accepted_rows = 0
    duplicate_protein_rows = duplicate_accession_rows = 0
    db = sqlite3.connect(":memory:")
    db.execute("PRAGMA synchronous=OFF")
    db.execute("CREATE TABLE proteins (id TEXT PRIMARY KEY)")
    db.execute("CREATE TABLE accessions (id TEXT PRIMARY KEY)")
    hits = parent / "results" / "hits_all.tsv"
    with hits.open(newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            model = row.get("family", "")
            if model not in rules:
                continue
            counts[f"{model}:raw_rows"] += 1
            try:
                ev, cov = float(row["E-value"]), float(row.get("cov", "0"))
            except (KeyError, ValueError):
                counts[f"{model}:malformed_rows"] += 1
                continue
            threshold, min_cov = rules[model]
            if cov < min_cov:
                counts[f"{model}:coverage_filtered_rows"] += 1
            if ev > threshold or cov < min_cov:
                continue
            accepted_rows += 1
            counts[f"{model}:accepted_rows"] += 1
            protein = row.get("protein", "")
            accession = protein.split("|", 1)[0]
            if protein and db.execute("INSERT OR IGNORE INTO proteins VALUES (?)", (protein,)).rowcount == 0:
                duplicate_protein_rows += 1
            if accession and db.execute("INSERT OR IGNORE INTO accessions VALUES (?)", (accession,)).rowcount == 0:
                duplicate_accession_rows += 1
            accepted_genomes[model].add(accession)
    model_names = list(rules)
    with (outdir / "model_counts.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["model", "raw_rows", "accepted_rows", "coverage_filtered_rows", "unique_genomes"])
        for model in model_names:
            writer.writerow([model, counts[f"{model}:raw_rows"], counts[f"{model}:accepted_rows"],
                             counts[f"{model}:coverage_filtered_rows"], len(accepted_genomes[model])])

    tier_table = tier / "results" / "tables" / "tier1_genome_family.tsv"
    tier_counts, combos, tier_union = _tier_stats(tier_table)
    with (outdir / "tier_summary.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["family", "tier1_genomes"])
        for family, number in sorted(tier_counts.items()):
            writer.writerow([family, number])
        writer.writerow(["strict_union", len(tier_union)])
    with (outdir / "tier_cooccurrence.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t"); writer.writerow(["family_set", "genomes"])
        writer.writerows(sorted(combos.items(), key=lambda item: (-item[1], item[0])))

    all_accepted_genomes = set().union(*accepted_genomes.values()) if accepted_genomes else set()
    overlap = {"status": "pending", "historical_genomes": None, "run13_genomes": len(all_accepted_genomes), "intersection": None}
    if historical and historical.is_file():
        old = set()
        with historical.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                value = line.rstrip("\n").split("\t", 1)[0]
                if value and value.lower() not in {"genome", "accession"}:
                    old.add(value)
        overlap = {"status": "verified", "historical_genomes": len(old), "run13_genomes": len(all_accepted_genomes),
                   "intersection": len(old & all_accepted_genomes)}

    input_paths = {
        "parent_hits": hits, "parent_scan_manifest": parent / "results" / "scan_manifest.json",
        "parent_input_contract": parent / "input_contract.json", "tier_manifest": tier / "results" / "tier_processing_manifest.json",
        "tier_input_contract": tier / "input_contract.json", "registry": registry,
        "tier_table": tier_table, "parallel_joblog": parent / "logs" / "parallel.joblog",
    }
    with (outdir / "input_hashes.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t"); writer.writerow(["name", "path", "status", "size", "sha256"])
        for name, path in input_paths.items():
            record = _file_record(path); writer.writerow([name, record["path"], record["status"], record["size"], record["sha256"]])

    job = _job_stats(parent)
    summary = {
        "status": "completed", "scope": "read-only QA; candidate-only interpretation",
        "parent_run": str(parent), "tier_run": str(tier), "accepted_rows": accepted_rows,
        "duplicate_protein_rows": duplicate_protein_rows, "duplicate_accession_rows": duplicate_accession_rows,
        "taxonomy": _taxonomy_stats(taxonomy, all_accepted_genomes), "jobs": job,
        "tier1": {"family_genomes": tier_counts, "strict_union": len(tier_union), "cooccurrence": dict(combos)},
        "historical_overlap": overlap,
        "input_hashes": {name: _file_record(path) for name, path in input_paths.items()},
    }
    (outdir / "qa_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    report = ["# Run 13 QA report", "", "Status: `completed`", "", "This is a read-only QA audit; all HMM/domain/tier evidence remains `candidate-only`.", "",
              f"- Accepted registry rows: `{accepted_rows}`", f"- Duplicate accepted protein rows: `{duplicate_protein_rows}`",
              f"- Duplicate accepted accession rows: `{duplicate_accession_rows}`",
              f"- Tasks completed/reused: `{job['manifest_task_completed']}/{job['manifest_task_total']}` / `{job['manifest_reused_tasks']}`",
              f"- Joblog executed/failed: `{job['joblog_tasks']}/{job['joblog_failed']}`",
              f"- Tier1 strict union: `{len(tier_union)}` genomes", "", "See `model_counts.tsv`, `tier_summary.tsv`, `tier_cooccurrence.tsv`, `input_hashes.tsv`, and `qa_summary.json` for machine-readable evidence."]
    (outdir / "QA_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--tier", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--taxonomy", type=Path)
    parser.add_argument("--historical", type=Path)
    args = parser.parse_args()
    audit(args.parent, args.tier, args.registry, args.outdir, args.taxonomy, args.historical)


if __name__ == "__main__":
    main()
