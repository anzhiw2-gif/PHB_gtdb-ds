#!/usr/bin/env python3
"""Prepare run-13 strict-core and broad-discovery inputs for tier processing."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


CORE_MODELS = {
    "ePhaZ_curated_core": "ePhaZ",
    "iPhaZ": "iPhaZ",
    "OH": "OH",
    "ArchPhaZ_hydrolase": "ArchPhaZ_hydrolase",
}
CORE_PRIORITY = {model: rank for rank, model in enumerate(CORE_MODELS)}
BROAD_MODEL = "ePhaZ_broad_discovery"


def _threshold(value: str) -> float:
    return float(value.replace("e-", "1e-"))


def _registry(path: Path) -> dict[str, tuple[float, float]]:
    with path.open(newline="") as handle:
        rows = {
            row["model"]: (_threshold(row["threshold"]), float(row["min_cov"]))
            for row in csv.DictReader(handle, delimiter="\t")
        }
    required = set(CORE_MODELS) | {BROAD_MODEL}
    missing = required - set(rows)
    if missing:
        raise ValueError(f"registry missing required models: {sorted(missing)}")
    return rows


def _accepted(row: dict[str, str], rules: dict[str, tuple[float, float]]) -> bool:
    model = row["family"]
    if model not in rules:
        return False
    threshold, min_cov = rules[model]
    return float(row["E-value"]) <= threshold and float(row["cov"]) >= min_cov


def _better(left: dict[str, str], right: dict[str, str]) -> bool:
    """Return whether left wins the core-family conflict for one protein."""
    left_rank = CORE_PRIORITY[left["family"]]
    right_rank = CORE_PRIORITY[right["family"]]
    return (left_rank, float(left["E-value"])) < (right_rank, float(right["E-value"]))


def prepare(hits: Path, registry: Path, outdir: Path) -> None:
    rules = _registry(registry)
    outdir.mkdir(parents=True, exist_ok=True)
    core: dict[str, dict[str, str]] = {}
    broad: dict[str, dict[str, str]] = {}
    with hits.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if not _accepted(row, rules):
                continue
            model = row["family"]
            protein = row["protein"]
            if model in CORE_MODELS:
                previous = core.get(protein)
                if previous is None or _better(row, previous):
                    core[protein] = row
            elif model == BROAD_MODEL:
                previous = broad.get(protein)
                if previous is None or float(row["E-value"]) < float(previous["E-value"]):
                    broad[protein] = row

    columns = ["family", "source_model", "genome", "protein", "E-value", "score", "cov"]
    with (outdir / "hits_filtered.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for protein, row in sorted(core.items()):
            writer.writerow({**row, "family": CORE_MODELS[row["family"]], "source_model": row["family"],
                             "genome": protein.split("|", 1)[0]})

    with (outdir / "unique_proteins.txt").open("w") as handle:
        for protein in sorted(core):
            handle.write(protein + "\n")

    with (outdir / "broad_discovery.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for protein, row in sorted(broad.items()):
            writer.writerow({**row, "source_model": BROAD_MODEL, "genome": protein.split("|", 1)[0]})

    (outdir / "summary.txt").write_text(
        f"strict_core_proteins={len(core)}\nbroad_discovery_proteins={len(broad)}\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hits", required=True, type=Path)
    parser.add_argument("--registry", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    prepare(args.hits, args.registry, args.outdir)


if __name__ == "__main__":
    main()
