#!/usr/bin/env python3
"""Summarize candidate placement in bounded PhaDED competition trees."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from Bio import Phylo


PANELS = ("extracellular_explicit", "phaC_vs_phaZ", "structure_anomaly", "ephaz_competition")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_membership(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    index = {}
    for row in rows:
        index[row["member_id"]] = row
        if row.get("source_accession"):
            index[row["source_accession"]] = row
    return index


def summarize_panel(tree_path: Path, membership_path: Path) -> list[dict[str, str]]:
    tree = Phylo.read(str(tree_path), "newick")
    members = read_membership(membership_path)
    terminals = {terminal.name: terminal for terminal in tree.get_terminals() if terminal.name}
    def row_for_leaf(name: str) -> dict[str, str] | None:
        if name in members:
            return members[name]
        parts = name.split("|", 3)
        if len(parts) != 4:
            return None
        if parts[0] == "candidate":
            return members.get(parts[2] + "|" + parts[3])
        if parts[0] == "reference_target":
            return members.get(parts[2])
        return None

    leaf_rows = {name: row_for_leaf(name) for name in terminals}
    refs = [name for name, row in leaf_rows.items() if row and row.get("role") == "reference_target"]
    candidates = [name for name, row in leaf_rows.items() if row and row.get("role") != "reference_target"]
    if not refs or not candidates:
        raise ValueError(f"panel lacks references or candidates: {tree_path}")
    output = []
    for candidate in candidates:
        nearest, distance = min(((ref, tree.distance(terminals[candidate], terminals[ref])) for ref in refs), key=lambda item: (item[1], item[0]))
        mrca = tree.common_ancestor(terminals[candidate], terminals[nearest])
        support = "" if mrca.confidence is None else f"{float(mrca.confidence):.1f}"
        row = dict(leaf_rows[candidate] or {})
        row.update({
            "tree_leaf_count": str(len(terminals)),
            "nearest_reference_member": nearest,
            "nearest_reference_accession": (leaf_rows[nearest] or {}).get("source_accession", ""),
            "nearest_reference_family": (leaf_rows[nearest] or {}).get("family", ""),
            "nearest_reference_distance": f"{distance:.6f}",
            "candidate_reference_mrca_support": support,
            "placement_interpretation": "candidate_only_phylogenetic_proximity",
        })
        output.append(row)
    return sorted(output, key=lambda row: (row.get("bucket", row.get("panel", "")), row["member_id"]))


def analyze(run_dir: Path) -> dict[str, object]:
    result_dir = run_dir / "results" / "phylogeny"
    rows: list[dict[str, str]] = []
    panel_stats = {}
    for panel in PANELS:
        tree = result_dir / "trees" / f"{panel}.treefile"
        membership = run_dir / "results" / "panels" / f"panel_{panel}.members.tsv"
        panel_rows = summarize_panel(tree, membership)
        rows.extend(panel_rows)
        panel_stats[panel] = {
            "candidate_count": len(panel_rows),
            "tree_leaf_count": int(panel_rows[0]["tree_leaf_count"]),
            "tree_sha256": sha256(tree),
            "membership_sha256": sha256(membership),
        }
    summary = result_dir / "phylogeny_review_summary.tsv"
    fields = sorted({key for row in rows for key in row})
    with summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    manifest = {
        "status": "completed_candidate_only",
        "panels": panel_stats,
        "summary": {"path": str(summary.resolve()), "size": summary.stat().st_size, "sha256": sha256(summary)},
        "phenotype_boundary": "Phylogenetic proximity and branch support are candidate homology evidence, not PHB/PHA degradation validation.",
        "leaf_count_note": "Leaf counts were recomputed from Newick terminal nodes; the wrapper's grep-based count overcounted internal support labels.",
    }
    manifest_path = result_dir / "phylogeny_review_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.run_dir)
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
