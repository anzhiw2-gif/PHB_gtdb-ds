#!/usr/bin/env python3
"""Validate Scheme A figure data, editable exports, and essential values."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.image as mpimg
import pandas as pd


STEMS = (
    "figure_0_seed_library",
    "figure_1_workflow_funnel",
    "figure_2_core_scale",
    "figure_3_phylum_distribution",
    "figure_4_ephaz_signal",
    "figure_5_neighborhood_context",
)
SUFFIXES = (".svg", ".pdf", ".tiff", ".png")
CORE_FAMILIES = {"ePhaZ", "iPhaZ", "OH", "ArchPhaZ_hydrolase"}
MARKERS = {"PhaC", "PhaE", "PhaJ", "BdhA", "phasin", "PHA_gran_rgn"}


class FigureQAError(ValueError):
    """Raised when a figure source or export does not meet its contract."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def check_svg_text(path: Path) -> None:
    if "<text" not in path.read_text(encoding="utf-8"):
        raise FigureQAError(f"SVG has no editable text nodes: {path}")


def check_required_exports(output_dir: Path, stem: str) -> list[Path]:
    files = []
    for suffix in SUFFIXES:
        path = output_dir / f"{stem}{suffix}"
        if not path.is_file() or path.stat().st_size == 0:
            raise FigureQAError(f"missing or empty export: {path}")
        files.append(path)
    return files


def read_tsv(data_dir: Path, name: str) -> pd.DataFrame:
    path = data_dir / name
    if not path.is_file() or path.stat().st_size == 0:
        raise FigureQAError(f"missing or empty source data: {path}")
    return pd.read_csv(path, sep="\t")


def check_source_invariants(data_dir: Path) -> dict[str, object]:
    union = read_tsv(data_dir, "figure2_union.tsv").iloc[0]
    if int(union["core_family_union_genomes"]) != 44_814:
        raise FigureQAError("Figure 2 union must equal 44,814 genomes")
    if int(union["gtdb_representatives"]) != 199_923:
        raise FigureQAError("Figure 2 denominator must equal 199,923 representatives")
    matrix = read_tsv(data_dir, "figure3_phylum_matrix.tsv")
    families = set(matrix["family"])
    if not families or not families.issubset(CORE_FAMILIES):
        raise FigureQAError("Figure 3 must contain only core candidate families")
    rates = read_tsv(data_dir, "figure5_neighborhood_rate.tsv")
    markers = set(rates["marker_family"])
    if markers != MARKERS:
        raise FigureQAError("Figure 5 marker set differs from the declared available marker set")
    return {
        "figure2_union": int(union["core_family_union_genomes"]),
        "figure2_denominator": int(union["gtdb_representatives"]),
        "figure3_families": sorted(families),
        "figure5_markers": sorted(markers),
    }


def check_provenance(data_dir: Path) -> None:
    path = data_dir / "provenance.json"
    if not path.is_file():
        raise FigureQAError("missing provenance.json")
    provenance = json.loads(path.read_text(encoding="utf-8"))
    for entry in provenance.get("outputs", {}).values():
        source_path = data_dir / entry["file"]
        if not source_path.is_file() or sha256(source_path) != entry["sha256"]:
            raise FigureQAError(f"source-data hash mismatch: {source_path}")


def run_checks(data_dir: Path, output_dir: Path) -> dict[str, object]:
    checks = check_source_invariants(data_dir)
    check_provenance(data_dir)
    exports: dict[str, dict[str, object]] = {}
    for stem in STEMS:
        files = check_required_exports(output_dir, stem)
        svg = next(path for path in files if path.suffix == ".svg")
        png = next(path for path in files if path.suffix == ".png")
        check_svg_text(svg)
        pixels = mpimg.imread(png)
        if pixels.shape[0] == 0 or pixels.shape[1] == 0:
            raise FigureQAError(f"PNG has zero dimensions: {png}")
        exports[stem] = {
            "files": [path.name for path in files],
            "png_pixels": [int(pixels.shape[1]), int(pixels.shape[0])],
            "editable_svg_text": True,
        }
    return {"source_invariants": checks, "exports": exports}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = {
        "status": "pass",
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        **run_checks(args.data_dir, args.output_dir),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"PASS: {len(STEMS)} figures; report: {args.report}")


if __name__ == "__main__":
    main()
