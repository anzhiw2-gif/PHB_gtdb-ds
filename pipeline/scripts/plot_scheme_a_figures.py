#!/usr/bin/env python3
"""Render the six approved Scheme A scientific figures from compact TSVs."""
from __future__ import annotations

import argparse
from pathlib import Path
from textwrap import fill

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import pandas as pd


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 7
plt.rcParams["axes.linewidth"] = 0.7
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
plt.rcParams["legend.frameon"] = False


WIDTH_IN = 183 / 25.4
COLORS = {
    "ePhaZ": "#0F4D92",
    "iPhaZ": "#42949E",
    "OH": "#B64342",
    "ArchPhaZ_hydrolase": "#9A4D8E",
    "ArchPhaZ_patatin": "#767676",
    "neutral": "#767676",
    "light": "#D8D8D8",
    "dark": "#272727",
}
CORE_FAMILIES = ["ePhaZ", "iPhaZ", "OH", "ArchPhaZ_hydrolase"]
MARKERS = ["PhaC", "PhaE", "PhaJ", "BdhA", "phasin", "PHA_gran_rgn"]
HEATMAP = LinearSegmentedColormap.from_list("support", ["#F2F2F2", "#B9DAD7", "#42949E", "#0F4D92"])


def read_source(data_dir: Path, name: str) -> pd.DataFrame:
    path = data_dir / name
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing or empty source data: {path}")
    return pd.read_csv(path, sep="\t")


def add_panel_label(ax, label: str) -> None:
    ax.text(-0.12, 1.04, label, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom")


def style_axis(ax, grid: str | None = "x") -> None:
    if grid:
        ax.grid(axis=grid, color="#E5E5E5", linewidth=0.6, zorder=0)
    ax.tick_params(length=2.5, width=0.7, pad=2)


def family_label(family: str) -> str:
    return {"ArchPhaZ_hydrolase": "ArchPhaZ\nhydrolase", "ArchPhaZ_patatin": "ArchPhaZ\npatatin"}.get(family, family)


def save_exports(fig, output_dir: Path, stem: str) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    base = output_dir / stem
    paths = [base.with_suffix(suffix) for suffix in (".svg", ".pdf", ".tiff", ".png")]
    fig.savefig(paths[0], bbox_inches="tight", facecolor="white")
    fig.savefig(paths[1], bbox_inches="tight", facecolor="white")
    fig.savefig(paths[2], dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(paths[3], dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return paths


def render_figure_0(data_dir: Path, output_dir: Path) -> list[Path]:
    family = read_source(data_dir, "figure0_seed_family.tsv")
    evidence = read_source(data_dir, "figure0_seed_evidence.tsv")
    split = read_source(data_dir, "figure0_seed_split.tsv")
    order = ["ePhaZ", "iPhaZ", "OH", "BdhA", "phasin"]
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH_IN, 2.55), gridspec_kw={"width_ratios": [1, 1, 1.25]})
    status = family.pivot(index="family", columns="review_status", values="seeds").fillna(0).reindex(order).fillna(0)
    y = np.arange(len(order))
    reviewed = status.get("reviewed", pd.Series(0, index=status.index))
    axes[0].barh(y, reviewed, color="#42949E", label="reviewed", zorder=2)
    axes[0].barh(y, status.get("unreviewed", pd.Series(0, index=status.index)), left=reviewed, color="#D8D8D8", label="unreviewed", zorder=2)
    axes[0].set(yticks=y, yticklabels=[family_label(x) for x in order], xlabel="Seed entries", title="Review status")
    axes[0].invert_yaxis(); axes[0].legend(loc="lower right", fontsize=6); style_axis(axes[0]); add_panel_label(axes[0], "a")

    allocation = split.pivot(index="family", columns="split", values="seeds").fillna(0).reindex(order).fillna(0)
    train = allocation.get("train", pd.Series(0, index=allocation.index))
    axes[1].barh(y, train, color="#0F4D92", label="train", zorder=2)
    axes[1].barh(y, allocation.get("validation", pd.Series(0, index=allocation.index)), left=train, color="#B9DAD7", label="validation", zorder=2)
    axes[1].set(yticks=y, yticklabels=[], xlabel="Seed entries", title="Allocation")
    axes[1].invert_yaxis(); axes[1].legend(loc="lower right", fontsize=6); style_axis(axes[1]); add_panel_label(axes[1], "b")

    evidence_pivot = evidence.pivot(index="family", columns="evidence_source", values="seeds").fillna(0).reindex(order).fillna(0)
    starts = np.zeros(len(order))
    evidence_colors = {"PMID and DOI": "#0F4D92", "PMID only": "#42949E", "DOI only": "#9A4D8E", "No linked ID": "#D8D8D8"}
    for key in ["PMID and DOI", "PMID only", "DOI only", "No linked ID"]:
        values = evidence_pivot.get(key, pd.Series(0, index=evidence_pivot.index)).to_numpy()
        axes[2].barh(y, values, left=starts, color=evidence_colors[key], label=key, zorder=2)
        starts += values
    axes[2].set(yticks=y, yticklabels=[], xlabel="Seed entries", title="Evidence-source record")
    axes[2].invert_yaxis(); axes[2].legend(loc="lower right", fontsize=5.6); style_axis(axes[2]); add_panel_label(axes[2], "c")
    fig.suptitle("Curated seed library composition and evidence provenance", x=0.5, y=1.04, fontsize=9, fontweight="bold")
    fig.subplots_adjust(wspace=0.25, top=0.78, bottom=0.22, left=0.10, right=0.99)
    return save_exports(fig, output_dir, "figure_0_seed_library")


def render_figure_1(data_dir: Path, output_dir: Path) -> list[Path]:
    funnel = read_source(data_dir, "figure1_funnel.tsv")
    genomes = read_source(data_dir, "figure1_genome_coverage.tsv")
    fig = plt.figure(figsize=(WIDTH_IN, 3.8))
    grid = fig.add_gridspec(1, 2, width_ratios=[1.15, 1], wspace=0.32)
    workflow = fig.add_subplot(grid[0])
    right = grid[1].subgridspec(2, 1, height_ratios=[0.78, 1.35], hspace=0.62)
    genome_ax, protein_ax = fig.add_subplot(right[0]), fig.add_subplot(right[1])
    workflow.axis("off")
    boxes = [
        "GTDB R232\nrepresentative genomes",
        "Protein prediction\nand QC",
        "Nine-family HMM\nscreening",
        "Family-aware coverage\nand arbitration",
        "Sequence validation\nand tier rescoring",
    ]
    y_positions = np.linspace(0.88, 0.12, len(boxes))
    for index, (text, y) in enumerate(zip(boxes, y_positions)):
        box = FancyBboxPatch((0.18, y - 0.065), 0.64, 0.13, boxstyle="round,pad=0.012,rounding_size=0.015", linewidth=0.8, edgecolor="#4D4D4D", facecolor="#F7F7F7")
        workflow.add_patch(box)
        workflow.text(0.5, y, text, ha="center", va="center", fontsize=7)
        if index < len(boxes) - 1:
            workflow.add_patch(FancyArrowPatch((0.5, y - 0.07), (0.5, y_positions[index + 1] + 0.07), arrowstyle="-|>", mutation_scale=10, linewidth=0.8, color="#4D4D4D"))
    workflow.text(0.5, 1.03, "Screening workflow", ha="center", va="bottom", fontsize=8, fontweight="bold", transform=workflow.transAxes)
    add_panel_label(workflow, "a")

    genome_values = genomes["value"].astype(float).to_numpy()
    genome_y = np.arange(len(genomes))[::-1]
    genome_ax.barh(genome_y, genome_values, color=["#D8D8D8", "#0F4D92"], zorder=2)
    for yi, value in zip(genome_y, genome_values):
        genome_ax.text(value, yi, f" {int(value):,}", va="center", fontsize=6.2)
    genome_ax.set(
        yticks=genome_y,
        yticklabels=[fill(stage.replace("Genomes with ", ""), 26) for stage in genomes["stage"]],
        xlabel="Genomes",
        title="Genome-level coverage",
    )
    coverage = genome_values[1] / genome_values[0]
    genome_ax.text(0.99, -0.48, f"{coverage * 100:.1f}% with >=1 tier1 core candidate", transform=genome_ax.transAxes, ha="right", fontsize=6.2, color="#0F4D92", fontweight="bold")
    genome_ax.tick_params(axis="y", labelsize=5.8)
    style_axis(genome_ax); add_panel_label(genome_ax, "b")

    values = funnel["value"].astype(float).to_numpy()
    log_values = np.log10(values)
    y = np.arange(len(funnel))[::-1]
    bar_colors = ["#D8D8D8", "#C3D3DE", "#90B5C8", "#42949E", "#0F4D92"]
    protein_ax.barh(y, log_values, color=bar_colors, zorder=2)
    for yi, value, row in zip(y, values, funnel.itertuples(index=False)):
        protein_ax.text(np.log10(value) + 0.08, yi, f"{int(value):,}\n{row.unit}", va="center", fontsize=5.8)
    protein_ax.set(yticks=y, yticklabels=[fill(label, 22) for label in funnel["stage"]], xlabel="log10(count)", title="Protein-level screening")
    protein_ax.set_xlim(0, max(log_values) + 0.9); protein_ax.tick_params(axis="y", labelsize=5.8)
    style_axis(protein_ax); add_panel_label(protein_ax, "c")
    fig.suptitle("Genome and protein-level screening of GTDB candidate homologs", x=0.5, y=1.02, fontsize=9, fontweight="bold")
    fig.subplots_adjust(top=0.84, bottom=0.13, left=0.08, right=0.98)
    return save_exports(fig, output_dir, "figure_1_workflow_funnel")


def render_figure_2(data_dir: Path, output_dir: Path) -> list[Path]:
    scale = read_source(data_dir, "figure2_core_scale.tsv")
    union = read_source(data_dir, "figure2_union.tsv").iloc[0]
    order = CORE_FAMILIES
    scale["family"] = pd.Categorical(scale["family"], categories=order, ordered=True)
    scale = scale.sort_values("family")
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH_IN, 2.75), gridspec_kw={"width_ratios": [1, 1, 0.86]})
    x = np.arange(len(scale))
    colors = [COLORS[family] for family in scale["family"]]
    for ax, column, title in zip(axes[:2], ["sequences", "genomes"], ["Tier1 sequences", "Tier1 genomes"]):
        values = scale[column].astype(int).to_numpy()
        ax.bar(x, values, color=colors, zorder=2)
        for xi, value in zip(x, values):
            ax.text(xi, value, f"{value:,}", ha="center", va="bottom", fontsize=6, rotation=90)
        ax.set(xticks=x, xticklabels=[family_label(f) for f in scale["family"]], ylabel="Count", title=title)
        ax.tick_params(axis="x", labelsize=6); style_axis(ax); add_panel_label(ax, "a" if column == "sequences" else "b")
    denominator, numerator = int(union["gtdb_representatives"]), int(union["core_family_union_genomes"])
    axes[2].barh([0], [denominator], color="#E0E0E0", height=0.35, zorder=1)
    axes[2].barh([0], [numerator], color="#0F4D92", height=0.35, zorder=2)
    axes[2].text(numerator / 2, 0, f"{numerator:,}", color="white", ha="center", va="center", fontsize=8, fontweight="bold")
    axes[2].text(denominator, 0.28, f"{denominator:,} GTDB representatives", ha="right", va="bottom", fontsize=6)
    axes[2].set(xlim=(0, denominator), yticks=[], xlabel="Genomes", title="Four-family genome union")
    axes[2].text(0.5, -0.36, f"{float(union['coverage_fraction']) * 100:.1f}% coverage", transform=axes[2].transAxes, ha="center", fontsize=7, color="#0F4D92", fontweight="bold")
    style_axis(axes[2]); add_panel_label(axes[2], "c")
    fig.suptitle("Core candidate-family scale and genome-union coverage", x=0.5, y=1.03, fontsize=9, fontweight="bold")
    fig.subplots_adjust(wspace=0.35, top=0.78, bottom=0.26, left=0.08, right=0.99)
    return save_exports(fig, output_dir, "figure_2_core_scale")


def render_figure_3(data_dir: Path, output_dir: Path) -> list[Path]:
    matrix = read_source(data_dir, "figure3_phylum_matrix.tsv")
    totals = read_source(data_dir, "figure3_phylum_totals.tsv")
    phyla = list(totals["phylum"])
    pivot = matrix.pivot(index="phylum", columns="family", values="genomes").reindex(index=phyla, columns=CORE_FAMILIES).fillna(0)
    fig = plt.figure(figsize=(WIDTH_IN, 4.15))
    grid = fig.add_gridspec(1, 2, width_ratios=[0.8, 1.25], wspace=0.35)
    bars, heat = fig.add_subplot(grid[0]), fig.add_subplot(grid[1])
    y = np.arange(len(totals))[::-1]
    bars.barh(y, totals["genomes"], color="#42949E", zorder=2)
    bars.set(yticks=y, yticklabels=totals["phylum"], xlabel="Core candidate genomes", title="Leading phyla")
    bars.tick_params(axis="y", labelsize=5.6); style_axis(bars); add_panel_label(bars, "a")
    values = np.log10(pivot.to_numpy(dtype=float) + 1)
    image = heat.imshow(values, cmap=HEATMAP, aspect="auto", vmin=0, vmax=max(1, values.max()))
    heat.set(xticks=np.arange(len(CORE_FAMILIES)), xticklabels=[family_label(x) for x in CORE_FAMILIES], yticks=np.arange(len(phyla)), yticklabels=phyla, title="Core-family matrix")
    heat.tick_params(axis="x", labelsize=5.5); heat.tick_params(axis="y", labelsize=5.6, length=0)
    for row in range(values.shape[0]):
        for col in range(values.shape[1]):
            count = int(pivot.iloc[row, col])
            if count:
                color = "white" if values[row, col] > values.max() * 0.55 else "#272727"
                heat.text(col, row, f"{count:,}", ha="center", va="center", fontsize=4.7, color=color)
    colorbar = fig.colorbar(image, ax=heat, fraction=0.04, pad=0.03)
    colorbar.set_label("log10(genomes + 1)", fontsize=6); colorbar.ax.tick_params(labelsize=5.5)
    add_panel_label(heat, "b")
    fig.suptitle("Phylum-level distribution of core candidate homologs", x=0.5, y=0.995, fontsize=9, fontweight="bold")
    fig.subplots_adjust(top=0.89, bottom=0.12, left=0.12, right=0.96)
    return save_exports(fig, output_dir, "figure_3_phylum_distribution")


def render_figure_4(data_dir: Path, output_dir: Path) -> list[Path]:
    classes = read_source(data_dir, "figure4_signal_classes.tsv")
    phyla = read_source(data_dir, "figure4_signal_phyla.tsv")
    order = ["SP", "LIPO", "TAT", "TATLIPO", "OTHER"]
    classes["signal_class"] = pd.Categorical(classes["signal_class"], categories=order, ordered=True)
    classes = classes.sort_values("signal_class")
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH_IN, 3.45), gridspec_kw={"width_ratios": [0.9, 1.1]})
    y = np.arange(len(classes))[::-1]
    signal_colors = {"SP": "#0F4D92", "LIPO": "#42949E", "TAT": "#9A4D8E", "TATLIPO": "#B64342", "OTHER": "#D8D8D8"}
    axes[0].barh(y, classes["sequences"], color=[signal_colors[x] for x in classes["signal_class"]], zorder=2)
    for yi, value in zip(y, classes["sequences"]):
        axes[0].text(value, yi, f" {int(value):,}", va="center", fontsize=6)
    axes[0].set(yticks=y, yticklabels=classes["signal_class"], xlabel="Tier1 ePhaZ sequences", title="SignalP prediction class")
    style_axis(axes[0]); add_panel_label(axes[0], "a")
    py = np.arange(len(phyla))[::-1]
    axes[1].barh(py, phyla["genomes"], color="#42949E", zorder=2)
    axes[1].set(yticks=py, yticklabels=phyla["phylum"], xlabel="Candidate genomes", title="Predicted signal-peptide candidates")
    axes[1].tick_params(axis="y", labelsize=5.7); style_axis(axes[1]); add_panel_label(axes[1], "b")
    fig.text(0.5, 0.01, "SignalP output is prediction evidence for tier1 ePhaZ candidates.", ha="center", fontsize=6, color="#4D4D4D")
    fig.suptitle("ePhaZ signal-peptide prediction and phylum representation", x=0.5, y=0.99, fontsize=9, fontweight="bold")
    fig.subplots_adjust(wspace=0.44, top=0.83, bottom=0.16, left=0.11, right=0.99)
    return save_exports(fig, output_dir, "figure_4_ephaz_signal")


def render_figure_5(data_dir: Path, output_dir: Path) -> list[Path]:
    rates = read_source(data_dir, "figure5_neighborhood_rate.tsv")
    patatin = read_source(data_dir, "figure5_patatin_context.tsv")
    core = rates.loc[rates["hit_family"].isin(CORE_FAMILIES)]
    pivot = core.pivot(index="hit_family", columns="marker_family", values="support_rate").reindex(index=CORE_FAMILIES, columns=MARKERS).fillna(0)
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH_IN, 3.55), gridspec_kw={"width_ratios": [1.2, 0.85]})
    values = pivot.to_numpy(dtype=float) * 100
    image = axes[0].imshow(values, cmap=HEATMAP, aspect="auto", vmin=0, vmax=max(1, values.max()))
    axes[0].set(xticks=np.arange(len(MARKERS)), xticklabels=MARKERS, yticks=np.arange(len(CORE_FAMILIES)), yticklabels=[family_label(x) for x in CORE_FAMILIES], title="Core candidate families")
    axes[0].tick_params(axis="x", labelrotation=35, labelsize=6); axes[0].tick_params(axis="y", labelsize=6, length=0)
    for row in range(values.shape[0]):
        for col in range(values.shape[1]):
            value = values[row, col]
            text_color = "white" if value > values.max() * 0.55 else "#272727"
            axes[0].text(col, row, f"{value:.1f}", ha="center", va="center", fontsize=5.3, color=text_color)
    colorbar = fig.colorbar(image, ax=axes[0], fraction=0.045, pad=0.03)
    colorbar.set_label("Support rate (%)", fontsize=6); colorbar.ax.tick_params(labelsize=5.5)
    add_panel_label(axes[0], "a")
    patatin = patatin.set_index("marker_family").reindex(MARKERS).reset_index()
    py = np.arange(len(patatin))[::-1]
    axes[1].barh(py, patatin["support_rate"] * 100, color="#767676", zorder=2)
    for yi, row in zip(py, patatin.itertuples(index=False)):
        axes[1].text(row.support_rate * 100, yi, f" {int(row.supported_loci):,} loci", va="center", fontsize=5.5)
    axes[1].set(yticks=py, yticklabels=patatin["marker_family"], xlabel="Support rate (%)", title="ArchPhaZ patatin candidates")
    style_axis(axes[1]); add_panel_label(axes[1], "b")
    fig.text(0.5, 0.01, "Genomic neighborhood support within +/-10 kb; denominator: analyzed candidate loci.", ha="center", fontsize=6, color="#4D4D4D")
    fig.suptitle("Current genomic neighborhood context of candidate loci", x=0.5, y=0.99, fontsize=9, fontweight="bold")
    fig.subplots_adjust(wspace=0.38, top=0.83, bottom=0.20, left=0.10, right=0.98)
    return save_exports(fig, output_dir, "figure_5_neighborhood_context")


def render_all(data_dir: Path, output_dir: Path) -> list[Path]:
    outputs: list[Path] = []
    for renderer in (render_figure_0, render_figure_1, render_figure_2, render_figure_3, render_figure_4, render_figure_5):
        outputs.extend(renderer(data_dir, output_dir))
    return outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    outputs = render_all(**vars(parse_args()))
    print(f"Wrote {len(outputs)} figure exports.")


if __name__ == "__main__":
    main()
