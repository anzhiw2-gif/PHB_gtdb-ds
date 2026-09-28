#!/usr/bin/env python3
"""Render the 2026-09-06 candidate-only PHB/PHA project figure suite."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
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
PALETTE = {
    "blue": "#0F4D92",
    "teal": "#42949E",
    "green": "#78A96B",
    "red": "#B64342",
    "violet": "#9A4D8E",
    "gold": "#C99000",
    "gray": "#767676",
    "light": "#D8D8D8",
    "dark": "#272727",
}
SUBTYPE_LABELS = {
    "ePhaZ_core_nonsecreted_evidence": "ePhaZ core nonsecreted",
    "ePhaZ_SCL_secreted_evidence": "ePhaZ SCL secreted",
    "ePhaZ_tier2_review_evidence": "ePhaZ tier2 review",
    "ePhaZ_competition_review_evidence": "ePhaZ competition review",
    "iPhaZ_intracellular_evidence": "iPhaZ intracellular",
    "iPhaZ_tier2_review_evidence": "iPhaZ tier2 review",
}
SUBTYPE_ORDER = list(SUBTYPE_LABELS)
LAYER_ORDER = [
    "PhaZh1_like_high",
    "PhaZh1_like_review",
    "archaeal_patatin_exploratory",
    "non_PhaZh1_patatin",
]
LAYER_LABELS = {
    "PhaZh1_like_high": "PhaZh1-like high",
    "PhaZh1_like_review": "PhaZh1-like review",
    "archaeal_patatin_exploratory": "patatin exploratory",
    "non_PhaZh1_patatin": "non-PhaZh1 patatin",
}
LAYER_COLORS = {
    "PhaZh1_like_high": PALETTE["blue"],
    "PhaZh1_like_review": PALETTE["gold"],
    "archaeal_patatin_exploratory": PALETTE["light"],
    "non_PhaZh1_patatin": PALETTE["red"],
}
CORE_FAMILY_ROWS = [
    ("ePhaZ curated core", 5646, 5080),
    ("iPhaZ", 32226, 25564),
    ("OH", 3570, 3446),
    ("ArchPhaZ hydrolase", 14571, 12469),
]


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing or empty source table: {path}")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def subtype_counts(frame: pd.DataFrame) -> pd.DataFrame:
    counts = frame["pha_subtype_evidence"].value_counts()
    rows = [
        {"source_category": key, "category": SUBTYPE_LABELS[key], "count": int(counts.get(key, 0))}
        for key in SUBTYPE_ORDER
        if counts.get(key, 0)
    ]
    return pd.DataFrame(rows)


def species_from_taxonomy(taxonomy: str) -> str:
    for token in taxonomy.split(";"):
        if token.startswith("s__") and token[3:].strip():
            return token[3:].strip()
    return "unclassified"


def validate_export_bundle(output_dir: Path, stems: list[str]) -> list[str]:
    return [
        f"{stem}{suffix}"
        for stem in stems
        for suffix in (".svg", ".pdf", ".tiff", ".png")
        if not (output_dir / f"{stem}{suffix}").is_file()
        or (output_dir / f"{stem}{suffix}").stat().st_size == 0
    ]


def save_exports(fig: plt.Figure, out_dir: Path, stem: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for suffix, dpi in ((".svg", None), (".pdf", None), (".tiff", 600), (".png", 300)):
        kwargs: dict[str, object] = {"bbox_inches": "tight", "facecolor": "white"}
        if dpi:
            kwargs["dpi"] = dpi
        fig.savefig(out_dir / f"{stem}{suffix}", **kwargs)
    plt.close(fig)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(0.01, 0.99, label, transform=ax.transAxes, fontweight="bold", fontsize=8, va="top", zorder=10)


def bar_label(value: int, ymax: float) -> tuple[str, str, str]:
    """Return text, placement and colour for a bar annotation.

    Large bars are labelled inside the bar so their vertical annotations cannot
    collide with the panel title; small bars remain outside for legibility.
    """
    text = f"{int(value):,}"
    if value >= ymax * 0.45:
        return text, "inside", "white"
    return text, "outside", PALETTE["dark"]


def style_axis(ax: plt.Axes, grid: str = "x") -> None:
    ax.grid(axis=grid, color="#E8E8E8", linewidth=0.6, zorder=0)
    ax.tick_params(length=2.5, width=0.7, pad=2)


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def resolve_input_tables(run_dir: Path) -> dict[str, Path]:
    """Resolve all renderer inputs from the immutable current-run snapshot."""
    input_dir = run_dir / "inputs"
    names = {
        "overall": "overall_pha_phb.tsv",
        "archaea": "archaeal_phazh1_full_audit.tsv",
        "domain": "domain_annotation.tsv",
        "calibration": "positive_negative_calibration.tsv",
        "manual": "cross_phylum_manual_review.tsv",
    }
    paths = {key: input_dir / name for key, name in names.items()}
    missing = [str(path) for path in paths.values() if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError("missing current-run input snapshot: " + ", ".join(missing))
    return paths


def marker_matrix(archaea: pd.DataFrame) -> pd.DataFrame:
    markers = ["BdhA", "PhaJ", "PhaC", "PhaE"]
    rows = []
    for layer in LAYER_ORDER:
        subset = archaea.loc[archaea["layer"] == layer]
        for marker in markers:
            count = sum(marker in value.split(";") for value in subset["nearby_markers"])
            rows.append({"layer": layer, "marker": marker, "count": count, "denominator": len(subset), "fraction": count / len(subset) if len(subset) else 0})
    return pd.DataFrame(rows)


def build_source_data(run_dir: Path, source_dir: Path) -> dict[str, pd.DataFrame]:
    inputs = resolve_input_tables(run_dir)
    overall = read_tsv(inputs["overall"])
    archaea = read_tsv(inputs["archaea"])
    domain = read_tsv(inputs["domain"])
    calibration = read_tsv(inputs["calibration"])
    manual = read_tsv(inputs["manual"])

    workflow = pd.DataFrame([
        {"step": "Formal HMM scan", "unit": "tasks", "count": 1000, "status": "completed", "note": "1000/1000; 0 failed"},
        {"step": "Strict tier1 core union", "unit": "genomes", "count": 38741, "status": "completed_candidate_only", "note": "four core families"},
        {"step": "PHA evidence layering", "unit": "proteins", "count": len(overall), "status": "completed_candidate_only", "note": "deduplicated ePhaZ/iPhaZ"},
        {"step": "Archaeal PhaZh1 audit", "unit": "records", "count": len(archaea), "status": "completed_candidate_only", "note": "full archaeal patatin audit"},
        {"step": "Cross-HMM domain annotation", "unit": "candidates", "count": len(domain), "status": "completed_candidate_only", "note": "197 high + 7 review"},
        {"step": "Manual cross-phylum review", "unit": "records", "count": len(manual), "status": "completed_candidate_only", "note": "targeted review subset"},
    ])
    subtype = subtype_counts(overall)
    core = pd.DataFrame(CORE_FAMILY_ROWS, columns=["family", "tier1_sequences", "tier1_genomes"])
    core_union = pd.DataFrame([{"metric": "strict_four_family_union", "genomes": 38741, "multi_family_genomes": 6578, "evidence": "tier1_computational_candidate"}])
    phylum = archaea.groupby(["phylum", "layer"], as_index=False).size().rename(columns={"size": "records"})
    phylum["layer"] = pd.Categorical(phylum["layer"], LAYER_ORDER, ordered=True)
    phylum = phylum.sort_values(["phylum", "layer"])
    audit = archaea.groupby("layer", as_index=False).size().rename(columns={"size": "records"})
    audit["layer"] = pd.Categorical(audit["layer"], LAYER_ORDER, ordered=True)
    audit = audit.sort_values("layer")
    motif = domain.groupby(["layer", "nterm_gxsxg_motif"], as_index=False).size().rename(columns={"size": "records"})
    high_review = domain.copy()
    high_review["species"] = high_review["taxonomy"].map(species_from_taxonomy)
    species = high_review.groupby("species", as_index=False).size().rename(columns={"size": "records"})
    species = species.sort_values(["records", "species"], ascending=[False, True]).head(15)
    calibration["bitscore"] = calibration["bitscore"].astype(float)
    marker = marker_matrix(archaea)
    data = {
        "figure_1_workflow_status.tsv": workflow,
        "figure_2_subtype_distribution.tsv": subtype,
        "figure_3_core_family_scale.tsv": core,
        "figure_3_core_union.tsv": core_union,
        "figure_4_archaeal_phylum_distribution.tsv": phylum,
        "figure_5_archaeal_audit_layers.tsv": audit,
        "figure_5_archaeal_motif.tsv": motif,
        "figure_6_top_species.tsv": species,
        "figure_6_manual_review.tsv": manual,
        "supplementary_figure_1_calibration.tsv": calibration,
        "supplementary_figure_2_neighborhood.tsv": marker,
    }
    for name, frame in data.items():
        write_tsv(frame, source_dir / name)
    return data


def figure_1(data: dict[str, pd.DataFrame], out: Path) -> None:
    table = data["figure_1_workflow_status.tsv"]
    fig = plt.figure(figsize=(WIDTH_IN, 4.0))
    grid = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.1], wspace=0.30)
    flow, counts = fig.add_subplot(grid[0]), fig.add_subplot(grid[1])
    flow.axis("off")
    y = np.linspace(0.88, 0.12, len(table))
    for index, (row, ypos) in enumerate(zip(table.itertuples(index=False), y)):
        color = PALETTE["blue"] if row.status == "completed" else "#DCEAF5"
        box = FancyBboxPatch((0.10, ypos - 0.052), 0.80, 0.104, boxstyle="round,pad=0.012,rounding_size=0.018", linewidth=0.8, edgecolor=PALETTE["gray"], facecolor=color)
        flow.add_patch(box)
        text_color = "white" if row.status == "completed" else PALETTE["dark"]
        flow.text(0.50, ypos + 0.013, row.step, ha="center", va="center", fontsize=7, color=text_color, fontweight="bold")
        flow.text(0.50, ypos - 0.018, row.note, ha="center", va="center", fontsize=5.5, color=text_color)
        if index < len(table) - 1:
            flow.add_patch(FancyArrowPatch((0.50, ypos - 0.065), (0.50, y[index + 1] + 0.065), arrowstyle="-|>", mutation_scale=10, linewidth=0.8, color=PALETTE["gray"]))
    flow.set_title("Auditable analysis workflow", fontsize=8, fontweight="bold", pad=12)
    panel_label(flow, "a")
    x = np.arange(len(table))
    colors = [PALETTE["blue"]] + [PALETTE["teal"]] * (len(table) - 1)
    counts.bar(x, np.log10(table["count"].astype(float)), color=colors, zorder=2)
    for xi, row in zip(x, table.itertuples(index=False)):
        counts.text(xi, np.log10(row.count), f"{row.count:,}\n{row.unit}", ha="center", va="bottom", fontsize=5.6, rotation=90)
    counts.set(xticks=x, xticklabels=["scan", "tier1", "layering", "archaea", "domain", "review"], ylabel="log10(record count)", title="Completed output scale")
    counts.tick_params(axis="x", labelsize=5.8); style_axis(counts); panel_label(counts, "b")
    fig.text(0.5, 0.01, "All stages shown here report candidate homology/function-potential evidence, not validated PHB degradation phenotypes.", ha="center", fontsize=6, color=PALETTE["gray"])
    fig.suptitle("PHB/PHA candidate-analysis workflow and completed outputs", y=0.99, fontsize=10, fontweight="bold")
    fig.subplots_adjust(top=0.84, bottom=0.16, left=0.08, right=0.98)
    save_exports(fig, out, "figure_1_workflow_and_status")


def figure_2(data: dict[str, pd.DataFrame], out: Path) -> None:
    table = data["figure_2_subtype_distribution.tsv"].sort_values("count")
    fig, ax = plt.subplots(figsize=(WIDTH_IN, 3.4))
    colors = [PALETTE["blue"], PALETTE["teal"], PALETTE["light"], PALETTE["gold"], PALETTE["green"], PALETTE["red"]]
    y = np.arange(len(table))
    ax.barh(y, table["count"], color=colors[:len(table)], zorder=2)
    for yi, row in zip(y, table.itertuples(index=False)):
        ax.text(row.count, yi, f" {row.count:,}", va="center", fontsize=6.5)
    ax.set(yticks=y, yticklabels=table["category"], xlabel="Deduplicated candidate proteins", title="Evidence-stratified PHA candidate distribution (n = 109,087)")
    ax.tick_params(axis="y", labelsize=6.5); style_axis(ax); panel_label(ax, "a")
    ax.text(0.99, -0.18, "Layers rank computational evidence; they are not experimental activity classes.", transform=ax.transAxes, ha="right", fontsize=6, color=PALETTE["gray"])
    fig.subplots_adjust(top=0.85, bottom=0.20, left=0.25, right=0.96)
    save_exports(fig, out, "figure_2_pha_subtype_distribution")


def figure_3(data: dict[str, pd.DataFrame], out: Path) -> None:
    core, union = data["figure_3_core_family_scale.tsv"], data["figure_3_core_union.tsv"].iloc[0]
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH_IN, 2.9), gridspec_kw={"width_ratios": [1, 1, 0.88]})
    x = np.arange(len(core))
    colors = [PALETTE["blue"], PALETTE["teal"], PALETTE["red"], PALETTE["violet"]]
    for ax, column, title, label in zip(axes[:2], ["tier1_sequences", "tier1_genomes"], ["Tier1 sequences", "Tier1 genomes"], ["a", "b"]):
        values = core[column].astype(int)
        ax.bar(x, values, color=colors, zorder=2)
        ymax = float(values.max())
        for xi, value in zip(x, values):
            text, placement, colour = bar_label(int(value), ymax)
            if placement == "inside":
                ax.text(xi, value * 0.70, text, ha="center", va="center", fontsize=5.6, rotation=90, color=colour, fontweight="bold")
            else:
                ax.text(xi, value, text, ha="center", va="bottom", fontsize=5.6, rotation=90, color=colour)
        ax.set(xticks=x, xticklabels=[name.replace(" ", "\n") for name in core["family"]], ylabel="Count", title=title)
        ax.tick_params(axis="x", labelsize=5.5); style_axis(ax); panel_label(ax, label)
    total = int(union.genomes)
    multi = int(union.multi_family_genomes)
    axes[2].barh([0], [total], color=PALETTE["blue"], height=0.38, zorder=2)
    axes[2].text(total / 2, 0, f"{total:,}\nstrict union", ha="center", va="center", color="white", fontsize=8, fontweight="bold")
    axes[2].text(total, 0.28, f"{multi:,} genomes with >=2 families", ha="right", fontsize=5.8)
    axes[2].set(xlim=(0, total * 1.08), yticks=[], xlabel="Genomes", title="Four-family union")
    style_axis(axes[2]); panel_label(axes[2], "c")
    fig.suptitle("Strict tier1 core-family candidate scale", y=1.02, fontsize=10, fontweight="bold")
    fig.text(0.5, 0.01, "Strict tier1 is a high-confidence computational candidate layer, not phenotype confirmation.", ha="center", fontsize=6, color=PALETTE["gray"])
    fig.subplots_adjust(top=0.80, bottom=0.22, left=0.08, right=0.98, wspace=0.42)
    save_exports(fig, out, "figure_3_core_family_scale")


def figure_4(data: dict[str, pd.DataFrame], out: Path) -> None:
    table = data["figure_4_archaeal_phylum_distribution.tsv"]
    pivot = table.pivot(index="phylum", columns="layer", values="records").reindex(columns=LAYER_ORDER, fill_value=0).fillna(0)
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH_IN, 3.2), gridspec_kw={"width_ratios": [1.25, 0.95]})
    y = np.arange(len(pivot))
    left = np.zeros(len(pivot))
    for layer in LAYER_ORDER:
        values = pivot[layer].to_numpy(dtype=float)
        axes[0].barh(y, values, left=left, color=LAYER_COLORS[layer], label=LAYER_LABELS[layer], zorder=2)
        left += values
    axes[0].set(yticks=y, yticklabels=pivot.index, xlabel="Archaeal audit records", title="Phylum distribution")
    axes[0].legend(fontsize=5.4, loc="lower right"); style_axis(axes[0]); panel_label(axes[0], "a")
    high_review = pivot[["PhaZh1_like_high", "PhaZh1_like_review"]]
    for index, layer in enumerate(high_review.columns):
        axes[1].bar(np.arange(len(high_review)) + (index - 0.5) * 0.32, high_review[layer], width=0.30, color=LAYER_COLORS[layer], label=LAYER_LABELS[layer], zorder=2)
    axes[1].set(xticks=np.arange(len(high_review)), xticklabels=high_review.index, ylabel="Records", title="Priority high/review candidates")
    axes[1].legend(fontsize=5.6); style_axis(axes[1]); panel_label(axes[1], "b")
    fig.suptitle("Archaeal PhaZh1-like candidate distribution by phylum", y=0.99, fontsize=10, fontweight="bold")
    fig.text(0.5, 0.01, "Taxonomy is retained from audit_04; high/review labels remain candidate-only review priorities.", ha="center", fontsize=6, color=PALETTE["gray"])
    fig.subplots_adjust(top=0.82, bottom=0.19, left=0.12, right=0.98, wspace=0.34)
    save_exports(fig, out, "figure_4_archaeal_phylum_distribution")


def figure_5(data: dict[str, pd.DataFrame], out: Path) -> None:
    layers, motif = data["figure_5_archaeal_audit_layers.tsv"], data["figure_5_archaeal_motif.tsv"]
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH_IN, 3.15), gridspec_kw={"width_ratios": [1.05, 0.95]})
    x = np.arange(len(layers))
    values = layers["records"].astype(int)
    axes[0].bar(x, values, color=[LAYER_COLORS[layer] for layer in layers["layer"]], zorder=2)
    for xi, value in zip(x, values):
        axes[0].text(xi, value, f"{value:,}", ha="center", va="bottom", fontsize=6, rotation=90)
    axes[0].set(xticks=x, xticklabels=[LAYER_LABELS[layer].replace(" ", "\n") for layer in layers["layer"]], ylabel="Archaeal records", title="Full PhaZh1-like audit (n = 2,307)")
    axes[0].tick_params(axis="x", labelsize=5.6); style_axis(axes[0]); panel_label(axes[0], "a")
    pivot = motif.pivot(index="layer", columns="nterm_gxsxg_motif", values="records").reindex(["PhaZh1_like_high", "PhaZh1_like_review"]).fillna(0)
    y = np.arange(len(pivot))
    present = pivot.get("present", pd.Series(0, index=pivot.index))
    axes[1].barh(y, present, color=PALETTE["green"], label="GxSxG present", zorder=2)
    axes[1].barh(y, pivot.get("absent", pd.Series(0, index=pivot.index)), left=present, color=PALETTE["light"], label="absent", zorder=2)
    axes[1].set(yticks=y, yticklabels=[LAYER_LABELS[layer] for layer in pivot.index], xlabel="Domain-annotated candidates", title="N-terminal motif audit")
    axes[1].legend(fontsize=5.8); style_axis(axes[1]); panel_label(axes[1], "b")
    fig.suptitle("Archaeal PhaZh1-like audit strata and N-terminal motif evidence", y=0.99, fontsize=10, fontweight="bold")
    fig.text(0.5, 0.01, "Motif presence is sequence-feature evidence only and does not establish enzymatic phenotype.", ha="center", fontsize=6, color=PALETTE["gray"])
    fig.subplots_adjust(top=0.82, bottom=0.20, left=0.10, right=0.98, wspace=0.34)
    save_exports(fig, out, "figure_5_archaeal_audit_and_motif")


def figure_6(data: dict[str, pd.DataFrame], out: Path) -> None:
    species, manual = data["figure_6_top_species.tsv"], data["figure_6_manual_review.tsv"].copy()
    manual["coverage"] = manual["coverage"].astype(float)
    manual["neglog_evalue"] = -np.log10(manual["evalue"].astype(float))
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH_IN, 4.0), gridspec_kw={"width_ratios": [1.05, 0.95]})
    species = species.sort_values("records")
    y = np.arange(len(species))
    axes[0].barh(y, species["records"], color=PALETTE["blue"], zorder=2)
    axes[0].set(yticks=y, yticklabels=species["species"], xlabel="High/review candidate records", title="Leading annotated species")
    axes[0].tick_params(axis="y", labelsize=5.2); style_axis(axes[0]); panel_label(axes[0], "a")
    for layer, color, marker in [("PhaZh1_like_high", PALETTE["blue"], "o"), ("PhaZh1_like_review", PALETTE["gold"], "s")]:
        subset = manual.loc[manual["layer"] == layer]
        axes[1].scatter(subset["coverage"], subset["neglog_evalue"], c=color, marker=marker, edgecolor="white", linewidth=0.5, s=40, label=LAYER_LABELS[layer], zorder=3)
    for row in manual.itertuples(index=False):
        axes[1].text(row.coverage + 0.003, row.neglog_evalue + 0.4, row.phylum, fontsize=4.8)
    axes[1].set(xlabel="Candidate-HMM coverage", ylabel="-log10(E-value)", title="Manual cross-phylum review (n = 9)")
    axes[1].legend(fontsize=5.6); style_axis(axes[1], grid="both"); panel_label(axes[1], "b")
    fig.suptitle("Species-level representation and targeted cross-phylum review", y=0.99, fontsize=10, fontweight="bold")
    fig.text(0.5, 0.01, "Species labels reflect taxonomy strings in audit_04; manual review records retain candidate-only disposition.", ha="center", fontsize=6, color=PALETTE["gray"])
    fig.subplots_adjust(top=0.83, bottom=0.16, left=0.18, right=0.98, wspace=0.38)
    save_exports(fig, out, "figure_6_species_and_manual_review")


def supplementary_1(data: dict[str, pd.DataFrame], out: Path) -> None:
    table = data["supplementary_figure_1_calibration.tsv"].copy()
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH_IN, 3.2), gridspec_kw={"width_ratios": [1.05, 0.95]})
    groups = ["positive_control", "negative_control"]
    colors = {"positive_control": PALETTE["green"], "negative_control": PALETTE["red"]}
    for index, group in enumerate(groups):
        subset = table.loc[table["panel"] == group]
        jitter = np.linspace(-0.10, 0.10, len(subset)) if len(subset) > 1 else np.array([0.0])
        axes[0].scatter(np.full(len(subset), index) + jitter, subset["bitscore"], s=35, c=colors[group], edgecolor="white", linewidth=0.5, zorder=3)
        axes[0].text(index, 0, f"{len(subset)} domain hits", ha="center", va="bottom", fontsize=6)
    axes[0].set(xticks=[0, 1], xticklabels=["Positive\ncontrol", "Negative\ncontrol"], ylabel="Domain bitscore", title="ArchPhaZ_patatin calibration")
    style_axis(axes[0]); panel_label(axes[0], "a")
    counts = table.groupby("panel").size().reindex(groups).reset_index(name="domain_hits")
    axes[1].bar([0, 1], counts["domain_hits"], color=[colors[x] for x in groups], zorder=2)
    for xi, value in enumerate(counts["domain_hits"]):
        axes[1].text(xi, value, str(value), ha="center", va="bottom", fontsize=7)
    axes[1].set(xticks=[0, 1], xticklabels=["Positive", "Negative"], ylabel="Reported domain hits", title="Panel-level result")
    style_axis(axes[1]); panel_label(axes[1], "b")
    fig.suptitle("Supplementary Figure 1 | Positive/negative HMM calibration", y=0.99, fontsize=10, fontweight="bold")
    fig.text(0.5, 0.01, "Negative-panel hits show that this registered model is not a zero-cross-reactivity discriminator; calibration does not establish phenotype.", ha="center", fontsize=6, color=PALETTE["gray"])
    fig.subplots_adjust(top=0.82, bottom=0.20, left=0.11, right=0.98, wspace=0.34)
    save_exports(fig, out, "supplementary_figure_1_hmm_calibration")


def supplementary_2(data: dict[str, pd.DataFrame], out: Path) -> None:
    table = data["supplementary_figure_2_neighborhood.tsv"]
    pivot = table.pivot(index="layer", columns="marker", values="fraction").reindex(index=LAYER_ORDER, columns=["BdhA", "PhaJ", "PhaC", "PhaE"]).fillna(0)
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH_IN, 3.35), gridspec_kw={"width_ratios": [1.05, 0.95]})
    image = axes[0].imshow(pivot.to_numpy() * 100, cmap="Blues", aspect="auto", vmin=0, vmax=100)
    axes[0].set(xticks=np.arange(len(pivot.columns)), xticklabels=pivot.columns, yticks=np.arange(len(pivot.index)), yticklabels=[LAYER_LABELS[layer] for layer in pivot.index], title="Marker occurrence by audit layer")
    axes[0].tick_params(axis="y", labelsize=5.8, length=0)
    for row in range(pivot.shape[0]):
        for col in range(pivot.shape[1]):
            value = pivot.iloc[row, col] * 100
            axes[0].text(col, row, f"{value:.1f}", ha="center", va="center", fontsize=5.5, color="white" if value > 50 else PALETTE["dark"])
    cbar = fig.colorbar(image, ax=axes[0], fraction=0.045, pad=0.03)
    cbar.set_label("Records with marker (%)", fontsize=6); cbar.ax.tick_params(labelsize=5.5)
    panel_label(axes[0], "a")
    total = table.groupby("marker", as_index=False).agg(count=("count", "sum"), denominator=("denominator", "sum"))
    total["fraction"] = total["count"] / total["denominator"]
    total = total.sort_values("fraction")
    y = np.arange(len(total))
    axes[1].barh(y, total["fraction"] * 100, color=PALETTE["teal"], zorder=2)
    axes[1].set(yticks=y, yticklabels=total["marker"], xlabel="All audit records with marker (%)", title="Overall marker context")
    style_axis(axes[1]); panel_label(axes[1], "b")
    fig.suptitle("Supplementary Figure 2 | Archaeal neighborhood-marker association", y=0.99, fontsize=10, fontweight="bold")
    fig.text(0.5, 0.01, "Marker co-occurrence is local genomic-context evidence, not proof of pathway membership, co-transcription, or PHB degradation.", ha="center", fontsize=6, color=PALETTE["gray"])
    fig.subplots_adjust(top=0.82, bottom=0.20, left=0.13, right=0.98, wspace=0.34)
    save_exports(fig, out, "supplementary_figure_2_neighborhood_markers")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args()
    source_dir = args.run_dir / "results/source_data"
    figures_dir = args.run_dir / "results/figures"
    data = build_source_data(args.run_dir, source_dir)
    figure_1(data, figures_dir)
    figure_2(data, figures_dir)
    figure_3(data, figures_dir)
    figure_4(data, figures_dir)
    figure_5(data, figures_dir)
    figure_6(data, figures_dir)
    supplementary_1(data, figures_dir)
    supplementary_2(data, figures_dir)
    stems = [
        "figure_1_workflow_and_status", "figure_2_pha_subtype_distribution", "figure_3_core_family_scale", "figure_4_archaeal_phylum_distribution",
        "figure_5_archaeal_audit_and_motif", "figure_6_species_and_manual_review", "supplementary_figure_1_hmm_calibration", "supplementary_figure_2_neighborhood_markers",
    ]
    missing = validate_export_bundle(figures_dir, stems)
    if missing:
        raise RuntimeError(f"missing figure exports: {', '.join(missing)}")
    manifest = {
        "run_id": args.run_dir.name,
        "status": "completed_candidate_only",
        "backend": "Python matplotlib",
        "figures": stems,
        "exports": [
            {"file": path.name, "size": path.stat().st_size, "sha256": sha256_file(path)}
            for path in sorted(figures_dir.iterdir())
            if path.is_file()
        ],
        "source_data_files": sorted(path.name for path in source_dir.glob("*.tsv")),
        "source_data_count": len(list(source_dir.glob("*.tsv"))),
        "phenotype_boundary": "HMM, domain, motif, neighborhood and taxonomy results are candidate homology/function-potential evidence, not validated PHB degradation phenotypes.",
    }
    (args.run_dir / "results/figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(stems)} figures and {len(manifest['source_data_files'])} source-data tables.")


if __name__ == "__main__":
    main()
