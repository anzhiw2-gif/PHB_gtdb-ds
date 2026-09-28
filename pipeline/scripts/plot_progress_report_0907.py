#!/usr/bin/env python3
"""Python-only Nature-style progress figures for the ePhaZ/iPhaZ main line."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd

MAIN_FAMILIES = ("ePhaZ", "iPhaZ")
WIDTH_IN = 183 / 25.4
COLORS = {"ePhaZ": "#0F4D92", "iPhaZ": "#42949E", "secreted": "#4C9F70", "review": "#C99000", "gray": "#777777", "dark": "#262626", "light": "#D9D9D9", "red": "#B64342"}
SUBTYPE_LABELS = {
    "ePhaZ_SCL_secreted_evidence": "ePhaZ curated secreted",
    "ePhaZ_core_nonsecreted_evidence": "ePhaZ curated nonsecreted",
    "ePhaZ_tier2_review_evidence": "ePhaZ tier2 review",
    "ePhaZ_competition_review_evidence": "ePhaZ competition review",
    "iPhaZ_intracellular_evidence": "iPhaZ tier1 intracellular",
    "iPhaZ_tier2_review_evidence": "iPhaZ tier2 review",
}


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def validate_main_family_scope(frame: pd.DataFrame) -> list[str]:
    return sorted(set(frame.get("family", [])) - set(MAIN_FAMILIES))


def signalp_scope_counts(frame: pd.DataFrame) -> dict[str, int]:
    scoped = frame[(frame.family == "ePhaZ") & (frame.source_layer == "ePhaZ_curated_core")]
    return {"input_count": int(len(scoped)), "excluded_ephaz_review": int(((frame.family == "ePhaZ") & (frame.source_layer != "ePhaZ_curated_core")).sum()), "excluded_iphaz": int((frame.family == "iPhaZ").sum())}


def initial_vs_tier1_counts(frame: pd.DataFrame) -> pd.DataFrame:
    """Summarize all initial candidates and the stricter tier1 subset by family."""
    scoped = frame[frame.family.isin(MAIN_FAMILIES)].copy()
    initial = scoped.groupby("family").size().reindex(MAIN_FAMILIES, fill_value=0)
    tier1_layers = {
        "ePhaZ": {"ePhaZ_curated_core"},
        "iPhaZ": {"tier1", "iPhaZ_tier1"},
    }
    tier1 = pd.Series(
        {
            family: int((scoped.family.eq(family) & scoped.source_layer.isin(layers)).sum())
            for family, layers in tier1_layers.items()
        }
    ).reindex(MAIN_FAMILIES, fill_value=0)
    return pd.DataFrame(
        {
            "family": list(MAIN_FAMILIES),
            "initial_candidates": initial.astype(int).to_numpy(),
            "tier1_candidates": tier1.astype(int).to_numpy(),
        }
    )


def genus_distribution_table(frame: pd.DataFrame) -> pd.DataFrame:
    """Count unique tier1 genomes by genus while retaining the family split."""
    scoped = frame[frame.family.isin(MAIN_FAMILIES)].copy()
    scoped = scoped.drop_duplicates(["genome", "family"])
    scoped["genus"] = scoped.gtdb_taxonomy.str.extract(r";g__([^;]+)")[0].fillna("unclassified")
    return (
        scoped.groupby(["genus", "family"], as_index=False)
        .genome.nunique()
        .rename(columns={"genome": "genomes"})
        .sort_values(["genus", "family"], kind="stable")
        .reset_index(drop=True)
    )


def report_figure_stems() -> list[str]:
    return ["figure_1_workflow_and_funnel", "figure_2_family_hierarchy", "figure_3_signalp_scope", "figure_4_phylum_distribution", "figure_5_genome_representation_and_boundary"]


def workflow_table() -> pd.DataFrame:
    return pd.DataFrame([
        {"step": "GTDB Release 11 R232", "tool": "GTDB taxonomy + metadata", "count": 199923, "scope": "main_line"},
        {"step": "Protein prediction and sequence QC", "tool": "Pyrodigal + FAA/GFF binding", "count": 199923, "scope": "main_line"},
        {"step": "HMM raw hits", "tool": "HMMER 3.4", "count": 6769772, "scope": "main_line"},
        {"step": "Profile search and tiering", "tool": "HMMER 3.4", "count": 6740900, "scope": "main_line"},
        {"step": "ePhaZ/iPhaZ evidence layering", "tool": "family-specific HMM + SignalP6 scope", "count": 109087, "scope": "main_line"},
        {"step": "Genome-level family representation", "tool": "taxonomy-bound de-duplication", "count": 76141, "scope": "main_line"},
    ])


def resolve_input_tables(run: Path) -> dict[str, Path]:
    paths = {k: run / "inputs" / n for k, n in {"overall": "overall_pha_phb.tsv", "taxonomy": "tier1_genome_family.tsv", "competition": "no_signal_competition.tsv"}.items()}
    missing = [str(p) for p in paths.values() if not p.is_file() or p.stat().st_size == 0]
    if missing:
        raise FileNotFoundError("missing inputs: " + ", ".join(missing))
    return paths


def style(ax, grid="x"):
    ax.grid(axis=grid, color="#E8E8E8", lw=.6, zorder=0)
    ax.tick_params(length=2.5, width=.7, pad=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def panel_label(ax, text):
    ax.text(.01, .99, text, transform=ax.transAxes, va="top", fontweight="bold", fontsize=8)


def figure4_axis_label() -> str:
    """Name Figure 4's unit without collapsing co-occurring families."""
    return "Tier1 genome-family assignments"


def format_figure4_genus_label(genus: str) -> str:
    """Make the GTDB placeholder explicit while preserving source identity."""
    return "Palsa-465 (GTDB placeholder)" if genus == "Palsa-465" else genus


def figure4_caption() -> str:
    return (
        "Counts are taxonomy-bound Tier1 genome-family assignments; families can "
        "co-occur in one genome. Taxonomy, HMM detection and annotation inclusion "
        "bias shape this candidate-only distribution and it is not an experimentally "
        "validated degrader list."
    )


def save_pub(fig, out: Path, stem: str):
    out.mkdir(parents=True, exist_ok=True)
    for ext, dpi in (("svg", None), ("pdf", None), ("tiff", 600), ("png", 300)):
        kwargs = {"facecolor": "white", "bbox_inches": "tight"}
        if dpi:
            kwargs["dpi"] = dpi
        fig.savefig(out / f"{stem}.{ext}", **kwargs)
    plt.close(fig)


def write_tsv(frame: pd.DataFrame, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def file_record(path: Path) -> dict[str, str | int]:
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def build_data(run: Path):
    paths = resolve_input_tables(run)
    overall, tax, comp = (read_tsv(paths[k]) for k in ("overall", "taxonomy", "competition"))
    contamination = validate_main_family_scope(overall)
    if contamination:
        raise ValueError(f"non-main families in overall table: {contamination}")
    main_tax = tax[tax.family.isin(MAIN_FAMILIES)].drop_duplicates(["genome", "family"])
    family = overall.groupby("family", as_index=False).agg(proteins=("accession", "nunique"), genomes=("genome", "nunique"))
    subtype = overall.groupby(["family", "pha_subtype_evidence"], as_index=False).size().rename(columns={"size": "proteins"})
    subtype["label"] = subtype.pha_subtype_evidence.map(SUBTYPE_LABELS)
    initial_tier1 = initial_vs_tier1_counts(overall)
    sig = overall[(overall.family == "ePhaZ") & (overall.source_layer == "ePhaZ_curated_core")].groupby("signal_type", as_index=False).size().rename(columns={"size": "proteins"})
    phylum = main_tax.groupby(["phylum", "family"], as_index=False).size().rename(columns={"size": "genomes"})
    genus_family = genus_distribution_table(main_tax)
    genus = genus_family.groupby("genus", as_index=False).genomes.sum().sort_values("genomes", ascending=False).head(15)
    gen = main_tax.assign(genus=main_tax.gtdb_taxonomy.str.extract(r";g__([^;]+)")[0].fillna("unclassified"))
    co = gen.groupby("genome").family.nunique().value_counts().rename_axis("families_per_genome").reset_index(name="genomes").sort_values("families_per_genome")
    comp_counts = comp["assignment"].value_counts().rename_axis("assignment").reset_index(name="proteins")
    return {"overall": overall, "family": family, "subtype": subtype, "initial_tier1": initial_tier1, "signal": sig, "phylum": phylum, "genus": genus, "genus_family": genus_family, "co": co, "competition": comp_counts, "workflow": workflow_table()}


def figure1(d, out):
    fig = plt.figure(figsize=(WIDTH_IN, 4.2)); gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.2], wspace=.28); ax = fig.add_subplot(gs[0]); funnel = fig.add_subplot(gs[1]); ax.axis("off")
    steps = [("1", "GTDB Release 11 R232", "199,923 representative genomes", "GTDB taxonomy + metadata"), ("2", "Protein prediction and QC", "FAA/GFF sequence binding", "Pyrodigal"), ("3", "HMM raw hits", "6,769,772 hit rows", "HMMER 3.4"), ("4", "HMMER profile search", "family-specific thresholds", "HMMER 3.4"), ("5", "ePhaZ / iPhaZ evidence layering", "109,087 deduplicated proteins", "SignalP6 only on ePhaZ curated-core")]
    y = np.linspace(.86, .18, len(steps))
    for i, (n, title, note, tool) in enumerate(steps):
        box = FancyBboxPatch((.08, y[i]-.065), .84, .13, boxstyle="round,pad=.012,rounding_size=.018", fc="#EAF2F8" if i < 4 else "#D8ECEA", ec="#666", lw=.8); ax.add_patch(box)
        ax.text(.14, y[i]+.03, n, ha="center", va="center", fontsize=8, fontweight="bold"); ax.text(.52, y[i]+.028, title, ha="center", va="center", fontsize=7, fontweight="bold"); ax.text(.52, y[i]-.003, note, ha="center", va="center", fontsize=5.5); ax.text(.52, y[i]-.033, tool, ha="center", va="center", fontsize=5, color="#555")
        if i < len(steps)-1: ax.add_patch(FancyArrowPatch((.5, y[i]-.08), (.5, y[i+1]+.08), arrowstyle="-|>", mutation_scale=9, lw=.8, color="#777"))
    ax.set_title("GTDB-to-candidate workflow", fontsize=8, fontweight="bold", pad=10); panel_label(ax, "a")
    vals = [199923, 6769772, 6740900, 70682, 38405, 76141]; names = ["GTDB R232", "HMM raw hits", "HMMER accepted", "ePhaZ proteins", "iPhaZ proteins", "genome×family"]; y2 = np.arange(len(vals)); funnel.barh(y2, vals, color=["#C4D3E2", "#8FAECC", "#6F95BC", COLORS["ePhaZ"], COLORS["iPhaZ"], "#6E6E6E"]); funnel.set_yticks(y2, names, fontsize=6); funnel.invert_yaxis(); funnel.set_xscale("log"); funnel.set_xlabel("Count (log10 scale)", fontsize=7); funnel.set_title("Main-line data funnel", fontsize=8, fontweight="bold"); style(funnel); panel_label(funnel, "b")
    for yi, v in zip(y2, vals): funnel.text(v * 1.08, yi, f"{v:,}", va="center", fontsize=6)
    fig.text(.5, .01, "Only ePhaZ/iPhaZ are included; computational evidence remains candidate-only.", ha="center", fontsize=6, color="#666"); fig.subplots_adjust(top=.84, bottom=.14, left=.06, right=.98); save_pub(fig, out, report_figure_stems()[0])


def figure2(d, out):
    counts = d["initial_tier1"]
    fig, ax = plt.subplots(figsize=(WIDTH_IN, 3.8)); x = np.arange(len(counts)); width = .34
    bars_initial = ax.bar(x - width / 2, counts.initial_candidates, width, color="#AFC2D2", label="Initial candidates", zorder=2)
    bars_tier1 = ax.bar(x + width / 2, counts.tier1_candidates, width, color=[COLORS["ePhaZ"], COLORS["iPhaZ"]], label="Tier1", zorder=2)
    ax.set_xticks(x, counts.family, fontsize=7); ax.set_ylabel("Candidate proteins"); ax.set_title("Initial candidates and tier1 subset", fontsize=9, fontweight="bold"); style(ax); panel_label(ax, "a"); ax.legend(fontsize=6, loc="upper right")
    for bars in (bars_initial, bars_tier1):
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.015, f"{int(bar.get_height()):,}", ha="center", va="bottom", fontsize=6)
    fig.text(.5, .035, "Initial candidates include tier1 and tier2 rows; tier1 is the stricter family-specific evidence layer.", ha="center", fontsize=6, color="#666"); fig.subplots_adjust(left=.12, right=.96, top=.85, bottom=.25); save_pub(fig, out, report_figure_stems()[1])


def figure3(d, out):
    sig = d["signal"].set_index("signal_type").reindex(["SP", "LIPO", "TAT", "TATLIPO", "OTHER"]).fillna(0); fig, ax = plt.subplots(figsize=(WIDTH_IN, 3.3)); vals = sig.proteins.astype(int); y = np.arange(len(vals)); ax.barh(y, vals, color=[COLORS["secreted"], "#6FA6C8", "#8C6BB1", "#B07AA1", COLORS["light"]], zorder=2); ax.set_yticks(y, sig.index); ax.invert_yaxis(); ax.set_xlabel("ePhaZ curated-core tier1 proteins"); ax.set_title("SignalP6 localization scope (input n = 38,692)", fontsize=9, fontweight="bold"); style(ax); panel_label(ax, "a")
    for yi, v in zip(y, vals): ax.text(v + 250, yi, f"{v:,}", va="center", fontsize=6)
    fig.text(.5, .035, "SignalP6 is an ePhaZ curated-core layer only; localization does not establish PHB substrate specificity.", ha="center", fontsize=6, color="#666"); fig.subplots_adjust(left=.18, right=.96, top=.85, bottom=.25); save_pub(fig, out, report_figure_stems()[2])


def figure4(d, out):
    ph = d["phylum"].pivot(index="phylum", columns="family", values="genomes").fillna(0); ph["total"] = ph.sum(axis=1); ph = ph.sort_values("total").tail(18)
    gen = d["genus_family"].copy(); gen["total"] = gen.groupby("genus").genomes.transform("sum"); top_genera = gen.groupby("genus").total.max().nlargest(15).sort_values().index; gen = gen[gen.genus.isin(top_genera)]; gp = gen.pivot(index="genus", columns="family", values="genomes").fillna(0).reindex(top_genera)
    fig, axs = plt.subplots(1, 2, figsize=(WIDTH_IN, 4.8), gridspec_kw={"width_ratios": [1.18, 1.0], "wspace": .42})
    y = np.arange(len(ph)); axs[0].barh(y, ph.get("ePhaZ", 0), color=COLORS["ePhaZ"], label="ePhaZ"); axs[0].barh(y, ph.get("iPhaZ", 0), left=ph.get("ePhaZ", 0), color=COLORS["iPhaZ"], label="iPhaZ"); axs[0].set_yticks(y, ph.index, fontsize=5.2); axs[0].set_xlabel(figure4_axis_label(), fontsize=7); axs[0].set_title("Phylum", fontsize=8, fontweight="bold"); axs[0].legend(fontsize=6, loc="lower right"); style(axs[0]); panel_label(axs[0], "a")
    y2 = np.arange(len(gp)); axs[1].barh(y2, gp.get("ePhaZ", 0), color=COLORS["ePhaZ"], label="ePhaZ"); axs[1].barh(y2, gp.get("iPhaZ", 0), left=gp.get("ePhaZ", 0), color=COLORS["iPhaZ"], label="iPhaZ"); axs[1].set_yticks(y2, [format_figure4_genus_label(str(name)) for name in gp.index], fontsize=5.2); axs[1].set_xlabel(figure4_axis_label(), fontsize=7); axs[1].set_title("Leading genera", fontsize=8, fontweight="bold"); style(axs[1]); panel_label(axs[1], "b")
    fig.suptitle("GTDB taxonomic distribution of the two main families", fontsize=9, fontweight="bold"); fig.text(.5, .01, "Counts are taxonomy-bound genome×family assignments; families can co-occur in one genome.", ha="center", fontsize=6, color="#666"); fig.subplots_adjust(left=.16, right=.98, top=.87, bottom=.15); save_pub(fig, out, report_figure_stems()[3])


def figure4_repaired(d, out):
    """Render the corrected Figure 4 without changing the frozen source data."""
    ph = d["phylum"].pivot(index="phylum", columns="family", values="genomes").fillna(0)
    ph["total"] = ph.sum(axis=1)
    ph = ph.sort_values("total").tail(18)
    gen = d["genus_family"].copy()
    gen["total"] = gen.groupby("genus").genomes.transform("sum")
    top_genera = gen.groupby("genus").total.max().nlargest(15).sort_values().index
    gen = gen[gen.genus.isin(top_genera)]
    gp = gen.pivot(index="genus", columns="family", values="genomes").fillna(0).reindex(top_genera)
    fig, axs = plt.subplots(1, 2, figsize=(WIDTH_IN, 4.8), gridspec_kw={"width_ratios": [1.18, 1.0], "wspace": .42})
    y = np.arange(len(ph))
    axs[0].barh(y, ph.get("ePhaZ", 0), color=COLORS["ePhaZ"], label="ePhaZ")
    axs[0].barh(y, ph.get("iPhaZ", 0), left=ph.get("ePhaZ", 0), color=COLORS["iPhaZ"], label="iPhaZ")
    axs[0].set_yticks(y, ph.index, fontsize=5.2)
    axs[0].set_xlabel(figure4_axis_label(), fontsize=7)
    axs[0].set_title("Phylum", fontsize=8, fontweight="bold")
    axs[0].legend(fontsize=6, loc="lower right")
    style(axs[0])
    panel_label(axs[0], "a")
    y2 = np.arange(len(gp))
    axs[1].barh(y2, gp.get("ePhaZ", 0), color=COLORS["ePhaZ"], label="ePhaZ")
    axs[1].barh(y2, gp.get("iPhaZ", 0), left=gp.get("ePhaZ", 0), color=COLORS["iPhaZ"], label="iPhaZ")
    axs[1].set_yticks(y2, [format_figure4_genus_label(str(name)) for name in gp.index], fontsize=5.2)
    axs[1].set_xlabel(figure4_axis_label(), fontsize=7)
    axs[1].set_title("Leading genera", fontsize=8, fontweight="bold")
    style(axs[1])
    panel_label(axs[1], "b")
    fig.suptitle("GTDB taxonomic distribution of the two main families", fontsize=9, fontweight="bold")
    fig.text(.5, .01, figure4_caption(), ha="center", fontsize=6, color="#666")
    fig.subplots_adjust(left=.16, right=.98, top=.87, bottom=.15)
    save_pub(fig, out, report_figure_stems()[3])


def figure5(d, out):
    fig, axs = plt.subplots(1, 3, figsize=(WIDTH_IN, 3.8), gridspec_kw={"width_ratios": [1.05, .8, 1.0]}); gen = d["genus"].sort_values("genomes"); y = np.arange(len(gen)); axs[0].barh(y, gen.genomes, color=COLORS["ePhaZ"]); axs[0].set_yticks(y, gen.genus, fontsize=5.5); axs[0].set_xlabel("Unique genomes"); axs[0].set_title("Leading genera", fontsize=8, fontweight="bold"); style(axs[0]); panel_label(axs[0], "a")
    co = d["co"]; axs[1].bar(co.families_per_genome, co.genomes, color=COLORS["iPhaZ"], width=.65); axs[1].set_xlabel("Families per genome"); axs[1].set_ylabel("Genomes"); axs[1].set_title("Family co-occurrence", fontsize=8, fontweight="bold"); style(axs[1]); panel_label(axs[1], "b")
    comp = d["competition"].set_index("assignment").reindex(["ePhaZ_like", "iPhaZ_like", "ambiguous"]).fillna(0).reset_index(); axs[2].barh(np.arange(len(comp)), comp.proteins, color=[COLORS["ePhaZ"], COLORS["iPhaZ"], COLORS["review"]]); axs[2].set_yticks(np.arange(len(comp)), ["ePhaZ-like", "iPhaZ-like", "ambiguous"], fontsize=6); axs[2].set_xlabel("Proteins"); axs[2].set_title("Competition audit", fontsize=8, fontweight="bold"); style(axs[2]); panel_label(axs[2], "c")
    for yi, value in enumerate(comp.proteins): axs[2].text(float(value) + max(float(comp.proteins.max()) * .015, 20), yi, f"{int(value):,}", va="center", fontsize=6)
    for ax in axs: ax.tick_params(labelsize=6)
    fig.suptitle("Genome representation and conservative interpretation boundary", fontsize=9, fontweight="bold"); fig.text(.5, .01, "Taxonomy, co-occurrence and competition are candidate-evidence summaries, not validated PHB degradation phenotypes.", ha="center", fontsize=6, color="#666"); fig.subplots_adjust(top=.82, bottom=.18, left=.12, right=.99, wspace=.42); save_pub(fig, out, report_figure_stems()[4])


def write_contract(run: Path, inputs: dict[str, Path], script: Path):
    contract = {"run_id": run.name, "purpose": "0907 progress report for the ePhaZ/iPhaZ main line", "status": "completed_candidate_only", "input_snapshot": {k: file_record(v) for k, v in sorted(inputs.items())}, "renderer": file_record(script), "environment": {"python_executable": sys.executable, "python_version": sys.version.split()[0], "pandas_version": pd.__version__, "matplotlib_version": matplotlib.__version__}, "scope": {"main_families": list(MAIN_FAMILIES), "archaeal_branch": "excluded from this figure set", "signalp": "ePhaZ curated-core tier1 only (38,692 proteins)"}, "phenotype_boundary": "HMM, SignalP, taxonomy, co-occurrence and competition evidence denotes candidate homology or functional potential, not validated PHB degradation phenotype.", "non_modification_statement": "Visualization only; no formal scan, registry, deploy directory, server state or historical run was modified."}
    (run / "input_contract.json").write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--run-dir", type=Path, required=True); args = ap.parse_args(); run = args.run_dir; out = run / "results" / "0907汇报"; src = run / "results" / "source_data"; inputs = resolve_input_tables(run); write_contract(run, inputs, Path(__file__).resolve()); d = build_data(run)
    for name, frame in {"figure_1_workflow.tsv": d["workflow"], "figure_2_hierarchy.tsv": d["initial_tier1"], "figure_3_signalp.tsv": d["signal"], "figure_4_phylum.tsv": d["phylum"], "figure_4_genus.tsv": d["genus_family"], "figure_5_genus.tsv": d["genus"], "figure_5_cooccurrence.tsv": d["co"], "figure_5_competition.tsv": d["competition"]}.items(): write_tsv(frame, src / name)
    for fn in (figure1, figure2, figure3, figure4_repaired, figure5): fn(d, out)
    files = sorted(p for p in out.iterdir() if p.is_file()); manifest = {"run_id": run.name, "status": "completed_candidate_only", "backend": "Python matplotlib", "main_families": list(MAIN_FAMILIES), "figure_stems": report_figure_stems(), "figure_count": 5, "supplementary_figures": [], "output_folder": "results/0907汇报", "exports": [{"file": p.name, "bytes": p.stat().st_size, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in files], "phenotype_boundary": "All displayed evidence remains candidate-only; no validated PHB phenotype is claimed."}; (run / "results" / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
