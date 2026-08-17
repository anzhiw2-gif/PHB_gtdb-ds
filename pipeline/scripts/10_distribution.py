#!/usr/bin/env python3
"""
10_distribution.py — PHB 降解基因的生态/分类学分布统计

输入:
  data/screen/genome_hits.tsv                       基因组×家族命中矩阵
  ~/GTDB/taxonomy/bac120_taxonomy_r232.tsv          分类学（只读）
  ~/GTDB/metadata/bac120_metadata_r232.tsv.gz       元数据（只读，含 isolation source）

输出:
  results/tables/phylum_family_distribution.tsv     门×家族
  results/tables/phylum_genome_counts.tsv           门总体计数
  results/tables/ecology_isolation_source.tsv       命中基因组的 isolation source 分布
  results/tables/ecology_family_isolation.tsv       家族×isolation source 共现
  results/figures/phylum_family_heatmap.png
  results/figures/family_genome_counts.png
  results/figures/ecology_isolation_source.png
"""
import argparse
import gzip
import os
import re
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# 将自由文本 isolation source 归并到粗粒度生态类别（小写关键词匹配）
ECOLOGY_RULES = [
    ("soil", re.compile(r"soil|sediment|rhizosphere|rhizoplane|compost|rhizobia|nodule|forest|peat|permafrost")),
    ("marine", re.compile(r"marine|sea|ocean|seawater|coastal|estuar|tidal|coral|sponge|deep.?sea")),
    ("freshwater", re.compile(r"freshwater|river|lake|pond|stream|aquifer|groundwater|spring|reservoir|wetland|sphagnum")),
    ("gut/host", re.compile(r"gut|feces|faeces|stool|intestin|rumen|termite|insect|host|human|clinical|patient|blood|urine|oral|skin|vagina|nasal|sputum|wound")),
    ("plant-associated", re.compile(r"plant|endophyte|phyllosphere|leaf|root|stem|seed|maize|rice|wheat|arabidopsis")),
    ("extreme", re.compile(r"hot spring|geothermal|hydrothermal|vent|hypersaline|salt|halo|alkali|soda|acid mine|volcan|desert|arid")),
    ("activated-sludge/bioreactor", re.compile(r"sludge|bioreactor|activated|wastewater|sewage|digester|ferment|industrial")),
    ("biofilm/misc", re.compile(r"biofilm|plaque|mat|ice|glacier|snow|air|dust")),
]


def categorize_isolation(src):
    if not src or str(src).lower() in ("nan", "none", "not available", "n/a", ""):
        return "unknown"
    s = str(src).lower()
    for cat, rx in ECOLOGY_RULES:
        if rx.search(s):
            return cat
    return "other"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", default="data/screen/genome_hits.tsv")
    ap.add_argument("--metadata", default=os.path.expanduser("~/GTDB/metadata/bac120_metadata_r232.tsv.gz"))
    ap.add_argument("--taxonomy", default=os.path.expanduser("~/GTDB/taxonomy/bac120_taxonomy_r232.tsv"))
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--top-phyla", type=int, default=25)
    ap.add_argument("--top-sources", type=int, default=20)
    args = ap.parse_args()
    os.makedirs(f"{args.outdir}/tables", exist_ok=True)
    os.makedirs(f"{args.outdir}/figures", exist_ok=True)

    hits = pd.read_csv(args.hits, sep="\t")
    print(f"genome-family pairs: {len(hits)}; genomes: {hits['genome'].nunique()}")

    # ---------- 分类学 ----------
    tax = {}
    if os.path.exists(args.taxonomy):
        with open(args.taxonomy) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    tax[p[0]] = p[1]
    hits["gtdb_taxonomy"] = hits["genome"].map(tax).fillna("unclassified")
    hits["phylum"] = hits["gtdb_taxonomy"].apply(
        lambda t: t.split(";")[1].replace("p__", "") if ";" in t and len(t.split(";")) > 1 else "unknown")
    hits["class"] = hits["gtdb_taxonomy"].apply(
        lambda t: t.split(";")[2].replace("c__", "") if ";" in t and len(t.split(";")) > 2 else "unknown")

    phylum_fam = hits.groupby(["phylum", "family"]).size().reset_index(name="genomes")
    phylum_fam.to_csv(f"{args.outdir}/tables/phylum_family_distribution.tsv", sep="\t", index=False)
    phylum_total = hits.groupby("phylum")["genome"].nunique().sort_values(ascending=False)
    phylum_total.to_csv(f"{args.outdir}/tables/phylum_genome_counts.tsv", sep="\t")

    # ---------- 生态元数据关联 ----------
    iso_map = {}
    if os.path.exists(args.metadata):
        try:
            meta = pd.read_csv(args.metadata, sep="\t", compression="gzip", usecols=lambda c: c in ("accession",) or any(
                k in c.lower() for k in ("isolation", "ecosystem", "habitat", "source", "environment")))
            acc_col = "accession" if "accession" in meta.columns else meta.columns[0]
            src_cols = [c for c in meta.columns if c != acc_col and any(
                k in c.lower() for k in ("isolation", "ecosystem", "habitat", "source", "environment"))]
            print(f"metadata 生态列: {src_cols}")
            if src_cols:
                # 优先 ncbi_isolation_source；否则取第一个可用列
                iso_col = "ncbi_isolation_source" if "ncbi_isolation_source" in meta.columns else src_cols[0]
                iso_map = dict(zip(meta[acc_col], meta[iso_col]))
        except Exception as e:
            print("metadata 读取失败:", e)
    else:
        print("[warn] metadata 不存在，跳过生态关联")

    if iso_map:
        hits["isolation_source"] = hits["genome"].map(iso_map).fillna("")
        hits["ecosystem_cat"] = hits["isolation_source"].apply(categorize_isolation)

        # 命中基因组（去重）的 isolation source 分布
        g = hits.drop_duplicates("genome")
        src_count = g["isolation_source"].replace("", "unknown").value_counts().head(args.top_sources)
        src_count.to_csv(f"{args.outdir}/tables/ecology_isolation_source.tsv", sep="\t")

        eco_count = g["ecosystem_cat"].value_counts()
        eco_count.to_csv(f"{args.outdir}/tables/ecology_ecosystem_cat.tsv", sep="\t")

        # 家族 × 生态类别共现
        fam_eco = hits.groupby(["family", "ecosystem_cat"])["genome"].nunique().reset_index(name="genomes")
        fam_eco.to_csv(f"{args.outdir}/tables/ecology_family_isolation.tsv", sep="\t", index=False)

        # 图：生态类别分布
        fig, ax = plt.subplots(figsize=(9, 5))
        eco_count.plot(kind="barh", ax=ax, color="#457b9d")
        ax.set_xlabel("genomes with hit")
        ax.set_title("PHB degradation gene hits across ecosystem categories")
        plt.tight_layout()
        plt.savefig(f"{args.outdir}/figures/ecology_isolation_source.png", dpi=150)
        plt.close()
        print("\n[生态类别分布]")
        print(eco_count.to_string())
    else:
        print("[warn] 无生态元数据可关联")

    # ---------- 图 1: 门×家族热图 ----------
    pivot = phylum_fam.pivot(index="phylum", columns="family", values="genomes").fillna(0)
    top = pivot.sum(axis=1).sort_values(ascending=False).head(args.top_phyla).index
    pivot = pivot.loc[top]
    fig, ax = plt.subplots(figsize=(10, 12))
    im = ax.imshow(pivot.values, cmap="YlGnBu", aspect="auto")
    ax.set_yticks(range(len(pivot.index)), pivot.index)
    ax.set_xticks(range(len(pivot.columns)), pivot.columns, rotation=45, ha="right")
    ax.set_title("PHB degradation gene families across GTDB phyla")
    plt.colorbar(im, label="genomes")
    plt.tight_layout()
    plt.savefig(f"{args.outdir}/figures/phylum_family_heatmap.png", dpi=150)
    plt.close()

    # ---------- 图 2: 家族检出基因组数 ----------
    fam_count = hits.groupby("family")["genome"].nunique().sort_values()
    fig, ax = plt.subplots(figsize=(8, 4))
    fam_count.plot(kind="barh", ax=ax, color="#2a9d8f")
    ax.set_xlabel("genomes with hit")
    ax.set_title("PHB degradation gene families: genome counts")
    plt.tight_layout()
    plt.savefig(f"{args.outdir}/figures/family_genome_counts.png", dpi=150)
    plt.close()

    print("\n[SUMMARY] 家族检出基因组数:")
    print(fam_count.to_string())


if __name__ == "__main__":
    main()
