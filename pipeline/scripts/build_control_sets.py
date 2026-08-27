#!/usr/bin/env python3
"""build_control_sets.py — 构建固定正负对照集（用于 HMM 阈值 / min-cov 校准）

正对照 = 已实验表征的 PHB/PHA 解聚酶（e-PhaZ / i-PhaZ / 寡聚体水解酶），
         排除尼龙水解酶同源物。
负对照 = 18 条明确非解聚酶同源物：
         13 条真核 BDH1/BDH2（酮体代谢，非 PHB 降解）+ 5 条尼龙水解酶 nylB/nylC。

依据：pipeline/seeds/seeds_annotation.md（同源物清单）+ seeds_manifest.tsv（query_group）。

输出（pipeline/seeds/controls/）：
  positive.faa   正对照序列
  negative.faa   负对照序列
  controls.tsv   accession / label / is_positive / family / reviewed / evidence

用法: python pipeline/scripts/build_control_sets.py
"""
import csv
import os

# 18 条非解聚酶同源物（来自 seeds_annotation.md）
NEGATIVE_ACC = {
    # 13 真核 BDH1/BDH2（酮体代谢）
    "P29147", "Q02337", "Q02338", "Q80XN0", "P86198", "Q5ZJZ5", "D4A1J4",
    "Q561X9", "Q8JZV9", "Q9BUT1", "C1C4R8", "Q3KPT7", "Q3T046",
    # 5 尼龙水解酶 nylB/nylC
    "Q79F77", "Q1EPR4", "Q1EPR5", "P07061", "P07062",
}

# 解聚酶家族（query_group 前缀）
DEPOL_PREFIXES = ("e-PhaZ", "i-PhaZ", "oligomer")

MANIFEST = "pipeline/seeds/seeds_manifest.tsv"
FASTA = "pipeline/seeds/seeds_curated.faa"
OUTDIR = "pipeline/seeds/controls"


def read_fasta(path):
    seqs = {}
    hdr = None
    buf = []
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if hdr is not None:
                seqs[hdr.split("|")[0]] = (hdr, "".join(buf))
            hdr = line[1:]
            buf = []
        else:
            buf.append(line.strip())
    if hdr is not None:
        seqs[hdr.split("|")[0]] = (hdr, "".join(buf))
    return seqs


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    fasta = read_fasta(FASTA)

    rows = list(csv.DictReader(open(MANIFEST, encoding="utf-8")))

    pos, neg = [], []
    for r in rows:
        acc = r["accession"]
        if acc not in fasta:
            print(f"[skip] {acc}: 序列缺失")
            continue
        qg = r.get("query_group", "")
        is_neg = acc in NEGATIVE_ACC
        is_depol = qg.startswith(DEPOL_PREFIXES)
        label = ("negative" if is_neg else ("positive" if is_depol else "excluded"))
        fam = r.get("protein_name", "")
        rec = (acc, qg, r.get("reviewed", ""), fam)
        if label == "positive":
            pos.append(rec)
        elif label == "negative":
            neg.append(rec)
        # excluded：BdhA 细菌 bdhA、phasin 等背景家族，不纳入对照

    with open(os.path.join(OUTDIR, "positive.faa"), "w", encoding="utf-8") as f:
        for acc, qg, rev, fam in pos:
            hdr, seq = fasta[acc]
            f.write(f">{hdr}\n{seq}\n")
    with open(os.path.join(OUTDIR, "negative.faa"), "w", encoding="utf-8") as f:
        for acc, qg, rev, fam in neg:
            hdr, seq = fasta[acc]
            f.write(f">{hdr}\n{seq}\n")

    with open(os.path.join(OUTDIR, "controls.tsv"), "w", encoding="utf-8") as f:
        f.write("accession\tlabel\tquery_group\treviewed\tprotein_name\n")
        for acc, qg, rev, fam in pos:
            f.write(f"{acc}\tpositive\t{qg}\t{rev}\t{fam}\n")
        for acc, qg, rev, fam in neg:
            f.write(f"{acc}\tnegative\t{qg}\t{rev}\t{fam}\n")

    print(f"正对照: {len(pos)} 条, 负对照: {len(neg)} 条")
    print(f"  positive.faa / negative.faa / controls.tsv -> {OUTDIR}")


if __name__ == "__main__":
    main()
