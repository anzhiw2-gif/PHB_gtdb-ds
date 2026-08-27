#!/usr/bin/env python3
"""extract_signalp_sp.py — 导出"高置信胞外 ePhaZ"子集（含信号肽的 tier1 序列）

按 SignalP6 结果，把 prediction != OTHER（即 Sec/SPI、Lipo/SPII、Tat、TatLipo、Pilin）
的 ePhaZ tier1 序列导出为单独 FASTA，header 追加 SignalP 类型，便于后续细分。

输出:
  data/screen/tiers/ePhaZ_tier1_signalpeptide.faa   含信号肽子集（21,856 条）
  results/tables/ePhaZ_signalp_subset.tsv           accession / type / length 明细
"""
import os
from collections import Counter

TIER = "data/screen/tiers/ePhaZ_tier1.faa"
PRED = "results/signalp/ePhaZ/prediction_results.txt"
OUT_FAA = "data/screen/tiers/ePhaZ_tier1_signalpeptide.faa"
OUT_TSV = "results/tables/ePhaZ_signalp_subset.tsv"


def read_faa(path):
    seqs = {}
    hdr = None
    buf = []
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if hdr is not None:
                seqs[hdr] = "".join(buf)
            hdr = line[1:].strip()
            buf = []
        else:
            buf.append(line.strip())
    if hdr is not None:
        seqs[hdr] = "".join(buf)
    return seqs


def read_pred(path):
    pred = {}
    for line in open(path, encoding="utf-8"):
        if line.startswith("#") or not line.strip():
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) >= 2:
            pred[c[0]] = c[1]
    return pred


def main():
    seqs = read_faa(TIER)
    pred = read_pred(PRED)

    sp_acc = [a for a, t in pred.items() if t != "OTHER"]
    types = Counter(pred[a] for a in sp_acc)

    # 写 FASTA（header 追加 SignalP 类型）
    n_written = 0
    with open(OUT_FAA, "w", encoding="utf-8") as f:
        for a in sp_acc:
            if a in seqs:
                f.write(f">{a}|{pred[a]}\n{seqs[a]}\n")
                n_written += 1

    # 写 TSV 明细
    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    with open(OUT_TSV, "w", encoding="utf-8") as f:
        f.write("accession\ttype\tlength\n")
        for a in sp_acc:
            if a in seqs:
                f.write(f"{a}\t{pred[a]}\t{len(seqs[a])}\n")

    print(f"含信号肽子集: {len(sp_acc)} 条（写入 {n_written} 条）")
    for t, c in types.most_common():
        print(f"  {t}: {c}")
    print(f"  FASTA -> {OUT_FAA}")
    print(f"  TSV   -> {OUT_TSV}")


if __name__ == "__main__":
    main()
