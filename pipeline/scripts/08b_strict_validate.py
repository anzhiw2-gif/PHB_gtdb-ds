#!/usr/bin/env python3
"""08b_strict_validate.py — 严格验证（文献四基序 + 疏水 x1）
对已命中的序列应用 PHA 解聚酶的四保守基序（LtPHBase/PDB 8DAJ, [T22]）：
  Ser 基序  IDXXXXYVXGLSXGG
  Asp 基序  GXXDYTV
  His 基序  GMXHXXPXXG
  oxyanion  HGCXQ
+ lipase box 疏水 x1（PhaDED [K09] 判别）
输出三级分类：
  tier1 = 四基序全中（高置信）
  tier2 = lipase box 疏水x1 + Ser-Asp-His 三联体（中置信）
  tier3 = 其余（低置信）
用法: python 08b_strict_validate.py --indir data/screen/family_seqs
"""
import argparse
import os
import re

# 四基序（LtPHBase 编号）
SER_MOTIF = re.compile(r"I[^A-Z]{0,1}D[^A-Z]{0,6}Y[^A-Z]{0,2}V[^A-Z]{0,2}G[^A-Z]{0,1}L[^A-Z]?S[^A-Z]?G{1,2}", re.I)
ASP_MOTIF = re.compile(r"G[^A-Z]{2}D[^A-Z]?Y[^A-Z]?T[^A-Z]?V", re.I)
HIS_MOTIF = re.compile(r"G[^A-Z]?M[^A-Z]?H[^A-Z]{2}P[^A-Z]{2}G", re.I)
OXY_MOTIF = re.compile(r"H[^A-Z]?G[^A-Z]?C[^A-Z]?Q", re.I)
LIPASE_BOX = re.compile(r"G([A-Z])S([A-Z])G", re.I)
HYDROPHOBIC = set("LIVMFWAY")


def strict_features(seq):
    f = {
        "ser_motif": bool(SER_MOTIF.search(seq)),
        "asp_motif": bool(ASP_MOTIF.search(seq)),
        "his_motif": bool(HIS_MOTIF.search(seq)),
        "oxy_motif": bool(OXY_MOTIF.search(seq)),
        "lipase_box": False, "x1_hydrophobic": False,
        "ser_pos": -1, "asp_pos": -1, "his_pos": -1,
    }
    m = LIPASE_BOX.search(seq)
    if m:
        f["lipase_box"] = True
        f["ser_pos"] = m.start() + 2
        f["x1_hydrophobic"] = m.group(1).upper() in HYDROPHOBIC
    # Ser-Asp-His 三联体
    if f["ser_pos"] >= 0:
        d = seq[f["ser_pos"] + 60:f["ser_pos"] + 240].find("D")
        if d >= 0:
            asp = f["ser_pos"] + 60 + d
            h = seq[asp + 30:asp + 150].find("H")
            if h >= 0:
                f["asp_pos"] = asp
                f["his_pos"] = asp + 30 + h
    return f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default="data/screen/family_seqs")
    args = ap.parse_args()

    for fam in ["ePhaZ", "iPhaZ", "OH", "ArchPhaZ_hydrolase"]:
        faa = os.path.join(args.indir, f"{fam}_validated.faa")
        if not os.path.exists(faa):
            continue
        seqs = {}
        hdr = None
        buf = []
        for line in open(faa):
            if line.startswith(">"):
                if hdr:
                    seqs[hdr] = "".join(buf)
                hdr = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
        if hdr:
            seqs[hdr] = "".join(buf)

        tiers = {"tier1": [], "tier2": [], "tier3": []}
        for hdr, seq in seqs.items():
            f = strict_features(seq)
            n_motifs = sum([f["ser_motif"], f["asp_motif"], f["his_motif"], f["oxy_motif"]])
            if n_motifs >= 4:
                tiers["tier1"].append(hdr)
            elif f["lipase_box"] and f["x1_hydrophobic"] and f["asp_pos"] >= 0 and f["his_pos"] >= 0:
                tiers["tier2"].append(hdr)
            else:
                tiers["tier3"].append(hdr)

        for tier, hdrs in tiers.items():
            out = os.path.join(args.indir, f"{fam}_validated_{tier}.faa")
            with open(out, "w") as fo:
                for h in hdrs:
                    fo.write(">" + h + "\n" + seqs[h] + "\n")
        print(f"{fam}: total={len(seqs)} tier1={len(tiers['tier1'])} "
              f"tier2={len(tiers['tier2'])} tier3={len(tiers['tier3'])}")


if __name__ == "__main__":
    main()
