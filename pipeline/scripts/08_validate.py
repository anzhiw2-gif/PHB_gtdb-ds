#!/usr/bin/env python3
"""
08_validate.py — 命中序列功能特征验证
  - lipase box (G-x-S-x-G) 扫描
  - 催化三联体 Ser-Asp-His / Ser-His-Asp 模式扫描（宽松）
  - ePhaZ 家族信号肽预测（调用 SignalP6，可选）
  - 输出: data/screen/validation.tsv + data/screen/family_seqs/*_validated.faa
用法: python 08_validate.py [--signalp 0|1] [--signalp-dir /path]
"""
import argparse
import os
import re
import subprocess
import sys

LIPASE_BOX = re.compile(r"G([A-Z])S([A-Z])G", re.I)  # 严格 G-x1-S-x2-G 五肽
HYDROPHOBIC = set("LIVMFWAY")  # 疏水 x1（PhaDED 判别关键 [K09]）
NAD_BINDING = re.compile(r"G[^A-Z]{0,3}G[^A-Z]{0,3}[GAS][^A-Z]{0,3}G", re.I)  # 短链脱氢酶 NAD 结合指纹

# 各家族验证规则：
#   ePhaZ/iPhaZ/OH: α/β 水解酶 → 催化三联体 + lipase box + 长度 200-900
#   BdhA: 短链脱氢酶 → NAD 结合基序 + 长度 200-500
#   phasin: 颗粒蛋白 → 长度 60-200（宽松）
#   patatin: 古菌 patatin 样解聚酶（PhaZh1 型）→ lipase box + Ser-Asp 催化
#     二元组（文献：Ser47+Asp195 必需）+ 长度 200-500；细菌命中需 PHA 上下文
FAMILY_RULES = {
    "ePhaZ": {"min_len": 200, "max_len": 900, "kind": "hydrolase"},
    "iPhaZ": {"min_len": 200, "max_len": 900, "kind": "hydrolase"},
    "OH":    {"min_len": 200, "max_len": 900, "kind": "hydrolase"},
    "BdhA":  {"min_len": 180, "max_len": 500, "kind": "dehydrogenase"},
    "phasin": {"min_len": 60, "max_len": 250, "kind": "phasin"},
    "patatin": {"min_len": 200, "max_len": 500, "kind": "patatin"},
    "ArchPhaZ_patatin": {"min_len": 200, "max_len": 500, "kind": "patatin"},
    "ArchPhaZ_hydrolase": {"min_len": 200, "max_len": 900, "kind": "hydrolase"},
    "PhaJ": {"min_len": 150, "max_len": 500, "kind": "hydratase"},
    "PhaC": {"min_len": 300, "max_len": 700, "kind": "synthase"},
}


def scan_features(seq: str, kind: str) -> dict:
    """按家族类型扫描功能特征。
    返回独立标志：ser_triad / cys_triad / patatin_dyad 互不覆盖。
    """
    feats = {"lipase_box": False, "x1_hydrophobic": False, "nad_binding": False,
             "ser_triad": False, "cys_triad": False, "patatin_dyad": False,
             "ser_pos": -1, "asp_pos": -1, "his_pos": -1, "cys_pos": -1}
    if kind == "dehydrogenase":
        feats["nad_binding"] = bool(NAD_BINDING.search(seq))
        return feats
    if kind in ("phasin", "hydratase", "synthase"):
        return feats

    # --- 严格 lipase box G-x1-S-x2-G + 疏水 x1 ---
    m = LIPASE_BOX.search(seq)
    ser_box = -1
    if m:
        feats["lipase_box"] = True
        ser_box = m.start() + 2
        feats["ser_pos"] = ser_box
        feats["x1_hydrophobic"] = m.group(1).upper() in HYDROPHOBIC

    if kind == "patatin":
        # patatin 催化 Ser-Asp 二元组（Ser 后 60-260aa 内 Asp）
        s = ser_box if ser_box >= 0 else seq.find("S", 30)
        if s >= 0:
            d = seq[s + 60:s + 260].find("D")
            if d >= 0:
                feats["asp_pos"] = s + 60 + d
                feats["patatin_dyad"] = True
        return feats

    # --- hydrolase：Ser 三联体（独立检测，不被 Cys 覆盖）---
    s = ser_box if ser_box >= 0 else seq.find("S", 50)
    if s >= 0:
        d = seq[s + 60:s + 240].find("D")
        if d >= 0:
            asp = s + 60 + d
            h = seq[asp + 30:asp + 150].find("H")
            if h >= 0:
                feats["ser_pos"] = s
                feats["asp_pos"] = asp
                feats["his_pos"] = asp + 30 + h
                feats["ser_triad"] = True

    # --- Cys 三联体（胞内 i-nPHASCL 无 lipase box，PhaDED 家族1）---
    c = seq.find("C", 100)
    if c >= 0:
        d = seq[c + 80:c + 300].find("D")
        if d >= 0:
            c_asp = c + 80 + d
            h = seq[c_asp + 20:c_asp + 120].find("H")
            if h >= 0:
                feats["cys_pos"] = c
                feats["cys_triad"] = True
    return feats


def run_signalp(faa: str, outdir: str, sig_dir: str) -> str:
    """SignalP6 批量预测，返回输出文件路径"""
    out = os.path.join(outdir, "signalp6")
    os.makedirs(out, exist_ok=True)
    cmd = ["bash", f"{sig_dir}/bin/signalp6", "-fasta", faa, "-organism", "bacteria",
           "-format", "txt", "-output_dir", out, "-mode", "fast"]
    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=3600)
        return os.path.join(out, "prediction_results.txt")
    except Exception as e:
        print(f"  [SignalP] 失败: {e}", file=sys.stderr)
        return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default="data/screen/family_seqs")
    ap.add_argument("--outdir", default="data/screen")
    ap.add_argument("--signalp", type=int, default=0)
    ap.add_argument("--signalp-dir", default="$HOME/software/signalp6")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    sig_dir = os.path.expandvars(args.signalp_dir)

    out_rows = []
    for fam in ["ePhaZ", "iPhaZ", "OH", "BdhA", "ArchPhaZ_patatin", "ArchPhaZ_hydrolase",
                "PhaJ", "phasin", "PhaC"]:
        faa = os.path.join(args.indir, f"{fam}.faa")
        if not os.path.exists(faa):
            continue
        seqs = []
        cur_hdr = None
        cur_seq = []
        with open(faa) as f:
            for line in f:
                line = line.rstrip()
                if line.startswith(">"):
                    if cur_hdr:
                        seqs.append((cur_hdr, "".join(cur_seq)))
                    cur_hdr = line[1:]
                    cur_seq = []
                else:
                    cur_seq.append(line)
            if cur_hdr:
                seqs.append((cur_hdr, "".join(cur_seq)))

        sig = {}
        if args.signalp and fam == "ePhaZ" and seqs:
            tmp = os.path.join(args.outdir, "family_seqs", f"{fam}_for_signalp.faa")
            with open(tmp, "w") as f:
                for i, (h, s) in enumerate(seqs):
                    f.write(f">s{i}\n{s}\n")
            res = run_signalp(tmp, args.outdir, sig_dir)
            if res and os.path.exists(res):
                for line in open(res):
                    p = line.split()
                    if len(p) >= 2 and p[0].startswith("s"):
                        idx = int(p[0][1:])
                        # SignalP6 输出列：ID, Prediction, ...（SP=分泌型）
                        sig[idx] = "SP" in line

        valid_faa = os.path.join(args.indir, f"{fam}_validated.faa")
        rules = FAMILY_RULES.get(fam, {"min_len": 100, "max_len": 2000, "kind": "hydrolase"})
        with open(valid_faa, "w") as fv:
            n_valid = 0
            for i, (hdr, seq) in enumerate(seqs):
                feats = scan_features(seq, rules["kind"])
                is_sp = sig.get(i, None)
                len_ok = rules["min_len"] <= len(seq) <= rules["max_len"]
                if rules["kind"] == "dehydrogenase":
                    high_conf = len_ok and feats["nad_binding"]
                elif rules["kind"] == "phasin":
                    high_conf = len_ok
                elif rules["kind"] == "patatin":
                    high_conf = len_ok and feats["patatin_dyad"]
                elif rules["kind"] in ("hydratase", "synthase"):
                    high_conf = len_ok
                elif fam == "iPhaZ":
                    # 胞内：Ser 型（lipase box + 疏水x1 + Ser三联体）或 Cys 型（无 lipase box）
                    high_conf = len_ok and (
                        (feats["lipase_box"] and feats["x1_hydrophobic"] and feats["ser_triad"])
                        or feats["cys_triad"])
                else:
                    # 胞外/经典：lipase box + 疏水x1 + Ser-Asp-His 三联体
                    high_conf = len_ok and feats["lipase_box"] and feats["x1_hydrophobic"] and feats["ser_triad"]
                row = {
                    "family": fam, "protein": hdr.split("|")[0], "genome": hdr.split("|")[0],
                    "length": len(seq),
                    "lipase_box": feats["lipase_box"],
                    "x1_hydrophobic": feats["x1_hydrophobic"],
                    "nad_binding": feats["nad_binding"],
                    "ser_triad": feats["ser_triad"],
                    "cys_triad": feats["cys_triad"],
                    "patatin_dyad": feats["patatin_dyad"],
                    "ser_pos": feats["ser_pos"], "asp_pos": feats["asp_pos"], "his_pos": feats["his_pos"],
                    "signal_peptide": "" if is_sp is None else ("yes" if is_sp else "no"),
                    "length_ok": len_ok,
                    "high_confidence": high_conf,
                }
                out_rows.append(row)
                if high_conf:
                    n_valid += 1
                    fv.write(f">{hdr}\n{seq}\n")
            print(f"  {fam}: {len(seqs)} hits, {n_valid} 高置信 -> {valid_faa}")

    with open(os.path.join(args.outdir, "validation.tsv"), "w") as f:
        cols = ["family", "protein", "genome", "length", "lipase_box", "x1_hydrophobic", "nad_binding",
                "ser_triad", "cys_triad", "patatin_dyad", "ser_pos", "asp_pos", "his_pos",
                "signal_peptide", "length_ok", "high_confidence"]
        f.write("\t".join(cols) + "\n")
        for r in out_rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"validation.tsv: {len(out_rows)} rows")


if __name__ == "__main__":
    main()
