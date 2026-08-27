#!/usr/bin/env python3
"""calibrate_hmms.py — 用固定正负对照集校准家族 HMM 的 E-value 阈值（含覆盖度）

对每个核心解聚酶家族 HMM，把"该家族正对照 + 全部负对照"跑 hmmsearch，
在阈值网格下计算 TP/FP/TN/FN、precision/recall/F1、MCC，并输出：
  - 每个命中序列的 best E-value 与 HMM 覆盖度（cov，来自 domtblout）
  - 推荐阈值（F1 最大；另给出 FPR=0 的最宽松阈值，供 --min-cov 校准参考）

输入: pipeline/seeds/controls/{positive,negative}.faa + controls.tsv
HMM: data/hmms/{ePhaZ,iPhaZ,OH}.hmm（curated 金标准 HMM）
输出: results/tables/calibration_summary.tsv + results/tables/calibration_hits.tsv

用法（服务器 T141）:
  ~/miniconda3/envs/phb_gtdb/bin/python scripts/calibrate_hmms.py
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys

# hmmsearch 可执行文件（优先 PATH，否则用 conda env 绝对路径）
HMMSEARCH = shutil.which("hmmsearch") or os.path.expanduser(
    "~/miniconda3/envs/phb_gtdb/bin/hmmsearch")

FAM_HMM = {
    "ePhaZ": "data/hmms/ePhaZ.hmm",
    "iPhaZ": "data/hmms/iPhaZ.hmm",
    "OH": "data/hmms/OH.hmm",
}
# query_group 前缀 → 家族
PREFIX_FAM = {"e-PhaZ": "ePhaZ", "i-PhaZ": "iPhaZ", "oligomer": "OH"}

# 服务器路径契约（cwd = 工作区根）：种子在 data/seeds/，本地 git 在 pipeline/seeds/
CONTROLS = "data/seeds/controls"
OUTDIR = "results/tables"
THRESHOLDS = [1e-2, 1e-3, 1e-5, 1e-8, 1e-10, 1e-15, 1e-20]


def fam_of(qg):
    for p, fam in PREFIX_FAM.items():
        if qg.startswith(p):
            return fam
    return None


def hmmsearch(hmm, faa, work):
    tbl = os.path.join(work, "hits.tbl")
    dom = os.path.join(work, "hits.dom")
    subprocess.run(
        [HMMSEARCH, "--tblout", tbl, "--domtblout", dom, "-E", "1e-2",
         "--cpu", "4", hmm, faa],
        check=True, capture_output=True, text=True)
    hits = {}
    for line in open(tbl):
        if line.startswith("#"):
            continue
        c = line.split()
        if len(c) >= 6:
            acc = c[0].split("|")[0]  # target 名取第一段 = accession
            hits[acc] = {"E": float(c[4])}
    # 覆盖度：max domain (hmm_to - hmm_from + 1) / qlen
    for line in open(dom):
        if line.startswith("#"):
            continue
        c = line.split()
        if len(c) >= 17:
            acc = c[0].split("|")[0]
            try:
                qlen = float(c[5])
                cov = (int(c[16]) - int(c[15]) + 1) / qlen
            except (ValueError, ZeroDivisionError):
                continue
            if acc in hits:
                hits[acc]["cov"] = max(hits[acc].get("cov", 0.0), min(1.0, cov))
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--controls", default=CONTROLS)
    ap.add_argument("--outdir", default=OUTDIR)
    args = ap.parse_args()
    controls = args.controls
    os.makedirs(args.outdir, exist_ok=True)

    # 读 controls.tsv
    pos_by_fam = {}
    neg_acc = []
    for r in csv.DictReader(open(os.path.join(controls, "controls.tsv"), encoding="utf-8"),
                             delimiter="\t"):
        if r["label"] == "positive":
            fam = fam_of(r["query_group"])
            if fam:
                pos_by_fam.setdefault(fam, []).append(r["accession"])
        else:
            neg_acc.append(r["accession"])

    # 读对照序列
    def read_faa(path):
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

    pos_seqs = read_faa(os.path.join(controls, "positive.faa"))
    neg_seqs = read_faa(os.path.join(controls, "negative.faa"))

    work = "/tmp/calibrate_work"
    os.makedirs(work, exist_ok=True)

    hit_rows = []
    summary_rows = []
    for fam, hmm in FAM_HMM.items():
        pos_acc = pos_by_fam.get(fam, [])
        if not pos_acc:
            print(f"[skip] {fam}: 无正对照")
            continue
        # 组合该家族正对照 + 全部负对照
        with open(os.path.join(work, f"{fam}_probe.faa"), "w", encoding="utf-8") as f:
            for acc in pos_acc:
                if acc in pos_seqs:
                    f.write(f">{pos_seqs[acc][0]}\n{pos_seqs[acc][1]}\n")
            for acc in neg_acc:
                if acc in neg_seqs:
                    f.write(f">{neg_seqs[acc][0]}\n{neg_seqs[acc][1]}\n")

        hits = hmmsearch(hmm, os.path.join(work, f"{fam}_probe.faa"), work)
        pos_set = set(pos_acc)
        for acc, v in hits.items():
            is_pos = acc in pos_set
            hit_rows.append({
                "family": fam, "accession": acc, "is_positive": is_pos,
                "E_value": f"{v['E']:.2e}", "cov": f"{v.get('cov', 0.0):.3f}",
            })

        # 阈值网格
        for t in THRESHOLDS:
            tp = sum(1 for a in pos_acc if a in hits and hits[a]["E"] <= t)
            fn = sum(1 for a in pos_acc if a not in hits or hits[a]["E"] > t)
            fp = sum(1 for a in neg_acc if a in hits and hits[a]["E"] <= t)
            tn = sum(1 for a in neg_acc if a not in hits or hits[a]["E"] > t)
            prec = tp / (tp + fp) if tp + fp else 0.0
            rec = tp / (tp + fn) if tp + fn else 0.0
            f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
            mcc_num = tp * tn - fp * fn
            mcc_den = ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5
            mcc = mcc_num / mcc_den if mcc_den else 0.0
            summary_rows.append({
                "family": fam, "threshold": f"{t:.0e}", "TP": tp, "FP": fp,
                "FN": fn, "TN": tn, "precision": round(prec, 3),
                "recall": round(rec, 3), "F1": round(f1, 3), "MCC": round(mcc, 3),
            })

    with open(os.path.join(args.outdir, "calibration_summary.tsv"), "w", encoding="utf-8") as f:
        cols = ["family", "threshold", "TP", "FP", "FN", "TN", "precision", "recall", "F1", "MCC"]
        f.write("\t".join(cols) + "\n")
        for r in summary_rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(os.path.join(args.outdir, "calibration_hits.tsv"), "w", encoding="utf-8") as f:
        cols = ["family", "accession", "is_positive", "E_value", "cov"]
        f.write("\t".join(cols) + "\n")
        for r in hit_rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"校准完成 -> {args.outdir}/calibration_summary.tsv + calibration_hits.tsv")
    # 打印每个家族 F1 最优阈值
    for fam in FAM_HMM:
        fam_rows = [r for r in summary_rows if r["family"] == fam]
        if not fam_rows:
            continue
        best = max(fam_rows, key=lambda r: r["F1"])
        fpr0 = [r for r in fam_rows if r["FP"] == 0]
        fpr0_best = max(fpr0, key=lambda r: r["recall"]) if fpr0 else None
        print(f"  {fam}: F1最优 阈值={best['threshold']} F1={best['F1']} "
              f"(TP={best['TP']}/FP={best['FP']}/FN={best['FN']})")
        if fpr0_best:
            print(f"       FPR=0 最宽松 阈值={fpr0_best['threshold']} recall={fpr0_best['recall']}")


if __name__ == "__main__":
    main()
