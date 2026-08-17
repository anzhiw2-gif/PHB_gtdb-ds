#!/usr/bin/env python3
"""06b_aggregate_hits.py — 稳健聚合 hmmsearch tbl 文件为 hits_all.tsv
（替代 bash 循环，容错空文件/异常行）
"""
import glob
import os
import sys

HMMOUT = "data/screen/hmmsearch"
OUT = "data/screen/hits_all.tsv"

tbl_files = sorted(glob.glob(os.path.join(HMMOUT, "*.tbl")))
print(f"tbl 文件数: {len(tbl_files)}")

with open(OUT, "w") as fo:
    fo.write("family\tshard\tprotein\ttacc\tE-value\tscore\tbias\tdomE\tqname\n")
    n = 0
    for f in tbl_files:
        base = os.path.basename(f)[:-4]  # 去掉 .tbl
        fam, sname = base.split("__", 1)
        with open(f, errors="replace") as fin:
            for line in fin:
                if line.startswith("#"):
                    continue
                p = line.rstrip("\n").split()
                if len(p) < 8:
                    continue
                # tblout 列: target(0) tacc(1) qname(2) qacc(3) E-value(4) score(5) bias(6) domE(7)
                fo.write(f"{fam}\t{sname}\t{p[0]}\t{p[1]}\t{p[4]}\t{p[5]}\t{p[6]}\t{p[7]}\t{p[2]}\n")
                n += 1
print(f"总命中: {n} 行 -> {OUT}")
