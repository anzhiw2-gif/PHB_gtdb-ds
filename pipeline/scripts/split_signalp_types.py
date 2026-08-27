#!/usr/bin/env python3
"""split_signalp_types.py — 把含信号肽的 ePhaZ 子集按 SignalP 类型拆分为独立 FASTA

输入: data/screen/tiers/ePhaZ_tier1_signalpeptide.faa（header 形如 >acc|TYPE）
输出: data/screen/tiers/ePhaZ_tier1_signalpeptide_{SP,LIPO,TAT}.faa
"""
import os
from collections import defaultdict

IN = "data/screen/tiers/ePhaZ_tier1_signalpeptide.faa"
OUTDIR = "data/screen/tiers"

# 读取并按类型分组
by_type = defaultdict(list)  # type -> list[(header, seq)]
hdr = None
buf = []
for line in open(IN, encoding="utf-8"):
    line = line.rstrip("\n")
    if line.startswith(">"):
        if hdr is not None:
            typ = hdr.rsplit("|", 1)[-1]
            by_type[typ].append((hdr, "".join(buf)))
        hdr = line[1:].strip()
        buf = []
    else:
        buf.append(line.strip())
if hdr is not None:
    typ = hdr.rsplit("|", 1)[-1]
    by_type[typ].append((hdr, "".join(buf)))

# 写各类型 FASTA
for typ, items in sorted(by_type.items(), key=lambda x: -len(x[1])):
    out = os.path.join(OUTDIR, f"ePhaZ_tier1_signalpeptide_{typ}.faa")
    with open(out, "w", encoding="utf-8") as f:
        for h, s in items:
            f.write(f">{h}\n{s}\n")
    print(f"{typ}: {len(items)} 条 -> {out}")

print(f"\n合计 {sum(len(v) for v in by_type.values())} 条")
