#!/usr/bin/env python3
"""rebuild_tier1_faa.py — 从 *_tier1.ids + *_validated.faa 重建 *_tier1.faa（幂等）

背景：tier1.faa 是派生文件（由 tier1.ids 从 validated.faa 提取序列）。
本脚本在 tier1.faa 被误写/损坏时从底层源重建，并验证序列数 == ids 行数。

用法（服务器 T141）:
  ~/miniconda3/envs/phb_gtdb/bin/python pipeline/scripts/rebuild_tier1_faa.py
"""
import os

TIER = "data/screen/tiers"
SEQ = "data/screen/family_seqs"
FAMS = ["ePhaZ", "iPhaZ", "OH", "ArchPhaZ_hydrolase", "ArchPhaZ_patatin"]


def extract(ids_path, validated_path, out_path):
    ids = set()
    with open(ids_path, encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if s:
                ids.add(s)
    n_out = 0
    hdr = None
    buf = []
    with open(validated_path, encoding="utf-8") as fin, \
         open(out_path, "w", encoding="utf-8") as fo:
        for line in fin:
            if line.startswith(">"):
                if hdr is not None and hdr in ids:
                    fo.write(">" + hdr + "\n" + "".join(buf) + "\n")
                    n_out += 1
                hdr = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
        if hdr is not None and hdr in ids:
            fo.write(">" + hdr + "\n" + "".join(buf) + "\n")
            n_out += 1
    return len(ids), n_out


def main():
    all_ok = True
    for fam in FAMS:
        ids_path = os.path.join(TIER, f"{fam}_tier1.ids")
        validated_path = os.path.join(SEQ, f"{fam}_validated.faa")
        out_path = os.path.join(TIER, f"{fam}_tier1.faa")
        if not os.path.exists(ids_path) or not os.path.exists(validated_path):
            print(f"[SKIP] {fam}: 缺少 ids 或 validated.faa")
            continue
        n_ids, n_out = extract(ids_path, validated_path, out_path)
        ok = (n_ids == n_out)
        all_ok = all_ok and ok
        print(f"{fam}: ids={n_ids} 重建序列={n_out} {'OK' if ok else 'MISMATCH!'}")
    print("RESULT:", "ALL_OK" if all_ok else "HAS_MISMATCH")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
