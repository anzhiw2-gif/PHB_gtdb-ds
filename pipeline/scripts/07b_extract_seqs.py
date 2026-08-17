#!/usr/bin/env python3
"""07b_extract_seqs.py — 高效提取命中蛋白序列
单遍流式扫描 shards，提取 unique_proteins.txt 中列出的蛋白序列，
按家族输出 FASTA。
用法: python 07b_extract_seqs.py --ids data/screen/unique_proteins.txt \
        --hits data/screen/hits_filtered.tsv --shards data/proteins/shards_filt \
        --outdir data/screen/family_seqs
"""
import argparse
import os
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", default="data/screen/unique_proteins.txt")
    ap.add_argument("--hits", default="data/screen/hits_filtered.tsv")
    ap.add_argument("--shards", default="data/proteins/shards_filt")
    ap.add_argument("--outdir", default="data/screen/family_seqs")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # 蛋白 ID → 家族（取过滤后命中的家族归属）
    prot2fam = {}
    with open(args.hits) as f:
        header = f.readline().rstrip("\n").split("\t")
        fi = header.index("family")
        pi = header.index("protein")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) > max(fi, pi):
                prot2fam[p[pi]] = p[fi]

    # 唯一蛋白集合（只保留有家族归属的）
    want = set()
    with open(args.ids) as f:
        for line in f:
            pid = line.strip()
            if pid and pid in prot2fam:
                want.add(pid)
    print(f"待提取蛋白: {len(want)}")

    # 流式扫描 shards
    found = defaultdict(list)  # fam -> list of (header, seq)
    nfound = 0
    import glob
    shards = sorted(glob.glob(os.path.join(args.shards, "shard_*.faa")))
    for si, sp in enumerate(shards):
        cur_name = None
        cur_seq = []
        cur_want = False
        with open(sp) as f:
            for line in f:
                if line.startswith(">"):
                    if cur_want and cur_name:
                        found[prot2fam[cur_name]].append((cur_name, "".join(cur_seq)))
                        nfound += 1
                    name = line[1:].split()[0]
                    cur_name = name
                    cur_want = name in want
                    cur_seq = []
                elif cur_want:
                    cur_seq.append(line.strip())
            if cur_want and cur_name:
                found[prot2fam[cur_name]].append((cur_name, "".join(cur_seq)))
                nfound += 1
        if (si + 1) % 10 == 0:
            print(f"  已扫描 {si+1}/{len(shards)} shards, 已提取 {nfound}")

    for fam, items in found.items():
        out = os.path.join(args.outdir, f"{fam}.faa")
        with open(out, "w") as f:
            for name, seq in items:
                f.write(f">{name}\n{seq}\n")
        print(f"  {fam}: {len(items)} seqs -> {out}")
    print(f"[DONE] 提取 {nfound} 条序列")


if __name__ == "__main__":
    main()
