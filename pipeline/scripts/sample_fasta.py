#!/usr/bin/env python3
"""sample_fasta.py — 从 FASTA 随机抽样（固定种子，可复现），并写出抽样名单。

用于系统发育建树的大族抽样（>2000 条），保证抽样可复现、可审计：
  - 输出抽样 FASTA（随机种子固定，顺序打乱后取前 N 条）
  - 输出抽样名单 .list（每行一个 header，即抽样序列的 ID）

用法: python sample_fasta.py INPUT OUTPUT N [--seed 42] [--list OUTPUT.list]
"""
import argparse
import random
import sys


def read_fasta(path):
    seqs = []
    hdr = None
    buf = []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if hdr is not None:
                seqs.append((hdr, "".join(buf)))
            hdr = line[1:].strip()
            buf = []
        else:
            buf.append(line.strip())
    if hdr is not None:
        seqs.append((hdr, "".join(buf)))
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("output")
    ap.add_argument("n", type=int)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--list", default=None)
    args = ap.parse_args()

    seqs = read_fasta(args.input)
    if len(seqs) <= args.n:
        # 无需抽样：直接复制
        with open(args.output, "w") as f:
            for h, s in seqs:
                f.write(f">{h}\n{s}\n")
        sampled = [h for h, _ in seqs]
        print(f"序列数 {len(seqs)} <= N={args.n}，无需抽样（全量复制）")
    else:
        rng = random.Random(args.seed)
        idx = list(range(len(seqs)))
        rng.shuffle(idx)
        sampled = [seqs[i][0] for i in idx[: args.n]]
        with open(args.output, "w") as f:
            for i in idx[: args.n]:
                f.write(f">{seqs[i][0]}\n{seqs[i][1]}\n")
        print(f"从 {len(seqs)} 条中抽样 {args.n} 条（seed={args.seed}）")

    if args.list:
        with open(args.list, "w") as f:
            for h in sampled:
                f.write(h + "\n")
        print(f"抽样名单 -> {args.list}（{len(sampled)} 条）")


if __name__ == "__main__":
    main()
