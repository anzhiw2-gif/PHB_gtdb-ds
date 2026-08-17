#!/usr/bin/env python3
"""检查种子 FASTA 质量：每条序列长度、家族分布、空序列"""
import sys
from collections import Counter, defaultdict

faa = sys.argv[1] if len(sys.argv) > 1 else "data/seeds/seeds_family.faa"
lengths = []
fam_lens = defaultdict(list)
cur_fam = "unknown"
cur_acc = ""
seqlen = 0
n = 0
empty = []

with open(faa) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if seqlen == 0 and cur_acc:
                empty.append(cur_acc)
            if cur_acc:
                lengths.append(seqlen)
                fam_lens[cur_fam].append(seqlen)
            parts = line[1:].split("|")
            cur_acc = parts[0]
            cur_fam = parts[1] if len(parts) > 1 else "unknown"
            seqlen = 0
            n += 1
        else:
            seqlen += len(line)
if cur_acc:
    lengths.append(seqlen)
    fam_lens[cur_fam].append(seqlen)

print(f"total seqs: {n}")
print(f"empty seqs: {len(empty)} {empty[:10]}")
if lengths:
    print(f"length: min={min(lengths)} max={max(lengths)} median={sorted(lengths)[len(lengths)//2]}")
    short = [l for l in lengths if l < 80]
    print(f"short(<80aa): {len(short)}")
print("\nper family:")
for fam, ls in sorted(fam_lens.items(), key=lambda x: -len(x[1])):
    print(f"  {fam}: {len(ls)} seqs, median_len={sorted(ls)[len(ls)//2] if ls else 0}")
