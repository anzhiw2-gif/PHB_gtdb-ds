#!/usr/bin/env python3
"""
03_prep_families.py — 将种子按干净家族整理为独立 FASTA
映射关系（manifest family -> 目标家族）：
  ePhaZ_EC3.1.1.75 / ePhaZ_EC3.1.1.76 / ePhaZ_pname / e-PhaZ_* / ePhaZ_ESTHER / PAZy_PHA -> ePhaZ
  iPhaZ / i-PhaZ_* / iPhaZ_PHAZ7 -> iPhaZ
  oligomer_hydrolase* -> OH
  3HB_dehydrogenase_* -> BdhA
  phasin -> phasin
输出：data/seeds/families/{ePhaZ,iPhaZ,OH,BdhA,phasin}.faa
"""
import argparse
import os
import sys

FAMILY_MAP = {
    "ePhaZ_EC3.1.1.75": "ePhaZ", "ePhaZ_EC3.1.1.76": "ePhaZ", "ePhaZ_pname": "ePhaZ",
    "ePhaZ_ESTHER": "ePhaZ", "PAZy_PHA": "ePhaZ", "ePhaZ": "ePhaZ",
    "iPhaZ": "iPhaZ", "iPhaZ_PHAZ7": "iPhaZ",
    "oligomer_hydrolase_EC3.1.1.22": "OH", "oligomer_hydrolase": "OH",
    "3HB_dehydrogenase_EC1.1.1.30": "BdhA",
    "phasin": "phasin",
}
# 带前缀的标签也按前缀映射
PREFIX_MAP = {
    "e-PhaZ": "ePhaZ", "i-PhaZ": "iPhaZ",
    "PAZy": "ePhaZ",
}


def target_family(tag: str) -> str:
    if tag in FAMILY_MAP:
        return FAMILY_MAP[tag]
    for prefix, fam in PREFIX_MAP.items():
        if tag.startswith(prefix):
            return fam
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="data/seeds/seeds_family.faa")
    ap.add_argument("--outdir", default="data/seeds/families")
    ap.add_argument("--min-len", type=int, default=80, help="最短序列长度")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    families = {}
    skipped = []
    cur_hdr = None
    cur_seq = []
    cur_fam = None
    with open(args.seeds) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur_hdr and cur_fam and cur_seq:
                    seq = "".join(cur_seq)
                    if len(seq) >= args.min_len:
                        families.setdefault(cur_fam, []).append((cur_hdr, seq))
                    else:
                        skipped.append((cur_hdr, len(seq)))
                parts = line[1:].split("|")
                cur_hdr = line[1:]
                cur_fam = target_family(parts[1] if len(parts) > 1 else "unknown")
                cur_seq = []
            else:
                cur_seq.append(line)
    if cur_hdr and cur_fam and cur_seq:
        seq = "".join(cur_seq)
        if len(seq) >= args.min_len:
            families.setdefault(cur_fam, []).append((cur_hdr, seq))
        else:
            skipped.append((cur_hdr, len(seq)))

    print("Family | seqs")
    for fam in ["ePhaZ", "iPhaZ", "OH", "BdhA", "phasin"]:
        items = families.get(fam, [])
        out = os.path.join(args.outdir, f"{fam}.faa")
        with open(out, "w") as fo:
            for hdr, seq in items:
                fo.write(f">{hdr}\n{seq}\n")
        print(f"  {fam}: {len(items)} -> {out}")
    unassigned = [t for t in [] if False]
    if skipped:
        print(f"skipped (short): {len(skipped)}", file=sys.stderr)
        for h, l in skipped[:10]:
            print(f"  {h}: {l}aa", file=sys.stderr)
    # 未映射到5家族的
    other = set()
    with open(args.seeds) as f:
        for line in f:
            if line.startswith(">"):
                parts = line[1:].split("|")
                tag = parts[1] if len(parts) > 1 else ""
                if target_family(tag) is None:
                    other.add(tag)
    if other:
        print("unmapped tags:", sorted(other))


if __name__ == "__main__":
    main()
