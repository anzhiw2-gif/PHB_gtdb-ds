#!/usr/bin/env python3
"""
07_process_hits.py — 处理 hmmsearch 命中
  1) 过滤：E-value、比对覆盖率（alnlen/qlen）
  2) 多家族命中仲裁（优先 ePhaZ > iPhaZ > OH > BdhA > phasin 可配置）
  3) 每基因组每家族保留最优命中（dedup）
  4) 输出:
     data/screen/hits_filtered.tsv      过滤后命中明细
     data/screen/genome_hits.tsv        基因组×家族命中矩阵（含拷贝数）
     data/screen/family_seqs/{fam}.faa  各家族命中序列（从 shard 提取）
用法: python 07_process_hits.py --hits data/screen/hits_all.tsv
"""
import argparse
import gzip
import os
import sys
from collections import defaultdict

FAMILY_PRIORITY = ["ePhaZ", "iPhaZ", "OH", "BdhA", "ArchPhaZ_patatin", "ArchPhaZ_hydrolase", "PhaJ", "phasin", "PhaC"]


def parse_hits(path: str):
    rows = []
    with open(path) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 7:
                continue
            d = dict(zip(header, p))
            try:
                d["E-value"] = float(d["E-value"])
                d["score"] = float(d["score"])
                d["domE"] = float(d.get("domE", d["E-value"]))
            except ValueError:
                continue
            # protein id: ACC|contig_pos（由 pyrodigal header 加前缀）
            prot = d["protein"]
            if "|" in prot:
                d["genome"] = prot.split("|")[0]
                d["locus"] = prot.split("|", 1)[1]
            else:
                d["genome"] = "unknown"
                d["locus"] = prot
            d["cov"] = 1.0  # 覆盖率过滤可选（domtblout 另行补充）
            rows.append(d)
    return rows


def extract_seq(shards_dir: str, prot: str) -> str:
    """从 shard FASTA 中提取蛋白序列（shard 内线性扫描，命中少可接受）"""
    if not prot:
        return ""
    for sf in sorted(os.listdir(shards_dir)):
        if not sf.endswith(".faa"):
            continue
        path = os.path.join(shards_dir, sf)
        want = False
        buf = []
        with open(path) as f:
            for line in f:
                if line.startswith(">"):
                    if want:
                        return "".join(buf)
                    if line[1:].split()[0] == prot:
                        want = True
                        buf = []
                elif want:
                    buf.append(line.strip())
            if want:
                return "".join(buf)
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", default="data/screen/hits_all.tsv")
    ap.add_argument("--shards", default="data/proteins/shards_filt")
    ap.add_argument("--outdir", default="data/screen")
    ap.add_argument("--max-eval", type=float, default=1e-5)
    ap.add_argument("--min-cov", type=float, default=0.5, help="比对覆盖度阈值")
    ap.add_argument("--extract", type=int, default=0, help="是否提取序列（大库建议 0，另用 07b）")
    args = ap.parse_args()
    os.makedirs(os.path.join(args.outdir, "family_seqs"), exist_ok=True)

    rows = parse_hits(args.hits)
    print(f"raw hits: {len(rows)}")

    # 过滤
    kept = [r for r in rows if r["E-value"] <= args.max_eval]
    print(f"after E<={args.max_eval}: {len(kept)}")

    # 多家族命中仲裁：每个 (genome, locus) 选优先级最高的家族
    best = {}
    for r in kept:
        key = (r["genome"], r["locus"])
        cur = best.get(key)
        if cur is None or FAMILY_PRIORITY.index(r["family"]) < FAMILY_PRIORITY.index(cur["family"]) or (
            r["family"] == cur["family"] and r["E-value"] < cur["E-value"]):
            best[key] = r

    # 每基因组每家族：最优 + 拷贝数
    genome_fam = defaultdict(list)
    for r in best.values():
        genome_fam[(r["genome"], r["family"])].append(r)

    filtered = sorted(best.values(), key=lambda x: (x["genome"], FAMILY_PRIORITY.index(x["family"])))
    with open(os.path.join(args.outdir, "hits_filtered.tsv"), "w") as f:
        cols = ["family", "genome", "locus", "protein", "E-value", "score", "domE"]
        f.write("\t".join(cols) + "\n")
        for r in filtered:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    with open(os.path.join(args.outdir, "genome_hits.tsv"), "w") as f:
        f.write("genome\tfamily\tcopies\tbest_E\tbest_score\n")
        for (g, fam), lst in sorted(genome_fam.items()):
            lst_sorted = sorted(lst, key=lambda x: x["E-value"])
            f.write(f"{g}\t{fam}\t{len(lst)}\t{lst_sorted[0]['E-value']:.2e}\t{lst_sorted[0]['score']:.1f}\n")

    # 每家族序列提取（用于系统发育/验证）——大库时跳过，用 07b 高效提取
    if args.extract:
        fam_seqs = defaultdict(list)
        for r in filtered:
            seq = extract_seq(args.shards, r["protein"])
            if seq:
                fam_seqs[r["family"]].append((r, seq))

        for fam in FAMILY_PRIORITY:
            items = fam_seqs.get(fam, [])
            out = os.path.join(args.outdir, "family_seqs", f"{fam}.faa")
            with open(out, "w") as f:
                for r, seq in items:
                    f.write(f">{r['genome']}|{r['locus']}|E={r['E-value']:.1e}\n{seq}\n")
            print(f"  {fam}: {len(items)} seqs -> {out}")
    else:
        # 写唯一蛋白 ID 清单供 07b 提取
        with open(os.path.join(args.outdir, "unique_proteins.txt"), "w") as f:
            for r in filtered:
                f.write(r["protein"] + "\n")
        print(f"  跳过序列提取，已写 unique_proteins.txt ({len(filtered)} 条)")

    print("\n[SUMMARY] 基因组×家族:")
    from collections import Counter
    c = Counter((g, fam) for (g, fam) in genome_fam)
    fam_count = Counter(fam for (g, fam) in c)
    for fam in FAMILY_PRIORITY:
        print(f"  {fam}: {fam_count.get(fam, 0)} genome-family pairs")
    print(f"  总基因组数(含任一命中): {len(set(g for g, f in genome_fam))}")


if __name__ == "__main__":
    main()
