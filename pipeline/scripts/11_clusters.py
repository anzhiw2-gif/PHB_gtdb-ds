#!/usr/bin/env python3
"""
11_clusters.py — 基因簇共定位分析（位点级 ±flank_kb 邻域）

补齐此前的占位实现：对每个含 PHB 降解基因命中的基因组，
  1. 重跑 Pyrodigal 获取全部 CDS 的基因组坐标（GFF + 蛋白 FASTA）
  2. 用标记家族 HMM（PhaC/PhaA/PhaB/PhaP/PhaR/PhaJ/BdhA/phasin）对
     该基因组全部蛋白做 hmmsearch，得到标记基因的位点
  3. 对每个命中位点，在同 contig 的 ±flank_kb 内查找标记基因
  4. 输出每个命中位点的邻域标记（cluster_context.tsv）与
     基因组级家族×标记共现矩阵（cluster_summary.tsv）

用途：
  - patatin 家族二次过滤（真 PhaZh1 型解聚酶须与 PhaC/PhaP 等颗粒基因共定位，
    无邻近 PhaC 的 patatin 命中多为广谱磷脂酶，应剔除）
  - 降解-动员通路共现（PhaZ + OH + BdhA + PhaJ 的基因簇背景）

用法（服务器 T141，conda activate phb_gtdb）:
  python scripts/11_clusters.py \
      --hits data/screen/hits_filtered.tsv \
      --marker-hmms data/hmms/v2 \
      --gtdb ~/GTDB/gtdb_genomes_reps_r232/database \
      --flank-kb 10 --threads 40 --max-genomes 0

说明：
  - 蛋白 header 形如 ACC|{contig}_{orf}（见 05_predict_proteins.sh），
    重跑 pyrodigal 时 GFF 的 CDS ID = {contig}_{orf}，可直接对齐命中位点。
  - GFF 缓存于 --workdir（可断点续跑）；pyrodigal 使用与 05 一致的 -p meta。
  - 标记家族可经 --marker-families 覆盖；默认取 hmm 目录中真实存在的子集。
"""
import argparse
import os
import subprocess
import sys
from collections import Counter, defaultdict

# 标记家族：PHA 合成/颗粒/动员相关（用于判断基因簇背景）
MARKER_FAMILIES = ["PhaC", "PhaA", "PhaB", "PhaP", "PhaR", "PhaJ", "BdhA", "phasin"]

# 关注的核心降解家族（对命中位点做邻域分析；默认取全部命中）
CORE_FAMILIES = ["ePhaZ", "iPhaZ", "OH", "ArchPhaZ_patatin", "ArchPhaZ_hydrolase"]


def parse_gff(gff_path):
    """返回 {cds_id: (contig, start, end, strand)}；cds_id 取 GFF 第 9 列 ID 属性。"""
    genes = {}
    if not os.path.exists(gff_path):
        return genes
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "CDS":
                continue
            attrs = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
            cid = attrs.get("ID") or attrs.get("locus_tag") or ""
            if not cid:
                continue
            genes[cid] = (p[0], int(p[3]), int(p[4]), p[6])
    return genes


def find_genome_path(acc, gtdb_root):
    """GTDB 路径: database/{GCA|GCF}/xx/yyy/zzz/{ACC}_genomic.fna.gz"""
    for prefix in ("GCA", "GCF"):
        cand = os.path.join(gtdb_root, prefix, acc[4:7], acc[7:10], acc[10:13],
                            f"{acc}_genomic.fna.gz")
        if os.path.exists(cand):
            return cand
    return None


def run(cmd, timeout=1200):
    """执行外部命令；失败抛异常（带命令回显，便于排错）。"""
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError(f"command failed ({r.returncode}): {' '.join(cmd)}\n{r.stderr[:800]}")
    return r


def annotate_markers(faa, marker_hmms, workdir, threads, evalue="1e-5"):
    """用标记家族 HMM 对基因组蛋白做 hmmsearch，返回 {cds_id: marker_family}。"""
    anno = {}
    for mh in sorted(marker_hmms):
        fam = os.path.basename(mh)[:-4]
        tbl = os.path.join(workdir, f"marker_{fam}.tbl")
        run(["hmmsearch", "--tblout", tbl, "-E", evalue, "--cpu",
             str(max(1, min(threads, 4))), mh, faa], timeout=600)
        with open(tbl) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                c = line.split()
                if len(c) >= 2:
                    # tblout 第 1 列 = target（蛋白 ID），只保留首个（最优）命中家族
                    anno.setdefault(c[0], fam)
    return anno


def read_hit_loci(hits_path):
    """读取命中位点：返回 [(genome, locus, family)]，family 过滤到核心降解家族。"""
    rows = []
    with open(hits_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        fi = header.index("family")
        gi = header.index("genome")
        li = header.index("locus")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(fi, gi, li):
                continue
            if p[fi] in CORE_FAMILIES:
                rows.append((p[gi], p[li], p[fi]))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", default="data/screen/hits_filtered.tsv")
    ap.add_argument("--marker-hmms", default="data/hmms/v2", help="标记家族 HMM 目录")
    ap.add_argument("--marker-families", default=",".join(MARKER_FAMILIES))
    ap.add_argument("--gtdb", default=os.path.expanduser("~/GTDB/gtdb_genomes_reps_r232/database"))
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--workdir", default="/tmp/cluster_work")
    ap.add_argument("--flank-kb", type=int, default=10)
    ap.add_argument("--threads", type=int, default=40)
    ap.add_argument("--max-genomes", type=int, default=0, help="0=全部命中基因组")
    args = ap.parse_args()

    os.makedirs(f"{args.outdir}/tables", exist_ok=True)
    os.makedirs(args.workdir, exist_ok=True)

    marker_fams = [x for x in args.marker_families.split(",") if x]
    marker_hmms = []
    for fam in marker_fams:
        p = os.path.join(args.marker_hmms, f"{fam}.hmm")
        if os.path.exists(p):
            marker_hmms.append(p)
        else:
            print(f"[warn] 标记 HMM 缺失，跳过: {fam}")
    print(f"标记家族: {[os.path.basename(x)[:-4] for x in marker_hmms]}")

    hits = read_hit_loci(args.hits)
    print(f"命中位点: {len(hits)}（核心降解家族）")

    # 按基因组分组
    by_genome = defaultdict(list)
    for g, l, fam in hits:
        by_genome[g].append((l, fam))
    genomes = sorted(by_genome)
    if args.max_genomes > 0:
        genomes = genomes[: args.max_genomes]
    print(f"待分析基因组: {len(genomes)}")

    flank = args.flank_kb * 1000
    ctx_rows = []
    cooccur = Counter()   # (hit_family, marker_family) -> 位点数
    for gi, g in enumerate(genomes):
        gpath = find_genome_path(g, args.gtdb)
        if not gpath:
            continue
        gff = os.path.join(args.workdir, f"{g}.gff")
        faa = os.path.join(args.workdir, f"{g}.faa")
        # 缓存：已生成则复用（断点续跑）
        if not (os.path.exists(gff) and os.path.exists(faa)):
            run(["pyrodigal", "-i", gpath, "-o", gff, "-f", "gff", "-a", faa, "-p", "meta"])
        genes = parse_gff(gff)
        if not genes:
            continue
        # 标记基因位点（locus -> marker_family）
        marker_loci = annotate_markers(faa, marker_hmms, args.workdir, args.threads)
        # 组织成 contig -> [(pos, marker_family)]
        contig_markers = defaultdict(list)
        for locus, mfam in marker_loci.items():
            if locus in genes:
                ctg, s, e, _ = genes[locus]
                contig_markers[ctg].append(((s + e) // 2, mfam))

        # 每个命中位点做 ±flank 邻域检索
        for locus, hfam in by_genome[g]:
            if locus not in genes:
                continue
            ctg, s, e, strand = genes[locus]
            mid = (s + e) // 2
            nearby = []
            for (mpos, mfam) in contig_markers.get(ctg, []):
                if abs(mpos - mid) <= flank:
                    nearby.append(mfam)
                    cooccur[(hfam, mfam)] += 1
            ctx_rows.append({
                "genome": g, "contig": ctg, "hit_locus": locus,
                "hit_family": hfam, "start": s, "end": e, "strand": strand,
                "nearby_markers": ";".join(sorted(set(nearby))),
            })
        if (gi + 1) % 100 == 0:
            print(f"  ... {gi + 1}/{len(genomes)}")

    # 写位点级邻域表
    out_ctx = f"{args.outdir}/tables/cluster_context.tsv"
    with open(out_ctx, "w") as f:
        f.write("genome\tcontig\thit_locus\thit_family\tstart\tend\tstrand\tnearby_markers\n")
        for r in ctx_rows:
            f.write("\t".join(str(r[k]) for k in
                     ["genome", "contig", "hit_locus", "hit_family",
                      "start", "end", "strand", "nearby_markers"]) + "\n")
    print(f"\n位点级邻域表: {out_ctx}（{len(ctx_rows)} 行）")

    # 写家族×标记共现矩阵
    out_sum = f"{args.outdir}/tables/cluster_summary.tsv"
    with open(out_sum, "w") as f:
        f.write("hit_family\tmarker_family\tcooccurring_loci\n")
        for (hf, mf), n in sorted(cooccur.items(), key=lambda x: -x[1]):
            f.write(f"{hf}\t{mf}\t{n}\n")

    print(f"家族×标记共现（前 20）:")
    for (hf, mf), n in cooccur.most_common(20):
        print(f"  {hf} + {mf}: {n} 位点")

    # 汇总：含 PhaC/PhaP 邻近的 patatin 位点数（patatin 二次过滤的直接依据）
    pat_total = sum(1 for r in ctx_rows if r["hit_family"] == "ArchPhaZ_patatin")
    pat_phac = sum(1 for r in ctx_rows if r["hit_family"] == "ArchPhaZ_patatin"
                   and ("PhaC" in r["nearby_markers"] or "PhaP" in r["nearby_markers"]))
    print(f"\n[patatin] 位点总数={pat_total}, 邻近 PhaC/PhaP={pat_phac}")


if __name__ == "__main__":
    main()
