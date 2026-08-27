#!/usr/bin/env python3
"""
11_clusters.py — 基因簇共定位分析（位点级 ±flank_kb 邻域）

对每个含 PHB 降解基因命中的基因组：
  1. 重跑 Pyrodigal 获取全部 CDS 的基因组坐标（GFF + 蛋白 FASTA）
  2. 用标记家族 HMM（PhaC/PhaE/PhaJ/BdhA/phasin/PHA_gran_rgn/PhaA/PhaB/PhaP/PhaR）
     对该基因组全部蛋白做 hmmsearch，得到标记基因的位点 + bit score（多 HMM 命中保留仲裁记录）
  3. 对每个命中位点，在同 contig 的 ±flank_kb 内查找标记基因，输出距离/方向/bit score
  4. 输出：
     - cluster_context.tsv   位点级：每个 (命中位点 × 邻近标记) 一行，含距离/方向/bit score/仲裁
     - cluster_summary.tsv   家族 × 标记：marker_hits(出现次数) + supporting_loci(唯一命中位点)
                             + supporting_genomes(唯一基因组)

fail-closed：--hits 必须是含 locus 列的位点级表（hits_filtered.tsv），
缺少 family/genome/locus 列时直接报错退出（不再静默跳过或记警告后继续）。

用法（服务器 T141，conda activate phb_gtdb）:
  python scripts/11_clusters.py \
      --hits data/screen/hits_filtered.tsv \
      --marker-hmms data/hmms/v2 \
      --gtdb ~/GTDB/gtdb_genomes_reps_r232/database \
      --flank-kb 10 --threads 40 --max-genomes 0
"""
import argparse
import os
import shutil
import subprocess
import sys
from collections import Counter, defaultdict

# 标记家族：PHA 合成/颗粒/动员相关（用于判断基因簇背景）
# PhaE = 古菌 PHA 合成酶亚基(PF09712, PHA_synth_III_E)；PHA_gran_rgn = 颗粒区蛋白(PF09650)
MARKER_FAMILIES = ["PhaC", "PhaE", "PhaA", "PhaB", "PhaP", "PhaR", "PhaJ", "BdhA", "phasin", "PHA_gran_rgn"]

# PHB 代谢上下文标记（用于 patatin 二次过滤的"局部邻域支持"判据，见 docs/STATUS.md）
PHB_CONTEXT_MARKERS = ["PhaC", "PhaE", "PhaJ", "BdhA", "phasin", "PHA_gran_rgn"]

# 关注的核心降解家族（对命中位点做邻域分析）
CORE_FAMILIES = ["ePhaZ", "iPhaZ", "OH", "ArchPhaZ_patatin", "ArchPhaZ_hydrolase"]


class ClusterInputError(RuntimeError):
    """Raised when a complete cluster analysis lacks required inputs."""


PYRODIGAL = shutil.which("pyrodigal") or os.path.expanduser(
    "~/miniconda3/envs/phb_gtdb/bin/pyrodigal"
)
HMMSEARCH = shutil.which("hmmsearch") or os.path.expanduser(
    "~/miniconda3/envs/phb_gtdb/bin/hmmsearch"
)


def resolve_marker_hmms(marker_dir, marker_families):
    """Return all declared marker HMMs, failing rather than silently changing the assay."""
    missing = []
    paths = []
    for family in marker_families:
        path = os.path.join(marker_dir, f"{family}.hmm")
        if not os.path.isfile(path):
            missing.append(family)
        else:
            paths.append(path)
    if missing:
        raise ClusterInputError(f"missing marker HMMs: {', '.join(missing)}")
    return paths


def write_cluster_summary(out_path, marker_hits, supporting_loci, supporting_genomes):
    """Write auditable occurrence, unique-locus, and unique-genome counts."""
    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write("hit_family\tmarker_family\tmarker_hits\tsupporting_loci\tsupporting_genomes\n")
        for key in sorted(marker_hits, key=lambda item: (-marker_hits[item], item)):
            hf, mf = key
            handle.write(
                f"{hf}\t{mf}\t{marker_hits[key]}\t"
                f"{len(supporting_loci.get(key, set()))}\t"
                f"{len(supporting_genomes.get(key, set()))}\n"
            )


def write_cluster_audit(out_path, rows):
    """Write one terminal analysis status for every requested hit locus."""
    columns = ("genome", "locus", "family", "status")
    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[column]) for column in columns) + "\n")


def write_cluster_genome_audit(out_path, rows):
    """Summarize terminal locus statuses per requested genome."""
    statuses = defaultdict(Counter)
    for row in rows:
        statuses[row["genome"]][row["status"]] += 1
    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write("genome\trequested_loci\tanalyzed_loci\tnot_analyzed_statuses\n")
        for genome in sorted(statuses):
            counts = statuses[genome]
            not_analyzed = ";".join(
                f"{status}:{count}"
                for status, count in sorted(counts.items())
                if status != "analyzed"
            )
            handle.write(
                f"{genome}\t{sum(counts.values())}\t{counts['analyzed']}\t{not_analyzed or 'none'}\n"
            )


def require_complete_locus_audit(rows):
    """Refuse a complete analysis when any requested locus was not analyzed."""
    incomplete = [row for row in rows if row["status"] != "analyzed"]
    if incomplete:
        counts = Counter(row["status"] for row in incomplete)
        reason_text = ", ".join(f"{reason}={count}" for reason, count in sorted(counts.items()))
        raise ClusterInputError(f"incomplete cluster analysis: {reason_text}")


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


def _hmmsearch_faa(faa, workdir):
    """Exclude HMMER-overlong targets as a tool-limit audit, not a biological no-hit."""
    filtered = os.path.join(workdir, os.path.basename(faa) + ".hmmsearch.faa")
    excluded = 0
    invalid = 0
    if os.path.getsize(faa) == 0:
        invalid = 1
    with open(faa) as source, open(filtered, "w") as target:
        header = None
        seq = []
        def flush():
            nonlocal excluded, invalid
            if header is None:
                return
            sequence = "".join(seq)
            if not header or header == ">" or not header[1:].strip() or not sequence:
                invalid += 1
            elif len(sequence) > 100000:
                excluded += 1
            else:
                target.write(header + "\n" + sequence + "\n")
        for line in source:
            if line.startswith(">"):
                flush()
                header = line.rstrip("\n")
                seq = []
            elif header is not None:
                seq.append(line.strip())
        flush()
    if invalid:
        with open(os.path.join(workdir, "invalid_fasta_records.tsv"), "a") as handle:
            handle.write(f"{os.path.basename(faa)}\t{invalid}\tinvalid_fasta_record\n")
    return filtered, excluded


def annotate_markers(faa, marker_hmms, workdir, threads, evalue="1e-5"):
    """用标记家族 HMM 对基因组蛋白做 hmmsearch。

    返回 {cds_id: [(marker_family, bitscore, evalue), ...]}，按 bitscore 降序。
    一个蛋白命中多个标记 HMM 时全部保留（用于多 HMM 命中仲裁记录）。
    """
    hits = defaultdict(list)
    scan_faa, excluded = _hmmsearch_faa(faa, workdir)
    if excluded:
        with open(os.path.join(workdir, "hmmsearch_overlong_exclusions.tsv"), "a") as handle:
            handle.write(f"{os.path.basename(faa)}\t{excluded}\thmmsearch_target_length_gt_100000\n")
    if os.path.getsize(scan_faa) == 0:
        return hits
    for mh in sorted(marker_hmms):
        fam = os.path.basename(mh)[:-4]
        tbl = os.path.join(workdir, f"marker_{fam}.tbl")
        if not os.path.isfile(HMMSEARCH) or not os.access(HMMSEARCH, os.X_OK):
            raise ClusterInputError(f"hmmsearch executable unavailable: {HMMSEARCH}")
        run([HMMSEARCH, "--tblout", tbl, "-E", evalue, "--cpu",
             str(max(1, min(threads, 4))), mh, scan_faa], timeout=600)
        with open(tbl) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                c = line.split()
                # tblout: 0 target,1 tacc,2 qname,3 qacc,4 E-value,5 score(bits),6 bias,7 domE
                if len(c) >= 6:
                    try:
                        hits[c[0]].append((fam, float(c[5]), float(c[4])))
                    except ValueError:
                        continue
    # 每个 cds_id 按 bitscore 降序（仲裁：bitscore 最高者为最佳标记家族）
    for cid in hits:
        hits[cid].sort(key=lambda x: -x[1])
    return hits


def read_hit_loci(hits_path, focus_fams):
    """读取命中位点：返回 [(genome, locus, family)]，family 过滤到关注家族。
    fail-closed：缺少 family/genome/locus 列时抛 ValueError（提示应传入 hits_filtered.tsv）。"""
    with open(hits_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for required in ("family", "genome", "locus"):
            if required not in header:
                raise ValueError(
                    f"[ERROR] --hits 缺少必需列 '{required}'（实际列: {header}）。"
                    f"基因簇分析需要位点级输入，请传入 data/screen/hits_filtered.tsv，"
                    f"而非 genome_hits.tsv（基因组×家族矩阵，无 locus）。")
        fi = header.index("family")
        gi = header.index("genome")
        li = header.index("locus")
        rows = []
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(fi, gi, li):
                continue
            if p[fi] in focus_fams:
                rows.append((p[gi], p[li], p[fi]))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", default="data/screen/hits_filtered.tsv")
    ap.add_argument("--hits-fasta", default="",
                    help="替代 --hits：从 FASTA 头（genome|locus）读命中位点，用于 tier1 序列集精确限定")
    ap.add_argument("--families", default=",".join(CORE_FAMILIES),
                    help="关注的核心降解家族（逗号分隔）")
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
    marker_hmms = resolve_marker_hmms(args.marker_hmms, marker_fams)
    print(f"标记家族: {[os.path.basename(x)[:-4] for x in marker_hmms]}")

    focus_fams = set(x for x in args.families.split(",") if x)
    if args.hits_fasta:
        # 从 FASTA 头读位点（>{genome}|{locus}），家族取 --families（单个）
        fam0 = next(iter(focus_fams)) if len(focus_fams) == 1 else "unknown"
        hits = []
        for line in open(args.hits_fasta):
            if line.startswith(">"):
                h = line[1:].strip().split()[0]
                parts = h.split("|")
                if len(parts) >= 2:
                    hits.append((parts[0], parts[1], fam0))
        print(f"命中位点(from fasta): {len(hits)}（关注家族: {sorted(focus_fams)}）")
    else:
        hits = read_hit_loci(args.hits, focus_fams)
        print(f"命中位点: {len(hits)}（关注家族: {sorted(focus_fams)}）")

    # 按基因组分组
    by_genome = defaultdict(list)
    for g, l, fam in hits:
        by_genome[g].append((l, fam))
    genomes = sorted(by_genome)
    sampled = args.max_genomes > 0 and args.max_genomes < len(genomes)
    selected_genomes = set(genomes)
    if args.max_genomes > 0:
        genomes = genomes[: args.max_genomes]
        selected_genomes = set(genomes)
    print(f"待分析基因组: {len(genomes)}" + ("（抽样）" if sampled else "（全部）"))

    flank = args.flank_kb * 1000
    ctx_rows = []
    # 家族×标记统计：marker_hits(出现次数) / supporting_loci(唯一命中位点) / supporting_genomes
    marker_hits = Counter()          # (hit_family, marker_family) -> 出现次数
    supporting_loci = defaultdict(set)   # (hit_family, marker_family) -> set(hit_locus)
    supporting_genomes = defaultdict(set)  # (hit_family, marker_family) -> set(genome)
    # patatin 局部邻域支持：位点/基因组（邻近任一 PHB 代谢标记）
    pat_context_loci = set()
    pat_context_genomes = set()
    pat_total_loci = 0
    pat_total_genomes = set()
    audit_rows = []
    for genome, locus, family in hits:
        if genome not in selected_genomes:
            audit_rows.append({
                "genome": genome, "locus": locus, "family": family,
                "status": "not_analyzed_sampled_out",
            })

    for gi, g in enumerate(genomes):
        gpath = find_genome_path(g, args.gtdb)
        if not gpath:
            audit_rows.extend(
                {"genome": g, "locus": locus, "family": family, "status": "missing_genome"}
                for locus, family in by_genome[g]
            )
            continue
        gff = os.path.join(args.workdir, f"{g}.gff")
        faa = os.path.join(args.workdir, f"{g}.faa")
        # 缓存：已生成则复用（断点续跑）
        if not (os.path.exists(gff) and os.path.exists(faa)):
            if not os.path.isfile(PYRODIGAL) or not os.access(PYRODIGAL, os.X_OK):
                raise ClusterInputError(f"pyrodigal executable unavailable: {PYRODIGAL}")
            run([PYRODIGAL, "-i", gpath, "-o", gff, "-f", "gff", "-a", faa, "-p", "meta"])
        genes = parse_gff(gff)
        if not genes:
            audit_rows.extend(
                {"genome": g, "locus": locus, "family": family, "status": "no_cds"}
                for locus, family in by_genome[g]
            )
            continue
        # 标记基因位点（cds_id -> 命中列表，含 bitscore + 仲裁）
        marker_anno = annotate_markers(faa, marker_hmms, args.workdir, args.threads)
        # 组织成 contig -> [(pos, cds_id, best_family, bitscore, arbitration)]
        contig_markers = defaultdict(list)
        for cid, matches in marker_anno.items():
            if cid in genes:
                ctg, s, e, _ = genes[cid]
                best_fam = matches[0][0]
                best_score = matches[0][1]
                arbitration = "|".join(f"{m[0]}:{m[1]:.1f}" for m in matches)
                contig_markers[ctg].append(((s + e) // 2, cid, best_fam, best_score, arbitration))

        # 每个命中位点做 ±flank 邻域检索
        for locus, hfam in by_genome[g]:
            if locus not in genes:
                audit_rows.append({
                    "genome": g, "locus": locus, "family": hfam, "status": "locus_not_found",
                })
                continue
            ctg, s, e, strand = genes[locus]
            mid = (s + e) // 2
            if hfam == "ArchPhaZ_patatin":
                pat_total_loci += 1
                pat_total_genomes.add(g)
            nearby_marker_fams = set()
            for (mpos, cid, mfam, mscore, arbitration) in contig_markers.get(ctg, []):
                if cid == locus:
                    continue  # 排除自身
                if abs(mpos - mid) <= flank:
                    signed = mpos - mid
                    direction = "downstream" if signed > 0 else ("upstream" if signed < 0 else "overlap")
                    nearby_marker_fams.add(mfam)
                    marker_hits[(hfam, mfam)] += 1
                    supporting_loci[(hfam, mfam)].add(locus)
                    supporting_genomes[(hfam, mfam)].add(g)
                    ctx_rows.append({
                        "genome": g, "contig": ctg, "hit_locus": locus,
                        "hit_family": hfam, "hit_start": s, "hit_end": e, "hit_strand": strand,
                        "marker_locus": cid, "marker_family": mfam,
                        "marker_bitscore": f"{mscore:.1f}",
                        "distance_bp": signed, "direction": direction,
                        "arbitration": arbitration,
                    })
            if hfam == "ArchPhaZ_patatin" and nearby_marker_fams & set(PHB_CONTEXT_MARKERS):
                pat_context_loci.add((g, locus))
                pat_context_genomes.add(g)
            audit_rows.append({"genome": g, "locus": locus, "family": hfam, "status": "analyzed"})
        if (gi + 1) % 100 == 0:
            print(f"  ... {gi + 1}/{len(genomes)}")

    # 写位点级邻域表
    out_ctx = f"{args.outdir}/tables/cluster_context.tsv"
    ctx_cols = ["genome", "contig", "hit_locus", "hit_family", "hit_start", "hit_end", "hit_strand",
                "marker_locus", "marker_family", "marker_bitscore", "distance_bp", "direction", "arbitration"]
    with open(out_ctx, "w") as f:
        f.write("\t".join(ctx_cols) + "\n")
        for r in ctx_rows:
            f.write("\t".join(str(r[k]) for k in ctx_cols) + "\n")
    print(f"\n位点级邻域表: {out_ctx}（{len(ctx_rows)} 行）")

    # 写家族×标记共现矩阵（区分"出现次数"与"唯一支持位点/基因组"）
    out_sum = f"{args.outdir}/tables/cluster_summary.tsv"
    write_cluster_summary(out_sum, marker_hits, supporting_loci, supporting_genomes)
    out_audit = f"{args.outdir}/tables/cluster_locus_audit.tsv"
    out_genome_audit = f"{args.outdir}/tables/cluster_genome_audit.tsv"
    write_cluster_audit(out_audit, audit_rows)
    write_cluster_genome_audit(out_genome_audit, audit_rows)
    print(f"locus audit: {out_audit}; genome audit: {out_genome_audit}")
    if not sampled:
        require_complete_locus_audit(audit_rows)

    print(f"\n家族×标记共现（前 20，按 marker_hits）:")
    for (hf, mf), n in marker_hits.most_common(20):
        print(f"  {hf} + {mf}: {n} 出现 / {len(supporting_loci[(hf, mf)])} 位点 / "
              f"{len(supporting_genomes[(hf, mf)])} 基因组")

    # 汇总：patatin 局部邻域支持（PHB 代谢上下文 = 任一 PHB_CONTEXT_MARKERS）
    print(f"\n[patatin] 位点总数={pat_total_loci}, 基因组总数={len(pat_total_genomes)}")
    print(f"[patatin] 局部邻域支持（邻近任一 PHB 代谢标记）: "
          f"位点={len(pat_context_loci)}, 基因组={len(pat_context_genomes)}")


if __name__ == "__main__":
    main()
