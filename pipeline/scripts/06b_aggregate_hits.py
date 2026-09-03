#!/usr/bin/env python3
"""06b_aggregate_hits.py — 聚合 hmmsearch tbl/dom 文件为 hits_all.tsv（含覆盖度）

读取 data/screen/hmmsearch/{family}__shard_XXXX.tbl 与同名的 .dom 文件：
  - tblout 提供命中行（E-value/score/bias/domE/qname）
  - domtblout 提供每个 domain 的 HMM 坐标（hmm_from/hmm_to）与 query 长度(qlen)，
    据此计算每 (family, shard, target) 的最大 HMM 覆盖度：
        cov = max_domain (hmm_to - hmm_from + 1) / qlen
    覆盖度用于 07_process_hits.py 的 --min-cov 过滤（见 docs/STATUS.md）。

输出: data/screen/hits_all.tsv
    列: family  shard  protein  tacc  E-value  score  bias  domE  qname  cov

fail-closed：任何 tbl 缺失、或 tbl 与 dom 行数不一致（dom 可能因 hmmsearch 失败而缺）
都会报错退出（exit 1），绝不静默产出不完整结果。
用法: python 06b_aggregate_hits.py [--hmmout data/screen/hmmsearch] [--out data/screen/hits_all.tsv]
"""
import argparse
import glob
import os
import sys


class AggregateError(RuntimeError):
    """Raised when HMMER output provenance is incomplete."""


def validate_pairs(tbl_files):
    """Validate that every tbl has a matching dom file, including empty tbl files."""
    missing = []
    malformed = []
    for tbl in tbl_files:
        base = os.path.basename(tbl)[:-4]
        if "__" not in base:
            malformed.append(base)
        dom = os.path.join(os.path.dirname(tbl), base + ".dom")
        if not os.path.isfile(dom):
            missing.append(base)
    if malformed:
        raise AggregateError(f"malformed tbl names: {', '.join(malformed[:5])}")
    if missing:
        raise AggregateError(f"missing dom files: {', '.join(missing[:5])}")
    return []

# domtblout 列索引（0-based，按空白切分）
# 0 target,1 tacc,2 tlen,3 qname,4 qacc,5 qlen,6 fullE,7 fullscore,8 fullbias,
# 9 dom#,10 ndom,11 cE,12 iE,13 domscore,14 dombias,15 hmm_from,16 hmm_to,
# 17 ali_from,18 ali_to,19 env_from,20 env_to,21 acc,22 desc
D_HMM_FROM = 15
D_HMM_TO = 16
D_QLEN = 5


def read_dom_cov(path: str) -> dict:
    """返回 {target: max_query_coverage}。查询覆盖度 = (hmm_to-hmm_from+1)/qlen。"""
    intervals = {}
    if not os.path.exists(path):
        return {}
    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split()
            if len(p) <= D_HMM_TO:
                continue
            target = p[0]
            try:
                qlen = float(p[D_QLEN])
                hfrom = int(p[D_HMM_FROM])
                hto = int(p[D_HMM_TO])
            except ValueError:
                continue
            if qlen <= 0:
                continue
            record = intervals.setdefault(target, {"qlen": qlen, "ranges": []})
            record["ranges"].append((min(hfrom, hto), max(hfrom, hto)))
    cov = {}
    for target, record in intervals.items():
        covered = 0
        end = 0
        for start, stop in sorted(record["ranges"]):
            if stop <= end:
                continue
            covered += stop - max(start, end + 1) + 1
            end = stop
        cov[target] = min(1.0, covered / record["qlen"])
    return cov


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hmmout", default="data/screen/hmmsearch")
    ap.add_argument("--out", default="data/screen/hits_all.tsv")
    args = ap.parse_args()

    tbl_files = sorted(glob.glob(os.path.join(args.hmmout, "*.tbl")))
    if not tbl_files:
        print(f"[ERROR] 未找到任何 .tbl 文件于 {args.hmmout}，无法聚合", file=sys.stderr)
        sys.exit(1)

    try:
        validate_pairs(tbl_files)
    except AggregateError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        sys.exit(1)

    missing_dom = []
    rows = 0
    with open(args.out, "w") as fo:
        fo.write("family\tshard\tprotein\ttacc\tE-value\tscore\tbias\tdomE\tqname\tcov\n")
        for f in tbl_files:
            base = os.path.basename(f)[:-4]  # 去掉 .tbl
            if "__" not in base:
                print(f"[ERROR] tbl 文件名缺少 __ 分隔符: {base}", file=sys.stderr)
                sys.exit(1)
            fam, sname = base.split("__", 1)
            dom_path = os.path.join(args.hmmout, base + ".dom")
            cov = read_dom_cov(dom_path)
            with open(f, errors="replace") as fin:
                for line in fin:
                    if line.startswith("#"):
                        continue
                    p = line.rstrip("\n").split()
                    if len(p) < 8:
                        continue
                    # tblout: 0 target,1 tacc,2 qname,3 qacc,4 E,5 score,6 bias,7 domE,8 domscore
                    c = cov.get(p[0], 0.0)
                    fo.write(f"{fam}\t{sname}\t{p[0]}\t{p[1]}\t{p[4]}\t{p[5]}\t{p[6]}\t{p[7]}\t{p[2]}\t{c:.6f}\n")
                    rows += 1

    print(f"聚合完成: {len(tbl_files)} 个 tbl, {rows} 行命中 -> {args.out}")


if __name__ == "__main__":
    main()
