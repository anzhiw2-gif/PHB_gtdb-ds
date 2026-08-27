#!/usr/bin/env python3
"""test_smoke.py — 小样本端到端冒烟测试（不依赖 GTDB 全库，无需 hmmsearch）

覆盖关键契约（对应审查阻断项/重要改进）：
  1. 06b_aggregate_hits.py：tbl/dom → hits_all.tsv 含 cov 列（=(hmm_to-hmm_from+1)/qlen）
  2. 07_process_hits.py：min-cov 过滤真正生效 + hits_filtered.tsv 含 locus
  3. 08_validate.py：NAD 基序正则（[^A-Z]→[A-Z] 修复后能匹配 GXGXXG）
  4. 11_clusters.py：read_hit_loci 缺 locus 列时 fail-closed（抛 ValueError）

用法: python pipeline/tests/test_smoke.py [--tmpdir ...]
退出码: 0=通过, 1=失败
"""
import argparse
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile

SCRIPTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts")


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def run_py(args, cwd):
    return subprocess.run([sys.executable, *args], cwd=cwd,
                          capture_output=True, text=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tmpdir", default=None)
    args = ap.parse_args()

    own_tmp = False
    if args.tmpdir:
        tmp = args.tmpdir
    else:
        # 默认在仓库内建临时目录（本地 sandbox 下系统 temp 可能不可写）
        tmp = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".smoke_tmp")
        shutil.rmtree(tmp, ignore_errors=True)
        os.makedirs(tmp, exist_ok=True)
        own_tmp = True
    os.makedirs(os.path.join(tmp, "hmmsearch"), exist_ok=True)
    os.makedirs(os.path.join(tmp, "shards"), exist_ok=True)

    # ---- 1. 伪造 hmmsearch 输出（tbl + dom），三个命中，覆盖度 1.0 / 0.5 / 0.2 ----
    tbl = os.path.join(tmp, "hmmsearch", "ePhaZ__shard_0001.tbl")
    dom = os.path.join(tmp, "hmmsearch", "ePhaZ__shard_0001.dom")
    with open(tbl, "w", encoding="utf-8") as f:
        f.write("# tblout header\n")
        for tgt, e in [("GCA_1|ctg_1", "1e-30"), ("GCA_1|ctg_2", "1e-25"), ("GCA_2|ctg_1", "1e-10")]:
            f.write(f"{tgt}  -  ePhaZ  -  {e}  120.0  0.1  {e}  120.0  0.1  "
                    f"1  1  0  0  1  1  1  1  -\n")
    with open(dom, "w", encoding="utf-8") as f:
        f.write("# domtblout header\n")
        for tgt, hto in [("GCA_1|ctg_1", 200), ("GCA_1|ctg_2", 100), ("GCA_2|ctg_1", 40)]:
            f.write(f"{tgt}  -  300  ePhaZ  -  200  1e-30  120.0  0.1  1  1  "
                    f"1e-30  1e-30  120.0  0.1  1  {hto}  1  {hto}  1  {hto}  0.99  -\n")

    # ---- 2. 06b 聚合，检查 cov 列 ----
    hits_all = os.path.join(tmp, "hits_all.tsv")
    r = run_py([os.path.join(SCRIPTS, "06b_aggregate_hits.py"),
                "--hmmout", os.path.join(tmp, "hmmsearch"), "--out", hits_all], tmp)
    assert r.returncode == 0, f"06b 失败: {r.stderr}"
    header = open(hits_all, encoding="utf-8").readline().rstrip("\n").split("\t")
    assert "cov" in header, f"hits_all.tsv 缺 cov 列: {header}"
    covs = {}
    for line in open(hits_all, encoding="utf-8"):
        if line.startswith("family"):
            continue
        p = line.rstrip("\n").split("\t")
        covs[p[2]] = float(p[9])
    assert abs(covs["GCA_1|ctg_1"] - 1.0) < 1e-6, covs
    assert abs(covs["GCA_1|ctg_2"] - 0.5) < 1e-6, covs
    assert abs(covs["GCA_2|ctg_1"] - 0.2) < 1e-6, covs
    print("[ok] 06b 覆盖度列正确:", covs)

    # ---- 3. 07 命中处理：min-cov 0.5 应丢弃 cov=0.2 的命中 ----
    r = run_py([os.path.join(SCRIPTS, "07_process_hits.py"),
                "--hits", hits_all, "--shards", os.path.join(tmp, "shards"),
                "--outdir", tmp, "--min-cov", "0.5"], tmp)
    assert r.returncode == 0, f"07 失败: {r.stderr}\n{r.stdout}"
    filt = os.path.join(tmp, "hits_filtered.tsv")
    assert os.path.exists(filt), "hits_filtered.tsv 未生成"
    filt_lines = [l for l in open(filt, encoding="utf-8") if not l.startswith("family")]
    assert len(filt_lines) == 2, f"min-cov 过滤后应剩 2 行，实际 {len(filt_lines)}"
    filt_header = open(filt, encoding="utf-8").readline().rstrip("\n").split("\t")
    assert "locus" in filt_header, f"hits_filtered.tsv 缺 locus 列: {filt_header}"
    gh = os.path.join(tmp, "genome_hits.tsv")
    gh_lines = [l for l in open(gh, encoding="utf-8") if not l.startswith("genome")]
    assert len(gh_lines) == 1, f"genome_hits 应 1 行（GCA_1 × ePhaZ），实际 {len(gh_lines)}"
    assert "GCA_1\tePhaZ\t2" in gh_lines[0], gh_lines[0]
    print("[ok] 07 min-cov 过滤 + locus 列正确")

    # ---- 4. 08_validate NAD 基序正则（GXGXXG 应匹配） ----
    v08 = load_module("v08", os.path.join(SCRIPTS, "08_validate.py"))
    assert v08.NAD_BINDING.search("GSGAAG"), "NAD 基序未匹配 GXGXXG"
    assert not v08.NAD_BINDING.search("AAAAAA"), "NAD 基序误匹配"
    print("[ok] 08 NAD 基序正则正确")

    # ---- 5. 11_clusters read_hit_loci 缺 locus 列 fail-closed ----
    v11 = load_module("v11", os.path.join(SCRIPTS, "11_clusters.py"))
    bad = os.path.join(tmp, "no_locus.tsv")
    with open(bad, "w", encoding="utf-8") as f:
        f.write("family\tgenome\n")
        f.write("ePhaZ\tGCA_1\n")
    try:
        v11.read_hit_loci(bad, {"ePhaZ"})
        raise AssertionError("缺 locus 列未抛异常（应 fail-closed）")
    except ValueError:
        print("[ok] 11 缺 locus 列 fail-closed")

    print(f"\nALL SMOKE TESTS PASSED (tmp={tmp})")
    if own_tmp:
        shutil.rmtree(tmp, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
