#!/bin/bash
# rerun_candidates.sh — 服务器专用：重建候选集（方案 A：OH 家族 min-cov 0.6）
#   重跑 06b(聚合cov) → 07(按家族min-cov) → 07b(提取) → 08(验证) → 08c(tier) → 09a(汇总)
# fail-closed：任一步失败即退出。硬编码工作区根（服务器扁平 scripts/ 结构）。
# 用法: nohup bash scripts/rerun_candidates.sh > results/logs/rerun_candidates.log 2>&1 &
set -euo pipefail

ROOT="/home/data/haoyu/PHB_gtdb-ds"
cd "$ROOT"
PY=~/miniconda3/envs/phb_gtdb/bin/python
LOG="results/logs/rerun_candidates.log"
mkdir -p results/logs

echo "[$(date)] === 重建候选集（方案 A：OH min-cov 0.6）==="

echo "[1/6] 06b 聚合 hits_all.tsv（含 cov 列）"
"$PY" scripts/06b_aggregate_hits.py --hmmout data/screen/hmmsearch --out data/screen/hits_all.tsv

echo "[2/6] 07 命中处理（--family-min-cov OH:0.6）"
"$PY" scripts/07_process_hits.py --hits data/screen/hits_all.tsv \
    --shards data/proteins/shards_filt --outdir data/screen --family-min-cov OH:0.6

echo "[3/6] 07b 提取序列"
"$PY" scripts/07b_extract_seqs.py --ids data/screen/unique_proteins.txt \
    --hits data/screen/hits_filtered.tsv --shards data/proteins/shards_filt \
    --outdir data/screen/family_seqs

echo "[4/6] 08 功能验证"
"$PY" scripts/08_validate.py --indir data/screen/family_seqs --outdir data/screen

echo "[5/6] 08c tier 重评分"
"$PY" scripts/08c_tier_rescore.py

echo "[6/6] 09a tier1 汇总"
"$PY" scripts/09a_tier1_summary.py

echo "[$(date)] === 重建候选集完成 ==="
