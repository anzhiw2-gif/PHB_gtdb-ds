#!/bin/bash
# resume_0611.sh — 从筛选阶段续跑（shards 已完整）
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
LOG="$ROOT/results/logs/resume_0611.log"
mkdir -p "$ROOT/results/logs"
THREADS_SCREEN=80
THREADS_PHYLO=40
log() { echo "[$(date '+%F %T')] $1" | tee -a "$LOG"; }

log "=== 续跑 06-11 启动 ==="

log "[06] 全库 HMM 筛选（9 家族 × 100 shard, $THREADS_SCREEN 线程）"
bash scripts/06_screen.sh --threads "$THREADS_SCREEN" --eval 1e-5 2>&1 | tee -a "$LOG" || { log "[ERROR] 06 失败"; exit 1; }
log "[06] 筛选完成: $(wc -l < data/screen/hits_all.tsv) 行"

log "[07] 命中处理"
source ~/miniconda3/etc/profile.d/conda.sh
conda run -n phb_gtdb python scripts/07_process_hits.py \
    --hits data/screen/hits_all.tsv --shards data/proteins/shards \
    --outdir data/screen 2>&1 | tee -a "$LOG" || { log "[ERROR] 07 失败"; exit 1; }

log "[08] 功能验证"
conda run -n phb_gtdb python scripts/08_validate.py \
    --indir data/screen/family_seqs --outdir data/screen 2>&1 | tee -a "$LOG" || log "[WARN] 08 部分失败"

log "[09] 系统发育"
bash scripts/09_phylogeny.sh --threads "$THREADS_PHYLO" 2>&1 | tee -a "$LOG" || log "[WARN] 09 部分失败"

log "[10] 生态分布"
conda run -n phb_gtdb python scripts/10_distribution.py \
    --hits data/screen/genome_hits.tsv 2>&1 | tee -a "$LOG" || log "[WARN] 10 部分失败"

log "[11] 基因簇共定位"
conda run -n phb_gtdb python scripts/11_clusters.py \
    --hits data/screen/genome_hits.tsv --max-genomes 2000 2>&1 | tee -a "$LOG" || log "[WARN] 11 部分失败"

log "=== 续跑 06-11 全部完成 ==="
