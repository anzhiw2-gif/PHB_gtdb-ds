#!/bin/bash
# =============================================================================
# run_pipeline.sh — 主编排：等待预测完成 → 全库筛选 → 处理 → 验证 → 建树 → 分布
# 用法: nohup bash run_pipeline.sh > results/logs/pipeline_master.log 2>&1 &
# 阶段:
#   [wait] 等待 05 预测完成（shard 生成）
#   [06]   hmmsearch 9 家族 × shard
#   [07]   命中处理（过滤/仲裁/序列提取）
#   [08]   功能验证（催化位点/SignalP）
#   [09]   各家族系统发育
#   [10]   生态分布统计
#   [11]   基因簇共定位
# =============================================================================
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"   # 项目根（scripts 的上一级）
cd "$ROOT"
LOG="$ROOT/results/logs/pipeline_master.log"
mkdir -p "$ROOT/results/logs"
THREADS_SCREEN=80
THREADS_PHYLO=40

log() { echo "[$(date '+%F %T')] $1" | tee -a "$LOG"; }

log "=== 主管线启动 ==="

# ---------- 等待预测完成 ----------
log "[wait] 等待蛋白预测完成（shard 合并稳定）"
prev=0
stable=0
for i in $(seq 1 480); do   # 最多等 80 小时（每 5 分钟检查）
    n=$(find data/proteins/shards -maxdepth 1 -name 'shard_*.faa' 2>/dev/null | wc -l)
    done_g=$(find data/proteins/per_genome -maxdepth 1 -name '*.faa.gz' 2>/dev/null | wc -l)
    if [ "$n" -gt 0 ] && [ "$n" -eq "$prev" ] && [ "$n" -ge 100 ]; then
        # shard 数稳定且 >=100 → 合并完成（05 脚本每 2000 基因组一个 shard）
        log "[wait] shard 数稳定在 $n（预测+合并完成）"
        stable=1
        break
    fi
    prev="$n"
    if [ $((i % 12)) -eq 0 ]; then
        log "[wait] 预测进度: $done_g / 199923, shard: $n"
    fi
    sleep 300
done

if [ "$stable" -eq 0 ]; then
    n=$(find data/proteins/shards -maxdepth 1 -name 'shard_*.faa' 2>/dev/null | wc -l)
    if [ "$n" -eq 0 ]; then
        log "[ERROR] 预测未完成，退出"
        exit 1
    fi
fi

# ---------- 06 全库筛选 ----------
log "[06] 启动全库 HMM 筛选（9 家族 × $n shard, $THREADS_SCREEN 线程）"
bash scripts/06_screen.sh --threads "$THREADS_SCREEN" --eval 1e-5 2>&1 | tee -a "$LOG" || { log "[ERROR] 06 失败"; exit 1; }
log "[06] 筛选完成: $(wc -l < data/screen/hits_all.tsv) 行"

# ---------- 07 命中处理 ----------
log "[07] 命中处理"
source ~/miniconda3/etc/profile.d/conda.sh
conda run -n phb_gtdb python scripts/07_process_hits.py \
    --hits data/screen/hits_all.tsv --shards data/proteins/shards \
    --outdir data/screen 2>&1 | tee -a "$LOG" || { log "[ERROR] 07 失败"; exit 1; }

# ---------- 08 功能验证 ----------
log "[08] 功能验证（催化位点/lipase box/长度）"
conda run -n phb_gtdb python scripts/08_validate.py \
    --indir data/screen/family_seqs --outdir data/screen 2>&1 | tee -a "$LOG" || log "[WARN] 08 部分失败"

# ---------- 09 系统发育 ----------
log "[09] 各家族系统发育（$THREADS_PHYLO 线程）"
bash scripts/09_phylogeny.sh --threads "$THREADS_PHYLO" 2>&1 | tee -a "$LOG" || log "[WARN] 09 部分失败"

# ---------- 10 生态分布 ----------
log "[10] 生态/分类学分布统计"
conda run -n phb_gtdb python scripts/10_distribution.py \
    --hits data/screen/genome_hits.tsv 2>&1 | tee -a "$LOG" || log "[WARN] 10 部分失败"

# ---------- 11 基因簇共定位 ----------
log "[11] 基因簇共定位（前 2000 个命中基因组）"
conda run -n phb_gtdb python scripts/11_clusters.py \
    --hits data/screen/genome_hits.tsv --max-genomes 2000 2>&1 | tee -a "$LOG" || log "[WARN] 11 部分失败"

log "=== 主管线全部完成 ==="
