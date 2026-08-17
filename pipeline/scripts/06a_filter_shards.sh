#!/bin/bash
# 06a_filter_shards.sh — 过滤 shards 中超长序列（>100K aa，规避 HMMER 硬限制）
# 输出到 data/proteins/shards_filt/（不改原 shards）
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
THREADS=30
mkdir -p data/proteins/shards_filt results/logs

echo "[$(date)] 过滤超长序列（>100K aa）..."
source ~/miniconda3/etc/profile.d/conda.sh

filter_one() {
    shard="$1"
    base=$(basename "$shard")
    out="data/proteins/shards_filt/$base"
    if [ -s "$out" ]; then return 0; fi
    conda run -n phb_gtdb python scripts/filter_long_seqs.py "$shard" "$out" 100000 \
        >> results/logs/filter_shards.log 2>&1
}
export -f filter_one
export ROOT

ls data/proteins/shards/shard_*.faa | sort | parallel -j "$THREADS" filter_one {} 
echo "[$(date)] 过滤完成"
echo "过滤后 shard 数: $(ls data/proteins/shards_filt/shard_*.faa 2>/dev/null | wc -l)"
du -sh data/proteins/shards_filt/ 2>/dev/null
tail -5 results/logs/filter_shards.log 2>/dev/null
