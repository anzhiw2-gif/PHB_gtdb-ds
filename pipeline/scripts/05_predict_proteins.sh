#!/bin/bash
# =============================================================================
# 05_predict_proteins.sh — GTDB R232 全基因组蛋白预测（Pyrodigal + GNU parallel）
#   - 只读使用 ~/GTDB/gtdb_genomes_reps_r232 基因组（.fna.gz）
#   - 输出到本工作区 data/proteins/per_genome/{ACC}.faa.gz（可断点续跑）
#   - 合并为 data/proteins/shards/shard_%04d.faa（供 hmmsearch 筛选）
# 用法: bash 05_predict_proteins.sh [--threads 80] [--genomes-per-shard 2000] [--dry-run]
# =============================================================================
set -euo pipefail

THREADS=80
GENOMES_PER_SHARD=2000
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
GTDB_DB="$HOME/GTDB/gtdb_genomes_reps_r232/database"
OUT="$ROOT/data/proteins"
PER_GENOME="$OUT/per_genome"
SHARDS="$OUT/shards"
LOG="$ROOT/results/logs"
GENOME_LIST="$OUT/genome_list.txt"
mkdir -p "$PER_GENOME" "$SHARDS" "$LOG"

DRY=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --genomes-per-shard) GENOMES_PER_SHARD="$2"; shift 2 ;;
        --dry-run) DRY=1; shift ;;
        *) echo "unknown: $1"; exit 1 ;;
    esac
done

echo "[$(date)] 生成基因组清单..."
find "$GTDB_DB" -name "*_genomic.fna.gz" | sort > "$GENOME_LIST"
TOTAL=$(wc -l < "$GENOME_LIST")
echo "  基因组总数: $TOTAL"

if [ "$DRY" -eq 1 ]; then
    echo "[dry-run] 预计运行: $TOTAL 基因组, $THREADS 并行, 约 $((TOTAL * 30 / THREADS / 3600)) 小时 (按 30s/基因组估算)"
    exit 0
fi

predict_one() {
    gz="$1"
    acc=$(basename "$gz" .gz)
    acc="${acc%_genomic.fna}"
    out="$PER_GENOME/${acc}.faa.gz"
    if [ -s "$out" ]; then
        return 0
    fi
    tmp="$PER_GENOME/${acc}.faa.tmp"
    # 蛋白预测（meta 模式，输出蛋白；header 加 ACC 前缀便于溯源）
    if pyrodigal -i "$gz" -a "$tmp" -p meta > /dev/null 2>&1; then
        # 重命名 header: 加 ACC| 前缀
        sed -i "s/^>/>${acc}|/" "$tmp"
        gzip -c "$tmp" > "$out"
        rm -f "$tmp"
    else
        echo "FAIL $gz" >> "$PER_GENOME/failed.log"
    fi
}
export -f predict_one
export PER_GENOME

echo "[$(date)] 开始蛋白预测（$THREADS 并行）..."
cat "$GENOME_LIST" | parallel -j "$THREADS" --progress predict_one {} 2> "$LOG/predict_progress.log" || true

DONE=$(find "$PER_GENOME" -name "*.faa.gz" ! -name "*.tmp*" | wc -l)
FAIL=$(wc -l < "$PER_GENOME/failed.log" 2>/dev/null || echo 0)
echo "[$(date)] 预测完成: $DONE 成功, $FAIL 失败"

echo "[$(date)] 合并分片..."
find "$PER_GENOME" -name "*.faa.gz" | sort > "$OUT/protein_files.txt"
: > "$OUT/shard_index.tsv"
i=0
s=0
while read -r pf; do
    if [ $((i % GENOMES_PER_SHARD)) -eq 0 ]; then
        s=$((s + 1))
        shard="$SHARDS/shard_$(printf '%04d' "$s").faa"
        : > "$shard"
    fi
    zcat "$pf" >> "$shard"
    echo -e "$(printf '%04d' "$s")\t$(basename "$pf" .faa.gz)" >> "$OUT/shard_index.tsv"
    i=$((i + 1))
done < "$OUT/protein_files.txt"
echo "[$(date)] 分片完成: $s 个 shard（每 shard $GENOMES_PER_SHARD 基因组）"
ls -la "$SHARDS" | head -5
echo "  shard 总数: $(ls "$SHARDS" | grep -c '^shard_')"
