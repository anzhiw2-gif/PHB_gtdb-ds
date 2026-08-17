#!/bin/bash
# =============================================================================
# 06_screen.sh — 用家族 HMM 对蛋白分片全库筛选（hmmsearch）+ 家族分类
#   - 输入: data/proteins/shards/shard_*.faa  （由 05 生成）
#   - HMM: data/hmms/v2/*.hmm（核心降解家族 + 辅助家族）
#   - 输出: data/screen/hmmsearch/{family}__shard_XXXX.tbl
#           data/screen/hits_all.tsv（含 family 标签 = 家族分类）
# 用法: bash 06_screen.sh [--threads 80] [--eval 1e-5]
#       [--families "ePhaZ iPhaZ OH BdhA ArchPhaZ_patatin ArchPhaZ_hydrolase"]
# =============================================================================
set -euo pipefail

# 激活 conda 环境（hmmsearch/parallel 依赖）
source ~/miniconda3/etc/profile.d/conda.sh
conda activate phb_gtdb

THREADS=80
EVAL=1e-5
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SHARDS="$ROOT/data/proteins/shards_filt"
HMM_DIR="$ROOT/data/hmms/v2"
SCREEN_DIR="$ROOT/data/screen"
HMMOUT="$SCREEN_DIR/hmmsearch"
LOG="$ROOT/results/logs"
mkdir -p "$HMMOUT" "$LOG"

# 核心降解家族（家族分类目标）+ 辅助（PhaJ/phasin/PhaC 供基因簇背景）
FAMILIES="ePhaZ iPhaZ OH BdhA ArchPhaZ_patatin ArchPhaZ_hydrolase"
AUX_FAMILIES="PhaJ phasin PhaC"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --eval) EVAL="$2"; shift 2 ;;
        --families) FAMILIES="$2"; shift 2 ;;
        *) echo "unknown: $1"; exit 1 ;;
    esac
done

mapfile -t SHARD_FILES < <(ls "$SHARDS"/shard_*.faa 2>/dev/null | sort)
echo "[$(date)] shard 数: ${#SHARD_FILES[@]}, 核心家族: $FAMILIES, 辅助: $AUX_FAMILIES, eval=$EVAL"

if [ "${#SHARD_FILES[@]}" -eq 0 ]; then
    echo "未找到 shard，请先运行 05_predict_proteins.sh"
    exit 1
fi

run_one() {
    hmm="$1"
    fam=$(basename "$hmm" .hmm)
    shard="$2"
    sname=$(basename "$shard" .faa)
    out="$HMMOUT/${fam}__${sname}.tbl"
    if [ -s "$out" ]; then
        return 0
    fi
    hmmsearch --tblout "$out" --domtblout "${out%.tbl}.dom" \
        -E "$EVAL" --cpu 1 "$hmm" "$shard" > /dev/null 2>&1 || true
}
export -f run_one
export HMMOUT EVAL

ALL_FAMILIES="$FAMILIES $AUX_FAMILIES"
echo "[$(date)] 开始 hmmsearch（全部家族: $ALL_FAMILIES）..."
for fam in $ALL_FAMILIES; do
    hmm="$HMM_DIR/${fam}.hmm"
    if [ ! -f "$hmm" ]; then echo "  $fam HMM 缺失，跳过"; continue; fi
    echo "  family: $fam (${#SHARD_FILES[@]} shards)"
    printf '%s\n' "${SHARD_FILES[@]}" | parallel -j "$THREADS" run_one "$hmm" {} 2> "$LOG/screen_${fam}.log" || true
done

echo "[$(date)] 汇总命中（含家族分类标签）..."
: > "$SCREEN_DIR/hits_all.tsv"
# tblout 列: $1 target, $2 tacc, $3 qname, $5 E-value, $6 score, $7 bias, $8 best-dom-E, $9 best-dom-score
echo -e "family\tshard\tprotein\ttacc\tE-value\tscore\tbias\tdomE\tqname" > "$SCREEN_DIR/hits_all.tsv"
for f in "$HMMOUT"/*.tbl; do
    [ -f "$f" ] || continue
    fam=$(basename "$f" .tbl); fam=${fam%%__*}
    sname=$(basename "$f" .tbl); sname=${sname#*__}
    grep -v "^#" "$f" | awk -v F="$fam" -v S="$sname" '{print F"\t"S"\t"$1"\t"$2"\t"$5"\t"$6"\t"$7"\t"$8"\t"$3}' \
        >> "$SCREEN_DIR/hits_all.tsv"
done
echo "[$(date)] 完成"
n=$(wc -l < "$SCREEN_DIR/hits_all.tsv")
echo "  总命中行数(含表头): $n"
echo "  按家族统计:"
cut -f1 "$SCREEN_DIR/hits_all.tsv" | tail -n +2 | sort | uniq -c
