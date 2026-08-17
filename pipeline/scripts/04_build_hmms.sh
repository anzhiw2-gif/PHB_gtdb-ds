#!/bin/bash
# =============================================================================
# 04_build_hmms.sh — 为每个 PHB 降解基因家族构建 HMM profile
#   CD-HIT 95% 去冗余 -> MAFFT 比对 -> hmmbuild
# 用法: bash 04_build_hmms.sh [--cdhit-id 0.95] [--threads 40]
# 输出: data/hmms/{ePhaZ,iPhaZ,OH,BdhA,phasin}.hmm (+ .faa.c95, *_aln.fasta)
# =============================================================================
set -euo pipefail

THREADS=40
CDHIT_ID=0.95
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SEED_DIR="$ROOT/data/seeds/families"
HMM_DIR="$ROOT/data/hmms"
ALN_DIR="$ROOT/data/alignments"
LOG_DIR="$ROOT/results/logs"
mkdir -p "$HMM_DIR" "$ALN_DIR" "$LOG_DIR"

echo "[$(date)] 开始构建 HMM profiles"
for fam in ePhaZ iPhaZ OH BdhA phasin; do
    FA="$SEED_DIR/$fam.faa"
    [ -f "$FA" ] || { echo "  $fam: 种子文件缺失，跳过"; continue; }
    n=$(grep -c '^>' "$FA")
    echo "=== $fam ($n 条种子) ==="
    if [ "$n" -lt 3 ]; then
        echo "  种子过少 (<3)，直接比对建模"
    fi
    # 1) CD-HIT 去冗余
    C95="$HMM_DIR/${fam}.c95.faa"
    cd-hit -i "$FA" -o "$C95" -c "$CDHIT_ID" -n 5 -T "$THREADS" -M 0 \
        > "$LOG_DIR/cdhit_${fam}.log" 2>&1 || { echo "  CD-HIT 失败"; continue; }
    n95=$(grep -c '^>' "$C95")
    echo "  CD-HIT c${CDHIT_ID}: $n -> $n95"
    # 2) MAFFT 比对
    ALN="$ALN_DIR/${fam}_aln.fasta"
    mafft --auto --thread "$THREADS" "$C95" > "$ALN" 2> "$LOG_DIR/mafft_${fam}.log" \
        || { echo "  MAFFT 失败"; continue; }
    echo "  比对完成: $ALN"
    # 3) hmmbuild
    HMM="$HMM_DIR/${fam}.hmm"
    hmmbuild --amino "$HMM" "$ALN" > "$LOG_DIR/hmmbuild_${fam}.log" 2>&1 \
        || { echo "  hmmbuild 失败"; continue; }
    # 模型统计
    ncons=$(grep -c '^CONS' "$HMM" || true)
    echo "  HMM: $HMM"
done

echo "[$(date)] 完成。HMM 文件："
ls -la "$HMM_DIR"/*.hmm 2>/dev/null || true
