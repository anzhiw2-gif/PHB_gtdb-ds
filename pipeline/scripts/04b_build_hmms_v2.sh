#!/bin/bash
# =============================================================================
# 04b_build_hmms_v2.sh — 基于全面种子库（v2）构建家族 HMM
#   ePhaZ/iPhaZ/OH/BdhA/PhaJ/PhaC/ArchPhaZ_patatin: CD-HIT c90 -> MAFFT -> hmmbuild
#   phasin: 用 Pfam PF09361 (phasin_2) 直接提取（种子过多且为辅助家族）
# 输出: data/hmms/v2/*.hmm
# 用法: bash 04b_build_hmms_v2.sh [--threads 40] [--cdhit-id 0.90]
# =============================================================================
set -euo pipefail

THREADS=40
CDHIT_ID=0.90
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"   # 仓库根（scripts 位于 pipeline/scripts/）
SEED_DIR="$ROOT/data/seeds/v2"
HMM_DIR="$ROOT/data/hmms/v2"
ALN_DIR="$ROOT/data/alignments/v2"
LOG="$ROOT/results/logs"
PFAM="$HOME/GTDB/pfam/Pfam-A.hmm"
mkdir -p "$HMM_DIR" "$ALN_DIR" "$LOG"

echo "[$(date)] 构建 v2 家族 HMM（CD-HIT c$CDHIT_ID）"
for fam in ePhaZ iPhaZ OH BdhA PhaJ PhaC ArchPhaZ_patatin; do
    FA="$SEED_DIR/$fam.faa"
    [ -f "$FA" ] || { echo "  $fam: 无种子，跳过"; continue; }
    n=$(grep -c '^>' "$FA")
    echo "=== $fam ($n) ==="
    C90="$HMM_DIR/${fam}.c90.faa"
    cd-hit -i "$FA" -o "$C90" -c "$CDHIT_ID" -n 5 -T "$THREADS" -M 0 \
        > "$LOG/cdhit_v2_${fam}.log" 2>&1 || { echo "  CD-HIT 失败"; continue; }
    n90=$(grep -c '^>' "$C90")
    echo "  c$CDHIT_ID: $n -> $n90"
    if [ "$n90" -lt 3 ]; then
        echo "  种子过少，跳过"; continue
    fi
    # 超大族抽样（>1500 用 c80 进一步去冗余）
    if [ "$n90" -gt 1500 ]; then
        C2="$HMM_DIR/${fam}.c80.faa"
        cd-hit -i "$C90" -o "$C2" -c 0.80 -n 5 -T "$THREADS" -M 0 \
            > "$LOG/cdhit_v2_${fam}_c80.log" 2>&1 || true
        C90="$C2"
        n90=$(grep -c '^>' "$C90")
        echo "  c80 二次去冗余: $n90"
    fi
    ALN="$ALN_DIR/${fam}_aln.fasta"
    mafft --auto --thread "$THREADS" "$C90" > "$ALN" 2> "$LOG/mafft_v2_${fam}.log" \
        || { echo "  MAFFT 失败"; continue; }
    HMM="$HMM_DIR/${fam}.hmm"
    hmmbuild --amino "$HMM" "$ALN" > "$LOG/hmmbuild_v2_${fam}.log" 2>&1 || { echo "  hmmbuild 失败"; continue; }
    echo "  -> $HMM"
done

# phasin: 用 Pfam PF09361（phasin_2）——若 Pfam-A 可用
echo "=== phasin (Pfam PF09361) ==="
if [ -f "$PFAM" ]; then
    hmmfetch "$PFAM" PF09361 > "$HMM_DIR/phasin.hmm" 2>/dev/null \
        && echo "  PF09361 提取成功 -> data/hmms/v2/phasin.hmm" \
        || { echo "  PF09361 提取失败，尝试按名称 phasin_2"; hmmfetch "$PFAM" phasin_2 > "$HMM_DIR/phasin.hmm" 2>/dev/null && echo "  ok by name"; }
else
    echo "  Pfam-A.hmm 不存在，跳过 phasin"
fi

echo "[$(date)] 完成。v2 HMM 列表:"
ls -la "$HMM_DIR"/*.hmm 2>/dev/null || true
