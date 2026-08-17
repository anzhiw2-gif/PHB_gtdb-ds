#!/bin/bash
# =============================================================================
# 09_phylogeny.sh — 各家族系统发育（MAFFT -> trimAl -> IQ-TREE2）
#   输入: data/screen/family_seqs/*_validated.faa（或 *.faa）
#   输出: results/trees/{fam}.treefile + 比对/修剪文件
# 用法: bash 09_phylogeny.sh [--threads 40] [--families "ePhaZ iPhaZ OH BdhA"]
# =============================================================================
set -euo pipefail

# 激活 conda 环境（mafft/trimal/iqtree/cd-hit 依赖）
source ~/miniconda3/etc/profile.d/conda.sh
conda activate phb_gtdb

THREADS=40
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SEQ_DIR="$ROOT/data/screen/family_seqs"
TREE_DIR="$ROOT/results/trees"
ALN_DIR="$ROOT/results/alignments"
LOG="$ROOT/results/logs"
mkdir -p "$TREE_DIR" "$ALN_DIR" "$LOG"

FAMILIES="ePhaZ iPhaZ OH BdhA ArchPhaZ_patatin ArchPhaZ_hydrolase"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --families) FAMILIES="$2"; shift 2 ;;
        *) echo "unknown: $1"; exit 1 ;;
    esac
done

for fam in $FAMILIES; do
    faa="$SEQ_DIR/${fam}_validated.faa"
    [ -f "$faa" ] || faa="$SEQ_DIR/${fam}.faa"
    [ -f "$faa" ] || { echo "  $fam: 无序列文件，跳过"; continue; }
    n=$(grep -c '^>' "$faa")
    echo "=== $fam ($n seqs) ==="
    if [ "$n" -lt 4 ]; then
        echo "  序列过少，跳过建树"
        continue
    fi
    # 大家族抽样（>2000 条时直接随机抽样，固定种子可复现）
    if [ "$n" -gt 2000 ]; then
        faa_sample="$SEQ_DIR/${fam}_sample2000.faa"
        python3 -c "
import random
random.seed(42)
seqs = {}
hdr = None; buf = []
for line in open('$faa'):
    if line.startswith('>'):
        if hdr: seqs[hdr] = ''.join(buf)
        hdr = line[1:].strip(); buf = []
    else: buf.append(line.strip())
if hdr: seqs[hdr] = ''.join(buf)
keys = list(seqs.keys())
random.shuffle(keys)
with open('$faa_sample', 'w') as f:
    for k in keys[:2000]:
        f.write('>' + k + '\n' + seqs[k] + '\n')
print('sampled 2000 from', len(keys))
"
        faa="$faa_sample"
        n=$(grep -c '^>' "$faa")
        echo "  随机抽样 2000 条用于建树"
    fi
    aln="$ALN_DIR/${fam}_aln.fasta"
    mafft --auto --thread "$THREADS" "$faa" > "$aln" 2> "$LOG/mafft_${fam}_phylo.log" || { echo "  MAFFT 失败"; continue; }
    trm="$ALN_DIR/${fam}_trim.fasta"
    trimal -in "$aln" -out "$trm" -automated1 2> "$LOG/trimal_${fam}.log" || trm="$aln"
    # IQ-TREE2（自动模型选择 + UFBoot，大族用较宽松设置）
    iqtree2 -s "$trm" -m LG+G4 -B 1000 -T "$THREADS" \
        --prefix "$TREE_DIR/${fam}" > "$LOG/iqtree_${fam}.log" 2>&1 \
        || iqtree -s "$trm" -m LG+G4 -bb 1000 -T "$THREADS" \
           --prefix "$TREE_DIR/${fam}" >> "$LOG/iqtree_${fam}.log" 2>&1 \
        || echo "  IQ-TREE 失败（两版均不可用）"
    echo "  树: $TREE_DIR/${fam}.treefile"
done
echo "[$(date)] 完成"
ls -la "$TREE_DIR"/*.treefile 2>/dev/null || true
