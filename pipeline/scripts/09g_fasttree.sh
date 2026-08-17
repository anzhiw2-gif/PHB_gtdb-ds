#!/bin/bash
# 09g_fasttree.sh — 用 FastTree 快速建树（大基因家族，近似 ML）
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate phb_gtdb
cd "$(dirname "$0")/.."
TIER="data/screen/tiers"
TREE="results/trees_tier1"
ALN="results/alignments_tier1"
LOG="results/logs"
mkdir -p "$TREE" "$ALN" "$LOG"

for fam in ePhaZ iPhaZ OH ArchPhaZ_hydrolase; do
  faa="$TIER/${fam}_tier1.faa"
  [ -f "$faa" ] || continue
  n=$(grep -c '^>' "$faa")
  echo "=== $fam ($n) ==="
  # 抽样到 1000（保证比对质量）
  if [ "$n" -gt 1000 ]; then
    python3 -c "
import random
random.seed(42)
seqs={}; hdr=None; buf=[]
for line in open('$faa'):
    if line.startswith('>'):
        if hdr: seqs[hdr]=''.join(buf)
        hdr=line[1:].strip(); buf=[]
    else: buf.append(line.strip())
if hdr: seqs[hdr]=''.join(buf)
keys=list(seqs.keys()); random.shuffle(keys)
with open('$ALN/${fam}_s1000.faa','w') as f:
    for k in keys[:1000]: f.write('>'+k+'\n'+seqs[k]+'\n')
"
    faa="$ALN/${fam}_s1000.faa"
  fi
  # MAFFT 比对（加 --maxiterate 提高质量）
  aln="$ALN/${fam}_aln.fasta"
  mafft --auto --maxiterate 1000 --thread 20 "$faa" > "$aln" 2> "$LOG/mafft_${fam}.log" || { echo "MAFFT 失败"; continue; }
  # FastTree 建树
  FastTree -lg -gamma -pseudo "$aln" > "$TREE/${fam}.fasttree.nwk" 2> "$LOG/fasttree_${fam}.log" || { echo "FastTree 失败"; continue; }
  echo "  树: $TREE/${fam}.fasttree.nwk"
done
echo "完成"
ls -la "$TREE"/*.nwk 2>/dev/null || true
