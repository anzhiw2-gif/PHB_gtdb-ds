#!/bin/bash
# 09b_tier1_phylogeny.sh — tier1 严格集建树（大族抽样 2000）
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate phb_gtdb
cd "$(dirname "$0")/.."
TIER="data/screen/tiers"
TREE="results/trees_tier1"
ALN="results/alignments_tier1"
LOG="results/logs"
mkdir -p "$TREE" "$ALN" "$LOG"
THREADS=30

for fam in ePhaZ iPhaZ OH ArchPhaZ_hydrolase; do
  faa="$TIER/${fam}_tier1.faa"
  [ -f "$faa" ] || continue
  n=$(grep -c '^>' "$faa")
  echo "=== $fam ($n) ==="
  if [ "$n" -gt 2000 ]; then
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
with open('$ALN/${fam}_sample2000.faa','w') as f:
    for k in keys[:2000]: f.write('>'+k+'\n'+seqs[k]+'\n')
print('sampled 2000 from', len(keys))
"
    faa="$ALN/${fam}_sample2000.faa"
  fi
  aln="$ALN/${fam}_aln.fasta"
  mafft --auto --thread $THREADS "$faa" > "$aln" 2> "$LOG/mafft_tier1_${fam}.log" || { echo "MAFFT 失败"; continue; }
  trm="$ALN/${fam}_trim.fasta"
  trimal -in "$aln" -out "$trm" -automated1 2>/dev/null || trm="$aln"
  iqtree -s "$trm" -m LG+G4 -bb 1000 -T $THREADS --prefix "$TREE/${fam}" > "$LOG/iqtree_tier1_${fam}.log" 2>&1 \
    || echo "IQ-TREE 失败"
  echo "  树: $TREE/${fam}.treefile"
done
echo "完成"
ls -la "$TREE"/*.treefile 2>/dev/null || true
