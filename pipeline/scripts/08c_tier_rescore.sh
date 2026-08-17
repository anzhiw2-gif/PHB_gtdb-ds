#!/bin/bash
# 08c_tier_rescore.sh — 用 curated 金标准 HMM 对命中做三级重评分
# tier1: curated HMM E<1e-20（严格）
# tier2: curated HMM E<1e-10（中等）
# tier3: 现有 validated（宽模型 + 通用验证）
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate phb_gtdb
cd "$(dirname "$0")/.."
SEQDIR="data/screen/family_seqs"
TIERDIR="data/screen/tiers"
mkdir -p "$TIERDIR"
LOG="results/logs/08c_tier.log"

# 核心家族 → curated HMM 映射
declare -A CURATED=(
  [ePhaZ]=data/hmms/ePhaZ.hmm
  [iPhaZ]=data/hmms/iPhaZ.hmm
  [OH]=data/hmms/OH.hmm
  [ArchPhaZ_patatin]=data/hmms/v2/ArchPhaZ_patatin.hmm
  [ArchPhaZ_hydrolase]=data/hmms/v2/ArchPhaZ_hydrolase.hmm
)

echo "=== 三级重评分 ===" | tee "$LOG"
for fam in ePhaZ iPhaZ OH ArchPhaZ_patatin ArchPhaZ_hydrolase; do
  faa="$SEQDIR/${fam}_validated.faa"
  [ -f "$faa" ] || continue
  hmm="${CURATED[$fam]}"
  n=$(grep -c '^>' "$faa")
  # tier2: E<1e-10
  hmmsearch --tblout "$TIERDIR/${fam}_tier2.tbl" -E 1e-10 --cpu 8 "$hmm" "$faa" > /dev/null 2>&1 || true
  n2=$(grep -vc '^#' "$TIERDIR/${fam}_tier2.tbl" 2>/dev/null || echo 0)
  # 提取 tier2 序列
  grep -v "^#" "$TIERDIR/${fam}_tier2.tbl" | awk '{print $1}' | sort -u > "$TIERDIR/${fam}_tier2.ids"
  python3 -c "
ids = set(open('$TIERDIR/${fam}_tier2.ids').read().split())
hdr=None; buf=[]
for line in open('$faa'):
    if line.startswith('>'):
        if hdr and hdr in ids: print('>'+hdr); print(''.join(buf))
        hdr=line[1:].strip(); buf=[]
    else: buf.append(line.strip())
if hdr and hdr in ids: print('>'+hdr); print(''.join(buf))
" > "$TIERDIR/${fam}_tier2.faa"
  # tier1: E<1e-20（对 tier2 子集再筛，更快）
  hmmsearch --tblout "$TIERDIR/${fam}_tier1.tbl" -E 1e-20 --cpu 8 "$hmm" "$TIERDIR/${fam}_tier2.faa" > /dev/null 2>&1 || true
  n1=$(grep -vc '^#' "$TIERDIR/${fam}_tier1.tbl" 2>/dev/null || echo 0)
  grep -v "^#" "$TIERDIR/${fam}_tier1.tbl" | awk '{print $1}' | sort -u > "$TIERDIR/${fam}_tier1.ids"
  python3 -c "
ids = set(open('$TIERDIR/${fam}_tier1.ids').read().split())
hdr=None; buf=[]
for line in open('$faa'):
    if line.startswith('>'):
        if hdr and hdr in ids: print('>'+hdr); print(''.join(buf))
        hdr=line[1:].strip(); buf=[]
    else: buf.append(line.strip())
if hdr and hdr in ids: print('>'+hdr); print(''.join(buf))
" > "$TIERDIR/${fam}_tier1.faa"
  n3=$n
  echo "$fam: tier1=$n1 tier2=$n2 tier3=$n3" | tee -a "$LOG"
done
