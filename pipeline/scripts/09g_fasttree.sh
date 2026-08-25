#!/bin/bash
# 09g_fasttree.sh — 用 FastTree 快速建树（大基因家族，近似 ML，抽样 1000）
# fail-closed：MAFFT/FastTree 任一失败立即退出。
# 树登记：追加写入 results/trees_tier1/tree_manifest.jsonl（由 09i_tree_manifest.py 固化）。
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate phb_gtdb

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RUN_ROOT="${PHB_RUN_ROOT:-$REPO_ROOT}"
cd "$RUN_ROOT"

TIER="data/screen/tiers"
TREE="results/trees_tier1"
ALN="results/alignments_tier1"
LOG="$RUN_ROOT/logs"
mkdir -p "$TREE" "$ALN" "$LOG"
SAMPLE_N=1000
SEED=42
MANIFEST_JSONL="$TREE/tree_manifest.jsonl"

for fam in ePhaZ iPhaZ OH ArchPhaZ_hydrolase; do
  faa="$TIER/${fam}_tier1.faa"
  [ -f "$faa" ] || { echo "  $fam: 无 tier1 序列，跳过"; continue; }
  n=$(grep -c '^>' "$faa")
  echo "=== $fam ($n) ==="
  sample_faa="$faa"
  sample_list=""
  if [ "$n" -gt "$SAMPLE_N" ]; then
    sample_faa="$ALN/${fam}_s${SAMPLE_N}.faa"
    sample_list="$ALN/${fam}_s${SAMPLE_N}.list"
    python "$SCRIPT_DIR/sample_fasta.py" "$faa" "$sample_faa" "$SAMPLE_N" \
        --seed "$SEED" --list "$sample_list"
  fi

  # MAFFT 比对（加 --maxiterate 提高质量）
  aln="$ALN/${fam}_aln.fasta"
  mafft --auto --maxiterate 1000 --thread 20 "$sample_faa" > "$aln" 2> "$LOG/mafft_${fam}.log" \
    || { echo "[ERROR] MAFFT 失败: $fam"; exit 1; }
  # FastTree 建树
  FastTree -lg -gamma -pseudo "$aln" > "$TREE/${fam}.fasttree.nwk" 2> "$LOG/fasttree_${fam}.log" \
    || { echo "[ERROR] FastTree 失败: $fam"; exit 1; }

  treefile="$TREE/${fam}.fasttree.nwk"
  n_leaves=$(grep -o '[A-Za-z][A-Za-z0-9_.|]*:' "$treefile" | wc -l)
  in_sha=$(sha256sum "$sample_faa" | awk '{print $1}')

  printf '{"family":"%s","tool":"FastTree (LG+gamma, 抽样 %d)","tree_file":"%s","n_leaves":%d,"input_faa":"%s","input_sha256":"%s","sample_list":"%s","sample_seed":%d,"sample_from":%d,"alignment":"%s","trim":"","mafft_log":"%s","trimal_log":"","tree_log":"%s"}\n' \
    "$fam" "$SAMPLE_N" "$treefile" "$n_leaves" "$sample_faa" "$in_sha" "$sample_list" "$SEED" "$n" \
    "$aln" "$LOG/mafft_${fam}.log" "$LOG/fasttree_${fam}.log" \
    >> "$MANIFEST_JSONL"
  echo "  树: $treefile（叶数≈$n_leaves）"
done
echo "完成"
ls -la "$TREE"/*.nwk 2>/dev/null || true
