#!/bin/bash
# 09b_tier1_phylogeny.sh — tier1 严格集建树（大族抽样 2000，seed=42）
# fail-closed：MAFFT/trimAl/IQ-TREE 任一失败立即退出（不再 continue/降级）。
# 每棵树登记：工具、输入 FASTA SHA-256、抽样名单、实际叶数、比对/修剪/建树日志路径，
# 追加写入 results/trees_tier1/tree_manifest.jsonl（最终由 09i_tree_manifest.py 固化）。
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
THREADS=30
SAMPLE_N=2000
SEED=42
MANIFEST_JSONL="$TREE/tree_manifest.jsonl"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        *) echo "unknown: $1"; exit 1 ;;
    esac
done

for fam in ePhaZ iPhaZ OH ArchPhaZ_hydrolase; do
  faa="$TIER/${fam}_tier1.faa"
  [ -f "$faa" ] || { echo "  $fam: 无 tier1 序列，跳过"; continue; }
  n=$(grep -c '^>' "$faa")
  echo "=== $fam ($n) ==="
  sample_faa="$faa"
  sample_list=""
  if [ "$n" -gt "$SAMPLE_N" ]; then
    sample_faa="$ALN/${fam}_sample${SAMPLE_N}.faa"
    sample_list="$ALN/${fam}_sample${SAMPLE_N}.list"
    python "$SCRIPT_DIR/sample_fasta.py" "$faa" "$sample_faa" "$SAMPLE_N" \
        --seed "$SEED" --list "$sample_list"
  fi

  aln="$ALN/${fam}_aln.fasta"
  mafft --auto --thread "$THREADS" "$sample_faa" > "$aln" 2> "$LOG/mafft_tier1_${fam}.log" \
    || { echo "[ERROR] MAFFT 失败: $fam"; exit 1; }
  trm="$ALN/${fam}_trim.fasta"
  trimal -in "$aln" -out "$trm" -automated1 > "$LOG/trimal_tier1_${fam}.log" 2>&1 \
    || { echo "[ERROR] trimAl 失败: $fam"; exit 1; }

  prefix="$TREE/${fam}"
  iqtree -s "$trm" -m LG+G4 -bb 1000 -T "$THREADS" --prefix "$prefix" \
      > "$LOG/iqtree_tier1_${fam}.log" 2>&1 \
    || { echo "[ERROR] IQ-TREE 失败: $fam"; exit 1; }

  treefile="$prefix.treefile"
  [ -f "$treefile" ] || { echo "[ERROR] 未生成 treefile: $fam"; exit 1; }
  n_leaves=$(grep -o '[A-Za-z0-9_.|]\+:' "$treefile" | wc -l)
  in_sha=$(sha256sum "$sample_faa" | awk '{print $1}')

  # 登记到 tree_manifest.jsonl
  printf '{"family":"%s","tool":"IQ-TREE2 (LG+G4, 1000 UFBoot)","tree_file":"%s","n_leaves":%d,"input_faa":"%s","input_sha256":"%s","sample_list":"%s","sample_seed":%d,"sample_from":%d,"alignment":"%s","trim":"%s","mafft_log":"%s","trimal_log":"%s","tree_log":"%s"}\n' \
    "$fam" "$treefile" "$n_leaves" "$sample_faa" "$in_sha" "$sample_list" "$SEED" "$n" \
    "$aln" "$trm" "$LOG/mafft_tier1_${fam}.log" "$LOG/trimal_tier1_${fam}.log" "$LOG/iqtree_tier1_${fam}.log" \
    >> "$MANIFEST_JSONL"
  echo "  树: $treefile（叶数≈$n_leaves）"
done

echo "完成"
ls -la "$TREE"/*.treefile 2>/dev/null || true
echo "树登记: $MANIFEST_JSONL"
