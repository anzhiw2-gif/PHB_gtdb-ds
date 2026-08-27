#!/bin/bash
# cdhit_tree.sh — 大族 CD-HIT 去冗余后建树（完整树的近似，替代全量 IQ-TREE）
#   输入 tier1.faa（如 ePhaZ 38,692 / iPhaZ 32,926），CD-HIT 去冗余 → MAFFT → IQ-TREE2。
# 用法: bash cdhit_tree.sh FAMILY [--cid 0.8] [--threads 40]
#   例: bash cdhit_tree.sh iPhaZ --cid 0.8 --threads 40
# 输出: results/trees_tier1/{fam}.cdhit.treefile + results/alignments_tier1/{fam}_cdhit_aln.fasta
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate phb_gtdb

# 服务器执行脚本：工作区根硬编码（本地 git 是 pipeline/scripts/ 结构，服务器是扁平 scripts/ 结构；
# 见 docs/STATUS.md 路径契约约定。本地运行的路径契约脚本用 SCRIPT_DIR/../..。）
ROOT="/home/data/haoyu/PHB_gtdb-ds"
cd "$ROOT"

FAM="$1"; shift || { echo "用法: bash cdhit_tree.sh FAMILY [--cid 0.8] [--threads 40]"; exit 1; }
CID=0.80
THREADS=40
while [[ $# -gt 0 ]]; do
    case "$1" in
        --cid) CID="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        *) echo "unknown: $1"; exit 1 ;;
    esac
done

TIER="data/screen/tiers"
ALN="results/alignments_tier1"
TREE="results/trees_tier1"
LOG="results/logs"
mkdir -p "$ALN" "$TREE" "$LOG"

FAA="$TIER/${FAM}_tier1.faa"
[ -f "$FAA" ] || { echo "[ERROR] 无 tier1 序列: $FAA"; exit 1; }
n=$(grep -c '^>' "$FAA")
echo "=== $FAM: $n 条 tier1，CD-HIT c=$CID 去冗余 ==="

C90="$ALN/${FAM}_cdhit.faa"
cd-hit -i "$FAA" -o "$C90" -c "$CID" -n 5 -T "$THREADS" -M 0 -d 0 \
    > "$LOG/cdhit_${FAM}.log" 2>&1 || { echo "[ERROR] cd-hit 失败"; exit 1; }
nc=$(grep -c '^>' "$C90")
echo "  cd-hit: $n -> $nc 条"

ALNF="$ALN/${FAM}_cdhit_aln.fasta"
mafft --auto --thread "$THREADS" "$C90" > "$ALNF" 2> "$LOG/mafft_cdhit_${FAM}.log" \
    || { echo "[ERROR] MAFFT 失败"; exit 1; }
echo "  MAFFT 完成: $ALNF"

iqtree -s "$ALNF" -m LG+G4 -bb 1000 -T "$THREADS" \
    --prefix "$TREE/${FAM}.cdhit" > "$LOG/iqtree_cdhit_${FAM}.log" 2>&1 \
    || { echo "[ERROR] IQ-TREE 失败"; exit 1; }

treefile="$TREE/${FAM}.cdhit.treefile"
n_leaves=$(grep -o '[A-Za-z][A-Za-z0-9_.|]*:' "$treefile" | wc -l)
echo "  树: $treefile（叶数≈$n_leaves）"

# 登记到 tree_manifest.jsonl（由 09i 固化；此处仅记录 cdhit 去冗余信息）
printf '{"family":"%s","tool":"IQ-TREE2 (CD-HIT c=%s 去冗余)","tree_file":"%s","n_leaves":%d,"note":"cdhit 去冗余后建树","cdhit_from":%d,"cdhit_to":%d}\n' \
    "$FAM" "$CID" "$treefile" "$n_leaves" "$n" "$nc" >> "$TREE/tree_manifest.jsonl"
echo "完成"
