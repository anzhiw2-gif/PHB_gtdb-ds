#!/bin/bash
# =============================================================================
# sync_from_server.sh — 从服务器 T141 同步结果与交付物到本地
#
# 同步内容:
#   results/tables/   results/figures/   results/trees/
#   data/hmms/        data/screen/tiers/   （HMM profiles + tier1 序列 = 交付物）
#
# 用法:
#   bash pipeline/sync_from_server.sh                # 实际同步
#   bash pipeline/sync_from_server.sh --dry-run      # 预览
#   SYNC_SERVER=user@host bash pipeline/sync_from_server.sh
#
# 说明: data/ 在 .gitignore 中不入 git；HMM 与 tier1 序列属发布物，
#       建议附为 GitHub Release / Zenodo DOI。
# =============================================================================
set -euo pipefail

: "${SYNC_SERVER:?set SYNC_SERVER to user@host}"
: "${SYNC_REMOTE:?set SYNC_REMOTE to the remote PHB_gtdb-ds path}"
SERVER="$SYNC_SERVER"
REMOTE="$SYNC_REMOTE"

DRY=""
for a in "$@"; do
    case "$a" in
        --dry-run|-n) DRY="--dry-run" ;;
        *) echo "unknown: $a"; exit 1 ;;
    esac
done

RSYNC=(rsync -avh --progress $DRY)

echo "[sync] 服务器: $SERVER  (remote=$REMOTE)"
echo "[sync] dry-run=${DRY:-no}"

mkdir -p results/tables results/figures results/trees
"${RSYNC[@]}" "$SERVER:$REMOTE/results/tables/"  ./results/tables/
"${RSYNC[@]}" "$SERVER:$REMOTE/results/figures/" ./results/figures/
"${RSYNC[@]}" "$SERVER:$REMOTE/results/trees/"   ./results/trees/

mkdir -p data/hmms data/screen/tiers
"${RSYNC[@]}" "$SERVER:$REMOTE/data/hmms/"         ./data/hmms/
"${RSYNC[@]}" "$SERVER:$REMOTE/data/screen/tiers/" ./data/screen/tiers/

echo "[sync] 完成。核对:"
echo "  results/tables:  $(ls results/tables 2>/dev/null | wc -l) 个文件"
echo "  results/figures: $(ls results/figures 2>/dev/null | wc -l) 个文件"
echo "  data/hmms:       $(ls data/hmms 2>/dev/null | wc -l) 个文件"
echo "  data/screen/tiers: $(ls data/screen/tiers 2>/dev/null | wc -l) 个文件"
