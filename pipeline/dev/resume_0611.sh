#!/bin/bash
# resume_0611.sh — 【已废弃】从筛选阶段续跑的旧编排。
# 原脚本把 pipeline/ 当项目根、且步骤 11 误传 genome_hits.tsv（无 locus），
# 已由 pipeline/scripts/run_pipeline.sh（fail-closed + run_manifest + 06a 接线）取代。
# 请改用: nohup bash pipeline/scripts/run_pipeline.sh > results/logs/pipeline_master.log 2>&1 &
set -uo pipefail
echo "[resume_0611.sh 已废弃] 请使用 pipeline/scripts/run_pipeline.sh（见 docs/STATUS.md §3）"
exit 1
