#!/bin/bash
# test_smoke.sh — 小样本端到端冒烟测试入口（无需 HMMER，纯 Python）
# 用法（服务器或本地均可）: bash pipeline/tests/test_smoke.sh
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$ROOT"
python "$SCRIPT_DIR/test_smoke.py"
