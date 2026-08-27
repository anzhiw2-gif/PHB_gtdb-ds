#!/bin/bash
# signalp_ePhaZ.sh — 对 ePhaZ tier1 序列跑 SignalP6，统计胞外(SP) vs 胞内(OTHER)
# 输出: results/signalp/ePhaZ/prediction_results.txt + 统计摘要
# 用法: bash signalp_ePhaZ.sh [MAX]   # MAX>0 时只取前 MAX 条（测速用）；0=全部
set -euo pipefail
cd /home/data/haoyu/PHB_gtdb-ds
SIG=~/miniconda3/envs/signalp6/bin/signalp6
PY=~/miniconda3/envs/phb_gtdb/bin/python
IN=data/screen/tiers/ePhaZ_tier1.faa
MAX=${1:-0}
OUT=results/signalp/ePhaZ
mkdir -p results/signalp

PROBE="$IN"
if [ "$MAX" -gt 0 ]; then
    PROBE=results/signalp/ePhaZ_probe.faa
    "$PY" - "$IN" "$PROBE" "$MAX" <<'PYEOF'
import sys
inp, outp, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
seqs = []
hdr = None
buf = []
for line in open(inp):
    line = line.rstrip("\n")
    if line.startswith(">"):
        if hdr is not None:
            seqs.append((hdr, "".join(buf)))
            if len(seqs) >= n:
                break
        hdr = line[1:]
        buf = []
    else:
        buf.append(line.strip())
if hdr is not None and len(seqs) < n:
    seqs.append((hdr, "".join(buf)))
with open(outp, "w") as fo:
    for h, s in seqs:
        fo.write(">" + h + "\n" + s + "\n")
print(f"probe: {len(seqs)} 条 -> {outp}")
PYEOF
fi

echo "[$(date)] SignalP6 开始: $PROBE"
"$SIG" --fastafile "$PROBE" --organism other --format txt --output_dir "$OUT" --mode fast
echo "[$(date)] SignalP6 完成"

# 统计摘要
"$PY" - "$OUT/prediction_results.txt" <<'PYEOF'
import sys
p = sys.argv[1]
n = sp = other = 0
for line in open(p):
    if line.startswith("#") or not line.strip():
        continue
    c = line.split()
    if len(c) < 3:
        continue
    n += 1
    if c[1] == "SP":
        sp += 1
    else:
        other += 1
print(f"总序列={n} 分泌型(SP/胞外)={sp} 非分泌(OTHER/胞内)={other} "
      f"SP比例={sp/n:.3f}" if n else "无结果")
PYEOF
