#!/bin/bash
# Build family alignments and trees inside the selected run directory.
set -euo pipefail

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate phb_gtdb

THREADS=40
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RUN_ROOT="${PHB_RUN_ROOT:-$REPO_ROOT}"
cd "$RUN_ROOT"
SEQ_DIR="$RUN_ROOT/data/screen/family_seqs"
TREE_DIR="$RUN_ROOT/results/trees"
ALN_DIR="$RUN_ROOT/results/alignments"
LOG="$RUN_ROOT/logs"
mkdir -p "$TREE_DIR" "$ALN_DIR" "$LOG"

FAMILIES="ePhaZ iPhaZ OH BdhA ArchPhaZ_patatin ArchPhaZ_hydrolase"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) [[ $# -ge 2 ]] || { echo "--threads requires a value" >&2; exit 1; }; THREADS="$2"; shift 2 ;;
        --families) [[ $# -ge 2 ]] || { echo "--families requires a value" >&2; exit 1; }; FAMILIES="$2"; shift 2 ;;
        *) echo "unknown argument: $1" >&2; exit 1 ;;
    esac
done

for fam in $FAMILIES; do
    faa="$SEQ_DIR/${fam}_validated.faa"
    [[ -f "$faa" ]] || faa="$SEQ_DIR/${fam}.faa"
    [[ -f "$faa" ]] || { echo "$fam: sequence file missing; skipped" >&2; continue; }
    n=$(grep -c '^>' "$faa")
    echo "=== $fam ($n sequences) ==="
    if [[ "$n" -lt 4 ]]; then
        echo "$fam: fewer than four sequences; skipped" >&2
        continue
    fi
    if [[ "$n" -gt 2000 ]]; then
        faa_sample="$SEQ_DIR/${fam}_sample2000.faa"
        python3 - "$faa" "$faa_sample" <<'PY'
import random
import sys

source, target = sys.argv[1:]
records = {}
header = None
sequence = []
with open(source, encoding="utf-8") as handle:
    for line in handle:
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(sequence)
            header = line[1:].strip()
            sequence = []
        else:
            sequence.append(line.strip())
if header is not None:
    records[header] = "".join(sequence)
keys = list(records)
random.Random(42).shuffle(keys)
with open(target, "w", encoding="utf-8", newline="\n") as handle:
    for key in keys[:2000]:
        handle.write(f">{key}\n{records[key]}\n")
PY
        faa="$faa_sample"
    fi
    aln="$ALN_DIR/${fam}_aln.fasta"
    if ! mafft --auto --thread "$THREADS" "$faa" > "$aln" 2> "$LOG/mafft_${fam}_phylo.log"; then
        echo "$fam: MAFFT failed" >&2
        continue
    fi
    trimmed="$ALN_DIR/${fam}_trim.fasta"
    trimal -in "$aln" -out "$trimmed" -automated1 2> "$LOG/trimal_${fam}.log" || trimmed="$aln"
    if ! iqtree2 -s "$trimmed" -m LG+G4 -B 1000 -T "$THREADS" \
        --prefix "$TREE_DIR/${fam}" > "$LOG/iqtree_${fam}.log" 2>&1; then
        iqtree -s "$trimmed" -m LG+G4 -bb 1000 -T "$THREADS" \
            --prefix "$TREE_DIR/${fam}" >> "$LOG/iqtree_${fam}.log" 2>&1 || {
            echo "$fam: IQ-TREE failed" >&2
            continue
        }
    fi
    echo "$fam: tree written to $TREE_DIR/${fam}.treefile"
done

ls -la "$TREE_DIR"/*.treefile 2>/dev/null || true
