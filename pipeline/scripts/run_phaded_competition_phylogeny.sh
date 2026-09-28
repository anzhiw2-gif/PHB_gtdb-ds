#!/usr/bin/env bash
# Run the bounded four-panel PhaDED competition phylogeny.
# This script is copied into a dated deploy/ directory before server execution.
set -euo pipefail

INPUT_DIR=""
OUTPUT_DIR=""
THREADS=40
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir) INPUT_DIR="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
[[ -n "$INPUT_DIR" && -n "$OUTPUT_DIR" ]] || { echo "--input-dir and --output-dir are required" >&2; exit 2; }
[[ "$THREADS" =~ ^[0-9]+$ && "$THREADS" -ge 1 && "$THREADS" -le 40 ]] || { echo "threads must be between 1 and 40" >&2; exit 2; }
command -v mafft >/dev/null || { echo "mafft not found" >&2; exit 127; }
command -v trimal >/dev/null || { echo "trimal not found" >&2; exit 127; }
command -v iqtree >/dev/null || { echo "iqtree not found" >&2; exit 127; }

mkdir -p "$OUTPUT_DIR/alignments" "$OUTPUT_DIR/trees" "$OUTPUT_DIR/logs"
manifest="$OUTPUT_DIR/phylogeny_manifest.tsv"
printf 'panel\tinput_fasta\tinput_sha256\tinput_leaves\talignment\talignment_sha256\ttrimmed_alignment\ttrimmed_sha256\ttree\ttree_sha256\tfinal_leaves\tthreads\tmafft_version\ttrimal_version\tiqtree_version\tstatus\n' > "$manifest"

sha256() { sha256sum "$1" | awk '{print $1}'; }
count_fasta() { grep -c '^>' "$1"; }
tool_version() { "$1" --version 2>&1 | head -n 1 | tr '\t' ' '; }

for panel in extracellular_explicit phaC_vs_phaZ structure_anomaly ephaz_competition; do
  input="$INPUT_DIR/panel_${panel}.phylo.faa"
  [[ -s "$input" ]] || { echo "missing panel input: $input" >&2; exit 1; }
  alignment="$OUTPUT_DIR/alignments/${panel}.aln.faa"
  trimmed="$OUTPUT_DIR/alignments/${panel}.trim.faa"
  prefix="$OUTPUT_DIR/trees/${panel}"
  mafft --auto --thread "$THREADS" "$input" > "$alignment" 2> "$OUTPUT_DIR/logs/${panel}.mafft.log"
  trimal -in "$alignment" -out "$trimmed" -automated1 > "$OUTPUT_DIR/logs/${panel}.trimal.log" 2>&1
  iqtree -s "$trimmed" -m LG+G4 -bb 1000 -T "$THREADS" --prefix "$prefix" > "$OUTPUT_DIR/logs/${panel}.iqtree.log" 2>&1
  tree="${prefix}.treefile"
  [[ -s "$tree" ]] || { echo "missing tree output: $tree" >&2; exit 1; }
  # Count terminal labels only; internal support/node labels also precede colons.
  final_leaves=$(grep -oE '(^|[,(])[^(),:;]+:' "$tree" | wc -l | tr -d ' ')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tcomplete\n' \
    "$panel" "$input" "$(sha256 "$input")" "$(count_fasta "$input")" \
    "$alignment" "$(sha256 "$alignment")" "$trimmed" "$(sha256 "$trimmed")" \
    "$tree" "$(sha256 "$tree")" "$final_leaves" "$THREADS" \
    "$(tool_version mafft)" "$(tool_version trimal)" "$(tool_version iqtree)" >> "$manifest"
done
echo "candidate-only phylogeny complete: $manifest"
