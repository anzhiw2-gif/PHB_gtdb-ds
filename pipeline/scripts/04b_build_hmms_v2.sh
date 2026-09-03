#!/bin/bash
# Build a complete, run-scoped HMM bundle from cleaned seeds.
# Layered ePhaZ example: --families ePhaZ_curated_core,ePhaZ_broad_discovery
set -Eeuo pipefail

THREADS=40
CDHIT_ID=0.90
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd -P)"
RUN_ROOT="${PHB_RUN_ROOT:-$REPO_ROOT}"
SEED_DIR="$RUN_ROOT/inputs/seeds_clean"
HMM_DIR="$RUN_ROOT/data/hmms/v2"
ALN_DIR="$RUN_ROOT/data/alignments/v2"
LOG_DIR="$RUN_ROOT/logs/hmm_build"
FAMILIES=(ePhaZ iPhaZ OH BdhA ArchPhaZ_patatin ArchPhaZ_hydrolase PhaJ PhaC phasin)
FAMILY_SELECTION="default"

parse_families() {
    local value="$1"
    local item
    local -a parsed
    IFS=',' read -r -a parsed <<< "$value"
    [[ "${#parsed[@]}" -gt 0 ]] || { echo "family list is empty" >&2; exit 2; }
    FAMILIES=()
    declare -A seen=()
    for item in "${parsed[@]}"; do
        item="${item//[[:space:]]/}"
        [[ "$item" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
            echo "invalid family name: $item" >&2; exit 2;
        }
        [[ -z "${seen[$item]+x}" ]] || { echo "duplicate family: $item" >&2; exit 2; }
        seen[$item]=1
        FAMILIES+=("$item")
    done
}

parse_family_file() {
    local list_path="$1"
    [[ -f "$list_path" && ! -L "$list_path" ]] || { echo "family list missing: $list_path" >&2; exit 1; }
    local joined=""
    local line
    while IFS= read -r line || [[ -n "$line" ]]; do
        line="${line%%#*}"
        line="${line//[[:space:]]/}"
        if [[ -n "$line" ]]; then
            [[ -z "$joined" ]] && joined="$line" || joined="$joined,$line"
        fi
    done < "$list_path"
    [[ -n "$joined" ]] || { echo "family list is empty: $list_path" >&2; exit 2; }
    parse_families "$joined"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --seed-dir) SEED_DIR="${2:?--seed-dir requires a path}"; shift 2 ;;
        --hmm-dir) HMM_DIR="${2:?--hmm-dir requires a path}"; shift 2 ;;
        --aln-dir) ALN_DIR="${2:?--aln-dir requires a path}"; shift 2 ;;
        --log-dir) LOG_DIR="${2:?--log-dir requires a path}"; shift 2 ;;
        --threads) THREADS="${2:?--threads requires a value}"; shift 2 ;;
        --cdhit-id) CDHIT_ID="${2:?--cdhit-id requires a value}"; shift 2 ;;
        --families)
            [[ "$FAMILY_SELECTION" == "default" ]] || { echo "--families/--family-list are mutually exclusive" >&2; exit 2; }
            parse_families "${2:?--families requires comma-separated names}"
            FAMILY_SELECTION="explicit"
            shift 2 ;;
        --family-list)
            [[ "$FAMILY_SELECTION" == "default" ]] || { echo "--families/--family-list are mutually exclusive" >&2; exit 2; }
            parse_family_file "${2:?--family-list requires a path}"
            FAMILY_SELECTION="explicit"
            shift 2 ;;
        --help|-h) echo "usage: $0 [--seed-dir DIR] [--hmm-dir DIR] [--aln-dir DIR] [--log-dir DIR] [--threads N] [--cdhit-id FLOAT] [--families NAME[,NAME...]] [--family-list FILE]"; exit 0 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "invalid --threads: $THREADS" >&2; exit 2; }
[[ -d "$SEED_DIR" ]] || { echo "seed directory missing: $SEED_DIR" >&2; exit 1; }
for fam in "${FAMILIES[@]}"; do
    fa="$SEED_DIR/$fam.faa"
    [[ -f "$fa" && ! -L "$fa" && -s "$fa" ]] || { echo "missing required seed FASTA: $fa" >&2; exit 1; }
done
mkdir -p "$HMM_DIR" "$ALN_DIR" "$LOG_DIR"

for fam in "${FAMILIES[@]}"; do
    fa="$SEED_DIR/$fam.faa"
    c90="$HMM_DIR/$fam.c90.faa"
    aln="$ALN_DIR/${fam}_aln.fasta"
    hmm="$HMM_DIR/$fam.hmm"
    echo "=== $fam ==="
    cd-hit -i "$fa" -o "$c90" -c "$CDHIT_ID" -n 5 -T "$THREADS" -M 0 \
        > "$LOG_DIR/cdhit_${fam}.log" 2>&1
    [[ -s "$c90" ]] || { echo "CD-HIT produced empty output: $c90" >&2; exit 1; }
    n90=$(grep -c '^>' "$c90")
    [[ "$n90" -ge 3 ]] || { echo "insufficient clustered seeds for $fam: $n90" >&2; exit 1; }
    if [[ "$n90" -gt 1500 ]]; then
        c80="$HMM_DIR/$fam.c80.faa"
        cd-hit -i "$c90" -o "$c80" -c 0.80 -n 5 -T "$THREADS" -M 0 \
            > "$LOG_DIR/cdhit_${fam}_c80.log" 2>&1
        [[ -s "$c80" ]] || { echo "CD-HIT c80 produced empty output: $c80" >&2; exit 1; }
        c90="$c80"
    fi
    mafft --auto --thread "$THREADS" "$c90" > "$aln" 2> "$LOG_DIR/mafft_${fam}.log"
    [[ -s "$aln" ]] || { echo "MAFFT produced empty output: $aln" >&2; exit 1; }
    hmmbuild --amino "$hmm" "$aln" > "$LOG_DIR/hmmbuild_${fam}.log" 2>&1
    [[ -s "$hmm" ]] || { echo "hmmbuild produced empty output: $hmm" >&2; exit 1; }
done

echo "HMM bundle complete: ${#FAMILIES[@]} families in $HMM_DIR"
