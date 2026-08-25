#!/bin/bash
# Build a complete family x shard HMMER matrix, validate provenance, then publish atomically.
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RUN_ROOT="${PHB_RUN_ROOT:-$REPO_ROOT}"
cd "$RUN_ROOT"

THREADS=70
EVAL=1e-5
SHARDS="$RUN_ROOT/data/proteins/shards_filt"
HMM_DIR="$RUN_ROOT/data/hmms/v2"
SCREEN_DIR="$RUN_ROOT/data/screen"
HMMOUT="$SCREEN_DIR/hmmsearch"
ARCHIVE_DIR="$SCREEN_DIR/archive"
BUILD_DIR="$SCREEN_DIR/hmmsearch.build.$$.tmp"
HITS_BUILD="$SCREEN_DIR/hits_all.tsv.build.$$.tmp"
MANIFEST_BUILD="$SCREEN_DIR/screen_manifest.json.build.$$.tmp"
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ).$$"
LOG="$RUN_ROOT/logs"
FAILED="$LOG/screen_failed.$RUN_ID.log"
mkdir -p "$BUILD_DIR" "$ARCHIVE_DIR" "$LOG"
cleanup() { rm -rf "$BUILD_DIR" "$HITS_BUILD" "$MANIFEST_BUILD"; }
trap cleanup EXIT

source ~/miniconda3/etc/profile.d/conda.sh
conda activate phb_gtdb

FAMILIES="ePhaZ iPhaZ OH BdhA ArchPhaZ_patatin ArchPhaZ_hydrolase"
AUX_FAMILIES="PhaJ phasin PhaC"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --eval) EVAL="$2"; shift 2 ;;
        --families) FAMILIES="$2"; shift 2 ;;
        *) echo "unknown: $1" >&2; exit 1 ;;
    esac
done
ALL_FAMILIES="$FAMILIES $AUX_FAMILIES"

mapfile -t SHARD_FILES < <(find "$SHARDS" -maxdepth 1 -type f -name 'shard_*.faa' | sort)
if [[ "${#SHARD_FILES[@]}" -eq 0 ]]; then
    echo "[ERROR] no filtered shards found: $SHARDS" >&2
    exit 1
fi
for shard in "${SHARD_FILES[@]}"; do
    if [[ ! -s "$shard" ]]; then
        echo "[ERROR] declared shard is missing or empty: $shard" >&2
        exit 1
    fi
done
for fam in $ALL_FAMILIES; do
    hmm="$HMM_DIR/$fam.hmm"
    if [ ! -s "$hmm" ]; then
        echo "[ERROR] $fam HMM 缺失" >&2
        exit 1
    fi
done
: > "$FAILED"

run_one() {
    local hmm="$1" fam="$2" shard="$3"
    local sname out rc
    sname=$(basename "$shard" .faa)
    out="$BUILD_DIR/${fam}__${sname}.tbl"
    if hmmsearch --tblout "$out" --domtblout "${out%.tbl}.dom" \
        -E "$EVAL" --cpu 1 "$hmm" "$shard" > /dev/null 2>&1; then
        return 0
    else
        rc=$?
        printf 'FAIL %s %s (rc=%s)\n' "$fam" "$sname" "$rc" >> "$FAILED"
        return "$rc"
    fi
}
export -f run_one
export BUILD_DIR EVAL FAILED

echo "[$(date)] screening ${#SHARD_FILES[@]} shards x families: $ALL_FAMILIES"
for fam in $ALL_FAMILIES; do
    hmm="$HMM_DIR/$fam.hmm"
    if printf '%s\n' "${SHARD_FILES[@]}" | parallel -j "$THREADS" run_one "$hmm" "$fam" {} 2> "$LOG/screen_${fam}.log"; then
        :
    else
        local_rc=$?
        echo "[ERROR] $fam parallel exit code $local_rc; refusing partial publication" >&2
        if [[ -s "$FAILED" ]]; then
            cat "$FAILED" >&2
        fi
        exit "$local_rc"
    fi
done
if [[ -s "$FAILED" ]]; then
    echo "[ERROR] HMMER tasks failed; refusing partial publication" >&2
    cat "$FAILED" >&2
    exit 1
fi

python "$SCRIPT_DIR/06_validate_screen_manifest.py" \
    --shard-dir "$SHARDS" --hmm-dir "$HMM_DIR" --hmmout "$BUILD_DIR" \
    --families "$ALL_FAMILIES" --eval "$EVAL" --out "$MANIFEST_BUILD"
for fam in $ALL_FAMILIES; do
    for shard in "${SHARD_FILES[@]}"; do
        stem=$(basename "$shard" .faa)
        tbl="$BUILD_DIR/${fam}__${stem}.tbl"
        dom="$BUILD_DIR/${fam}__${stem}.dom"
        if [[ ! -s "$tbl" || ! -s "$dom" ]]; then
            echo "[ERROR] missing or empty HMMER output: $fam x $stem" >&2
            exit 1
        fi
    done
done
python "$SCRIPT_DIR/06b_aggregate_hits.py" --hmmout "$BUILD_DIR" --out "$HITS_BUILD"
[[ -s "$HITS_BUILD" ]] || { echo "[ERROR] aggregate output is empty" >&2; exit 1; }

if [[ -e "$HMMOUT" ]]; then mv "$HMMOUT" "$ARCHIVE_DIR/hmmsearch.$RUN_ID"; fi
if [[ -e "$SCREEN_DIR/hits_all.tsv" ]]; then mv "$SCREEN_DIR/hits_all.tsv" "$ARCHIVE_DIR/hits_all.tsv.$RUN_ID"; fi
if [[ -e "$SCREEN_DIR/screen_manifest.json" ]]; then mv "$SCREEN_DIR/screen_manifest.json" "$ARCHIVE_DIR/screen_manifest.json.$RUN_ID"; fi
mv "$BUILD_DIR" "$HMMOUT"
mv "$HITS_BUILD" "$SCREEN_DIR/hits_all.tsv"
mv "$MANIFEST_BUILD" "$SCREEN_DIR/screen_manifest.json"
trap - EXIT
echo "[$(date)] screen manifest verified and published: $SCREEN_DIR/screen_manifest.json"
