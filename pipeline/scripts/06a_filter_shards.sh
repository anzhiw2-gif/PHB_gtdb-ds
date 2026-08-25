#!/bin/bash
# Build filtered shards in a fresh directory and publish only a complete, validated run.
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RUN_ROOT="${PHB_RUN_ROOT:-$REPO_ROOT}"
cd "$RUN_ROOT"

THREADS=30
MAX_AA=100000
SOURCE_DIR="$RUN_ROOT/data/proteins/shards"
FILTERED_DIR="$RUN_ROOT/data/proteins/shards_filt"
ARCHIVE_DIR="$RUN_ROOT/data/proteins/archive"
BUILD_DIR="$RUN_ROOT/data/proteins/shards_filt.build.$$.tmp"
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ).$$"
FAILED="$RUN_ROOT/logs/filter_shards_failed.log"
mkdir -p "$BUILD_DIR" "$ARCHIVE_DIR" "$RUN_ROOT/logs"
trap 'rm -rf "$BUILD_DIR"' EXIT

source ~/miniconda3/etc/profile.d/conda.sh
: > "$FAILED"

filter_one() {
    local source="$1"
    local name
    name=$(basename "$source")
    if ! conda run -n phb_gtdb python "$SCRIPT_DIR/filter_long_seqs.py" "$source" "$BUILD_DIR/$name" "$MAX_AA" >> "$RUN_ROOT/logs/filter_shards.log" 2>&1; then
        printf 'FAIL %s\n' "$name" >> "$FAILED"
        return 1
    fi
}
export -f filter_one
export BUILD_DIR FAILED MAX_AA RUN_ROOT SCRIPT_DIR

mapfile -t SOURCE_SHARDS < <(find "$SOURCE_DIR" -maxdepth 1 -type f -name 'shard_*.faa' | sort)
if [ "${#SOURCE_SHARDS[@]}" -eq 0 ]; then
    echo "[ERROR] no source shards found: $SOURCE_DIR" >&2
    exit 1
fi
printf '%s\n' "${SOURCE_SHARDS[@]}" | parallel -j "$THREADS" filter_one {}
if [ -s "$FAILED" ]; then
    echo "[ERROR] one or more filtering tasks failed" >&2
    cat "$FAILED" >&2
    exit 1
fi

python "$SCRIPT_DIR/06a_validate_filter_manifest.py" --source-dir "$SOURCE_DIR" --filtered-dir "$BUILD_DIR" --max-aa "$MAX_AA" --filter-script "$SCRIPT_DIR/filter_long_seqs.py" --out "$BUILD_DIR/filter_manifest.json"
if [ -e "$FILTERED_DIR" ]; then
    mv "$FILTERED_DIR" "$ARCHIVE_DIR/shards_filt.$RUN_ID"
fi
mv "$BUILD_DIR" "$FILTERED_DIR"
trap - EXIT
echo "filter manifest verified: $FILTERED_DIR/filter_manifest.json"
