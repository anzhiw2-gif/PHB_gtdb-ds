#!/bin/bash
# Predict GTDB proteins and emit an auditable, fail-closed shard manifest.
set -Eeuo pipefail

THREADS=70
GENOMES_PER_SHARD=2000
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RUN_ROOT="${PHB_RUN_ROOT:-$REPO_ROOT}"
cd "$RUN_ROOT"
GTDB_DB="$HOME/GTDB/gtdb_genomes_reps_r232/database"
OUT="$RUN_ROOT/data/proteins"
PER_GENOME="$OUT/per_genome"
SHARDS="$OUT/shards"
LOG="$RUN_ROOT/logs"
GENOME_LIST="$OUT/genome_list.txt"
GENOME_RECORDS="$OUT/genome_records.tsv"
PROTEIN_FILES="$OUT/protein_files.txt"
SHARD_INDEX="$OUT/shard_index.tsv"
MANIFEST="$OUT/prediction_manifest.json"
FAILED_LOG="$PER_GENOME/failed.log"
mkdir -p "$PER_GENOME" "$SHARDS" "$LOG"

DRY=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --genomes-per-shard) GENOMES_PER_SHARD="$2"; shift 2 ;;
        --dry-run) DRY=1; shift ;;
        *) echo "unknown argument: $1" >&2; exit 1 ;;
    esac
done
if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ && "$GENOMES_PER_SHARD" =~ ^[1-9][0-9]*$ ]]; then
    echo "[ERROR] --threads and --genomes-per-shard must be positive integers" >&2
    exit 1
fi

find "$GTDB_DB" -type f -name "*_genomic.fna.gz" | sort > "$GENOME_LIST"
TOTAL=$(wc -l < "$GENOME_LIST")
if [ "$TOTAL" -le 0 ]; then
    echo "[ERROR] no GTDB genomes found under: $GTDB_DB" >&2
    exit 1
fi
awk -F/ '{name=$NF; sub(/_genomic\.fna\.gz$/, "", name); print name "\t" $0}' "$GENOME_LIST" > "$GENOME_RECORDS"
if [ "$(cut -f1 "$GENOME_RECORDS" | sort -u | wc -l)" -ne "$TOTAL" ]; then
    echo "[ERROR] duplicate genome accessions in prediction input" >&2
    exit 1
fi
EXPECTED_SHARDS=$(( (TOTAL + GENOMES_PER_SHARD - 1) / GENOMES_PER_SHARD ))
echo "[$(date)] input genomes: $TOTAL; expected shards: $EXPECTED_SHARDS"
if [ "$DRY" -eq 1 ]; then
    echo "[dry-run] no prediction files were written"
    exit 0
fi

predict_one() {
    local genome_path="$1" accession output temporary
    accession=$(basename "$genome_path" _genomic.fna.gz)
    output="$PER_GENOME/$accession.faa.gz"
    if [ -s "$output" ]; then
        return 0
    fi
    temporary="$PER_GENOME/$accession.faa.tmp"
    rm -f "$temporary"
    if pyrodigal -i "$genome_path" -a "$temporary" -p meta > /dev/null 2>&1 && [ -s "$temporary" ]; then
        sed -i "s/^>/>$accession|/" "$temporary"
        gzip -c "$temporary" > "$output"
        rm -f "$temporary"
        return 0
    fi
    rm -f "$temporary"
    printf '%s\n' "$genome_path" >> "$FAILED_LOG"
    return 1
}
export -f predict_one
export PER_GENOME FAILED_LOG

: > "$FAILED_LOG"
parallel -j "$THREADS" --progress predict_one {} :::: "$GENOME_LIST" 2> "$LOG/predict_progress.log"

: > "$PROTEIN_FILES"
DONE=0
while IFS=$'\t' read -r accession genome_path; do
    protein="$PER_GENOME/$accession.faa.gz"
    if [ ! -s "$protein" ]; then
        echo "[ERROR] missing predicted protein: $accession ($genome_path)" >&2
        exit 1
    fi
    printf '%s\n' "$protein" >> "$PROTEIN_FILES"
    DONE=$((DONE + 1))
done < "$GENOME_RECORDS"
FAIL=$(wc -l < "$FAILED_LOG")
if [ "$DONE" -ne "$TOTAL" ] || [ "$FAIL" -ne 0 ]; then
    echo "[ERROR] incomplete prediction: predicted=$DONE/$TOTAL failed=$FAIL" >&2
    exit 1
fi

EXPECTED_NAMES=$(seq -f 'shard_%04g.faa' 1 "$EXPECTED_SHARDS")
if find "$SHARDS" -maxdepth 1 -type f -name 'shard_*.faa' -printf '%f\n' | sort | comm -23 - <(printf '%s\n' "$EXPECTED_NAMES" | sort) | grep -q .; then
    echo "[ERROR] unexpected old shards present; isolate them before this run" >&2
    exit 1
fi

BUILD_DIR="$OUT/shards.build.$$.tmp"
mkdir "$BUILD_DIR"
trap 'rm -rf "$BUILD_DIR"' EXIT
: > "$SHARD_INDEX"
index=0
shard_number=0
while read -r protein; do
    if [ $((index % GENOMES_PER_SHARD)) -eq 0 ]; then
        shard_number=$((shard_number + 1))
        shard_path="$BUILD_DIR/shard_$(printf '%04d' "$shard_number").faa"
        : > "$shard_path"
    fi
    zcat "$protein" >> "$shard_path"
    printf '%04d\t%s\n' "$shard_number" "$(basename "$protein" .faa.gz)" >> "$SHARD_INDEX"
    index=$((index + 1))
done < "$PROTEIN_FILES"
if [ "$index" -ne "$TOTAL" ] || [ "$shard_number" -ne "$EXPECTED_SHARDS" ]; then
    echo "[ERROR] shard build count mismatch" >&2
    exit 1
fi
for name in $EXPECTED_NAMES; do
    if [ ! -s "$BUILD_DIR/$name" ]; then
        echo "[ERROR] missing or empty shard: $name" >&2
        exit 1
    fi
    mv -f "$BUILD_DIR/$name" "$SHARDS/$name"
done

python - "$TOTAL" "$EXPECTED_SHARDS" "$GENOMES_PER_SHARD" "$GENOME_LIST" "$SHARD_INDEX" "$SHARDS" "$MANIFEST" <<'PY'
import datetime as dt
import hashlib
import json
import pathlib
import sys

total, expected_count, per_shard, genome_list, shard_index, shard_dir, manifest_path = sys.argv[1:]
total, expected_count, per_shard = map(int, (total, expected_count, per_shard))

def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

counts = {}
with open(shard_index, encoding="utf-8") as handle:
    for line in handle:
        shard, _ = line.rstrip("\n").split("\t", 1)
        counts[shard] = counts.get(shard, 0) + 1
names = [f"shard_{number:04d}.faa" for number in range(1, expected_count + 1)]
records = []
for name in names:
    path = pathlib.Path(shard_dir, name)
    records.append({"name": name, "sha256": sha256(path), "genomes": counts.get(path.stem, 0)})
payload = {
    "schema_version": 1,
    "created_utc": dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
    "total_genomes": total,
    "predicted_genomes": total,
    "failed_genomes": 0,
    "genomes_per_shard": per_shard,
    "genome_list_sha256": sha256(genome_list),
    "expected_shards": expected_count,
    "expected_shard_names": names,
    "shards": records,
}
with open(manifest_path, "w", encoding="utf-8", newline="\n") as handle:
    json.dump(payload, handle, ensure_ascii=False, indent=2)
    handle.write("\n")
PY
python "$SCRIPT_DIR/05_validate_prediction_manifest.py" --manifest "$MANIFEST" --shard-dir "$SHARDS"
echo "[$(date)] prediction manifest verified: $MANIFEST"
