#!/bin/bash
# Build a complete family x shard HMMER matrix, validate provenance, then publish atomically.
#
# Task 9 — HMMER shard database size (Z)
# --------------------------------------
# HMMER E-values are E = Z * P, so E scales linearly with the number of target
# sequences Z the search is told it examined.  If a sharded search omits -Z,
# HMMER uses its default Z = <sequences in this shard>, and every shard then
# reports E-values on its own private scale.  Those numbers are not comparable
# across shards and must never be presented as one full-library scale.
#
# This entrypoint therefore REFUSES TO RUN without --database-size-z:
#   * it exits 1 when the flag is absent (no silent unscaled shard search);
#   * it exits 1 when the value is empty, zero, negative, or non-numeric;
#   * every HMMER call passes that same -Z "$DATABASE_SIZE_Z";
#   * the value, its stated basis, and a run contract file are recorded under
#     $SCREEN_DIR so a later reader can prove the scale;
#   * the manifest step records database_size_Z, database_size_basis, every
#     shard's sequence count, the total, and the command template.
# --domZ is never set: no downstream consumer rescales domain E-values, and
# adding a second scale would be undocumented (see the Task 9 plan step).
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RUN_ROOT="${PHB_RUN_ROOT:-$REPO_ROOT}"
cd "$RUN_ROOT"

THREADS=40
EVAL=1e-5
DATABASE_SIZE_Z="${PHB_DATABASE_SIZE_Z:-}"
DATABASE_SIZE_BASIS=""
SHARDS="$RUN_ROOT/data/proteins/shards_filt"
HMM_DIR="$RUN_ROOT/data/hmms/v2"
SCREEN_DIR="$RUN_ROOT/data/screen"
HMMOUT="$SCREEN_DIR/hmmsearch"
ARCHIVE_DIR="$SCREEN_DIR/archive"
BUILD_DIR="$SCREEN_DIR/hmmsearch.build.$$.tmp"
HITS_BUILD="$SCREEN_DIR/hits_all.tsv.build.$$.tmp"
MANIFEST_BUILD="$SCREEN_DIR/screen_manifest.json.build.$$.tmp"
RUN_CONTRACT_BUILD="$SCREEN_DIR/screen_run_contract.json.build.$$.tmp"
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ).$$"
LOG="$RUN_ROOT/logs"
FAILED="$LOG/screen_failed.$RUN_ID.log"
mkdir -p "$BUILD_DIR" "$ARCHIVE_DIR" "$LOG"
cleanup() { rm -rf "$BUILD_DIR" "$HITS_BUILD" "$MANIFEST_BUILD" "$RUN_CONTRACT_BUILD"; }
trap cleanup EXIT

source ~/miniconda3/etc/profile.d/conda.sh
conda activate phb_gtdb

FAMILIES="ePhaZ iPhaZ OH BdhA ArchPhaZ_patatin ArchPhaZ_hydrolase"
AUX_FAMILIES="PhaJ phasin PhaC"
SERVER=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --eval) EVAL="$2"; shift 2 ;;
        --families) FAMILIES="$2"; shift 2 ;;
        --database-size-z) DATABASE_SIZE_Z="$2"; shift 2 ;;
        --database-size-basis) DATABASE_SIZE_BASIS="$2"; shift 2 ;;
        --server) SERVER=1; shift ;;
        *) echo "unknown: $1" >&2; exit 1 ;;
    esac
done
ALL_FAMILIES="$FAMILIES $AUX_FAMILIES"
# Full-database scale is fail-closed.  A screening run that cannot state the
# complete target sequence count Z has no way to prove that the shard
# E-values it produces share one scale, so it must not run at all.
if [ -z "$DATABASE_SIZE_Z" ]; then
    echo "[ERROR] --database-size-z is required: it must be the COMPLETE library" >&2
    echo "        target sequence count Z (all shards, not this shard).  Without" >&2
    echo "        it every shard would fall back to HMMER's default Z = shard" >&2
    echo "        sequence count, producing incomparable per-shard E-values." >&2
    echo "        Measure or obtain the full-library count first, then rerun with" >&2
    echo "        --database-size-z <positive integer> --database-size-basis <source>." >&2
    echo "        This refit is fail-closed by design (Task 9): there is no default." >&2
    exit 1
fi
if ! [[ "$DATABASE_SIZE_Z" =~ ^[1-9][0-9]*$ ]]; then
    echo "[ERROR] --database-size-z must be a positive integer (no sign, spaces, or decimals): got '$DATABASE_SIZE_Z'" >&2
    exit 1
fi
if [ -z "$DATABASE_SIZE_BASIS" ]; then
    echo "[WARN] --database-size-basis not supplied; recording it as 'unspecified_unverified'" >&2
    DATABASE_SIZE_BASIS="unspecified_unverified"
fi
# Shared-server governance: 40 is a ceiling, never a launch value.  With
# --server active the requested count must fit inside the freshly measured
# limit min(40, nproc - /proc/loadavg 1m - 10); otherwise abort instead of
# overloading the host.
if [ "$SERVER" -eq 1 ]; then
    RESOURCE_RECORD="$(python "$SCRIPT_DIR/server_resources.py")" || {
        echo "[ERROR] cannot measure server capacity; refusing to launch" >&2; exit 1; }
    MEASURED_LIMIT="$(printf '%s' "$RESOURCE_RECORD" | sed -n 's/.*"threads": \([0-9]\+\).*/\1/p')"
    if ! [[ "$MEASURED_LIMIT" =~ ^[1-9][0-9]*$ ]]; then
        echo "[ERROR] server_resources.py returned no usable limit: $RESOURCE_RECORD" >&2
        exit 1
    fi
    if [ "$THREADS" -gt "$MEASURED_LIMIT" ]; then
        echo "[ERROR] requested $THREADS threads exceeds the measured limit ($MEASURED_LIMIT): $RESOURCE_RECORD" >&2
        exit 1
    fi
    echo "[$(date)] server capacity: $RESOURCE_RECORD"
fi

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
# Run contract: bind this screening run to the one full-database Z it used, on
# the same "missing evidence is never invented" convention as the run manifests
# (records the exact command template and the stated basis of the number).
python - "$RUN_CONTRACT_BUILD" "$RUN_ID" "$DATABASE_SIZE_Z" "$DATABASE_SIZE_BASIS" "$EVAL" "$THREADS" "$SHARDS" <<'PY'
import json, sys
out, run_id, size_z, basis, evalue, threads, shards = sys.argv[1:]
template = ('hmmsearch --tblout "$out" --domtblout "${out%.tbl}.dom" '
            '-Z "$DATABASE_SIZE_Z" -E "$EVAL" --cpu "$HMM_CPU" "$hmm" "$shard"')
payload = {
    "schema_version": 1,
    "run_id": run_id,
    "database_size_Z": int(size_z),
    "database_size_basis": basis,
    "database_size_basis_status": "unspecified_unverified" if basis == "unspecified_unverified" else "stated",
    "evalue": evalue,
    "threads": int(threads),
    "shard_dir": shards,
    "hmmsearch_command_template": template,
    "domz_used": False,
    "note": ("every shard is searched with the same full-library -Z; validation "
             "fails if any task contradicts this scale"),
}
with open(out, "w", encoding="utf-8", newline="\n") as handle:
    json.dump(payload, handle, indent=2, ensure_ascii=False)
    handle.write("\n")
PY
[[ -s "$RUN_CONTRACT_BUILD" ]] || { echo "[ERROR] run contract was not written" >&2; exit 1; }

run_one() {
    local hmm="$1" fam="$2" shard="$3"
    local sname out rc
    sname=$(basename "$shard" .faa)
    out="$BUILD_DIR/${fam}__${sname}.tbl"
    # -Z is the FULL-library database size, identical for every shard of every
    # family, so all shard E-values share one scale.  Do not add --domZ.
    if hmmsearch --tblout "$out" --domtblout "${out%.tbl}.dom" \
        -Z "$DATABASE_SIZE_Z" -E "$EVAL" --cpu 1 "$hmm" "$shard" > /dev/null 2>&1; then
        return 0
    else
        rc=$?
        printf 'FAIL %s %s (rc=%s)\n' "$fam" "$sname" "$rc" >> "$FAILED"
        return "$rc"
    fi
}
export -f run_one
export BUILD_DIR EVAL FAILED DATABASE_SIZE_Z

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
    --families "$ALL_FAMILIES" --eval "$EVAL" --out "$MANIFEST_BUILD" \
    --database-size-z "$DATABASE_SIZE_Z" \
    --database-size-basis "$DATABASE_SIZE_BASIS"
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
if [[ -e "$SCREEN_DIR/screen_run_contract.json" ]]; then mv "$SCREEN_DIR/screen_run_contract.json" "$ARCHIVE_DIR/screen_run_contract.json.$RUN_ID"; fi
mv "$BUILD_DIR" "$HMMOUT"
mv "$HITS_BUILD" "$SCREEN_DIR/hits_all.tsv"
mv "$MANIFEST_BUILD" "$SCREEN_DIR/screen_manifest.json"
mv "$RUN_CONTRACT_BUILD" "$SCREEN_DIR/screen_run_contract.json"
trap - EXIT
echo "[$(date)] screen manifest verified and published: $SCREEN_DIR/screen_manifest.json"
echo "[$(date)] database scale bound: Z=$DATABASE_SIZE_Z ($DATABASE_SIZE_BASIS)"
