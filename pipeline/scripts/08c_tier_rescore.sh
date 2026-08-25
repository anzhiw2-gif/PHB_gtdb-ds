#!/bin/bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RUN_ROOT="${PHB_RUN_ROOT:-$REPO_ROOT}"
# ROOT is a deprecated runtime alias kept for existing HMM path assertions.
ROOT="$RUN_ROOT"
cd "$RUN_ROOT"
SEQDIR="$RUN_ROOT/data/screen/family_seqs"
TIERDIR="$RUN_ROOT/data/screen/tiers"
BUILD_DIR="$RUN_ROOT/data/screen/tiers.build.$$.tmp"
ARCHIVE_DIR="$RUN_ROOT/data/screen/archive"
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ).$$"
LOG="$RUN_ROOT/logs/08c_tier.log"
mkdir -p "$BUILD_DIR" "$ARCHIVE_DIR" "$RUN_ROOT/logs"
cleanup() { rm -rf "$BUILD_DIR"; }
trap cleanup EXIT

fail() {
  local message="$1"
  printf '[ERROR] %s\n' "$message" | tee -a "$LOG" >&2
  exit 1
}

require_nonempty() {
  local label="$1" path="$2"
  if [[ ! -s "$path" ]]; then
    fail "$label missing or empty: $path"
  fi
}

source ~/miniconda3/etc/profile.d/conda.sh
conda activate phb_gtdb
declare -A CURATED=(  [ePhaZ]=data/hmms/ePhaZ.hmm  [iPhaZ]=data/hmms/iPhaZ.hmm  [OH]=data/hmms/OH.hmm  [ArchPhaZ_patatin]=data/hmms/v2/ArchPhaZ_patatin.hmm  [ArchPhaZ_hydrolase]=data/hmms/v2/ArchPhaZ_hydrolase.hmm)
FAMILIES="ePhaZ iPhaZ OH ArchPhaZ_patatin ArchPhaZ_hydrolase"
for fam in $FAMILIES; do
  faa="$SEQDIR/${fam}_validated.faa"
  hmm="${CURATED[$fam]}"
  if [ ! -s "$faa" ]; then
    fail "$fam validated input missing or empty: $faa"
  fi
  if [ ! -s "$ROOT/$hmm" ]; then
    fail "$fam curated HMM missing or empty: $ROOT/$hmm"
  fi
  tier2_tbl="$BUILD_DIR/${fam}_tier2.tbl"
  tier2_ids="$BUILD_DIR/${fam}_tier2.ids"
  tier2_faa="$BUILD_DIR/${fam}_tier2.faa"
  tier1_tbl="$BUILD_DIR/${fam}_tier1.tbl"
  tier1_ids="$BUILD_DIR/${fam}_tier1.ids"
  tier1_faa="$BUILD_DIR/${fam}_tier1.faa"
  if hmmsearch --tblout "$tier2_tbl" -E 1e-10 --cpu 8 "$ROOT/$hmm" "$faa" > /dev/null 2>&1; then
    :
  else
    rc=$?
    fail "$fam tier2 hmmsearch failed (exit $rc)"
  fi
  require_nonempty "$fam tier2 tblout" "$tier2_tbl"
  awk '!/^#/ {print $1}' "$tier2_tbl" | sort -u > "$tier2_ids"
  require_nonempty "$fam tier2 hit IDs" "$tier2_ids"
  if python "$SCRIPT_DIR/08c_tier_rescore.py" --extract "$tier2_ids" "$faa" --output "$tier2_faa"; then
    :
  else
    rc=$?
    fail "$fam tier2 sequence extraction failed (exit $rc)"
  fi
  require_nonempty "$fam tier2 FASTA" "$tier2_faa"
  if hmmsearch --tblout "$tier1_tbl" -E 1e-20 --cpu 8 "$ROOT/$hmm" "$tier2_faa" > /dev/null 2>&1; then
    :
  else
    rc=$?
    fail "$fam tier1 hmmsearch failed (exit $rc)"
  fi
  require_nonempty "$fam tier1 tblout" "$tier1_tbl"
  awk '!/^#/ {print $1}' "$tier1_tbl" | sort -u > "$tier1_ids"
  require_nonempty "$fam tier1 hit IDs" "$tier1_ids"
  if python "$SCRIPT_DIR/08c_tier_rescore.py" --extract "$tier1_ids" "$faa" --output "$tier1_faa"; then
    :
  else
    rc=$?
    fail "$fam tier1 sequence extraction failed (exit $rc)"
  fi
  require_nonempty "$fam tier1 FASTA" "$tier1_faa"
done
if python "$SCRIPT_DIR/08c_tier_rescore.py" --validate-build "$BUILD_DIR" --families "$FAMILIES"; then
  :
else
  rc=$?
  fail "tier build validation failed (exit $rc)"
fi
if [ -e "$TIERDIR" ]; then mv "$TIERDIR" "$ARCHIVE_DIR/tiers.$RUN_ID"; fi
mv "$BUILD_DIR" "$TIERDIR"
trap - EXIT
echo "tier rescore verified and published: $TIERDIR"
