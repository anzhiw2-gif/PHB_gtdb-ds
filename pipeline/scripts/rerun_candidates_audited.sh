#!/bin/bash
# T141 scheme-A rerun in an isolated run directory. Existing results are read-only.
set -Eeuo pipefail

SOURCE_ROOT="/home/data/haoyu/PHB_gtdb-ds"
RUN_ROOT="${1:?usage: rerun_candidates_audited.sh RUN_ROOT [--preflight-only]}"
MODE="${2:-run}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python"
FAMILIES="ePhaZ iPhaZ OH BdhA ArchPhaZ_patatin ArchPhaZ_hydrolase PhaJ phasin PhaC"
CORE_FAMILIES="ePhaZ iPhaZ OH ArchPhaZ_patatin ArchPhaZ_hydrolase"
LOG="$RUN_ROOT/results/logs/rerun_candidates.log"
MANIFEST_JSONL="$RUN_ROOT/results/run_manifest.jsonl"

mkdir -p "$RUN_ROOT/data" "$RUN_ROOT/results/logs" "$RUN_ROOT/results/tables" "$RUN_ROOT/inputs"
link_input() {
    local target="$1" link="$2"
    if [ -L "$link" ]; then
        [ "$(readlink "$link")" = "$target" ] || { echo "[ERROR] link target mismatch: $link" >&2; exit 1; }
    elif [ -e "$link" ]; then
        echo "[ERROR] expected symlink path is occupied: $link" >&2
        exit 1
    else
        mkdir -p "$(dirname "$link")"
        ln -s "$target" "$link"
    fi
}
link_input "$SOURCE_ROOT/data/proteins/shards_filt" "$RUN_ROOT/data/proteins/shards_filt"
link_input "$SOURCE_ROOT/data/screen/hmmsearch" "$RUN_ROOT/data/screen/hmmsearch"
link_input "$SOURCE_ROOT/data/hmms" "$RUN_ROOT/data/hmms"
cd "$RUN_ROOT"

source /home/data/haoyu/miniconda3/etc/profile.d/conda.sh
conda activate phb_gtdb

record_step() {
    local name="$1" note="$2" rc="$3"
    printf '{"step":"%s","exit_code":%d,"started":"%s","ended":"%s","note":"%s"}\n' \
        "$name" "$rc" "$STEP_STARTED" "$(date -u +%FT%TZ)" "$note" >> "$MANIFEST_JSONL"
}

run_step() {
    local name="$1" note="$2"
    shift 2
    STEP_STARTED=$(date -u +%FT%TZ)
    echo "[$STEP_STARTED] [$name] $note" | tee -a "$LOG"
    set +e
    "$@" 2>&1 | tee -a "$LOG"
    local rc=${PIPESTATUS[0]}
    set -e
    record_step "$name" "$note" "$rc"
    if [ "$rc" -ne 0 ]; then
        echo "[ERROR] step $name failed: $rc" | tee -a "$LOG" >&2
        exit "$rc"
    fi
}

preflight() {
    test -x "$PY"
    test -d data/proteins/shards_filt
    test -d data/screen/hmmsearch
    for fam in $FAMILIES; do
        test -s "data/hmms/v2/${fam}.hmm"
    done
    for fam in $CORE_FAMILIES; do
        test -s "$SOURCE_ROOT/data/screen/family_seqs/${fam}_validated.faa"
    done
    if [ -s inputs/screen_manifest.json ]; then
        "$PY" "$SCRIPT_DIR/06_validate_screen_manifest.py" \
            --validate inputs/screen_manifest.json
    else
        "$PY" "$SCRIPT_DIR/06_validate_screen_manifest.py" \
            --shard-dir data/proteins/shards_filt \
            --hmm-dir data/hmms/v2 \
            --hmmout data/screen/hmmsearch \
            --families "$FAMILIES" --eval 1e-5 \
            --out inputs/screen_manifest.json
    fi
}

echo "[$(date -u +%FT%TZ)] scheme-A run root: $RUN_ROOT" | tee "$LOG"
if [ "$MODE" = "--resume-after-validation" ]; then
    test -s inputs/screen_manifest.json
    "$PY" "$SCRIPT_DIR/06_validate_screen_manifest.py" --validate inputs/screen_manifest.json
    test -s data/screen/hits_all.tsv
    test -s data/screen/hits_filtered.tsv
    test -s data/screen/unique_proteins.txt
    test -s data/screen/genome_hits.tsv
    test -s data/screen/validation.tsv
    for fam in $FAMILIES; do
        test -s "data/screen/family_seqs/${fam}_validated.faa"
    done
    echo "resume-after-validation: validated existing upstream outputs" | tee -a "$LOG"
elif [ "$MODE" = "--resume-after-process" ]; then
    test -s inputs/screen_manifest.json
    "$PY" "$SCRIPT_DIR/06_validate_screen_manifest.py" --validate inputs/screen_manifest.json
    test -s data/screen/hits_all.tsv
    test -s data/screen/hits_filtered.tsv
    test -s data/screen/unique_proteins.txt
    test -s data/screen/genome_hits.tsv
    echo "resume-after-process: validated existing pre-process outputs" | tee -a "$LOG"
else
    run_step preflight "validate immutable HMMER matrix and required inputs" preflight
    if [ "$MODE" = "--preflight-only" ]; then
        echo "preflight-only complete" | tee -a "$LOG"
        exit 0
    fi

    run_step aggregate_hits "rebuild hits_all.tsv with domain coverage" \
        "$PY" "$SCRIPT_DIR/06b_aggregate_hits.py" \
        --hmmout data/screen/hmmsearch --out data/screen/hits_all.tsv
    run_step process_hits "apply scheme-A OH min-cov 0.6 arbitration" \
        "$PY" "$SCRIPT_DIR/07_process_hits.py" \
        --hits data/screen/hits_all.tsv --shards data/proteins/shards_filt \
        --outdir data/screen --family-min-cov OH:0.6
fi
if [ "$MODE" != "--resume-after-validation" ]; then
    run_step extract_sequences "extract selected protein sequences" \
        "$PY" "$SCRIPT_DIR/07b_extract_seqs.py" \
        --ids data/screen/unique_proteins.txt --hits data/screen/hits_filtered.tsv \
        --shards data/proteins/shards_filt --outdir data/screen/family_seqs
    run_step validate_sequences "apply sequence-level validation" \
        "$PY" "$SCRIPT_DIR/08_validate.py" --indir data/screen/family_seqs --outdir data/screen
fi
run_step tier_rescore "curated HMM tier1/tier2 rescoring" \
    "$PY" "$SCRIPT_DIR/08c_tier_rescore.py"
run_step validate_tier_outputs "verify tier ids and FASTA counts" \
    "$PY" "$SCRIPT_DIR/08c_tier_rescore.py" --validate-build data/screen/tiers --families "$CORE_FAMILIES"
run_step tier_summary "derive tier1 genome-family and phylum tables" \
    "$PY" "$SCRIPT_DIR/09a_tier1_summary.py"
run_step distribution "recompute taxonomy and ecology tables" \
    "$PY" "$SCRIPT_DIR/10_distribution.py" \
    --hits data/screen/genome_hits.tsv \
    --taxonomy /home/data/haoyu/GTDB/taxonomy/bac120_taxonomy_r232.tsv \
    --metadata /home/data/haoyu/GTDB/metadata/bac120_metadata_r232.tsv.gz \
    --outdir results
run_step clusters "recompute locus-level neighborhoods" \
    "$PY" "$SCRIPT_DIR/11_clusters.py" \
    --hits data/screen/hits_filtered.tsv --marker-hmms data/hmms/v2 \
    --gtdb /home/data/haoyu/GTDB/gtdb_genomes_reps_r232/database \
    --outdir results --workdir data/cluster_work --flank-kb 10 --threads 40 --max-genomes 0

STEP_STARTED=$(date -u +%FT%TZ)
run_step finalize_manifest "freeze run inputs, steps, and outputs" \
    "$PY" "$SCRIPT_DIR/run_manifest.py" finalize \
    --jsonl "$MANIFEST_JSONL" --out results/run_manifest.json \
    --final-step-name finalize_manifest \
    --final-step-note "freeze run inputs, steps, and outputs" \
    --final-step-started "$STEP_STARTED" \
    --inputs inputs/screen_manifest.json \
             data/hmms/v2/ePhaZ.hmm data/hmms/v2/iPhaZ.hmm data/hmms/v2/OH.hmm \
             data/hmms/v2/ArchPhaZ_hydrolase.hmm \
    --outputs data/screen/hits_all.tsv data/screen/hits_filtered.tsv data/screen/genome_hits.tsv \
              data/screen/tiers/ePhaZ_tier1.faa data/screen/tiers/iPhaZ_tier1.faa \
              data/screen/tiers/OH_tier1.faa data/screen/tiers/ArchPhaZ_hydrolase_tier1.faa \
              results/tables/tier1_genome_family.tsv results/tables/tier1_phylum_distribution.tsv \
              results/tables/phylum_family_distribution.tsv results/tables/cluster_context.tsv \
              results/tables/cluster_summary.tsv results/tables/cluster_locus_audit.tsv \
              results/tables/cluster_genome_audit.tsv
echo "[$(date -u +%FT%TZ)] scheme-A run complete: $RUN_ROOT" | tee -a "$LOG"
