#!/bin/bash
# Full pipeline orchestrator. Each stage must prove its declared handoff.
set -Eeuo pipefail

CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"
if [[ -f "$CONDA_SH" ]]; then
    # Make `conda run` available on hosts where conda is a shell function.
    source "$CONDA_SH"
fi
conda activate phb_gtdb
PYTHON_RUN=(conda run -n phb_gtdb python)

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd -P)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd -P)"
# ROOT is retained as a source-tree compatibility alias.  Runtime data and
# results must use RUN_ROOT so a run can never overwrite the repository.
ROOT="$REPO_ROOT"

THREADS_PREDICT=70
THREADS_SCREEN=70
THREADS_PHYLO=40
RUN_PHYLOGENY=0
LEGACY_ROOT=0
REQUESTED_RUN_ID=""
REQUESTED_RUN_ROOT=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-phylogeny) RUN_PHYLOGENY=1; shift ;;
        --run-id)
            [[ $# -ge 2 ]] || { echo "--run-id requires a value" >&2; exit 1; }
            REQUESTED_RUN_ID="$2"; shift 2 ;;
        --run-root|--run-dir)
            [[ $# -ge 2 ]] || { echo "--run-root requires a path" >&2; exit 1; }
            REQUESTED_RUN_ROOT="$2"; shift 2 ;;
        --legacy-root|--legacy-root-results) LEGACY_ROOT=1; shift ;;
        --help|-h)
            echo "usage: $0 [--run-id ID | --run-root PATH] [--run-phylogeny] [--legacy-root-results]"
            echo "default: execute in a new $REPO_ROOT/runs/<UTC-id> directory"
            echo "--legacy-root-results: explicitly reuse repository data/results (not isolated)"
            exit 0 ;;
        *) echo "unknown option: $1" >&2; exit 1 ;;
    esac
done

if [[ -n "$REQUESTED_RUN_ID" && -n "$REQUESTED_RUN_ROOT" ]]; then
    echo "--run-id and --run-root are mutually exclusive" >&2
    exit 1
fi
if [[ "$LEGACY_ROOT" -eq 1 && ( -n "$REQUESTED_RUN_ID" || -n "$REQUESTED_RUN_ROOT" ) ]]; then
    echo "--legacy-root cannot be combined with --run-id or --run-root" >&2
    exit 1
fi
if [[ -n "$REQUESTED_RUN_ID" && ! "$REQUESTED_RUN_ID" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]]; then
    echo "invalid --run-id (use letters, digits, '.', '_' or '-'): $REQUESTED_RUN_ID" >&2
    exit 1
fi
if [[ -n "$REQUESTED_RUN_ID" && ${#REQUESTED_RUN_ID} -gt 64 ]]; then
    echo "invalid --run-id (maximum length is 64): $REQUESTED_RUN_ID" >&2
    exit 1
fi
if [[ -n "$REQUESTED_RUN_ID" && "$REQUESTED_RUN_ID" == *..* ]]; then
    echo "invalid --run-id (path traversal marker '..' is forbidden): $REQUESTED_RUN_ID" >&2
    exit 1
fi

if [[ "$LEGACY_ROOT" -eq 1 ]]; then
    RUN_ID="legacy"
    RUN_ROOT="$REPO_ROOT"
else
    GIT_SHORT="$(cd "$REPO_ROOT" && git rev-parse --short HEAD 2>/dev/null || echo nogit)"
    GENERATED_RUN_ID="${REQUESTED_RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)_${GIT_SHORT}.$$}"
    RUN_ID="$GENERATED_RUN_ID"
    RUNS_ROOT="$REPO_ROOT/runs"
    if [[ -L "$RUNS_ROOT" ]]; then
        echo "refusing to use symlinked runs directory: $RUNS_ROOT" >&2
        exit 1
    fi
    mkdir -p "$RUNS_ROOT"
    RUNS_ROOT="$(cd "$RUNS_ROOT" && pwd -P)"
    if [[ -n "$REQUESTED_RUN_ROOT" ]]; then
        RUN_ROOT="$REQUESTED_RUN_ROOT"
        if [[ "$RUN_ROOT" != /* ]]; then
            case "$RUN_ROOT" in
                ..|../*|*/../*|*/..)
                    echo "relative --run-root must stay within the repository: $RUN_ROOT" >&2
                    exit 1 ;;
            esac
            RUN_ROOT="$REPO_ROOT/${RUN_ROOT#./}"
        fi
        case "$RUN_ROOT" in
            "$RUNS_ROOT"/*) : ;;
            *) echo "isolated --run-root must be inside $RUNS_ROOT: $RUN_ROOT" >&2; exit 1 ;;
        esac
        RELATIVE_RUN_ROOT="${RUN_ROOT#"$RUNS_ROOT"/}"
        if [[ "$RELATIVE_RUN_ROOT" == */* || "$RELATIVE_RUN_ROOT" == *..* ]]; then
            echo "--run-root must name one new child directly under $RUNS_ROOT" >&2
            exit 1
        fi
        RUN_ID="$RELATIVE_RUN_ROOT"
    else
        RUN_ROOT="$RUNS_ROOT/$RUN_ID"
    fi
    if [[ -e "$RUN_ROOT" || -L "$RUN_ROOT" ]]; then
        echo "run root already exists; refusing to overwrite: $RUN_ROOT" >&2
        exit 1
    fi
    mkdir "$RUN_ROOT" || {
        echo "could not atomically create run root (it may already exist): $RUN_ROOT" >&2
        exit 1
    }
    mkdir -p "$RUN_ROOT/data" "$RUN_ROOT/inputs" "$RUN_ROOT/logs" "$RUN_ROOT/results"
    ln -s ../logs "$RUN_ROOT/results/logs"
    HMM_SOURCE="$REPO_ROOT/data/hmms"
    HMM_LINK="$RUN_ROOT/data/hmms"
    if [[ ! -d "$HMM_SOURCE" ]]; then
        echo "required HMM input directory is missing: $HMM_SOURCE" >&2
        exit 1
    fi
    ln -s "$HMM_SOURCE" "$HMM_LINK"
fi

if [[ "$LEGACY_ROOT" -eq 1 ]]; then
    mkdir -p "$RUN_ROOT/inputs" "$RUN_ROOT/logs" "$RUN_ROOT/results/logs"
fi
export PHB_RUN_ROOT="$RUN_ROOT"
export PHB_RUN_ID="$RUN_ID"
cd "$RUN_ROOT"

LOG="$RUN_ROOT/logs/pipeline_master.log"
MANIFEST_JSONL="$RUN_ROOT/results/run_manifest.jsonl"
mkdir -p "$RUN_ROOT/logs" "$RUN_ROOT/results"
: > "$MANIFEST_JSONL"
printf 'run_id=%s\nrun_root=%s\nrepo_root=%s\nlegacy_root=%s\n' \
    "$RUN_ID" "$RUN_ROOT" "$REPO_ROOT" "$LEGACY_ROOT" > "$RUN_ROOT/run_context.env"

# Record GTDB, parameter, environment, and HMM inputs before execution.
# Missing GTDB files remain explicitly pending; no placeholder hash is made.
conda run -n phb_gtdb python - "$SCRIPT_DIR" "$RUN_ROOT" "$REPO_ROOT" <<'PY'
import importlib.util
import os
import pathlib
import sys

script_dir, run_root, repo_root = map(pathlib.Path, sys.argv[1:])
spec = importlib.util.spec_from_file_location("run_context", script_dir / "run_context.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
hmm_root = run_root / "data" / "hmms"
hmm_source = hmm_root.resolve()
hmm_paths = {
    f"hmm:{relative}": hmm_root / relative
    for path in sorted(hmm_source.rglob("*.hmm"))
    for relative in (path.relative_to(hmm_source),)
}
module.write_input_contract(
    run_root,
    gtdb_inputs={
        "taxonomy": pathlib.Path(os.path.expanduser("~/GTDB/taxonomy/bac120_taxonomy_r232.tsv")),
        "metadata": pathlib.Path(os.path.expanduser("~/GTDB/metadata/bac120_metadata_r232.tsv.gz")),
        "tree": pathlib.Path(os.path.expanduser("~/GTDB/GTDB_tree/bac120_r232.tree")),
    },
    inputs={
        "environment": repo_root / "environment.yml",
        "parameters": repo_root / "pipeline" / "config" / "params.yaml",
        **hmm_paths,
    },
    run_id=os.environ.get("PHB_RUN_ID", run_root.name),
)
PY
STEP_NAME=""
STEP_STARTED=""
STEP_NOTE=""
STEP_COMMAND=()

log() {
    echo "[$(date '+%F %T')] $1" | tee -a "$LOG"
}

step_begin() {
    STEP_NAME="$1"
    STEP_NOTE="$2"
    STEP_STARTED=$(date -u +%FT%TZ)
}

step_end() {
    local rc="$1"
    local ended
    ended=$(date -u +%FT%TZ)
    "${PYTHON_RUN[@]}" - "$MANIFEST_JSONL" "$STEP_NAME" "$rc" "$STEP_STARTED" "$ended" "$STEP_NOTE" \
        __COMMAND__ "${STEP_COMMAND[@]}" <<'PY'
import json
import sys

manifest_path, step, exit_code, started, ended, note = sys.argv[1:7]
separator = sys.argv.index("__COMMAND__", 7)
command = sys.argv[separator + 1:]
record = {
    "step": step,
    "exit_code": int(exit_code),
    "started": started,
    "ended": ended,
    "note": note,
    "command": command,
}
with open(manifest_path, "a", encoding="utf-8") as handle:
    handle.write(json.dumps(record, ensure_ascii=False, separators=(",", ":")) + "\n")
PY
    if [ "$rc" -ne 0 ]; then
        log "[ERROR] step $STEP_NAME failed (exit $rc); stopping"
        exit "$rc"
    fi
}

run_step() {
    local name="$1"
    local note="$2"
    shift 2
    STEP_COMMAND=("$@")
    log "[$name] $note"
    step_begin "$name" "$note"
    set +e
    "$@" 2>&1 | tee -a "$LOG"
    local -a pipe_status=("${PIPESTATUS[@]}")
    local rc="${pipe_status[0]}"
    if [ "$rc" -eq 0 ] && [ "${pipe_status[1]}" -ne 0 ]; then
        rc="${pipe_status[1]}"
    fi
    set -e
    step_end "$rc"
}

log "=== pipeline start (fail-closed) ==="
log "repository: $REPO_ROOT; run_root: $RUN_ROOT; git: $(cd "$REPO_ROOT" && git rev-parse --short HEAD 2>/dev/null || echo unknown)"

# Stage 05 is reusable only after the manifest proves the current GTDB input
# has one non-empty prediction per genome and a complete, hashed shard set.
run_step 05_predict_proteins "predict proteins and build verified shards" \
    bash "$SCRIPT_DIR/05_predict_proteins.sh" --threads "$THREADS_PREDICT"
run_step 05_validate_prediction_manifest "verify prediction input and shard hashes" \
    "${PYTHON_RUN[@]}" "$SCRIPT_DIR/05_validate_prediction_manifest.py" \
        --manifest data/proteins/prediction_manifest.json \
        --shard-dir data/proteins/shards

run_step 06a_filter_shards "filter sequences above the HMMER length limit" \
    bash "$SCRIPT_DIR/06a_filter_shards.sh"
run_step 06_screen "screen all filtered shards with declared HMMs" \
    bash "$SCRIPT_DIR/06_screen.sh" --threads "$THREADS_SCREEN" --eval 1e-5
run_step 07_process_hits "filter, arbitrate, and deduplicate HMM hits" \
    conda run -n phb_gtdb python "$SCRIPT_DIR/07_process_hits.py" \
        --hits data/screen/hits_all.tsv --shards data/proteins/shards_filt \
        --outdir data/screen
run_step 07b_extract_sequences "extract every selected protein sequence" \
    conda run -n phb_gtdb python "$SCRIPT_DIR/07b_extract_seqs.py" \
        --ids data/screen/unique_proteins.txt \
        --hits data/screen/hits_filtered.tsv \
        --shards data/proteins/shards_filt \
        --outdir data/screen/family_seqs
run_step 08_validate "apply sequence-level validation rules" \
    conda run -n phb_gtdb python "$SCRIPT_DIR/08_validate.py" \
        --indir data/screen/family_seqs --outdir data/screen
run_step 08c_tier_rescore "rescore validated sequences into evidence tiers" \
    bash "$SCRIPT_DIR/08c_tier_rescore.sh"
run_step 09a_tier1_summary "derive tier1 genome-family and phylum tables" \
    conda run -n phb_gtdb python "$SCRIPT_DIR/09a_tier1_summary.py"
if [ "$RUN_PHYLOGENY" -eq 1 ]; then
    run_step 09b_tier1_phylogeny "build registered tier1 phylogenetic trees" \
        bash "$SCRIPT_DIR/09b_tier1_phylogeny.sh" --threads "$THREADS_PHYLO"
else
    log "[09b_tier1_phylogeny] skip paused phylogeny (use --run-phylogeny to opt in)"
    step_begin 09b_tier1_phylogeny "skip paused phylogeny"
    STEP_COMMAND=("skip" "paused phylogeny")
    step_end 0
fi
if [ "$RUN_PHYLOGENY" -eq 1 ]; then
    run_step 09i_tree_manifest "register trees against current tier1 inputs" \
        conda run -n phb_gtdb python "$SCRIPT_DIR/09i_tree_manifest.py"
else
    log "[09i_tree_manifest] skip paused phylogeny registration"
    step_begin 09i_tree_manifest "skip paused phylogeny registration"
    STEP_COMMAND=("skip" "paused phylogeny registration")
    step_end 0
fi
run_step 10_distribution "compute taxonomy and ecology summaries" \
    conda run -n phb_gtdb python "$SCRIPT_DIR/10_distribution.py" \
        --hits data/screen/genome_hits.tsv
run_step 11_clusters "compute locus-level cluster context" \
    conda run -n phb_gtdb python "$SCRIPT_DIR/11_clusters.py" \
        --hits data/screen/hits_filtered.tsv

FINAL_OUTPUTS=(
    data/screen/hits_all.tsv
    data/screen/hits_filtered.tsv
    results/tables/tier1_genome_family.tsv
    results/tables/cluster_context.tsv
    results/tables/cluster_summary.tsv
    results/tables/cluster_locus_audit.tsv
    results/tables/cluster_genome_audit.tsv
)
HMM_INPUTS=(
    data/hmms/ePhaZ.hmm
    data/hmms/iPhaZ.hmm
    data/hmms/OH.hmm
    data/hmms/v2/ePhaZ.hmm
    data/hmms/v2/iPhaZ.hmm
    data/hmms/v2/OH.hmm
    data/hmms/v2/BdhA.hmm
    data/hmms/v2/ArchPhaZ_patatin.hmm
    data/hmms/v2/ArchPhaZ_hydrolase.hmm
    data/hmms/v2/PhaJ.hmm
    data/hmms/v2/phasin.hmm
    data/hmms/v2/PhaC.hmm
)
GTDB_INPUTS=(
    "$HOME/GTDB/taxonomy/bac120_taxonomy_r232.tsv"
    "$HOME/GTDB/metadata/bac120_metadata_r232.tsv.gz"
    "$HOME/GTDB/GTDB_tree/bac120_r232.tree"
)
if [ "$RUN_PHYLOGENY" -eq 1 ]; then
    FINAL_OUTPUTS+=(results/trees_tier1/tree_manifest.tsv)
fi

run_step finalize_manifest "freeze required inputs, outputs, and step results" \
    conda run -n phb_gtdb python "$SCRIPT_DIR/run_manifest.py" finalize \
        --jsonl "$MANIFEST_JSONL" \
        --out "$RUN_ROOT/results/run_manifest.json" \
        --inputs data/proteins/prediction_manifest.json \
                 "${HMM_INPUTS[@]}" \
        --outputs "${FINAL_OUTPUTS[@]}" \
        --gtdb-inputs "${GTDB_INPUTS[@]}" \
        --hmm-inputs "${HMM_INPUTS[@]}" \
        --source-files "$SCRIPT_DIR"/*.py "$SCRIPT_DIR"/*.sh \
                      "$ROOT/environment.yml" "$ROOT/pipeline/config/params.yaml" \
        --input-contract "$RUN_ROOT/input_contract.json" \
        --run-id "$RUN_ID" \
        --run-root "$RUN_ROOT" \
        --allow-pending-gtdb \
        --strict-provenance \
        --final-step-command conda run -n phb_gtdb python "$SCRIPT_DIR/run_manifest.py" finalize

# The wrapper appends the successful finalize step to JSONL after the first
# write. Re-freeze once so the JSON and JSONL contain the same step records.
if ! conda run -n phb_gtdb python "$SCRIPT_DIR/run_manifest.py" finalize \
    --jsonl "$MANIFEST_JSONL" \
    --out "$RUN_ROOT/results/run_manifest.json" \
    --inputs data/proteins/prediction_manifest.json \
             "${HMM_INPUTS[@]}" \
    --outputs "${FINAL_OUTPUTS[@]}" \
    --gtdb-inputs "${GTDB_INPUTS[@]}" \
    --hmm-inputs "${HMM_INPUTS[@]}" \
    --source-files "$SCRIPT_DIR"/*.py "$SCRIPT_DIR"/*.sh \
                  "$ROOT/environment.yml" "$ROOT/pipeline/config/params.yaml" \
    --input-contract "$RUN_ROOT/input_contract.json" \
    --run-id "$RUN_ID" \
    --run-root "$RUN_ROOT" \
    --allow-pending-gtdb \
    --strict-provenance; then
    log "[ERROR] final manifest re-freeze failed"
    exit 1
fi

log "=== pipeline complete: every recorded step exited zero ==="
