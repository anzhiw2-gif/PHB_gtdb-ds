#!/usr/bin/env bash
# Resume a frozen formal scan in a new dated run without overwriting its parent.
set -Eeuo pipefail

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate phb_gtdb

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd -P)"
if [[ -d "$SCRIPT_DIR/../../pipeline" ]]; then
    REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd -P)"
else
    REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd -P)"
fi
REGISTRY="${FORMAL_SCAN_REGISTRY:-$SCRIPT_DIR/../config/formal_scan_models.tsv}"
RUN_ID=""; SHARD_DIR=""; SPLIT_ROOT=""; SEEDCLEAN_ROOT=""; FROZEN_MODEL_ROOT=""; RESUME_FROM=""
THREADS=60

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-id) RUN_ID="${2:?--run-id requires a value}"; shift 2 ;;
        --shard-dir) SHARD_DIR="${2:?--shard-dir requires a path}"; shift 2 ;;
        --split-root) SPLIT_ROOT="${2:?--split-root requires a path}"; shift 2 ;;
        --seedclean-root) SEEDCLEAN_ROOT="${2:?--seedclean-root requires a path}"; shift 2 ;;
        --frozen-model-root) FROZEN_MODEL_ROOT="${2:?--frozen-model-root requires a path}"; shift 2 ;;
        --resume-from) RESUME_FROM="${2:?--resume-from requires a path}"; shift 2 ;;
        --threads) THREADS="${2:?--threads requires a value}"; shift 2 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

[[ "$RUN_ID" =~ ^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$ && "$RUN_ID" != *..* ]] || exit 2
# THREADS -le 60 is the hard concurrency ceiling.
[[ "$THREADS" =~ ^[1-9][0-9]*$ && "$THREADS" -le 60 ]] || { echo "--threads must be 1..60" >&2; exit 2; }
[[ -s "$REGISTRY" && -d "$SHARD_DIR" && -d "$SPLIT_ROOT" && -d "$SEEDCLEAN_ROOT" && -d "$FROZEN_MODEL_ROOT" ]] || exit 1
[[ -d "$RESUME_FROM" ]] || { echo "--resume-from must be an existing parent run" >&2; exit 1; }

RUN_ROOT="$REPO_ROOT/runs/$RUN_ID"
[[ ! -e "$RUN_ROOT" && ! -L "$RUN_ROOT" ]] || { echo "run already exists: $RUN_ROOT" >&2; exit 1; }
mkdir -p "$RUN_ROOT/logs/task_errors" "$RUN_ROOT/inputs/hmms" "$RUN_ROOT/inputs/scan_shards" \
    "$RUN_ROOT/results/hmmsearch.build" "$RUN_ROOT/results/task_status"

for script in formal_frozen_screen_parallel.sh monitor_formal_scan.sh 06b_aggregate_hits.py filter_hmmsearch_shard.py run_context.py; do
    cp "$SCRIPT_DIR/$script" "$RUN_ROOT/inputs/$script"
done
cp "$REGISTRY" "$RUN_ROOT/inputs/formal_scan_models.tsv"
printf '%s\n' "$RESUME_FROM" > "$RUN_ROOT/inputs/resume_parent.txt"

python - "$SCRIPT_DIR/run_context.py" "$RUN_ROOT" "$RUN_ID" "$REGISTRY" "$SCRIPT_DIR/formal_frozen_screen_parallel.sh" "$SCRIPT_DIR/06b_aggregate_hits.py" "$SCRIPT_DIR/filter_hmmsearch_shard.py" "$SCRIPT_DIR/run_context.py" <<'PY'
import importlib.util, pathlib, sys
script, run, run_id, registry = map(pathlib.Path, sys.argv[1:5])
spec = importlib.util.spec_from_file_location("run_context", script)
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
module.write_input_contract(run, run_id=str(run_id), gtdb_inputs={
    "taxonomy": pathlib.Path.home()/"GTDB/taxonomy/bac120_taxonomy_r232.tsv",
    "metadata": pathlib.Path.home()/"GTDB/metadata/bac120_metadata_r232.tsv.gz",
    "tree": pathlib.Path.home()/"GTDB/GTDB_tree/bac120_r232.tree",
}, inputs={
    "registry": registry, "formal_frozen_screen_parallel": pathlib.Path(sys.argv[5]),
    "aggregate_hits": pathlib.Path(sys.argv[6]), "filter_hmmsearch_shard": pathlib.Path(sys.argv[7]),
    "run_context": pathlib.Path(sys.argv[8]),
})
PY

models="$(awk 'NR > 1 && $0 !~ /^[[:space:]]*$/ { n++ } END { print n+0 }' "$REGISTRY")"
[[ "$models" -eq 10 ]] || { echo "expected 10 frozen models" >&2; exit 1; }
while IFS=$'\t' read -r model source threshold min_cov group; do
    [[ "$model" == model ]] && continue
    [[ "$threshold" == e-5 ]] || exit 1
    case "$source" in
        split_run) root="$SPLIT_ROOT/data/hmms/v2" ;;
        seedclean_run) root="$SEEDCLEAN_ROOT/data/hmms/v2" ;;
        frozen_data_root) root="$FROZEN_MODEL_ROOT" ;;
        *) exit 1 ;;
    esac
    [[ -s "$root/$model.hmm" ]] || { echo "missing HMM: $model" >&2; exit 1; }
    cp "$root/$model.hmm" "$RUN_ROOT/inputs/hmms/$model.hmm"
done < "$REGISTRY"

find "$SHARD_DIR" -maxdepth 1 -type f -name 'shard_*.faa' | sort > "$RUN_ROOT/inputs/shard_paths.txt"
shards="$(wc -l < "$RUN_ROOT/inputs/shard_paths.txt")"
[[ "$shards" -eq 100 ]] || { echo "expected 100 shards" >&2; exit 1; }

EXCLUSIONS="$RUN_ROOT/results/overlength_exclusions.tsv"
if [[ -s "$RESUME_FROM/results/overlength_exclusions.tsv" ]]; then
    cp "$RESUME_FROM/results/overlength_exclusions.tsv" "$EXCLUSIONS"
else
    printf 'source_shard\taccession\tlength_aa\treason\n' > "$EXCLUSIONS"
fi
while read -r shard; do
    stem="$(basename "$shard" .faa)"
    filtered="$RUN_ROOT/inputs/scan_shards/$stem.faa"
    parent="$RESUME_FROM/inputs/scan_shards/$stem.faa"
    if [[ -s "$parent" ]]; then cp "$parent" "$filtered"; else python "$RUN_ROOT/inputs/filter_hmmsearch_shard.py" "$shard" "$filtered" "$EXCLUSIONS"; fi
    [[ -s "$filtered" ]] || { echo "empty filtered shard: $stem" >&2; exit 1; }
done < "$RUN_ROOT/inputs/shard_paths.txt"

TASKS="$RUN_ROOT/inputs/tasks.tsv"; BUILD="$RUN_ROOT/results/hmmsearch.build"; STATUS="$RUN_ROOT/results/task_status"
printf 'model\tstem\thmm\tfiltered\n' > "$TASKS"; reused=0
while IFS=$'\t' read -r model source threshold min_cov group; do
    [[ "$model" == model ]] && continue
    while read -r shard; do
        stem="$(basename "$shard" .faa)"; base="$model"__"$stem"
        if [[ -s "$RESUME_FROM/results/hmmsearch.build/$base.tbl" && -s "$RESUME_FROM/results/hmmsearch.build/$base.dom" ]]; then
            cp "$RESUME_FROM/results/hmmsearch.build/$base.tbl" "$BUILD/$base.tbl"
            cp "$RESUME_FROM/results/hmmsearch.build/$base.dom" "$BUILD/$base.dom"
            : > "$STATUS/$base.ok"; reused=$((reused+1))
        else
            printf '%s\t%s\t%s\t%s\n' "$model" "$stem" "$RUN_ROOT/inputs/hmms/$model.hmm" "$RUN_ROOT/inputs/scan_shards/$stem.faa" >> "$TASKS"
        fi
    done < "$RUN_ROOT/inputs/shard_paths.txt"
done < "$REGISTRY"

run_task() {
    local model="$1" stem="$2" hmm="$3" filtered="$4" base="$1"__"$2"
    if hmmsearch --tblout "$BUILD/$base.tbl" --domtblout "$BUILD/$base.dom" -E 1e-5 --cpu 1 "$hmm" "$filtered" > /dev/null 2> "$RUN_ROOT/logs/task_errors/$base.stderr"; then
        [[ -e "$BUILD/$base.tbl" && -e "$BUILD/$base.dom" ]] || return 1
        : > "$STATUS/$base.ok"
    else
        printf '%s\t%s\thmmsearch_failed\n' "$model" "$stem" > "$RUN_ROOT/logs/task_errors/$base.tsv"
        return 1
    fi
}
export -f run_task; export BUILD STATUS RUN_ROOT
pending="$(($(wc -l < "$TASKS") - 1))"
if (( pending > 0 )); then
    tail -n +2 "$TASKS" | parallel -j "$THREADS" --halt soon,fail=1 --joblog "$RUN_ROOT/logs/parallel.joblog" --colsep '\t' run_task
fi

FAILED="$RUN_ROOT/logs/failed_tasks.tsv"; printf 'model\tshard\treason\n' > "$FAILED"
for file in "$RUN_ROOT/logs/task_errors"/*.tsv; do [[ -e "$file" ]] && cat "$file" >> "$FAILED"; done
completed="$(find "$STATUS" -maxdepth 1 -type f -name '*.ok' | wc -l)"; expected=$((models * shards))
[[ "$completed" -eq "$expected" && "$(wc -l < "$FAILED")" -eq 1 ]] || { echo "incomplete: $completed/$expected" >&2; exit 1; }

python "$RUN_ROOT/inputs/06b_aggregate_hits.py" --hmmout "$BUILD" --out "$RUN_ROOT/results/hits_all.tsv"
[[ -s "$RUN_ROOT/results/hits_all.tsv" ]] || exit 1
mv "$BUILD" "$RUN_ROOT/results/hmmsearch"
PHB_PARENT_RUN="$RESUME_FROM" PHB_REUSED="$reused" PHB_THREADS="$THREADS" python - "$RUN_ROOT/results/scan_manifest.json" "$RUN_ROOT" <<'PY'
import datetime, hashlib, json, os, pathlib, sys
out, run = map(pathlib.Path, sys.argv[1:])
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20), b""): h.update(b)
    return h.hexdigest()
models=[]
for row in (run/"inputs/formal_scan_models.tsv").read_text().splitlines()[1:]:
    if row.strip():
        model, source, threshold, min_cov, group=row.split("\t"); hmm=run/"inputs/hmms"/(model+".hmm")
        models.append({"model":model,"hmm_source":source,"threshold":threshold,"min_cov":float(min_cov),"report_group":group,"sha256":sha(hmm),"bytes":hmm.stat().st_size})
out.write_text(json.dumps({"schema_version":2,"status":"completed","run_id":run.name,"created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"parent_run":os.environ["PHB_PARENT_RUN"],"reused_tasks":int(os.environ["PHB_REUSED"]),"threads":int(os.environ["PHB_THREADS"]),"task_total":len(models)*100,"task_completed":len(list((run/"results/task_status").glob("*.ok"))),"registry_sha256":sha(run/"inputs/formal_scan_models.tsv"),"models":models,"overlength_exclusions_sha256":sha(run/"results/overlength_exclusions.tsv"),"hits_all_sha256":sha(run/"results/hits_all.tsv")},ensure_ascii=False,indent=2)+"\n")
PY
echo "complete: run=$RUN_ROOT threads=$THREADS reused=$reused"
