#!/usr/bin/env bash
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
RUN_ID=""; SHARD_DIR=""; SPLIT_ROOT=""; SEEDCLEAN_ROOT=""; FROZEN_MODEL_ROOT=""; THREADS=70; PREFLIGHT_ONLY=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-id) RUN_ID="$2"; shift 2 ;;
    --shard-dir) SHARD_DIR="$2"; shift 2 ;;
    --split-root) SPLIT_ROOT="$2"; shift 2 ;;
    --seedclean-root) SEEDCLEAN_ROOT="$2"; shift 2 ;;
    --frozen-model-root) FROZEN_MODEL_ROOT="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --preflight-only) PREFLIGHT_ONLY=true; shift ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done
[[ "$RUN_ID" =~ ^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$ && "$RUN_ID" != *..* ]] || exit 2
[[ "$THREADS" =~ ^[1-9][0-9]*$ && -s "$REGISTRY" && -d "$SHARD_DIR" && -d "$SPLIT_ROOT" && -d "$SEEDCLEAN_ROOT" && -d "$FROZEN_MODEL_ROOT" ]] || exit 1
RUN_ROOT="$REPO_ROOT/runs/$RUN_ID"
[[ ! -e "$RUN_ROOT" && ! -L "$RUN_ROOT" ]] || exit 1
mkdir "$RUN_ROOT"; mkdir -p "$RUN_ROOT/logs" "$RUN_ROOT/inputs/hmms" "$RUN_ROOT/results"
cp "$REGISTRY" "$RUN_ROOT/inputs/formal_scan_models.tsv"
cp "$SCRIPT_DIR/formal_frozen_screen.sh" "$RUN_ROOT/inputs/formal_frozen_screen.sh"
cp "$SCRIPT_DIR/06b_aggregate_hits.py" "$RUN_ROOT/inputs/06b_aggregate_hits.py"
cp "$SCRIPT_DIR/filter_hmmsearch_shard.py" "$RUN_ROOT/inputs/filter_hmmsearch_shard.py"
cp "$SCRIPT_DIR/run_context.py" "$RUN_ROOT/inputs/run_context.py"
python - "$SCRIPT_DIR/run_context.py" "$RUN_ROOT" "$RUN_ID" "$REGISTRY" "$SCRIPT_DIR/formal_frozen_screen.sh" "$SCRIPT_DIR/06b_aggregate_hits.py" "$SCRIPT_DIR/filter_hmmsearch_shard.py" "$SCRIPT_DIR/run_context.py" <<'PY'
import importlib.util, pathlib, sys
script = pathlib.Path(sys.argv[1])
run_root = pathlib.Path(sys.argv[2])
run_id = sys.argv[3]
registry = pathlib.Path(sys.argv[4])
spec = importlib.util.spec_from_file_location("run_context", script)
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
module.write_input_contract(
    run_root,
    run_id=run_id,
    gtdb_inputs={
        "taxonomy": pathlib.Path.home() / "GTDB/taxonomy/bac120_taxonomy_r232.tsv",
        "metadata": pathlib.Path.home() / "GTDB/metadata/bac120_metadata_r232.tsv.gz",
        "tree": pathlib.Path.home() / "GTDB/GTDB_tree/bac120_r232.tree",
    },
    inputs={
        "registry": registry,
        "formal_frozen_screen": pathlib.Path(sys.argv[5]),
        "aggregate_hits": pathlib.Path(sys.argv[6]),
        "filter_hmmsearch_shard": pathlib.Path(sys.argv[7]),
        "run_context": pathlib.Path(sys.argv[8]),
    },
)
PY
expected_count=$(awk 'BEGIN { count=0 } NR > 1 && $0 !~ /^[[:space:]]*$/ { count++ } END { print count }' "$REGISTRY")
[[ "$expected_count" -ge 1 ]] || exit 1
count=0
while IFS=$'\t' read -r model source threshold min_cov report_group; do
  [[ "$model" == "model" ]] && continue
  [[ "$threshold" == "e-5" ]] || exit 1
  case "$source" in
    split_run) base="$SPLIT_ROOT/data/hmms/v2" ;;
    seedclean_run) base="$SEEDCLEAN_ROOT/data/hmms/v2" ;;
    frozen_data_root) base="$FROZEN_MODEL_ROOT" ;;
    *) echo "unsupported registry source: $source" >&2; exit 1 ;;
  esac
  [[ -d "$base" ]] || { echo "model source directory missing: $base" >&2; exit 1; }
  [[ -s "$base/$model.hmm" ]] || exit 1
  cp "$base/$model.hmm" "$RUN_ROOT/inputs/hmms/$model.hmm"
  count=$((count+1))
done < "$REGISTRY"
[[ "$count" -eq "$expected_count" ]] || exit 1
find "$SHARD_DIR" -maxdepth 1 -type f -name 'shard_*.faa' | sort > "$RUN_ROOT/inputs/shard_paths.txt"
[[ "$(wc -l < "$RUN_ROOT/inputs/shard_paths.txt")" -eq 100 ]] || exit 1
while read -r shard; do [[ -s "$shard" ]] || exit 1; done < "$RUN_ROOT/inputs/shard_paths.txt"
python - "$RUN_ROOT/inputs/shards.tsv" "$SHARD_DIR" <<'PY'
import hashlib,pathlib,sys
out,root=map(pathlib.Path,sys.argv[1:])
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for b in iter(lambda:f.read(1<<20),b""): h.update(b)
 return h.hexdigest()
with out.open("w",encoding="utf-8") as h:
 h.write("name\tpath\tbytes\tsha256\n")
 for p in sorted(root.glob("shard_*.faa")): h.write(f"{p.name}\t{p}\t{p.stat().st_size}\t{sha(p)}\n")
PY
if [[ "$PREFLIGHT_ONLY" == true ]]; then
  python - "$RUN_ROOT/results/preflight_manifest.json" "$RUN_ROOT" "$SHARD_DIR" <<'PY'
import datetime as dt, hashlib, json, pathlib, subprocess, sys
out, run, shard_dir = map(pathlib.Path, sys.argv[1:])
def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()
registry = run / "inputs" / "formal_scan_models.tsv"
models = []
for line in registry.read_text(encoding="utf-8").splitlines()[1:]:
    if not line.strip():
        continue
    model, source, threshold, min_cov, report_group = line.split("\t")
    hmm = run / "inputs" / "hmms" / f"{model}.hmm"
    models.append({
        "model": model,
        "hmm_source": source,
        "threshold": threshold,
        "min_cov": float(min_cov),
        "report_group": report_group,
        "bytes": hmm.stat().st_size,
        "sha256": sha(hmm),
    })
tool = subprocess.run(["hmmsearch", "-h"], capture_output=True, text=True, check=True)
version = (tool.stdout or tool.stderr).splitlines()[0]
payload = {
    "schema_version": 1,
    "status": "planned_not_run",
    "run_id": run.name,
    "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
    "registry_sha256": sha(registry),
    "source_shard_dir": str(shard_dir),
    "shard_count": len(list(shard_dir.glob("shard_*.faa"))),
    "shards_tsv_sha256": sha(run / "inputs" / "shards.tsv"),
    "models": models,
    "hmmsearch_version": version,
    "hmmsearch_outputs_created": False,
}
out.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
PY
  [[ ! -e "$RUN_ROOT/results/hmmsearch.build" && ! -e "$RUN_ROOT/results/hmmsearch" ]] || exit 1
  echo "formal frozen screen preflight complete: $RUN_ROOT"
  exit 0
fi
BUILD="$RUN_ROOT/results/hmmsearch.build"; FAILED="$RUN_ROOT/logs/failed_tasks.tsv"; : > "$FAILED"
mkdir -p "$BUILD" "$RUN_ROOT/inputs/scan_shards"
EXCLUSIONS="$RUN_ROOT/results/overlength_exclusions.tsv"
printf 'source_shard\taccession\tlength_aa\treason\n' > "$EXCLUSIONS"
while read -r shard; do
  stem="$(basename "$shard" .faa)"; filtered="$RUN_ROOT/inputs/scan_shards/$stem.faa"
  python "$RUN_ROOT/inputs/filter_hmmsearch_shard.py" "$shard" "$filtered" "$EXCLUSIONS"
  while IFS=$'\t' read -r model source threshold min_cov report_group; do
    [[ "$model" == "model" ]] && continue
    hmm="$RUN_ROOT/inputs/hmms/$model.hmm"; out="$BUILD/$model"__"$stem"
    hmmsearch --tblout "$out.tbl" --domtblout "$out.dom" -E 1e-5 --cpu 1 "$hmm" "$filtered" > /dev/null 2>>"$RUN_ROOT/logs/$model.stderr" || printf '%s\t%s\n' "$model" "$stem" >> "$FAILED"
  done < "$REGISTRY"
done < "$RUN_ROOT/inputs/shard_paths.txt"
[[ ! -s "$FAILED" ]] || exit 1
python "$RUN_ROOT/inputs/06b_aggregate_hits.py" --hmmout "$BUILD" --out "$RUN_ROOT/results/hits_all.tsv"
[[ -s "$RUN_ROOT/results/hits_all.tsv" ]] || exit 1
mv "$BUILD" "$RUN_ROOT/results/hmmsearch"
python - "$RUN_ROOT/results/scan_manifest.json" "$RUN_ROOT" "$REGISTRY" "$SHARD_DIR" <<'PY'
import datetime as dt,hashlib,json,pathlib,subprocess,sys
out,run,registry,shards=map(pathlib.Path,sys.argv[1:])
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for b in iter(lambda:f.read(1<<20),b""): h.update(b)
 return h.hexdigest()
models=[]
for line in registry.read_text(encoding="utf-8").splitlines()[1:]:
 if line.strip():
  model,source,threshold,min_cov,group=line.split("\t"); p=run/"inputs"/"hmms"/(model+".hmm")
  models.append({"model":model,"threshold":threshold,"min_cov":float(min_cov),"report_group":group,"sha256":sha(p),"bytes":p.stat().st_size})
payload={"schema_version":1,"status":"completed","run_id":run.name,"created_utc":dt.datetime.now(dt.timezone.utc).isoformat(),"registry_sha256":sha(run/"inputs"/"formal_scan_models.tsv"),"source_shard_dir":str(shards),"shard_count":len(list(pathlib.Path(shards).glob("shard_*.faa"))),"models":models,"overlength_exclusions_sha256":sha(run/"results"/"overlength_exclusions.tsv"),"hits_all_sha256":sha(run/"results"/"hits_all.tsv"),"hmmsearch_version":subprocess.run(["hmmsearch","-h"],capture_output=True,text=True).stdout.splitlines()[0]}
out.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
PY
echo "formal frozen screen complete: $RUN_ROOT"
