#!/usr/bin/env bash
# Build a strict tier1-core result from the completed formal scan 13.
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd -P)"
DEPLOY_DIR="$(cd "$SCRIPT_DIR/.." && pwd -P)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd -P)"
RUN_ID=""
PARENT_RUN="$PROJECT_ROOT/runs/20260901_formal_frozen_scan_13"
HMM_CPU=60
PYTHON="$HOME/miniconda3/envs/phb_gtdb/bin/python"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-id) RUN_ID="${2:?--run-id requires a value}"; shift 2 ;;
        --parent-run) PARENT_RUN="${2:?--parent-run requires a path}"; shift 2 ;;
        --hmm-cpu) HMM_CPU="${2:?--hmm-cpu requires a value}"; shift 2 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

[[ "$RUN_ID" =~ ^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$ && "$RUN_ID" != *..* ]] || {
    echo "--run-id must be a safe dated identifier" >&2; exit 2;
}
[[ "$HMM_CPU" =~ ^[1-9][0-9]*$ && "$HMM_CPU" -le 60 ]] || {
    echo "--hmm-cpu must be 1..60" >&2; exit 2;
}
[[ -x "$PYTHON" ]] || { echo "missing phb_gtdb Python: $PYTHON" >&2; exit 1; }
[[ -s "$PARENT_RUN/results/hits_all.tsv" && -s "$PARENT_RUN/input_contract.json" ]] || {
    echo "parent run lacks required accepted hit evidence" >&2; exit 1;
}

"$PYTHON" - "$SCRIPT_DIR/run_context.py" "$PROJECT_ROOT" "$RUN_ID" <<'PY'
import importlib.util, pathlib, sys
script, root, run_id = map(pathlib.Path, sys.argv[1:])
spec = importlib.util.spec_from_file_location("run_context", script)
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
module.create_run_layout(root, str(run_id))
PY

RUN_ROOT="$PROJECT_ROOT/runs/$RUN_ID"
BUILD="$RUN_ROOT/results/tier_processing.build"
mkdir -p "$RUN_ROOT/inputs/hmms" "$BUILD"

for name in formal_scan13_tier_processing.sh prepare_formal_scan13_tier.py parallel_extract_sequences.py 07b_extract_seqs.py 08_validate.py 08c_tier_rescore.py run_context.py; do
    cp "$SCRIPT_DIR/$name" "$RUN_ROOT/inputs/$name"
done
cp "$PARENT_RUN/inputs/formal_scan_models.tsv" "$RUN_ROOT/inputs/formal_scan_models.tsv"
cp "$PARENT_RUN/results/hits_all.tsv" "$RUN_ROOT/inputs/parent_hits_all.tsv"
cp "$PARENT_RUN/results/scan_manifest.json" "$RUN_ROOT/inputs/parent_scan_manifest.json"
for hmm in ePhaZ_curated_core iPhaZ OH ArchPhaZ_hydrolase; do
    cp "$PARENT_RUN/inputs/hmms/$hmm.hmm" "$RUN_ROOT/inputs/hmms/$hmm.hmm"
done
printf '%s\n' "$PARENT_RUN" > "$RUN_ROOT/inputs/parent_run.txt"

"$PYTHON" - "$RUN_ROOT/inputs/run_context.py" "$RUN_ROOT" "$RUN_ID" "$PARENT_RUN" <<'PY'
import importlib.util, pathlib, sys
script, run, run_id, parent = map(pathlib.Path, sys.argv[1:])
spec = importlib.util.spec_from_file_location("run_context", script)
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
module.write_input_contract(run, run_id=str(run_id), gtdb_inputs={
    "taxonomy": pathlib.Path.home()/"GTDB/taxonomy/bac120_taxonomy_r232.tsv",
    "metadata": pathlib.Path.home()/"GTDB/metadata/bac120_metadata_r232.tsv.gz",
    "tree": pathlib.Path.home()/"GTDB/GTDB_tree/bac120_r232.tree",
}, inputs={
    "parent_hits_all": run/"inputs/parent_hits_all.tsv",
    "parent_scan_manifest": run/"inputs/parent_scan_manifest.json",
    "registry": run/"inputs/formal_scan_models.tsv",
    "tier_driver": run/"inputs/formal_scan13_tier_processing.sh",
    "tier_prepare": run/"inputs/prepare_formal_scan13_tier.py",
    "extract_sequences": run/"inputs/07b_extract_seqs.py",
    "validate_sequences": run/"inputs/08_validate.py",
    "tier_rescore": run/"inputs/08c_tier_rescore.py",
    "run_context": run/"inputs/run_context.py",
    "ePhaZ_curated_core_hmm": run/"inputs/hmms/ePhaZ_curated_core.hmm",
    "iPhaZ_hmm": run/"inputs/hmms/iPhaZ.hmm",
    "OH_hmm": run/"inputs/hmms/OH.hmm",
    "ArchPhaZ_hydrolase_hmm": run/"inputs/hmms/ArchPhaZ_hydrolase.hmm",
})
PY

"$PYTHON" "$RUN_ROOT/inputs/prepare_formal_scan13_tier.py" \
    --hits "$RUN_ROOT/inputs/parent_hits_all.tsv" \
    --registry "$RUN_ROOT/inputs/formal_scan_models.tsv" \
    --outdir "$BUILD/screen" | tee "$RUN_ROOT/logs/01_prepare.log"
[[ -s "$BUILD/screen/broad_discovery.tsv" ]] || { echo "missing broad_discovery.tsv" >&2; exit 1; }

mkdir -p "$BUILD/shard_extract"
find "$PARENT_RUN/inputs/scan_shards" -maxdepth 1 -type f -name 'shard_*.faa' | sort | \
    parallel -j 20 --halt soon,fail=1 --joblog "$RUN_ROOT/logs/extract.joblog" \
    "$PYTHON" "$RUN_ROOT/inputs/parallel_extract_sequences.py" --shard {} \
    --ids "$BUILD/screen/unique_proteins.txt" --out "$BUILD/shard_extract/{/.}.faa" \
    > "$RUN_ROOT/logs/02_extract.log"

"$PYTHON" - "$BUILD/screen/hits_filtered.tsv" "$BUILD/shard_extract" "$BUILD/family_seqs" <<'PY'
import csv, pathlib, sys
hits, shard_dir, out_dir = map(pathlib.Path, sys.argv[1:])
out_dir.mkdir(parents=True, exist_ok=True)
protein_family = {}
with hits.open() as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        protein_family[row["protein"]] = row["family"]
handles = {}
try:
    for shard in sorted(shard_dir.glob("shard_*.faa")):
        current = None; seq = []
        def flush():
            if current in protein_family:
                fam = protein_family[current]
                if fam not in handles:
                    handles[fam] = (out_dir / f"{fam}.faa").open("w")
                handles[fam].write(f">{current}\n{''.join(seq)}\n")
        for line in shard.open():
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush(); current = line[1:].split()[0]; seq = []
            else:
                seq.append(line.strip())
        flush()
finally:
    for handle in handles.values(): handle.close()
PY

expected="$(wc -l < "$BUILD/screen/unique_proteins.txt")"
observed="$(grep -h -c '^>' "$BUILD/family_seqs"/*.faa | awk -F: '{s+=$NF} END{print s+0}')"
[[ "$expected" -eq "$observed" ]] || { echo "extraction mismatch: $observed/$expected" >&2; exit 1; }

"$PYTHON" "$RUN_ROOT/inputs/08_validate.py" \
    --indir "$BUILD/family_seqs" --outdir "$BUILD/validation" --signalp 0 \
    | tee "$RUN_ROOT/logs/03_validate.log"

mkdir -p "$BUILD/data/screen" "$BUILD/data/hmms/v2"
ln -s "$BUILD/family_seqs" "$BUILD/data/screen/family_seqs"
cp "$RUN_ROOT/inputs/hmms/ePhaZ_curated_core.hmm" "$BUILD/data/hmms/ePhaZ.hmm"
cp "$RUN_ROOT/inputs/hmms/iPhaZ.hmm" "$BUILD/data/hmms/iPhaZ.hmm"
cp "$RUN_ROOT/inputs/hmms/OH.hmm" "$BUILD/data/hmms/OH.hmm"
cp "$RUN_ROOT/inputs/hmms/ArchPhaZ_hydrolase.hmm" "$BUILD/data/hmms/v2/ArchPhaZ_hydrolase.hmm"
(
    cd "$BUILD"
    "$PYTHON" "$RUN_ROOT/inputs/08c_tier_rescore.py" --cpu "$HMM_CPU"
) | tee "$RUN_ROOT/logs/04_tier_rescore.log"

"$PYTHON" "$RUN_ROOT/inputs/08c_tier_rescore.py" --validate-build "$BUILD/data/screen/tiers" \
    --families "ePhaZ iPhaZ OH ArchPhaZ_hydrolase"

mv "$BUILD" "$RUN_ROOT/results/tier_processing"
"$PYTHON" - "$RUN_ROOT" "$PARENT_RUN" "$HMM_CPU" <<'PY'
import datetime, hashlib, json, pathlib, sys
run, parent = map(pathlib.Path, sys.argv[1:3])
cpu = int(sys.argv[3])
def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""): h.update(block)
    return h.hexdigest()
out = run/"results"/"tier_processing_manifest.json"
tree = run/"results"/"tier_processing"
out.write_text(json.dumps({
    "schema_version": 1, "status": "completed", "run_id": run.name,
    "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "parent_run": str(parent), "hmm_cpu": cpu,
    "input_contract_sha256": sha(run/"input_contract.json"),
    "tier_processing_summary_sha256": sha(tree/"screen"/"summary.txt"),
    "tier1_counts": {p.stem.replace("_tier1", ""): sum(1 for line in p.open() if line.startswith(">"))
                     for p in (tree/"data"/"screen"/"tiers").glob("*_tier1.faa")},
}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
PY
echo "completed tier processing: $RUN_ROOT"
