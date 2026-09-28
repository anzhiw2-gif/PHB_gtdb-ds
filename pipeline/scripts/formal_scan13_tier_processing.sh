#!/usr/bin/env bash
# Build a strict tier1-core result from the completed formal scan 13.
#
# Task F4 — HMMER shard database size (Z), same contract as 06_screen.sh
# ----------------------------------------------------------------------
# This entrypoint re-runs HMMER inside 08c_tier_rescore.py, and HMMER E-values
# are E = Z * P for the number of target sequences Z it is told about.  Its
# inputs are extracted subsets, so with no -Z every E-value would silently use
# that subset's own sequence count.  Therefore:
#   * --database-size-z is REQUIRED; without it the script exits 1 before any
#     run directory, copy, or HMMER call exists (no derivation, no default);
#   * the value must be a positive integer (empty/zero/negative/non-numeric all
#     exit non-zero with a clear message);
#   * that one Z is handed to 08c_tier_rescore.py, which routes every
#     hmmsearch through hmmer_command.build_hmmsearch_command (-Z always
#     present, --domZ never emitted);
#   * results/tier_processing_manifest.json records database_size_Z, its stated
#     basis, the parent library's per-shard sequence counts, their total, the
#     exact command template, and domz_used=false.  Missing scale evidence makes
#     the manifest invalid instead of being defaulted.
# Per-task HMMER runs single-threaded (HMM_CPU=1) so outer parallelism owns
# concurrency; --hmm-cpu may still raise it explicitly up to the 60 ceiling.
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd -P)"
DEPLOY_DIR="$(cd "$SCRIPT_DIR/.." && pwd -P)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd -P)"
RUN_ID=""
PARENT_RUN="$PROJECT_ROOT/runs/20260901_formal_frozen_scan_13"
HMM_CPU=1
DATABASE_SIZE_Z="${PHB_DATABASE_SIZE_Z:-}"
DATABASE_SIZE_BASIS="${PHB_DATABASE_SIZE_BASIS:-}"
PYTHON="$HOME/miniconda3/envs/phb_gtdb/bin/python"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-id) RUN_ID="${2:?--run-id requires a value}"; shift 2 ;;
        --parent-run) PARENT_RUN="${2:?--parent-run requires a path}"; shift 2 ;;
        --hmm-cpu) HMM_CPU="${2:?--hmm-cpu requires a value}"; shift 2 ;;
        --database-size-z) DATABASE_SIZE_Z="${2:?--database-size-z requires a value}"; shift 2 ;;
        --database-size-basis) DATABASE_SIZE_BASIS="${2:?--database-size-basis requires a value}"; shift 2 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

# Full-database scale is fail-closed.  This check runs before the run layout is
# created and before any evidence is copied, so a refused invocation leaves no
# half-built run directory behind.
if [ -z "$DATABASE_SIZE_Z" ]; then
    echo "[ERROR] --database-size-z is required: it must be the COMPLETE library" >&2
    echo "        target sequence count Z (the whole parent scan, not one shard or" >&2
    echo "        this extracted subset).  Without it every rescore call would fall" >&2
    echo "        back to HMMER's default Z = that subset's sequence count, so the" >&2
    echo "        tier E-values would not share one scale." >&2
    echo "        Measure or obtain the full-library count first, then rerun with" >&2
    echo "        --database-size-z <positive integer> --database-size-basis <source>." >&2
    echo "        This entrypoint is fail-closed by design (Task F4): there is no default." >&2
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

for name in formal_scan13_tier_processing.sh prepare_formal_scan13_tier.py parallel_extract_sequences.py 07b_extract_seqs.py 08_validate.py 08c_tier_rescore.py hmmer_command.py run_context.py; do
    cp "$SCRIPT_DIR/$name" "$RUN_ROOT/inputs/$name"
done
cp "$PARENT_RUN/inputs/formal_scan_models.tsv" "$RUN_ROOT/inputs/formal_scan_models.tsv"
cp "$PARENT_RUN/results/hits_all.tsv" "$RUN_ROOT/inputs/parent_hits_all.tsv"
cp "$PARENT_RUN/results/scan_manifest.json" "$RUN_ROOT/inputs/parent_scan_manifest.json"
for hmm in ePhaZ_curated_core iPhaZ OH ArchPhaZ_hydrolase; do
    cp "$PARENT_RUN/inputs/hmms/$hmm.hmm" "$RUN_ROOT/inputs/hmms/$hmm.hmm"
done
printf '%s\n' "$PARENT_RUN" > "$RUN_ROOT/inputs/parent_run.txt"

"$PYTHON" - "$RUN_ROOT/inputs/run_context.py" "$RUN_ROOT" "$RUN_ID" "$PARENT_RUN" "$DATABASE_SIZE_Z" "$DATABASE_SIZE_BASIS" <<'PY'
import importlib.util, json, pathlib, sys
script, run, run_id, parent = map(pathlib.Path, sys.argv[:4])
database_size_z, database_size_basis = sys.argv[4], sys.argv[5]
spec = importlib.util.spec_from_file_location("run_context", script)
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
contract = module.write_input_contract(run, run_id=str(run_id), gtdb_inputs={
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
    "hmmer_command": run/"inputs/hmmer_command.py",
    "run_context": run/"inputs/run_context.py",
    "ePhaZ_curated_core_hmm": run/"inputs/hmms/ePhaZ_curated_core.hmm",
    "iPhaZ_hmm": run/"inputs/hmms/iPhaZ.hmm",
    "OH_hmm": run/"inputs/hmms/OH.hmm",
    "ArchPhaZ_hydrolase_hmm": run/"inputs/hmms/ArchPhaZ_hydrolase.hmm",
})
# The full-database scale is part of the run contract, not a side note: fold it
# into the contract file that carries the SHA-256 of every bound input, so the
# one Z this rerun used is covered by the manifest hash that follows.
contract["database_size_Z"] = int(database_size_z)
contract["database_size_basis"] = database_size_basis
(run/"input_contract.json").write_text(
    json.dumps(contract, ensure_ascii=False, indent=2) + "\n",
    encoding="utf-8",
)
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
    "$PYTHON" "$RUN_ROOT/inputs/08c_tier_rescore.py" --database-size-z "$DATABASE_SIZE_Z" --database-size-basis "$DATABASE_SIZE_BASIS" --cpu "$HMM_CPU"
) | tee "$RUN_ROOT/logs/04_tier_rescore.log"

"$PYTHON" "$RUN_ROOT/inputs/08c_tier_rescore.py" --validate-build "$BUILD/data/screen/tiers" \
    --families "ePhaZ iPhaZ OH ArchPhaZ_hydrolase"

mv "$BUILD" "$RUN_ROOT/results/tier_processing"
"$PYTHON" - "$RUN_ROOT" "$PARENT_RUN" "$HMM_CPU" "$DATABASE_SIZE_Z" "$DATABASE_SIZE_BASIS" <<'PY'
import datetime, hashlib, json, pathlib, sys

# The full-library command shape this rerun used.  Kept verbatim identical to
# pipeline/scripts/hmmer_command.py::COMMAND_TEMPLATE (and 06_screen.sh) so a
# later reader can see the exact flags that produced the rescore E-values.
command_template = ('hmmsearch --tblout "$out" --domtblout "${out%.tbl}.dom" '
                    '-Z "$DATABASE_SIZE_Z" -E "$EVAL" --cpu "$HMM_CPU" "$hmm" "$shard"')
# Fail-closed scale fields: these are the fields a tier manifest must carry.  A
# missing or unusable one makes the manifest invalid rather than defaulted.
required_scale_fields = (
    "database_size_Z",
    "database_size_basis",
    "shards",
    "shard_sequence_total",
    "hmmsearch_command_template",
    "domz_used",
)


def count_fasta_records(path):
    """Count FASTA records by header lines (the number HMMER uses as Z)."""
    count = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""): h.update(block)
    return h.hexdigest()


def require_database_size_scale(shard_dir, declared_z, basis):
    """Prove the manifest's scale from the parent library's real shard counts.

    Every rescore call used ``declared_z`` for the *whole* parent scan, so the
    measured per-shard record counts must sum to exactly that number.  Anything
    else — no shards, absent Z, blank basis, or a total that contradicts the
    declared scale — exits non-zero instead of writing a manifest whose scale
    cannot be checked.
    """
    if declared_z is None or (isinstance(declared_z, str) and not declared_z.strip()):
        raise SystemExit("[ERROR] database_size_Z is missing: refusing to write a "
                         "tier manifest with no full-library scale")
    if isinstance(declared_z, str):
        if not declared_z.strip().isdigit() or int(declared_z) < 1:
            raise SystemExit(f"[ERROR] database_size_Z must be a positive integer, "
                             f"got {declared_z!r}")
        declared_z = int(declared_z)
    if isinstance(declared_z, bool) or not isinstance(declared_z, int) or declared_z < 1:
        raise SystemExit(f"[ERROR] database_size_Z must be a positive integer, "
                         f"got {declared_z!r}")
    if not isinstance(basis, str) or not basis.strip():
        raise SystemExit("[ERROR] database_size_basis is missing: refusing to write a "
                         "tier manifest that cannot say where Z came from")
    shard_paths = sorted(pathlib.Path(shard_dir).glob("shard_*.faa"))
    if not shard_paths:
        raise SystemExit(f"[ERROR] no shard_*.faa under {shard_dir}: the per-shard "
                         "sequence counts are required to prove the scale")
    shards = []
    for path in shard_paths:
        count = count_fasta_records(path)
        if count < 1:
            raise SystemExit(f"[ERROR] shard {path.name} has no FASTA records")
        shards.append({"name": path.name, "sha256": sha(path),
                       "sequence_count": count})
    total = sum(item["sequence_count"] for item in shards)
    if total != declared_z:
        raise SystemExit(
            f"[ERROR] shard_sequence_total {total} does not match the declared "
            f"database_size_Z {declared_z}: the stated scale contradicts the "
            "measured per-shard counts")
    return shards, total


def validate_tier_manifest(manifest):
    """Refuse to publish a tier manifest without the full scale evidence."""
    missing = [field for field in required_scale_fields if field not in manifest]
    if missing:
        raise SystemExit("[ERROR] tier manifest is missing the database-size "
                         f"fields: {', '.join(missing)}")
    if manifest["domz_used"] is not False:
        raise SystemExit("[ERROR] domz_used must be false: no downstream consumer "
                         "rescales domain E-values")
    if "-Z" not in manifest["hmmsearch_command_template"]:
        raise SystemExit("[ERROR] hmmsearch_command_template must pin -Z")
    if "--domZ" in manifest["hmmsearch_command_template"]:
        raise SystemExit("[ERROR] hmmsearch_command_template must not set --domZ")
    if manifest["shard_sequence_total"] != manifest["database_size_Z"]:
        raise SystemExit("[ERROR] shard_sequence_total does not match database_size_Z")
    return manifest


def build_tier_processing_manifest(run, parent, cpu, declared_z, basis):
    """Build, prove, and return the tier manifest — or exit non-zero.

    The parent scan-13 library shards are the target space every rescore call's
    Z describes, so their measured record counts are the evidence for the scale.
    Nothing is written until the manifest validates.
    """
    shards, shard_sequence_total = require_database_size_scale(
        str(parent/"inputs"/"scan_shards"), declared_z, basis)
    tree = run/"results"/"tier_processing"
    manifest = {
        "schema_version": 1, "status": "completed", "run_id": run.name,
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "parent_run": str(parent), "hmm_cpu": cpu,
        "database_size_Z": int(declared_z),
        "database_size_basis": basis,
        "database_size_source": str(parent/"inputs"/"scan_shards"),
        "shards": shards,
        "shard_sequence_total": shard_sequence_total,
        "hmmsearch_command_template": command_template,
        "domz_used": False,
        "input_contract_sha256": sha(run/"input_contract.json"),
        "hmmer_command_sha256": sha(run/"inputs"/"hmmer_command.py"),
        "tier_processing_summary_sha256": sha(tree/"screen"/"summary.txt"),
        "tier1_counts": {p.stem.replace("_tier1", ""): sum(1 for line in p.open() if line.startswith(">"))
                         for p in (tree/"data"/"screen"/"tiers").glob("*_tier1.faa")},
    }
    return validate_tier_manifest(manifest)


def main(argv):
    run, parent = map(pathlib.Path, argv[1:3])
    cpu = int(argv[3])
    # A missing Z or basis reaches the gate as an empty string and is refused
    # there; there is no default and no derivation.
    manifest = build_tier_processing_manifest(run, parent, cpu, argv[4], argv[5])
    out = run/"results"/"tier_processing_manifest.json"
    out.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
                   encoding="utf-8")
    return manifest


# The ``__name__`` guard matters: this program is also executed as a library by
# pipeline/tests/test_hmmer_database_size.py, which must be able to call the
# scale gate without the module-level argparse-style unpacking running.
if __name__ == "__main__":
    main(sys.argv)
PY
echo "completed tier processing: $RUN_ROOT"
echo "database scale bound: Z=$DATABASE_SIZE_Z ($DATABASE_SIZE_BASIS), domZ=false"
