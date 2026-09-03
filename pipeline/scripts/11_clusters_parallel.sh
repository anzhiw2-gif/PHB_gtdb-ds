#!/bin/bash
set -Eeuo pipefail
RUN_ROOT="${1:?usage: 11_clusters_parallel.sh RUN_ROOT [BATCHES]}"
BATCHES="${2:-80}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PHB_PYTHON:-python}"
GTDB_ROOT="${PHB_GTDB_ROOT:?set PHB_GTDB_ROOT}"
cd "$RUN_ROOT"
MARKERS="PhaC,PhaE,PhaJ,BdhA,phasin,PHA_gran_rgn"
BATCH_ROOT="$RUN_ROOT/data/cluster_batches"
mkdir -p "$BATCH_ROOT"
if [ ! -s "$BATCH_ROOT/batches.tsv" ]; then
  "$PY" "$SCRIPT_DIR/11_split_hits.py" --hits data/screen/hits_filtered.tsv --outdir "$BATCH_ROOT" --batches "$BATCHES"
fi
while IFS=$'\t' read -r batch genomes loci; do
  [ "$batch" = "batch" ] && continue
  root="$BATCH_ROOT/$batch"
  mkdir -p "$root/results/logs" "$root/results/tables" "$root/data/work"
  if [ ! -s "$root/rc" ] || ! grep -qx '0' "$root/rc"; then
    rm -f "$root/rc"
    setsid -f bash -c "$PY $SCRIPT_DIR/11_clusters.py --hits $root/hits.tsv --marker-hmms data/hmms/v2 --marker-families $MARKERS --gtdb $GTDB_ROOT --outdir $root/results --workdir $root/data/work --flank-kb 10 --threads 1 --max-genomes 0 > $root/results/logs/cluster.log 2>&1; printf '%s\\n' \$? > $root/rc" </dev/null
  fi
done < "$BATCH_ROOT/batches.tsv"
while true; do
  done_count=$(find "$BATCH_ROOT" -mindepth 2 -maxdepth 2 -name rc | wc -l)
  [ "$done_count" -ge "$BATCHES" ] && break
  sleep 30
done
if find "$BATCH_ROOT" -mindepth 2 -maxdepth 2 -name rc -exec grep -L '^0$' {} + | grep -q .; then
  echo "one or more cluster batches failed" >&2
  exit 1
fi
"$PY" "$SCRIPT_DIR/11_merge_cluster_outputs.py" --batch-root "$BATCH_ROOT" --outdir results --marker-families "$MARKERS"
