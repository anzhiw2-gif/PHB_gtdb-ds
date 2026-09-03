#!/usr/bin/env bash
# Read-only progress monitor for a dated formal scan run.
set -Eeuo pipefail

RUN_DIR=""
INTERVAL=30
ONCE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-dir) RUN_DIR="${2:?--run-dir requires a path}"; shift 2 ;;
        --interval) INTERVAL="${2:?--interval requires seconds}"; shift 2 ;;
        --once) ONCE=1; shift ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

[[ -n "$RUN_DIR" && -d "$RUN_DIR" ]] || { echo "--run-dir must be an existing directory" >&2; exit 1; }
[[ "$INTERVAL" =~ ^[1-9][0-9]*$ ]] || { echo "--interval must be a positive integer" >&2; exit 2; }
RUN_DIR="$(cd "$RUN_DIR" && pwd -P)"

snapshot() {
    local registry="$RUN_DIR/inputs/formal_scan_models.tsv"
    local shards="$RUN_DIR/inputs/shard_paths.txt"
    local status_dir="$RUN_DIR/results/task_status"
    local errors="$RUN_DIR/logs/failed_tasks.tsv"
    local tasks="$RUN_DIR/inputs/tasks.tsv"
    local expected=0 completed=0 failed=0 tbl=0 dom=0 pending=0
    local latest="unknown" manifest="missing" hits="missing" parent="unknown"

    if [[ -s "$registry" && -s "$shards" ]]; then
        expected=$(( $(awk 'NR > 1 && $0 !~ /^[[:space:]]*$/ { n++ } END { print n+0 }' "$registry") * $(wc -l < "$shards") ))
    fi
    [[ -d "$status_dir" ]] && completed="$(find "$status_dir" -maxdepth 1 -type f -name '*.ok' | wc -l)"
    [[ -f "$errors" ]] && failed="$(awk 'NR > 1 && $0 !~ /^[[:space:]]*$/ { n++ } END { print n+0 }' "$errors")"
    [[ -f "$tasks" ]] && pending="$(awk 'NR > 1 && $0 !~ /^[[:space:]]*$/ { n++ } END { print n+0 }' "$tasks")"
    [[ -d "$RUN_DIR/results/hmmsearch.build" ]] && tbl="$(find "$RUN_DIR/results/hmmsearch.build" -maxdepth 1 -type f -name '*.tbl' -size +0c | wc -l)"
    [[ -d "$RUN_DIR/results/hmmsearch.build" ]] && dom="$(find "$RUN_DIR/results/hmmsearch.build" -maxdepth 1 -type f -name '*.dom' -size +0c | wc -l)"
    [[ -s "$RUN_DIR/results/scan_manifest.json" ]] && manifest="present"
    [[ -s "$RUN_DIR/results/hits_all.tsv" ]] && hits="present"
    [[ -s "$RUN_DIR/inputs/resume_parent.txt" ]] && parent="$(tr '\n' ' ' < "$RUN_DIR/inputs/resume_parent.txt")"
    latest="$(find "$RUN_DIR" -type f -printf '%T@ %p\n' 2>/dev/null | sort -nr | head -n 1 | cut -d' ' -f2- || true)"

    printf '\n[%s] run=%s\n' "$(date -Is)" "$(basename "$RUN_DIR")"
    printf 'progress: completed=%s expected=%s pending_new=%s failed=%s\n' "$completed" "$expected" "$pending" "$failed"
    printf 'outputs: nonempty_tbl=%s nonempty_dom=%s hits_all=%s manifest=%s\n' "$tbl" "$dom" "$hits" "$manifest"
    printf 'parent: %s\n' "$parent"
    printf 'latest_file: %s\n' "$latest"
    printf 'disk: %s\n' "$(du -sh "$RUN_DIR" 2>/dev/null | awk '{print $1}')"
    printf 'hmmsearch_processes:\n'
    pgrep -af hmmsearch || printf '  none\n'
}

if (( ONCE )); then
    snapshot
else
    while true; do
        snapshot
        sleep "$INTERVAL"
    done
fi
