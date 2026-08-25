#!/usr/bin/env python3
"""Create and validate the filtered-protein shard provenance contract."""
import argparse
import hashlib
import json
import os
import sys


class FilterManifestError(RuntimeError):
    """Raised when filtered shards cannot be bound to their source shards."""


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_manifest(manifest):
    required = ("source_shards", "filtered_shards", "max_aa")
    missing = [field for field in required if field not in manifest]
    if missing:
        raise FilterManifestError(f"missing manifest fields: {', '.join(missing)}")
    if not isinstance(manifest["max_aa"], int) or manifest["max_aa"] <= 0:
        raise FilterManifestError("max_aa must be a positive integer")
    source = manifest["source_shards"]
    filtered = manifest["filtered_shards"]
    if not isinstance(source, list) or not source:
        raise FilterManifestError("source_shards must be a non-empty list")
    if not isinstance(filtered, list) or len(filtered) != len(source):
        raise FilterManifestError("filtered_shards must cover every source shard exactly once")
    source_by_name = {item.get("name"): item for item in source if isinstance(item, dict)}
    filtered_by_name = {item.get("name"): item for item in filtered if isinstance(item, dict)}
    if len(source_by_name) != len(source) or any(not name or not item.get("sha256") for name, item in source_by_name.items()):
        raise FilterManifestError("source shard names and hashes must be unique and non-empty")
    if set(source_by_name) != set(filtered_by_name):
        raise FilterManifestError("filtered shard names do not match source shard names")
    for name, item in filtered_by_name.items():
        if not item.get("sha256") or item.get("source_sha256") != source_by_name[name]["sha256"]:
            raise FilterManifestError(f"invalid filtered shard provenance: {name}")
    return manifest


def build_manifest(source_dir, filtered_dir, max_aa, filter_script):
    source_paths = sorted(
        os.path.join(source_dir, name) for name in os.listdir(source_dir)
        if name.startswith("shard_") and name.endswith(".faa")
    )
    if not source_paths:
        raise FilterManifestError(f"no source shards found: {source_dir}")
    source_shards = []
    filtered_shards = []
    for source_path in source_paths:
        name = os.path.basename(source_path)
        filtered_path = os.path.join(filtered_dir, name)
        if not os.path.isfile(filtered_path):
            raise FilterManifestError(f"missing filtered shard: {filtered_path}")
        source_hash = sha256(source_path)
        source_shards.append({"name": name, "sha256": source_hash})
        filtered_shards.append({
            "name": name,
            "source_sha256": source_hash,
            "sha256": sha256(filtered_path),
        })
    manifest = {
        "schema_version": 1,
        "max_aa": max_aa,
        "filter_script_sha256": sha256(filter_script),
        "source_shards": source_shards,
        "filtered_shards": filtered_shards,
    }
    validate_manifest(manifest)
    return manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", default="data/proteins/shards")
    parser.add_argument("--filtered-dir", default="data/proteins/shards_filt")
    parser.add_argument("--max-aa", type=int, required=True)
    parser.add_argument("--filter-script", required=True)
    parser.add_argument("--out", default="data/proteins/filter_manifest.json")
    args = parser.parse_args()
    manifest = build_manifest(args.source_dir, args.filtered_dir, args.max_aa, args.filter_script)
    with open(args.out, "w", encoding="utf-8", newline="\n") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    print(f"filter manifest verified: {len(manifest['source_shards'])} shards")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, FilterManifestError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
