#!/usr/bin/env python3
"""Create and validate a complete HMM-family by shard screening contract."""
import argparse
import hashlib
import json
import os
import sys


class ScreenManifestError(RuntimeError):
    """Raised when HMMER outputs do not form a complete, bound task matrix."""


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_manifest(manifest):
    required = ("families", "shards", "tasks")
    missing = [field for field in required if field not in manifest]
    if missing:
        raise ScreenManifestError(f"missing manifest fields: {', '.join(missing)}")
    families = manifest["families"]
    shards = manifest["shards"]
    tasks = manifest["tasks"]
    if not isinstance(families, list) or not families or len(set(families)) != len(families):
        raise ScreenManifestError("families must be a unique non-empty list")
    if not isinstance(shards, list) or not shards:
        raise ScreenManifestError("shards must be a non-empty list")
    shard_hashes = {item.get("name"): item.get("sha256") for item in shards if isinstance(item, dict)}
    if len(shard_hashes) != len(shards) or any(not name or not digest for name, digest in shard_hashes.items()):
        raise ScreenManifestError("shard names and hashes must be unique and non-empty")
    if not isinstance(tasks, list):
        raise ScreenManifestError("tasks must be a list")
    expected = {(family, name) for family in families for name in shard_hashes}
    observed = set()
    for task in tasks:
        if not isinstance(task, dict):
            raise ScreenManifestError("task record must be an object")
        key = (task.get("family"), task.get("shard"))
        if key in observed or key not in expected:
            raise ScreenManifestError(f"invalid or duplicate task: {key}")
        observed.add(key)
        if task.get("input_sha256") != shard_hashes[key[1]]:
            raise ScreenManifestError(f"input hash does not match declared shard: {key}")
        if not all(task.get(field) for field in ("hmm_sha256", "tbl_sha256", "dom_sha256", "evalue")):
            raise ScreenManifestError(f"missing task provenance: {key}")
    if observed != expected:
        raise ScreenManifestError(f"incomplete family-shard matrix: expected={len(expected)}, observed={len(observed)}")
    return manifest


def build_manifest(shard_dir, hmm_dir, hmmout, families, evalue):
    shard_paths = sorted(
        os.path.join(shard_dir, name) for name in os.listdir(shard_dir)
        if name.startswith("shard_") and name.endswith(".faa")
    )
    if not shard_paths:
        raise ScreenManifestError(f"no shards found: {shard_dir}")
    shards = [{"name": os.path.basename(path), "sha256": sha256(path)} for path in shard_paths]
    tasks = []
    for family in families:
        hmm_path = os.path.join(hmm_dir, f"{family}.hmm")
        if not os.path.isfile(hmm_path):
            raise ScreenManifestError(f"missing HMM: {hmm_path}")
        hmm_hash = sha256(hmm_path)
        for shard in shards:
            stem = os.path.splitext(shard["name"])[0]
            tbl_path = os.path.join(hmmout, f"{family}__{stem}.tbl")
            dom_path = os.path.join(hmmout, f"{family}__{stem}.dom")
            if not os.path.isfile(tbl_path) or not os.path.isfile(dom_path):
                raise ScreenManifestError(f"missing HMMER output: {family} x {shard['name']}")
            tasks.append({
                "family": family,
                "shard": shard["name"],
                "input_sha256": shard["sha256"],
                "hmm_sha256": hmm_hash,
                "tbl_sha256": sha256(tbl_path),
                "dom_sha256": sha256(dom_path),
                "evalue": str(evalue),
            })
    manifest = {"schema_version": 1, "families": families, "shards": shards, "tasks": tasks}
    validate_manifest(manifest)
    return manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--validate", help="validate an existing JSON manifest without rescanning input files")
    parser.add_argument("--shard-dir", default="data/proteins/shards_filt")
    parser.add_argument("--hmm-dir", default="data/hmms/v2")
    parser.add_argument("--hmmout", default="data/screen/hmmsearch")
    parser.add_argument("--families", default="", help="space-separated HMM family names")
    parser.add_argument("--eval", default="")
    parser.add_argument("--out", default="data/screen/screen_manifest.json")
    args = parser.parse_args()
    if args.validate:
        with open(args.validate, encoding="utf-8-sig") as handle:
            manifest = json.load(handle)
        validate_manifest(manifest)
        print(f"screen manifest already verified: {len(manifest['tasks'])} family-shard tasks")
        return
    if not args.families or not args.eval:
        parser.error("--families and --eval are required when building a manifest")
    families = args.families.split()
    manifest = build_manifest(args.shard_dir, args.hmm_dir, args.hmmout, families, args.eval)
    with open(args.out, "w", encoding="utf-8", newline="\n") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    print(f"screen manifest verified: {len(manifest['tasks'])} family-shard tasks")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, ScreenManifestError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
