#!/usr/bin/env python3
"""Validate the complete protein-prediction contract before downstream screening."""
import argparse
import hashlib
import json
import os
import sys


class PredictionManifestError(RuntimeError):
    """Raised when predicted genomes or shards do not cover the declared GTDB input."""


def validate_manifest(manifest):
    required = ("total_genomes", "predicted_genomes", "failed_genomes", "expected_shards", "expected_shard_names", "shards")
    missing = [key for key in required if key not in manifest]
    if missing:
        raise PredictionManifestError(f"missing manifest fields: {', '.join(missing)}")
    total = manifest["total_genomes"]
    predicted = manifest["predicted_genomes"]
    failed = manifest["failed_genomes"]
    expected_shards = manifest["expected_shards"]
    expected_shard_names = manifest["expected_shard_names"]
    shards = manifest["shards"]
    if not isinstance(total, int) or total <= 0:
        raise PredictionManifestError(f"invalid total_genomes: {total!r}")
    if predicted != total or failed != 0:
        raise PredictionManifestError(
            f"incomplete prediction: total={total}, predicted={predicted}, failed={failed}"
        )
    if not isinstance(expected_shards, int) or expected_shards <= 0:
        raise PredictionManifestError(f"invalid expected_shards: {expected_shards!r}")
    if (
        not isinstance(expected_shard_names, list)
        or len(expected_shard_names) != expected_shards
        or any(not isinstance(name, str) or not name for name in expected_shard_names)
        or len(set(expected_shard_names)) != len(expected_shard_names)
    ):
        raise PredictionManifestError("expected_shard_names must be a unique non-empty list")
    if not isinstance(shards, list) or len(shards) != expected_shards:
        raise PredictionManifestError(
            f"incomplete shard set: expected={expected_shards}, observed={len(shards) if isinstance(shards, list) else 'invalid'}"
        )
    names = [item.get("name") for item in shards if isinstance(item, dict)]
    if len(names) != len(shards) or len(set(names)) != len(names) or any(not name for name in names):
        raise PredictionManifestError("shard names must be present and unique")
    if set(names) != set(expected_shard_names):
        raise PredictionManifestError("observed shard names do not match expected_shard_names")
    for item in shards:
        if not item.get("sha256") or not isinstance(item.get("genomes"), int) or item["genomes"] <= 0:
            raise PredictionManifestError(f"invalid shard record: {item!r}")
    if sum(item["genomes"] for item in shards) != total:
        raise PredictionManifestError("shard genome counts do not sum to total_genomes")
    return manifest


def validate_shard_files(manifest, shard_dir):
    for item in manifest["shards"]:
        path = os.path.join(shard_dir, item["name"])
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            raise PredictionManifestError(f"missing or empty shard: {path}")
        digest = hashlib.sha256()
        with open(path, "rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        if digest.hexdigest() != item["sha256"]:
            raise PredictionManifestError(f"shard hash mismatch: {path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="data/proteins/prediction_manifest.json")
    parser.add_argument("--shard-dir", default="data/proteins/shards")
    args = parser.parse_args()
    with open(args.manifest, encoding="utf-8") as handle:
        manifest = json.load(handle)
    validate_manifest(manifest)
    validate_shard_files(manifest, args.shard_dir)
    print(f"prediction manifest verified: {manifest['total_genomes']} genomes, {manifest['expected_shards']} shards")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, PredictionManifestError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
