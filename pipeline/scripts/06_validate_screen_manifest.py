#!/usr/bin/env python3
"""Create and validate a complete HMM-family by shard screening contract.

Task 9 — the database-size (``Z``) scale is part of the contract
----------------------------------------------------------------
HMMER E-values are ``E = Z * P``, so a sharded search that omits ``-Z`` reports
every shard's E-values on that shard's *own* sequence count.  Such numbers are
not comparable across shards and cannot be presented as one full-library scale.

This validator therefore binds the scale into the manifest:

* ``database_size_Z`` — the one full-library ``Z`` passed to every shard call;
* ``database_size_basis`` — where that number came from (never inferred here);
* ``shard_sequence_total`` — the sum of the per-shard counts;
* per-shard ``sequence_count`` — the real FASTA record count of each shard;
* ``hmmsearch_command_template`` — the exact command shape used;
* per-task ``database_size_Z`` — the value each task actually ran with.

Validation fails when shard runs used different ``Z`` values, when a shard's
sequence count contradicts the declared total, and when ``database_size_Z`` is
missing: it is **never** defaulted to a shard-local count.
"""
import argparse
import hashlib
import json
import os
import sys


class ScreenManifestError(RuntimeError):
    """Raised when HMMER outputs do not form a complete, bound task matrix."""

# The exact command shape the screening entrypoint must use.  Kept in sync with
# pipeline/scripts/hmmer_command.py::COMMAND_TEMPLATE and asserted by
# pipeline/tests/test_hmmer_database_size.py.  ``--domZ`` is deliberately
# absent: no downstream consumer rescales domain E-values.
COMMAND_TEMPLATE = (
    'hmmsearch --tblout "$out" --domtblout "${out%.tbl}.dom" '
    '-Z "$DATABASE_SIZE_Z" -E "$EVAL" --cpu "$HMM_CPU" "$hmm" "$shard"'
)
DATABASE_SIZE_FIELDS = (
    "database_size_Z",
    "database_size_basis",
    "shard_sequence_total",
    "hmmsearch_command_template",
)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def count_fasta_records(path):
    """Count FASTA records by header lines (the number HMMER uses as ``Z``)."""
    count = 0
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def validate_database_size_scale(manifest):
    """Prove that every shard task shares one full-database ``Z``.

    Fails closed: a missing ``database_size_Z`` is an error, never a default.
    """
    missing = [field for field in DATABASE_SIZE_FIELDS if field not in manifest]
    if missing:
        raise ScreenManifestError(
            "missing database-size fields: "
            + ", ".join(missing)
            + " (refusing to default Z to a shard-local sequence count)"
        )
    declared = manifest["database_size_Z"]
    if isinstance(declared, bool) or not isinstance(declared, int) or declared < 1:
        raise ScreenManifestError(
            f"database_size_Z must be a positive integer, got {declared!r}"
        )
    basis = manifest["database_size_basis"]
    if not isinstance(basis, str) or not basis.strip():
        raise ScreenManifestError("database_size_basis must be a non-empty string")
    template = manifest["hmmsearch_command_template"]
    if not isinstance(template, str) or not template.strip():
        raise ScreenManifestError(
            "hmmsearch_command_template must be a non-empty string"
        )
    if "-Z" not in template:
        raise ScreenManifestError("hmmsearch_command_template must pin -Z")
    if "--domZ" in template:
        raise ScreenManifestError(
            "hmmsearch_command_template must not set --domZ "
            "(no downstream consumer rescales domain E-values)"
        )
    total = manifest["shard_sequence_total"]
    if isinstance(total, bool) or not isinstance(total, int) or total < 1:
        raise ScreenManifestError(
            f"shard_sequence_total must be a positive integer, got {total!r}"
        )
    observed_total = 0
    for shard in manifest["shards"]:
        count = shard.get("sequence_count")
        if isinstance(count, bool) or not isinstance(count, int) or count < 1:
            raise ScreenManifestError(
                f"shard sequence_count must be a positive integer: "
                f"{shard.get('name')!r} -> {count!r}"
            )
        observed_total += count
    if observed_total != total:
        raise ScreenManifestError(
            f"shard sequence_count total {observed_total} does not match "
            f"database_size_Z {total} / shard_sequence_total {total}"
        )
    for task in manifest["tasks"]:
        task_z = task.get("database_size_Z")
        if task_z != declared:
            raise ScreenManifestError(
                "shard runs used different database sizes: expected "
                f"database_size_Z={declared}, got {task_z!r} for "
                f"({task.get('family')}, {task.get('shard')})"
            )
    return manifest


def validate_manifest(manifest, *, require_database_size=True):
    """Validate a family x shard screening contract.

    ``require_database_size=False`` is the explicit legacy escape hatch for
    manifests written before the database-size scale was recorded; callers must
    opt in, so a missing scale can never be silently accepted by default.
    """
    if not isinstance(manifest, dict):
        raise ScreenManifestError("manifest must be a JSON object")
    required = ("families", "shards", "tasks")
    if require_database_size:
        required = required + DATABASE_SIZE_FIELDS
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
    if require_database_size:
        validate_database_size_scale(manifest)
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


def build_manifest(shard_dir, hmm_dir, hmmout, families, evalue, database_size_z=None, database_size_basis=None):
    """Build and self-validate the screening contract, scale included.

    ``database_size_z`` is required: a screening run that does not know the full
    library size cannot prove that its shard E-values share one scale.
    """
    if database_size_z is None:
        raise ScreenManifestError(
            "--database-size-z is required: the full-library target count Z must be "
            "recorded, never defaulted to a shard-local sequence count"
        )
    if isinstance(database_size_z, bool) or not isinstance(database_size_z, int) or database_size_z < 1:
        raise ScreenManifestError(
            f"database_size_Z must be a positive integer, got {database_size_z!r}"
        )
    if not isinstance(database_size_basis, str) or not database_size_basis.strip():
        raise ScreenManifestError(
            "database_size_basis is required and must be a non-empty string "
            "(where the full-library Z came from)"
        )
    shard_paths = sorted(
        os.path.join(shard_dir, name) for name in os.listdir(shard_dir)
        if name.startswith("shard_") and name.endswith(".faa")
    )
    if not shard_paths:
        raise ScreenManifestError(f"no shards found: {shard_dir}")
    shards = [
        {
            "name": os.path.basename(path),
            "sha256": sha256(path),
            "sequence_count": count_fasta_records(path),
        }
        for path in shard_paths
    ]
    shard_sequence_total = sum(shard["sequence_count"] for shard in shards)
    if shard_sequence_total != database_size_z:
        raise ScreenManifestError(
            f"shard sequence total {shard_sequence_total} does not match the declared "
            f"full-database database_size_Z {database_size_z}"
        )
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
                # Every task carries the same full-library scale; validate below.
                "database_size_Z": database_size_z,
            })
    manifest = {
        "schema_version": 1,
        "families": families,
        "shards": shards,
        "tasks": tasks,
        "database_size_Z": database_size_z,
        "database_size_basis": database_size_basis,
        "shard_sequence_total": shard_sequence_total,
        "hmmsearch_command_template": COMMAND_TEMPLATE,
    }
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
    parser.add_argument(
        "--database-size-z",
        default=None,
        help="full-library target sequence count Z passed to every shard's hmmsearch -Z",
    )
    parser.add_argument(
        "--database-size-basis",
        default=None,
        help="where the full-library Z came from (required with --database-size-z)",
    )
    parser.add_argument(
        "--allow-legacy-manifest",
        action="store_true",
        help="accept a manifest written before the database-size scale was recorded",
    )
    parser.add_argument("--out", default="data/screen/screen_manifest.json")
    args = parser.parse_args()
    if args.database_size_z is not None:
        if not args.database_size_z.isdigit() or int(args.database_size_z) < 1:
            parser.error("--database-size-z must be a positive integer")
        args.database_size_z = int(args.database_size_z)
    if args.validate:
        with open(args.validate, encoding="utf-8-sig") as handle:
            manifest = json.load(handle)
        try:
            validate_manifest(
                manifest, require_database_size=not args.allow_legacy_manifest
            )
        except ScreenManifestError as exc:
            if args.allow_legacy_manifest:
                raise
            raise ScreenManifestError(
                f"{exc}; pass --allow-legacy-manifest only to read a historical "
                "manifest that predates the database-size scale"
            ) from exc
        if args.database_size_z is not None and manifest.get("database_size_Z") != args.database_size_z:
            raise ScreenManifestError(
                f"manifest database_size_Z {manifest.get('database_size_Z')!r} contradicts "
                f"the supplied --database-size-z {args.database_size_z}"
            )
        print(f"screen manifest already verified: {len(manifest['tasks'])} family-shard tasks")
        return
    if not args.families or not args.eval:
        parser.error("--families and --eval are required when building a manifest")
    if args.database_size_z is None:
        parser.error(
            "--database-size-z is required when building a manifest: the full-library "
            "target count Z must be recorded, never defaulted to a shard-local sequence count"
        )
    if not args.database_size_basis:
        parser.error("--database-size-basis is required with --database-size-z")
    families = args.families.split()
    manifest = build_manifest(
        args.shard_dir,
        args.hmm_dir,
        args.hmmout,
        families,
        args.eval,
        database_size_z=args.database_size_z,
        database_size_basis=args.database_size_basis,
    )
    with open(args.out, "w", encoding="utf-8", newline="\n") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    print(
        f"screen manifest verified: {len(manifest['tasks'])} family-shard tasks on "
        f"one full-database Z={manifest['database_size_Z']} "
        f"({manifest['database_size_basis']})"
    )


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, ScreenManifestError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
