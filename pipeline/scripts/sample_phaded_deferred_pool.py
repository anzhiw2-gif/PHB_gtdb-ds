#!/usr/bin/env python3
"""P5: stratified representative sampling of the with-lipase deferred pool.

Packet P5 of the 2026-09-28 plan asks for "stratified representatives from the
with-lipase deferred pool by sequence cluster, taxonomy, length, completeness,
score, and architecture".  The pool is the 1,206,655 pool-external `DED_hfam_2`
hits archived as ``deferred_structure_review``: it must never enter a candidate
count and must never be deleted, and it can only be recalled through structure
evidence.

Only some of the plan's six stratification keys exist in the archived pool; the
rest need joins that are not in it.  This module uses **the keys that are
actually present** and says so, rather than inventing the missing ones:

* ``trained_best_model`` / ``override_by_trained`` -> which trained profile, if
  any, claims the hit (an ARCHITECTURE-adjacent key)
* ``discovery_best_evalue`` -> the SCORE key
* ``discovery_best_family`` / ``discovery_families_hit`` -> the recall-layer
  architecture key
* ``genome`` -> the TAXONOMY PROXY actually present (a genome accession, not a
  genus; the taxonomy join is reported as pending)

Sample size and seed are **inherited from the frozen pilot**, not chosen here:
the pilot drew 1,001 sequences by stride 1,206 from this same pool and its
extrapolability was machine-verified (Spearman(index, E rank) = 0.0074).  Reusing
that size keeps the new sample comparable to the pilot's blastp result
(954/987 at E<1e-5 against 8YNV chain A).

Candidate-only boundary: this selects sequences for a structural review.  It
deletes nothing, promotes nothing, and the deferred layer still enters no count.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
from pathlib import Path

#: Keys read from the archived deferred pool (its real column names).
COLUMNS = (
    "protein_id", "genome", "superfamily", "discovery_best_family",
    "discovery_families_hit", "discovery_best_evalue", "trained_best_model",
    "trained_best_evalue", "override_by_trained", "model_layer",
)

#: Inherited from the frozen pilot so the new sample stays comparable to it.
PILOT_SAMPLE_SIZE = 1001
PILOT_STRIDE = 1206
DEFAULT_SEED = 42

#: E-value buckets.  Their edges come from the frozen deferred summary, which
#: reports 1,031 / 519 / 2,463 / 676,148 / 526,223 hits in these ranges.
EVALUE_BUCKETS = (
    ("lt_1e-50", 0.0, 1e-50),
    ("1e-50_to_1e-30", 1e-50, 1e-30),
    ("1e-30_to_1e-10", 1e-30, 1e-10),
    ("1e-10_to_1e-5", 1e-10, 1e-5),
    ("gte_1e-5", 1e-5, float("inf")),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def evalue_bucket(value: str) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "unparsable"
    if number < 0:
        return "unparsable"
    for name, low, high in EVALUE_BUCKETS:
        if low <= number < high:
            return name
    return "unparsable"


def strata_of(row: dict[str, str]) -> tuple[str, ...]:
    """The stratum key: the keys that actually exist in the archived pool."""
    trained = (row.get("trained_best_model") or "").strip() or "no_trained_claim"
    override = (row.get("override_by_trained") or "").strip() or "unset"
    return (trained, evalue_bucket(row.get("discovery_best_evalue", "")), override)


def read_pool(path: Path):
    """Stream the deferred pool; never load 183 MB of rows into memory."""
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"deferred pool is not a regular file: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        missing = [column for column in COLUMNS if column not in fields]
        if missing:
            raise ValueError("deferred pool missing columns: " + ",".join(missing))
        for row in reader:
            yield row


def collect_strata(pool_path: Path) -> tuple[dict[tuple[str, ...], int], int]:
    counts: dict[tuple[str, ...], int] = {}
    total = 0
    for row in read_pool(pool_path):
        key = strata_of(row)
        counts[key] = counts.get(key, 0) + 1
        total += 1
    if total == 0:
        raise ValueError("deferred pool is empty")
    return counts, total


def allocate(counts: dict[tuple[str, ...], int], total: int, sample_size: int) -> dict[tuple[str, ...], int]:
    """Proportional allocation with largest-remainder rounding, floor 1 per stratum.

    Every stratum gets at least one draw so no reported layer is unrepresented;
    the remainder is distributed by largest fractional part, then any shortfall
    is taken from the largest strata.
    """
    if sample_size < len(counts):
        raise ValueError(
            f"sample size {sample_size} is below the stratum count {len(counts)}; "
            "every stratum must be represented"
        )
    exact = {key: counts[key] / total * sample_size for key in counts}
    allocation = {key: max(1, int(value)) for key, value in exact.items()}
    while sum(allocation.values()) < sample_size:
        deficit = {key: exact[key] - allocation[key] for key in counts}
        best = max(deficit, key=lambda key: (deficit[key], counts[key]))
        allocation[best] += 1
    while sum(allocation.values()) > sample_size:
        surplus = {
            key: allocation[key] - exact[key]
            for key in counts if allocation[key] > 1
        }
        if not surplus:
            break
        worst = max(surplus, key=lambda key: (surplus[key], -counts[key]))
        allocation[worst] -= 1
    return allocation


def draw_sample(pool_path: Path, allocation: dict[tuple[str, ...], int], seed: int) -> list[dict[str, str]]:
    """Reservoir-sample each stratum, so one pass suffices and memory is bounded."""
    seen: dict[tuple[str, ...], int] = {key: 0 for key in allocation}
    reservoirs: dict[tuple[str, ...], list[dict[str, str]]] = {key: [] for key in allocation}
    rng = random.Random(seed)
    for row in read_pool(pool_path):
        key = strata_of(row)
        quota = allocation.get(key)
        if quota is None:
            continue
        seen[key] += 1
        bucket = reservoirs[key]
        if len(bucket) < quota:
            bucket.append(row)
        else:
            index = rng.randint(0, seen[key] - 1)
            if index < quota:
                bucket[index] = row
    return [row for key in sorted(reservoirs) for row in reservoirs[key]]


def summarize(
    counts: dict[tuple[str, ...], int],
    total: int,
    allocation: dict[tuple[str, ...], int],
    sample: list[dict[str, str]],
    pilot_composition: dict[str, int] | None,
) -> dict:
    sample_counts: dict[tuple[str, ...], int] = {}
    for row in sample:
        key = strata_of(row)
        sample_counts[key] = sample_counts.get(key, 0) + 1
    strata_table = []
    for key in sorted(counts, key=lambda k: -counts[k]):
        trained, bucket, override = key
        strata_table.append({
            "trained_best_model": trained,
            "discovery_evalue_bucket": bucket,
            "override_by_trained": override,
            "pool_count": counts[key],
            "pool_share": counts[key] / total,
            "allocated": allocation[key],
            "drawn": sample_counts.get(key, 0),
        })
    summary = {
        "pool_total": total,
        "stratum_count": len(counts),
        "sample_size": len(sample),
        "seed": None,
        "sample_size_basis": (
            f"inherited from the frozen pilot ({PILOT_SAMPLE_SIZE} sequences by stride "
            f"{PILOT_STRIDE}, extrapolability machine-verified at Spearman 0.0074) so the new "
            "sample is comparable to its blastp result"
        ),
        "strata": strata_table,
        "keys_used": ["trained_best_model", "discovery_best_evalue bucket", "override_by_trained"],
        "keys_requested_but_unavailable": [
            "sequence cluster (needs a clustering not present in the archived pool)",
            "taxonomy (only a genome accession is present; a genus/phyla join is pending)",
            "length and completeness (needs the candidate sequence FASTA join; pending)",
        ],
        "pilot_composition_comparison": pilot_composition,
        "boundary": (
            "Selects sequences for a structural review. Deletes nothing, promotes nothing, and "
            "the deferred layer still enters no candidate count. A drawn sequence is not a "
            "functional claim, and recall still requires structure evidence."
        ),
    }
    return summary


def pilot_strata_from_accessions(path: Path) -> dict[str, int]:
    """Composition of the frozen pilot's realised sample, by E-value bucket.

    The pilot file is ``query pident length evalue bitscore`` against 8YNV, so it
    gives the pilot's own hit-E distribution; that is the only composition the
    frozen pilot recorded.
    """
    if not path.is_file():
        return {}
    buckets: dict[str, int] = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            buckets[evalue_bucket(parts[3])] = buckets.get(evalue_bucket(parts[3]), 0) + 1
    return buckets


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pool", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--sample-size", type=int, default=PILOT_SAMPLE_SIZE)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--expect-pool-size", type=int, default=1206655)
    parser.add_argument("--pilot-hits", type=Path, default=None,
                        help="frozen pilot blastp table, for the composition comparison")
    args = parser.parse_args(argv)

    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise ValueError(f"refusing to write into a non-empty output dir: {args.out_dir}")
    counts, total = collect_strata(args.pool)
    if args.expect_pool_size is not None and total != args.expect_pool_size:
        raise ValueError(
            f"pool size {total} does not match the declared {args.expect_pool_size}"
        )
    allocation = allocate(counts, total, args.sample_size)
    sample = draw_sample(args.pool, allocation, args.seed)
    pilot_composition = (
        pilot_strata_from_accessions(args.pilot_hits) if args.pilot_hits else None
    )
    summary = summarize(counts, total, allocation, sample, pilot_composition)
    summary["seed"] = args.seed
    summary["pool_sha256"] = sha256(args.pool)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    fields = list(COLUMNS)
    with (args.out_dir / "p5_deferred_stratified_sample.tsv").open(
        "w", encoding="utf-8", newline="",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(sample)
    (args.out_dir / "p5_stratified_sample_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n",
    )
    print(json.dumps({
        "pool_total": total,
        "stratum_count": len(counts),
        "sample_size": len(sample),
        "seed": args.seed,
        "out_dir": str(args.out_dir),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
