#!/usr/bin/env python3
"""Build the pool-external profile-score layer with honest labelling.

The pool-external candidates carry only profile scores (no architecture,
localization, or InterPro evidence).  Per §11.9.4.2/§11.9.4.3 of the project
handoff they must:

* never be called "high confidence" — this layer is ``score_only``;
* carry an explicit E-value tier so weak "kept-entirely" specific-family hits
  (E >= 1e-10) are visibly marked rather than sitting next to strong hits;
* reproduce the published 35,558 pool-external anchors: specific families
  (Cys / dPHAMCL / PhaZ7-like / periplasmic) kept whole, broad families
  (type 1 / type 2 / nPHAMCL) kept only when effective best E < 1e-30, and
  with-lipase excluded (it is already archived as the deferred
  structural-validation tier).

The per-protein superfamily assignment reuses the frozen, tested rule in
``build_phaded_with_lipase_deferred_tier.py`` (trained profile overrides the
discovery assignment iff its best E-value is STRICTLY smaller; otherwise the
discovery assignment wins).
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
from collections import defaultdict
from pathlib import Path

STRONG_E = 1e-30
MID_E = 1e-10

WITH_LIPASE = "intracellular nPHASCL with lipase box"

# Specific families: kept whole (no E-value floor), but their weak hits are
# tier-labelled so they cannot be mistaken for strong hits.
SPECIFIC_KEEP_ALL = {
    "intracellular nPHASCL without lipase box",      # Cys
    "extracellular dPHAMCL",
    "extracellular native-SCL/PhaZ7-like",
    "periplasmic PHA depolymerases",
}

# Broad families: kept only when the effective best E-value is < 1e-30.
BROAD_STRONG_ONLY = {
    "extracellular dPHASCL type 1",
    "extracellular dPHASCL type 2",
    "intracellular nPHAMCL",
}

EVIDENCE_LAYER = "score_only"
OUTPUT_NAME = "pool_external_profile_hits.tsv"


def _load_deferred():
    """Load the frozen, tested with-lipase deferred-tier module by path."""
    here = Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location(
        "build_phaded_with_lipase_deferred_tier",
        here / "build_phaded_with_lipase_deferred_tier.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def score_tier(evalue: float) -> str:
    if evalue < STRONG_E:
        return "strong"
    if evalue < MID_E:
        return "mid"
    return "weak"


def _parse_e(s: str) -> float:
    if not s:
        return 1.0
    try:
        return float(s)
    except ValueError:
        return 1.0


def effective_evalue(row: dict[str, str]) -> float:
    """Winning-assignment best E-value (trained when it overrode, else discovery)."""
    if row.get("override_by_trained") == "1":
        return _parse_e(row.get("trained_best_evalue", ""))
    return _parse_e(row.get("discovery_best_evalue", ""))


def best_evalue_str(row: dict[str, str]) -> str:
    if row.get("override_by_trained") == "1":
        return row.get("trained_best_evalue", "") or ""
    return row.get("discovery_best_evalue", "") or ""


def select_and_tier(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Keep/tier each pool-external protein; exclude the deferred with-lipase tier."""
    out: list[dict[str, str]] = []
    for row in rows:
        sf = row["superfamily"]
        if sf == WITH_LIPASE:
            continue
        eff = effective_evalue(row)
        if sf in BROAD_STRONG_ONLY and eff >= STRONG_E:
            continue
        out.append({
            "protein_id": row["protein_id"],
            "genome": row.get("genome", ""),
            "superfamily": sf,
            "best_evalue": best_evalue_str(row),
            "score_tier": score_tier(eff),
            "evidence_layer": EVIDENCE_LAYER,
            "architecture_status": "pending",
            "localization_status": "pending",
            "interpro_status": "pending",
            "discovery_best_family": row.get("discovery_best_family", ""),
            "discovery_best_evalue": row.get("discovery_best_evalue", ""),
            "trained_best_model": row.get("trained_best_model", ""),
            "trained_best_evalue": row.get("trained_best_evalue", ""),
            "override_by_trained": row.get("override_by_trained", ""),
            "model_layer": row.get("model_layer", ""),
        })
    return out


def write_layer(rows: list[dict[str, str]], out_dir: Path) -> dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / OUTPUT_NAME
    fieldnames = [
        "protein_id", "genome", "superfamily", "best_evalue", "score_tier",
        "evidence_layer", "architecture_status", "localization_status",
        "interpro_status", "discovery_best_family", "discovery_best_evalue",
        "trained_best_model", "trained_best_evalue", "override_by_trained",
        "model_layer",
    ]
    ordered = sorted(rows, key=lambda r: r["protein_id"])
    with out_path.open("w", encoding="utf-8", newline="\n") as handle:
        w = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(ordered)

    per_sf: dict[str, int] = defaultdict(int)
    per_tier: dict[str, int] = defaultdict(int)
    for r in rows:
        per_sf[r["superfamily"]] += 1
        per_tier[r["score_tier"]] += 1

    summary = {
        "layer": EVIDENCE_LAYER,
        "honesty_note": (
            "score-only profile layer: NOT high confidence; architecture / "
            "localization / InterPro evidence is pending for every row. "
            "Specific families (Cys / dPHAMCL / PhaZ7-like / periplasmic) are "
            "kept whole and their weak hits are tier-labelled; broad families "
            "(type 1 / type 2 / nPHAMCL) are kept only when best E-value < 1e-30."
        ),
        "rule": (
            "trained profile overrides discovery assignment iff trained best "
            "E-value < discovery best E-value (strict); otherwise discovery wins"
        ),
        "total": len(rows),
        "per_superfamily": dict(sorted(per_sf.items())),
        "per_score_tier": dict(sorted(per_tier.items())),
        "output": str(out_path),
    }
    (out_dir / "pool_external_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--discovery", type=Path, required=True,
                        help="shard_pool_external.tsv (discovery-layer shard merge)")
    parser.add_argument("--trained", type=Path, required=True,
                        help="trained_pool_external.tsv (trained-profile scan)")
    parser.add_argument("--profile-manifest", type=Path, required=True,
                        help="profile_manifest.tsv (frozen profile registry)")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)

    deferred = _load_deferred()
    fam2sf = deferred.read_family_superfamily(args.profile_manifest)
    model2sf = deferred.read_trained_model_superfamily(args.profile_manifest)
    trained = deferred.load_trained_best(args.trained)
    rows, _ = deferred.classify_pool_external(
        args.discovery, trained, fam2sf, model2sf
    )
    selected = select_and_tier(rows)
    summary = write_layer(selected, args.out_dir)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
