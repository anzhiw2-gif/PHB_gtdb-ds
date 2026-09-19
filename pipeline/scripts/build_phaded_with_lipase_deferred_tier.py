#!/usr/bin/env python3
"""Build the with-lipase deferred structural-validation tier.

The pool-external high-confidence filter (2026-09-19) excluded the
``intracellular nPHASCL with lipase box`` superfamily from the 35,558
high-confidence pool-external candidates: 1,206,655 pool-external hits were
dropped as noise.  Per project decision those hits are NOT discarded — they are
preserved as a *deferred tier* for later structural validation (Foldseek
against PDB 8YNV / predicted-structure screening).  This script rebuilds that
tier, bit-for-bit, from the two authoritative scoring tables and the frozen
profile manifest.

Assignment rule (reproduces the published per-superfamily anchor counts,
validated 2026-09-19):

* discovery layer: superfamily = phaded_superfamily of the best family hit in
  ``shard_pool_external.tsv`` (33 discovery HMMs, recall-only, uncalibrated);
* trained overlay: a trained profile (``model_status == trained`` in the
  profile manifest) overrides the discovery assignment ONLY when its best
  E-value is STRICTLY smaller than the discovery best E-value;
* every pool-external protein is assigned exactly once.

The deferred tier is recall-only evidence.  It never enters the
high-confidence count, never changes a family/subtype call, and is labelled
``discovery_hmm_uncalibrated`` / ``deferred_structural_validation``.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

WITH_LIPASE = "intracellular nPHASCL with lipase box"
MODEL_LAYER = "discovery_hmm_uncalibrated"
DEFERRED_REASON = "deferred_structural_validation"

# Mirror of SUPERFAMILY_CRITERIA keys in filter_phaded_high_confidence.py
# (frozen 2026-09-19).  Used only to reproduce that filter's assignment
# precedence for the pool-internal deferred extraction.
SUPERFAMILY_NAMES = {
    "extracellular dPHAMCL",
    "extracellular dPHASCL type 1",
    "extracellular dPHASCL type 2",
    "extracellular native-SCL/PhaZ7-like",
    "intracellular nPHAMCL",
    "intracellular nPHASCL with lipase box",
    "intracellular nPHASCL without lipase box",
    "periplasmic PHA depolymerases",
}
TRAINED_MATRIX_EVIDENCE = {"profile_trained_hit", "profile_ambiguous_family"}


def read_family_superfamily(manifest: Path) -> dict[str, str]:
    """phaded_family_id -> phaded_superfamily from the profile manifest.

    Superfamily-level profiles have an empty family id and are skipped.
    """
    fam2sf: dict[str, str] = {}
    with manifest.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            fam = (row.get("phaded_family_id") or "").strip()
            sf = (row.get("phaded_superfamily") or "").strip()
            if fam and sf:
                fam2sf[fam] = sf
    return fam2sf


def read_trained_model_superfamily(manifest: Path) -> dict[str, str]:
    """profile_id -> phaded_superfamily for trained profiles only."""
    model2sf: dict[str, str] = {}
    with manifest.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (row.get("model_status") or "").strip() != "trained":
                continue
            pid = (row.get("profile_id") or "").strip()
            sf = (row.get("phaded_superfamily") or "").strip()
            if pid and sf:
                model2sf[pid] = sf
    return model2sf


def load_trained_best(path: Path) -> dict[str, tuple[str, float]]:
    """protein_id -> (best_model, best_evalue) from the trained-profile scan."""
    trained: dict[str, tuple[str, float]] = {}
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pid = (row.get("protein_id") or "").strip()
            if not pid:
                continue
            model = (row.get("best_model") or "").strip()
            try:
                e = float(row.get("best_evalue") or "1")
            except ValueError:
                e = 1.0
            trained[pid] = (model, e)
    return trained


def genome_of(protein_id: str) -> str:
    return protein_id.split("|", 1)[0]


def classify_pool_external(
    discovery: Path,
    trained: dict[str, tuple[str, float]],
    fam2sf: dict[str, str],
    model2sf: dict[str, str],
) -> tuple[list[dict[str, str]], dict[str, int]]:
    """Assign each pool-external protein once; return rows + per-superfamily counts.

    Discovery rows stream from ``discovery`` (1.42M rows); trained-only proteins
    are appended afterwards.  Unknown trained model stems are a hard error
    (fail-closed provenance, never silently dropped).
    """
    rows: list[dict[str, str]] = []
    counts: dict[str, int] = defaultdict(int)
    seen_trained: set[str] = set()

    with discovery.open(encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pid = (row.get("protein_id") or "").strip()
            if not pid:
                continue
            best_family = (row.get("best_family") or "").strip()
            families_hit = (row.get("families_hit") or "").strip()
            try:
                disc_e = float(row.get("best_evalue") or "1")
            except ValueError:
                disc_e = 1.0
            sf = fam2sf.get(best_family, "")
            override = "0"
            trained_model = ""
            trained_e_str = ""
            if pid in trained:
                seen_trained.add(pid)
                model, t_e = trained[pid]
                if model not in model2sf:
                    raise ValueError(
                        f"trained model {model!r} for {pid!r} is not a trained "
                        f"profile in the manifest; refusing to guess its superfamily"
                    )
                if t_e < disc_e:  # strict: trained wins only when stronger
                    sf = model2sf[model]
                    override = "1"
                    trained_model = model
                    trained_e_str = "%.3g" % t_e
            if not sf:
                raise ValueError(
                    f"family {best_family!r} for {pid!r} has no superfamily mapping"
                )
            counts[sf] += 1
            rows.append({
                "protein_id": pid,
                "genome": genome_of(pid),
                "superfamily": sf,
                "discovery_best_family": best_family,
                "discovery_families_hit": families_hit,
                "discovery_best_evalue": "%.3g" % disc_e,
                "trained_best_model": trained_model,
                "trained_best_evalue": trained_e_str,
                "override_by_trained": override,
                "model_layer": MODEL_LAYER,
            })

    # Trained-only proteins (no discovery hit at all).
    for pid, (model, t_e) in sorted(trained.items()):
        if pid in seen_trained:
            continue
        if model not in model2sf:
            raise ValueError(
                f"trained model {model!r} for {pid!r} is not a trained profile "
                f"in the manifest; refusing to guess its superfamily"
            )
        sf = model2sf[model]
        counts[sf] += 1
        rows.append({
            "protein_id": pid,
            "genome": genome_of(pid),
            "superfamily": sf,
            "discovery_best_family": "",
            "discovery_families_hit": "",
            "discovery_best_evalue": "",
            "trained_best_model": model,
            "trained_best_evalue": "%.3g" % t_e,
            "override_by_trained": "1",
            "model_layer": "trained_profile",
        })

    return rows, dict(counts)


def write_deferred_tier(rows: list[dict[str, str]], out_dir: Path) -> tuple[Path, int]:
    """Write the pool-external with-lipase deferred TSV; return (path, n_rows)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "pool_external_with_lipase_deferred.tsv"
    fieldnames = [
        "protein_id", "genome", "superfamily",
        "discovery_best_family", "discovery_families_hit", "discovery_best_evalue",
        "trained_best_model", "trained_best_evalue", "override_by_trained",
        "model_layer",
    ]
    deferred = [r for r in rows if r["superfamily"] == WITH_LIPASE]
    with out_path.open("w", encoding="utf-8", newline="\n") as handle:
        w = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(deferred, key=lambda r: r["protein_id"]))
    return out_path, len(deferred)


def extract_pool_internal_with_lipase(
    classification: Path, candidate_faa: Path
) -> list[dict[str, str]]:
    """Pool-internal with-lipase rows: unique claim only, intersected with the
    frozen candidate universe (FASTA headers)."""
    universe: set[str] = set()
    with candidate_faa.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                universe.add(line[1:].split()[0])
    out: list[dict[str, str]] = []
    with classification.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pid = (row.get("protein_id") or "").strip()
            claim = (row.get("superfamily_claim") or "").strip()
            if pid in universe and claim == WITH_LIPASE:
                out.append({
                    "protein_id": pid,
                    "genome": genome_of(pid),
                    "superfamily": claim,
                    "superfamily_confidence": row.get("superfamily_confidence", ""),
                    "best_evalue": row.get("best_evalue", ""),
                    "n_families_hit": row.get("n_families_hit", ""),
                    "families_hit": row.get("families_hit", ""),
                    "model_layer": MODEL_LAYER,
                })
    return sorted(out, key=lambda r: r["protein_id"])


def extract_pool_internal_filter_consistent(
    classification: Path, matrix: Path, confounders: Path
) -> list[dict[str, str]]:
    """Pool-internal with-lipase rows using the SAME assignment precedence as
    filter_phaded_high_confidence.py: matrix branch (superfamily in the frozen
    criteria AND profile evidence in {trained_hit, ambiguous_family}) first,
    then the discovery-unique fallback.  Rows are NOT screened by the common /
    architecture criteria — the deferred tier keeps them all for later
    structural validation.  Reproduces the published ``total: 13``."""
    matrix_cols: dict[str, dict[str, str]] = {}
    with matrix.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            acc = (row.get("accession") or "").strip()
            if acc:
                matrix_cols[acc] = row
    cls_rows: dict[str, dict[str, str]] = {}
    with classification.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pid = (row.get("protein_id") or "").strip()
            if pid:
                cls_rows[pid] = row
    scope = set(matrix_cols) | set(cls_rows)

    out: list[dict[str, str]] = []
    for acc in sorted(scope):
        mc = matrix_cols.get(acc, {})
        cl = cls_rows.get(acc, {})
        sf: str | None = None
        evidence = ""
        matrix_sf = (mc.get("phaded_superfamily_best") or "").strip()
        matrix_pev = (mc.get("profile_evidence_status") or "").strip()
        if matrix_sf in SUPERFAMILY_NAMES and matrix_pev in TRAINED_MATRIX_EVIDENCE:
            sf = matrix_sf
            evidence = matrix_pev
        else:
            claim = (cl.get("superfamily_claim") or "").strip()
            confidence = (cl.get("superfamily_confidence") or "").strip()
            if confidence == "unique" and claim in SUPERFAMILY_NAMES:
                sf = claim
                evidence = "discovery_unique"
        if sf == WITH_LIPASE:
            out.append({
                "protein_id": acc,
                "genome": mc.get("genome", "") or genome_of(acc),
                "superfamily": sf,
                "assignment_evidence": evidence,
                "matrix_superfamily_best": matrix_sf,
                "classification_claim": (cl.get("superfamily_claim") or "").strip(),
                "best_evalue": cl.get("best_evalue", ""),
                "n_families_hit": cl.get("n_families_hit", ""),
                "families_hit": cl.get("families_hit", ""),
                "model_layer": MODEL_LAYER,
            })
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--discovery", type=Path, required=True,
                        help="shard_pool_external.tsv (discovery-layer shard merge)")
    parser.add_argument("--trained", type=Path, required=True,
                        help="trained_pool_external.tsv (trained-profile scan)")
    parser.add_argument("--profile-manifest", type=Path, required=True,
                        help="profile_manifest.tsv (frozen profile registry)")
    parser.add_argument("--pool-classification", type=Path, required=True,
                        help="superfamily_classification.tsv (pool-internal)")
    parser.add_argument("--candidate-faa", type=Path, required=True,
                        help="candidate_union.faa (frozen 109,087-protein universe)")
    parser.add_argument("--matrix", type=Path, default=None,
                        help="optional phaded_full_library_foldseek_subtype_matrix.tsv "
                             "(enables the filter-consistent 13-row extraction)")
    parser.add_argument("--confounders", type=Path, default=None,
                        help="optional confounder_candidates.tsv (with --matrix)")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)

    fam2sf = read_family_superfamily(args.profile_manifest)
    model2sf = read_trained_model_superfamily(args.profile_manifest)
    trained = load_trained_best(args.trained)
    rows, counts = classify_pool_external(args.discovery, trained, fam2sf, model2sf)

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    deferred_path, n_deferred = write_deferred_tier(rows, out_dir)

    internal = extract_pool_internal_with_lipase(
        args.pool_classification, args.candidate_faa
    )
    internal_path = out_dir / "pool_internal_with_lipase_discovery_unique.tsv"
    internal_fields = [
        "protein_id", "genome", "superfamily", "superfamily_confidence",
        "best_evalue", "n_families_hit", "families_hit", "model_layer",
    ]
    with internal_path.open("w", encoding="utf-8", newline="\n") as handle:
        w = csv.DictWriter(handle, fieldnames=internal_fields, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(internal)

    filter_consistent = []
    filter_consistent_path = ""
    if args.matrix is not None and args.confounders is not None:
        filter_consistent = extract_pool_internal_filter_consistent(
            args.pool_classification, args.matrix, args.confounders
        )
        filter_consistent_path = str(
            out_dir / "pool_internal_with_lipase_filter_consistent.tsv"
        )
        fc_fields = [
            "protein_id", "genome", "superfamily", "assignment_evidence",
            "matrix_superfamily_best", "classification_claim", "best_evalue",
            "n_families_hit", "families_hit", "model_layer",
        ]
        with open(filter_consistent_path, "w", encoding="utf-8", newline="\n") as handle:
            w = csv.DictWriter(handle, fieldnames=fc_fields, delimiter="\t",
                               lineterminator="\n", extrasaction="ignore")
            w.writeheader()
            w.writerows(filter_consistent)

    summary = {
        "rule": (
            "trained profile overrides discovery assignment iff trained best "
            "E-value < discovery best E-value (strict); otherwise discovery wins"
        ),
        "deferred_reason": DEFERRED_REASON,
        "with_lipase_pool_external_deferred": n_deferred,
        "with_lipase_pool_internal_discovery_unique": len(internal),
        "with_lipase_pool_internal_filter_consistent": len(filter_consistent),
        "per_superfamily_final_counts": dict(sorted(counts.items())),
        "outputs": {
            "pool_external_deferred": str(deferred_path),
            "pool_internal_discovery_unique": str(internal_path),
            "pool_internal_filter_consistent": filter_consistent_path,
        },
    }
    summary_path = out_dir / "deferred_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
