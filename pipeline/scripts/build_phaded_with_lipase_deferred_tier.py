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
* validated overlay: a profile whose ``model_layer`` is in
  :data:`TRAINED_OVERLAY_MODEL_LAYERS` (``sequence_family_hmm_validated`` or
  ``calibrated_candidate_model``) overrides the discovery assignment ONLY when
  its best E-value is STRICTLY smaller than the discovery best E-value;
* every pool-external protein is assigned exactly once.

Model-layer semantics (F3)
--------------------------
The overlay is keyed on ``model_layer``, **never** on the legacy ``model_status``
column: ``model_status`` still reads ``trained`` for discovery-layer HMMs, so
keying on it would let a recall-only model reassign a pool-external protein out
of the deferred tier — i.e. delete candidates from the layer that must never be
deleted.  A ``discovery_hmm_uncalibrated`` profile is therefore never an overlay
model; when it would have been the strongest "trained" hit, the fact is recorded
in ``discovery_layer_overlay_suppressed`` for audit instead of being applied.

Legacy fallback: a manifest that predates the ``model_layer`` column cannot
express a validated sequence-family layer.  The absent layer is read as
``reference_query_only`` (non-discriminating), so the overlay is empty, a warning
is printed and the fallback is written into ``deferred_summary.json`` — a
``model_status == "trained"`` row is never assumed validated.

The deferred tier is recall-only evidence.  It never enters the
high-confidence count, never changes a family/subtype call, and is labelled
``discovery_hmm_uncalibrated`` / ``deferred_structural_validation``.  It is never
deleted: no layer may upgrade it out of the tier or downgrade it away.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

WITH_LIPASE = "intracellular nPHASCL with lipase box"
MODEL_LAYER = "discovery_hmm_uncalibrated"
DEFERRED_REASON = "deferred_structural_validation"

#: The four governed model layers (AGENTS.md / evidence-model redesign Task 5).
REFERENCE_QUERY_ONLY = "reference_query_only"
DISCOVERY_LAYER = "discovery_hmm_uncalibrated"
SEQUENCE_FAMILY_LAYER = "sequence_family_hmm_validated"
CALIBRATED_LAYER = "calibrated_candidate_model"
MODEL_LAYERS = frozenset({
    REFERENCE_QUERY_ONLY, DISCOVERY_LAYER, SEQUENCE_FAMILY_LAYER, CALIBRATED_LAYER,
})
#: The only layers that may override the discovery-layer assignment.  A
#: discovery-layer or HMM-less reference profile can never be an overlay model.
TRAINED_OVERLAY_MODEL_LAYERS = frozenset({SEQUENCE_FAMILY_LAYER, CALIBRATED_LAYER})
assert TRAINED_OVERLAY_MODEL_LAYERS <= MODEL_LAYERS, "overlay layers must be model layers"
assert DISCOVERY_LAYER not in TRAINED_OVERLAY_MODEL_LAYERS, (
    "the discovery layer is recall-only and must never override an assignment"
)
assert REFERENCE_QUERY_ONLY not in TRAINED_OVERLAY_MODEL_LAYERS, (
    "a reference_query_only profile has no HMM and must never override an assignment"
)

#: Layer assumed when a manifest has no ``model_layer`` column at all.
ABSENT_LAYER_FALLBACK = REFERENCE_QUERY_ONLY
ABSENT_LAYER_NOTE = (
    "profile_manifest.tsv has no model_layer column: every profile is treated as "
    f"'{ABSENT_LAYER_FALLBACK}' (non-discriminating), so the validated overlay is empty "
    "and the deferred tier keeps the discovery-layer assignment. A legacy model_status "
    "value of 'trained' is never assumed validated; re-emit the manifest with model_layer."
)
LAYER_SOURCE_COLUMN = "model_layer_column"
LAYER_SOURCE_ABSENT = "absent_model_layer_column"

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


def resolve_model_layer(row: dict[str, str]) -> str:
    """Return the row's declared ``model_layer``, fail-closed on anything else.

    ``model_status`` is deliberately not a fallback (it reads ``trained`` for
    discovery-layer HMMs too).  A row without a declared layer resolves to
    :data:`ABSENT_LAYER_FALLBACK`, which is non-discriminating.
    """
    declared = (row.get("model_layer") or "").strip()
    if not declared:
        return ABSENT_LAYER_FALLBACK
    if declared not in MODEL_LAYERS:
        profile_id = (row.get("profile_id") or "").strip() or "<unknown profile>"
        raise ValueError(
            f"profile {profile_id!r} declares model_layer={declared!r}, which is not a "
            f"recognized model layer; expected one of {sorted(MODEL_LAYERS)}"
        )
    return declared


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


def read_manifest_summary(manifest: Path) -> dict[str, object]:
    """Where the model layer came from, the per-layer counts and the profile ids.

    ``profile_ids`` makes the unknown-model provenance check below exact: a scan
    that names a model the manifest does not contain is fail-closed, while a
    scan that names a *recall-only* profile present in the manifest is merely
    suppressed (it cannot move a candidate).
    """
    with manifest.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        raw_rows = list(reader)
    column_present = "model_layer" in fieldnames
    counts: dict[str, int] = defaultdict(int)
    profile_ids: set[str] = set()
    for raw in raw_rows:
        counts[resolve_model_layer(raw)] += 1
        pid = (raw.get("profile_id") or "").strip()
        if pid:
            profile_ids.add(pid)
    return {
        "model_layer_column_present": column_present,
        "layer_source": LAYER_SOURCE_COLUMN if column_present else LAYER_SOURCE_ABSENT,
        "fallback_layer": "" if column_present else ABSENT_LAYER_FALLBACK,
        "fallback_note": "" if column_present else ABSENT_LAYER_NOTE,
        "model_layer_counts": dict(sorted(counts.items())),
        "profile_ids": profile_ids,
    }


def read_model_layer_report(manifest: Path) -> dict[str, object]:
    """Alias kept for callers that only want the report (no profile id set)."""
    report = read_manifest_summary(manifest)
    report.pop("profile_ids", None)
    return report


def read_trained_model_superfamily(manifest: Path) -> dict[str, str]:
    """profile_id -> phaded_superfamily for the validated overlay layers only.

    ``model_layer`` decides; the legacy ``model_status`` column is ignored, so a
    discovery-layer HMM whose ``model_status`` still says ``trained`` is NOT an
    overlay model.  A manifest without the column falls back to
    :data:`ABSENT_LAYER_FALLBACK` (non-discriminating) and yields an empty
    overlay — loudly, never silently assuming "trained => validated".
    """
    model2sf: dict[str, str] = {}
    with manifest.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if resolve_model_layer(row) not in TRAINED_OVERLAY_MODEL_LAYERS:
                continue
            pid = (row.get("profile_id") or "").strip()
            sf = (row.get("phaded_superfamily") or "").strip()
            if pid and sf:
                model2sf[pid] = sf
    return model2sf


def read_profile_id_superfamily(manifest: Path) -> dict[str, str]:
    """profile_id -> phaded_superfamily for EVERY declared profile (audit only).

    Used solely so a protein that is seen only by a recall-only profile keeps that
    profile's superfamily label instead of being dropped: recalled hits may add a
    label, they may never delete a candidate.  This mapping is never used to
    decide whether a candidate is kept, filtered or demoted.
    """
    profile2sf: dict[str, str] = {}
    with manifest.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            resolve_model_layer(row)  # fail-closed on unrecognized layer values
            pid = (row.get("profile_id") or "").strip()
            sf = (row.get("phaded_superfamily") or "").strip()
            if pid and sf:
                profile2sf[pid] = sf
    return profile2sf


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
    known_profile_ids: "set[str]" = frozenset(),
    profile2sf: "dict[str, str] | None" = None,
) -> tuple[list[dict[str, str]], dict[str, int]]:
    """Assign each pool-external protein once; return rows + per-superfamily counts.

    Discovery rows stream from ``discovery`` (1.42M rows); trained-only proteins
    are appended afterwards.  Only validated overlay layers may override the
    discovery assignment; a discovery-layer model that would have been stronger
    is recorded as suppressed, not applied.  A model stem absent from the
    manifest (``known_profile_ids``) is a hard error: fail-closed provenance,
    never silently dropped.
    """
    rows: list[dict[str, str]] = []
    counts: dict[str, int] = defaultdict(int)
    seen_trained: set[str] = set()
    if known_profile_ids:
        known = set(known_profile_ids)
    else:
        # Backwards-compatible inference for callers that only read the manifest
        # in fragments; the CLI always passes the full profile-id set.
        known = set(model2sf) | set(fam2sf) | {model for model, _ in trained.values()}

    def _check_known_stem(model: str, pid: str) -> None:
        if model not in known:
            raise ValueError(
                f"trained model {model!r} for {pid!r} is not a profile in the "
                f"manifest; refusing to guess its superfamily"
            )

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
            suppressed = ""
            if pid in trained:
                seen_trained.add(pid)
                model, t_e = trained[pid]
                _check_known_stem(model, pid)
                if model in model2sf:
                    if t_e < disc_e:  # strict: validated overlay wins only when stronger
                        sf = model2sf[model]
                        override = "1"
                        trained_model = model
                        trained_e_str = "%.3g" % t_e
                else:
                    # Present in the manifest but not a validated overlay layer
                    # (recall-only): it cannot move the candidate.  Recorded for
                    # audit instead of being applied.
                    suppressed = model
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
                # Audit only: a non-validated (recall-only) model was the best hit
                # but could not override, so the candidate stays where it was.
                "discovery_layer_overlay_suppressed": suppressed,
                "model_layer": MODEL_LAYER,
            })

    # Trained-only proteins (no discovery hit at all).
    for pid, (model, t_e) in sorted(trained.items()):
        if pid in seen_trained:
            continue
        _check_known_stem(model, pid)
        if model in model2sf:
            sf = model2sf[model]
            layer = SEQUENCE_FAMILY_LAYER
            override = "1"
            suppressed = ""
        else:
            # Best hit is a non-validated (recall-only) profile: it cannot move a
            # protein anywhere.  There is no discovery row to fall back to, so the
            # protein keeps the superfamily declared by that profile — recalled
            # hits may add a label, never delete a candidate.
            sf = (profile2sf or {}).get(model, "") or fam2sf.get(model, "")
            layer = MODEL_LAYER
            override = "0"
            suppressed = model
        if not sf:
            raise ValueError(
                f"model {model!r} for {pid!r} has no superfamily mapping"
            )
        counts[sf] += 1
        rows.append({
            "protein_id": pid,
            "genome": genome_of(pid),
            "superfamily": sf,
            "discovery_best_family": "",
            "discovery_families_hit": "",
            "discovery_best_evalue": "",
            "trained_best_model": model if override == "1" else "",
            "trained_best_evalue": "%.3g" % t_e,
            "override_by_trained": override,
            "discovery_layer_overlay_suppressed": suppressed,
            "model_layer": layer,
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
        "discovery_layer_overlay_suppressed", "model_layer",
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

    layer_report = read_manifest_summary(args.profile_manifest)
    if not layer_report["model_layer_column_present"]:
        print(f"[warn] {ABSENT_LAYER_NOTE}", file=sys.stderr)
    fam2sf = read_family_superfamily(args.profile_manifest)
    profile2sf = read_profile_id_superfamily(args.profile_manifest)
    model2sf = read_trained_model_superfamily(args.profile_manifest)
    trained = load_trained_best(args.trained)
    rows, counts = classify_pool_external(
        args.discovery, trained, fam2sf, model2sf,
        known_profile_ids=layer_report["profile_ids"],
        profile2sf=profile2sf,
    )
    suppressed = sum(
        1 for row in rows if row.get("discovery_layer_overlay_suppressed")
    )

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
            "a validated overlay profile (model_layer in trained_overlay_model_layers) "
            "overrides the discovery assignment iff its best E-value < discovery best "
            "E-value (strict); otherwise discovery wins. model_status is never a decision "
            "input."
        ),
        "decision_field": "model_layer",
        "legacy_model_status_is_a_decision_input": False,
        "model_layer": MODEL_LAYER,
        "trained_overlay_model_layers": sorted(TRAINED_OVERLAY_MODEL_LAYERS),
        "model_layer_column_present": layer_report["model_layer_column_present"],
        "layer_source": layer_report["layer_source"],
        "fallback_layer": layer_report["fallback_layer"],
        "fallback_note": layer_report["fallback_note"],
        "model_layer_counts": layer_report["model_layer_counts"],
        "trained_overlay_profiles": sorted(model2sf),
        "discovery_layer_overlay_suppressed": suppressed,
        "deferred_reason": DEFERRED_REASON,
        "deferred_tier_enters_confidence_counts": False,
        "deferred_tier_is_never_deleted": True,
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
