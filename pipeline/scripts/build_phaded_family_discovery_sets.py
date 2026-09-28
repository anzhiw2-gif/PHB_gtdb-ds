#!/usr/bin/env python3
"""Prepare per-family discovery-layer training sets with experimental-positive provenance.

Reads the frozen PhaDED reference ledger (723 rows) + reference FASTA, groups by
phaded_family_id, and emits, for each reference-only family profile, a training
FASTA plus a provenance manifest marking every training sequence as
experimental_positive vs annotation_only (with PMID/DOI for positives).

This builder is **recall-only and purely additive**: it constructs *training
sets*.  It never filters, drops, demotes or excludes a candidate, never assigns a
family call, and its output can never delete a row from any candidate table (the
discovery layer only recalls — AGENTS.md, 78.69% sensitivity on the 563 x1
confounder hits, no discriminating power).

Which families get a discovery set (F3)
---------------------------------------
The selector is keyed on ``model_layer``, never on the legacy ``model_status``
column (``model_status`` still reads ``trained`` for discovery-layer HMMs):

* ``reference_query_only`` — no HMM exists yet: build a discovery set
  (:data:`DISCOVERY_SET_MODEL_LAYERS`);
* ``discovery_hmm_uncalibrated`` — a recall-only HMM already exists: no new set;
* ``sequence_family_hmm_validated`` / ``calibrated_candidate_model`` — the family
  already has a discriminating model
  (:data:`EXCLUDED_DISCOVERY_SET_MODEL_LAYERS`): no discovery set, and a
  recall-only model must never be presented as that family's model.

Legacy fallback: a manifest that predates the ``model_layer`` column is read as
``reference_query_only`` (non-discriminating) so historical manifests keep
parsing and no row is silently treated as validated; the fallback is stated in
the summary output (``layer_source`` / ``fallback_note``).  ``model_layer`` and
``functional_calibration_status`` are independent and neither is inferred from
the other.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

REFERENCE_ONLY_FAMILIES_EXPECTED = 33
TRAINED_FAMILIES = {"DED_hfam_4", "DED_hfam_52", "DED_hfam_55", "DED_hfam_70", "DED_hfam_8"}

#: The four governed model layers (AGENTS.md, redesign Task 5).
REFERENCE_QUERY_ONLY = "reference_query_only"
DISCOVERY_LAYER = "discovery_hmm_uncalibrated"
SEQUENCE_FAMILY_LAYER = "sequence_family_hmm_validated"
CALIBRATED_LAYER = "calibrated_candidate_model"
MODEL_LAYERS = frozenset({
    REFERENCE_QUERY_ONLY, DISCOVERY_LAYER, SEQUENCE_FAMILY_LAYER, CALIBRATED_LAYER,
})

#: Layers for which a discovery set is built: no HMM exists at all.
DISCOVERY_SET_MODEL_LAYERS = frozenset({REFERENCE_QUERY_ONLY})
#: Layers that already carry a discriminating family model: no discovery set.
EXCLUDED_DISCOVERY_SET_MODEL_LAYERS = frozenset({SEQUENCE_FAMILY_LAYER, CALIBRATED_LAYER})
assert DISCOVERY_SET_MODEL_LAYERS | EXCLUDED_DISCOVERY_SET_MODEL_LAYERS | {DISCOVERY_LAYER} \
    == MODEL_LAYERS, "the discovery-set classification must cover the model-layer vocabulary"
assert not (DISCOVERY_SET_MODEL_LAYERS & EXCLUDED_DISCOVERY_SET_MODEL_LAYERS), (
    "a layer cannot both need and be excluded from a discovery set"
)
assert DISCOVERY_LAYER not in EXCLUDED_DISCOVERY_SET_MODEL_LAYERS, (
    "a recall-only discovery HMM is not a validated family model"
)

#: Layer assumed when a manifest has no ``model_layer`` column at all.
ABSENT_LAYER_FALLBACK = REFERENCE_QUERY_ONLY
ABSENT_LAYER_NOTE = (
    "profile_manifest.tsv has no model_layer column: every profile is treated as "
    f"'{ABSENT_LAYER_FALLBACK}' (non-discriminating). A legacy model_status value of "
    "'trained' is never assumed validated; re-emit the manifest with model_layer so "
    "already validated families are excluded from discovery-set construction."
)
LAYER_SOURCE_COLUMN = "model_layer_column"
LAYER_SOURCE_ABSENT = "absent_model_layer_column"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_model_layer(row: dict[str, str]) -> str:
    """Return the row's declared ``model_layer``, fail-closed on anything else.

    ``model_status`` is deliberately not a fallback: it reads ``trained`` for
    discovery-layer HMMs too.  A row without a declared layer resolves to
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


def read_profile_manifest(manifest: Path) -> "tuple[list[dict[str, str]], dict[str, object]]":
    """Read the profile manifest, resolving ``model_layer`` for every row.

    Rows gain a resolved ``model_layer`` plus ``layer_source`` and a verbatim
    ``functional_calibration_status`` (``unknown`` when absent — never inferred
    from the layer).  The report documents where the layer came from.
    """
    with manifest.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        raw_rows = list(reader)
    column_present = "model_layer" in fieldnames
    layer_source = LAYER_SOURCE_COLUMN if column_present else LAYER_SOURCE_ABSENT
    rows: list[dict[str, str]] = []
    counts: dict[str, int] = {layer: 0 for layer in sorted(MODEL_LAYERS)}
    for raw in raw_rows:
        layer = resolve_model_layer(raw)
        row = dict(raw)
        row["model_layer"] = layer
        row["layer_source"] = layer_source
        row["functional_calibration_status"] = (
            (raw.get("functional_calibration_status") or "").strip() or "unknown"
        )
        counts[layer] += 1
        rows.append(row)
    return rows, {
        "model_layer_column_present": column_present,
        "layer_source": layer_source,
        "fallback_layer": "" if column_present else ABSENT_LAYER_FALLBACK,
        "fallback_note": "" if column_present else ABSENT_LAYER_NOTE,
        "model_layer_counts": {layer: n for layer, n in counts.items() if n},
    }


def select_discovery_set_families(rows: "list[dict[str, str]]") -> "set[str]":
    """Family ids that need a discovery set: the no-HMM layer, family kind only.

    This selects *what to build*, never what to exclude from a candidate table.
    """
    families: set[str] = set()
    for row in rows:
        if row["model_layer"] not in DISCOVERY_SET_MODEL_LAYERS:
            continue
        if (row.get("profile_kind") or "").strip() != "family":
            continue
        family = (row.get("phaded_family_id") or "").strip()
        if family:
            families.add(family)
    return families


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    current: list[str] = []
    current_id = ""
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(current)
                header = line[1:].split()[0]
                # Reference FASTA headers are "reference_id|accession"; key by
                # the accession (last |-delimited token) to match the ledger.
                current_id = header.split("|")[-1]
                current = []
            else:
                current.append(line.strip())
        if current_id:
            seqs[current_id] = "".join(current)
    return seqs


def build(
    ledger: Path, fasta: Path, profile_manifest: Path, out_dir: Path
) -> dict[str, object]:
    # 1) Families that need a discovery set (authoritative list from the profile
    #    manifest, keyed on model_layer: reference_query_only == no HMM yet).
    profile_rows, layer_report = read_profile_manifest(profile_manifest)
    reference_only_families = select_discovery_set_families(profile_rows)
    if len(reference_only_families) != REFERENCE_ONLY_FAMILIES_EXPECTED:
        raise ValueError(
            f"expected {REFERENCE_ONLY_FAMILIES_EXPECTED} reference-only families, "
            f"got {len(reference_only_families)}"
        )

    # 2) Ledger rows grouped by family.
    by_family: dict[str, list[dict[str, str]]] = defaultdict(list)
    with ledger.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            fid = row.get("phaded_family_id", "").strip()
            if fid:
                by_family[fid].append(row)

    # 3) Reference FASTA sequences keyed by accession.
    seqs = read_fasta(fasta)

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "training").mkdir(exist_ok=True)

    manifest_rows: list[dict[str, str]] = []
    for fid in sorted(reference_only_families):
        rows = by_family[fid]
        if not rows:
            raise ValueError(f"no ledger rows for family {fid}")
        n_total = len(rows)
        positives = [r for r in rows if r.get("evidence_status", "") == "experimental_positive"]
        annotations = [r for r in rows if r.get("evidence_status", "") == "annotation_only"]
        accessions = [r["accession"] for r in rows]
        missing = [a for a in accessions if a not in seqs]
        if missing:
            raise ValueError(f"{fid}: {len(missing)} accessions missing from FASTA: {missing[:5]}")

        faa = out_dir / "training" / f"{fid}.faa"
        with faa.open("w", encoding="utf-8", newline="\n") as handle:
            for r in rows:
                handle.write(f">{r['accession']} {fid} {r['evidence_status']}\n")
                handle.write(seqs[r["accession"]] + "\n")

        manifest_rows.append(
            {
                "phaded_family_id": fid,
                "n_sequences": str(n_total),
                "n_experimental_positive": str(len(positives)),
                "n_annotation_only": str(len(annotations)),
                "experimental_positive_accessions": ";".join(
                    sorted(r["accession"] for r in positives)
                ),
                "experimental_positive_pmids": ";".join(
                    sorted({r.get("pmid", "") for r in positives if r.get("pmid", "")})
                ),
                "experimental_positive_dois": ";".join(
                    sorted({r.get("primary_doi", "") for r in positives if r.get("primary_doi", "")})
                ),
                "training_faa": str(faa.name),
                "training_faa_sha256": sha256(faa),
                # Layer bookkeeping: the training set exists to build a
                # recall-only discovery model, never a validated one.
                "model_layer": DISCOVERY_LAYER,
                "functional_calibration_status": "not_function_calibrated",
                "recall_only": "true",
            }
        )

    manifest_path = out_dir / "family_training_provenance.tsv"
    fieldnames = [
        "phaded_family_id", "n_sequences", "n_experimental_positive",
        "n_annotation_only", "experimental_positive_accessions",
        "experimental_positive_pmids", "experimental_positive_dois",
        "training_faa", "training_faa_sha256",
        "model_layer", "functional_calibration_status", "recall_only",
    ]
    with manifest_path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest_rows)

    return {
        "reference_only_families": len(reference_only_families),
        "training_fasta_files": len(manifest_rows),
        "provenance_manifest": str(manifest_path),
        "excluded_trained_families": sorted(TRAINED_FAMILIES),
        "decision_field": "model_layer",
        "legacy_model_status_is_a_decision_input": False,
        "discovery_set_model_layers": sorted(DISCOVERY_SET_MODEL_LAYERS),
        "excluded_discovery_set_model_layers": sorted(EXCLUDED_DISCOVERY_SET_MODEL_LAYERS),
        "output_model_layer": DISCOVERY_LAYER,
        "recall_only": True,
        "filters_candidates": False,
        "deletes_candidates": False,
        "model_layer_column_present": layer_report["model_layer_column_present"],
        "layer_source": layer_report["layer_source"],
        "fallback_layer": layer_report["fallback_layer"],
        "fallback_note": layer_report["fallback_note"],
        "model_layer_counts": layer_report["model_layer_counts"],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--profile-manifest", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    print(json.dumps(build(args.ledger, args.fasta, args.profile_manifest, args.out_dir), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
