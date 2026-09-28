#!/usr/bin/env python3
"""Build the **v2 calibration readiness** consumed by the promotion run (F15).

The frozen amendment readiness
(``runs/20260915_phaded_reference_panel_acquisition_02/results/amendment/profile_calibration_readiness.tsv``,
SHA-256 ``5a88a4c71e53f4b72e8473c4e5e79b6ff348f66f728095abf99301c58038c7d6``) predates the
2026-09-28 evidence-model redesign: its positive slot was the *historical row count*, so
five families that had already been curated away from the gate
(``DED_hfam_4/8/52/55/70``) are not in it at all, and the four trained superfamilies are
missing too.  This builder produces the input the redesigned family gate actually reads,
without touching any frozen artifact:

1. **append-only over the frozen readiness.**  Every frozen row is copied verbatim, in the
   frozen row order, and is a tab-field prefix of the emitted line.  No frozen column value
   is recomputed -- not even when the frozen value disagrees with today's manifest (the
   frozen value is authoritative for frozen rows).  Two qualification columns are appended:

   * ``qualified_e2_e3_positive_count`` -- unique ``independence_group`` count of curated
     ledger rows whose ``experimental_evidence_grade`` is E2/E3 **and** whose
     ``functional_calibration_eligible`` is true;
   * ``independent_genus_count`` -- distinct genera (first token of ``organism``) among
     those groups.

   A row count of ungraded historical positives can no longer carry the gate, and a family
   with no qualifying row fails closed to ``0``.

2. **a closed profile frame.**  The frame is the frozen profile manifest (38 families +
   8 superfamilies).  Reconciliation runs in both directions and fails closed: a family in
   the curated ledger that has no profile, a profile whose family is absent from the ledger,
   a profile that is neither frozen nor resolvable to cited trained evidence, a profile that
   is both, or a duplicate profile id all raise :class:`ReadinessReconciliationError`.
   Nothing is ever dropped silently.

3. **materialisation is cited, never invented.**  Rows the frozen readiness does not cover
   are materialised from the frozen ``profile_manifest.tsv`` (identity, kind, superfamily,
   ``training_count``) plus:

   * ``heldout_validation.tsv`` -> ``heldout_positive`` (``pass == "pass"``);
   * ``challenge_validation.tsv`` -> ``external_bound_challenge`` (rows bound **and**
     tested) and ``challenge_status`` (``sufficient`` only when every bound challenge is
     rejected; otherwise ``not_sufficient``, never laundered);
   * the cited trained-evidence table -> ``external_bound_negative`` and
     ``family_resolved_negative_status`` (the only values no frozen table resolves
     mechanically);
   * the frozen superfamily-row convention (all legs ``0``, both statuses ``missing``,
     decision ``planned_not_run_superfamily_requires_family_resolution``) for superfamily
     profiles, whose resolution gate needs a member-family rollup and an axis-anchor audit
     that this readiness deliberately does not carry.

   ``external_bound_positive`` is ``0`` for materialised rows, following the frozen
   amendment convention (``independent_positive_records = 0``); ``existing_experimental_positive``
   keeps the manifest's historical ``training_count`` for audit only -- it no longer feeds
   any gate slot.

4. **explicit promotion scope.**  ``promotion_scope`` is ``authorized_promotion_frame`` for
   every row except the profiles named by ``--defer-profile``, which are marked
   ``deferred_design_spec_special_family`` together with the mandatory ``--defer-reason``.
   The promotion runner enforces the same scope, so a design-spec deferral cannot be
   dropped silently between the input and the promotion.

Outputs (all deterministic; no timestamp is embedded, the run log carries the clock):
``calibration_readiness_v2.tsv``, ``family_qualification_mapping.tsv`` (per family: counts,
contributing ``independence_group`` values, genera, and the ``file::column`` source) and
``readiness_reconciliation.json``.

Boundary: this is a governance input.  A gate pass is a *candidate* statement about
reference sequence evidence; no column here is a PHB/PHA phenotype claim.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]

#: Frozen readiness schema (order and names are part of the frozen contract).
FROZEN_FIELDS = (
    "profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id",
    "existing_experimental_positive", "external_bound_positive",
    "external_bound_negative", "external_bound_challenge", "heldout_positive",
    "family_resolved_negative_status", "challenge_status", "calibration_decision",
    "new_family_call_created",
)
QUALIFICATION_FIELDS = ("qualified_e2_e3_positive_count", "independent_genus_count")
PROVENANCE_FIELDS = (
    "qualification_source", "readiness_row_origin", "readiness_notes", "promotion_scope",
)
OUTPUT_FIELDS = FROZEN_FIELDS + QUALIFICATION_FIELDS + PROVENANCE_FIELDS

MAPPING_FIELDS = (
    "phaded_family_id", "profile_id", "phaded_superfamily",
    "qualified_e2_e3_positive_count", "independent_genus_count",
    "required_qualified_e2_e3_positive_count", "required_independent_genus_count",
    "qualifying_ledger_row_count", "contributing_independence_groups", "contributing_genera",
    "source_file", "source_column", "source_file_sha256",
)

LEDGER_REQUIRED = (
    "phaded_family_id", "independence_group", "experimental_evidence_grade",
    "functional_calibration_eligible", "organism",
)
MANIFEST_REQUIRED = (
    "profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id",
    "training_count", "model_status",
)
TRAINED_EVIDENCE_REQUIRED = (
    "profile_id", "external_bound_negative", "external_bound_negative_accessions",
    "family_resolved_negative_status", "calibration_decision", "evidence_citation",
)
HELDOUT_REQUIRED = ("family", "heldout", "pass")
CHALLENGE_REQUIRED = ("family", "challenge", "rejected")

#: Gate thresholds this input must feed (governance constants, never relaxed here).
REQUIRED_QUALIFIED_POSITIVES = 3
REQUIRED_INDEPENDENT_GENERA = 3

QUALIFYING_GRADES = ("E2", "E3")
KNOWN_GRADES = ("E3", "E2", "E1", "A", "pending")
TRUE_TOKENS = ("true", "1", "yes")
FALSE_TOKENS = ("false", "0", "no")
FAMILY_ID_PATTERN = re.compile(r"^DED_hfam_[0-9]+$")

FROZEN_ROW_ORIGIN = "frozen_readiness_verbatim"
MATERIALISED_FAMILY_ORIGIN = "materialised_frozen_manifest_and_trained_evidence"
MATERIALISED_SUPERFAMILY_ORIGIN = "materialised_frozen_manifest_superfamily_resolution_not_run"
SCOPE_AUTHORIZED = "authorized_promotion_frame"
SCOPE_DEFERRED = "deferred_design_spec_special_family"
SCOPES = (SCOPE_AUTHORIZED, SCOPE_DEFERRED)
SUPERFAMILY_DECISION = "planned_not_run_superfamily_requires_family_resolution"
SUPPORTED_DECISIONS = {
    "calibrated_candidate_model", "candidate_gate_passed_not_promoted",
    "reference_only_insufficient_panel", "planned_not_run_superfamily_requires_family_resolution",
    "planned_not_run",
}
POSITIVE_LEG_NOT_ASSESSED = "missing"

QUALIFICATION_COLUMNS = (
    "experimental_evidence_grade", "functional_calibration_eligible",
    "independence_group", "organism",
)
EXTERNAL_POSITIVE_BASIS = (
    "frozen amendment convention: independent_positive_records = 0 "
    "(runs/20260915_phaded_reference_panel_acquisition_02/results/amendment/"
    "panel_acquisition_amendment_report.json)"
)
SUPERFAMILY_LEG_BASIS = (
    "frozen superfamily-row convention: superfamily-resolution legs are not assessed in "
    "this readiness (they need a member-family rollup plus an axis-anchor audit)"
)

OUTPUT_READINESS = "calibration_readiness_v2.tsv"
OUTPUT_MAPPING = "family_qualification_mapping.tsv"
OUTPUT_RECONCILIATION = "readiness_reconciliation.json"


class ReadinessReconciliationError(RuntimeError):
    """Raised when the readiness frame cannot be closed; always fails closed."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO)).replace("\\", "/")
    except ValueError:
        return str(path.resolve()).replace("\\", "/")


def _describe(path: Path) -> dict[str, object]:
    return {"path": _display_path(path), "size": path.stat().st_size, "sha256": _sha256(path)}


def _read_tsv(path: Path, required: tuple[str, ...], what: str) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.is_symlink():
        raise ReadinessReconciliationError(f"{what} is not a regular file: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        missing = [field for field in required if field not in fieldnames]
        if missing:
            raise ReadinessReconciliationError(
                f"{what} is missing columns: {','.join(missing)} ({path})")
        rows = [{key: (value or "") for key, value in row.items() if key is not None}
                for row in reader]
    if not rows:
        raise ReadinessReconciliationError(f"{what} has no data rows: {path}")
    return fieldnames, rows


def _require_int(value: str, what: str) -> int:
    text = (value or "").strip()
    try:
        number = int(text)
    except ValueError as error:
        raise ReadinessReconciliationError(f"{what} is not an integer: {text!r}") from error
    if number < 0:
        raise ReadinessReconciliationError(f"{what} is negative: {text!r}")
    return number


def _bool_token(value: str, what: str) -> bool:
    text = (value or "").strip().lower()
    if text in TRUE_TOKENS:
        return True
    if text in FALSE_TOKENS:
        return False
    raise ReadinessReconciliationError(f"{what} is not a boolean token: {value!r}")


def _family_short(family_id: str) -> str:
    return family_id[len("DED_"):] if family_id.startswith("DED_") else family_id


def read_frozen_readiness(path: Path) -> list[dict[str, str]]:
    """Read the frozen readiness, requiring its exact schema (append-only safety)."""
    if not path.is_file() or path.is_symlink():
        raise ReadinessReconciliationError(f"frozen readiness is not a regular file: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        rows = [{key: (value or "") for key, value in row.items() if key is not None}
                for row in reader]
    if tuple(fieldnames) != FROZEN_FIELDS:
        raise ReadinessReconciliationError(
            "frozen readiness schema changed; refusing to treat it as the frozen contract: "
            f"expected {list(FROZEN_FIELDS)}, found {fieldnames}")
    if not rows:
        raise ReadinessReconciliationError(f"frozen readiness has no data rows: {path}")
    seen: set[str] = set()
    for row in rows:
        profile_id = row["profile_id"].strip()
        if not profile_id:
            raise ReadinessReconciliationError("frozen readiness has an empty profile_id")
        if profile_id in seen:
            raise ReadinessReconciliationError(f"frozen readiness has a duplicate profile_id: {profile_id}")
        seen.add(profile_id)
        for field in ("existing_experimental_positive", "external_bound_positive",
                      "external_bound_negative", "external_bound_challenge", "heldout_positive"):
            _require_int(row[field], f"{profile_id}: {field}")
    return rows


def qualify_ledger(ledger_rows: list[dict[str, str]], ledger_path: Path,
                   what: str = "curated ledger") -> dict[str, dict[str, object]]:
    """Per-family qualified E2/E3 evidence, keyed by ``phaded_family_id``."""
    per_family: dict[str, dict[str, object]] = {}
    for index, row in enumerate(ledger_rows):
        context = f"{what} row {index + 2}"
        family = row["phaded_family_id"].strip()
        if not family:
            raise ReadinessReconciliationError(f"{context}: empty phaded_family_id")
        grade = row["experimental_evidence_grade"].strip()
        if grade not in KNOWN_GRADES:
            raise ReadinessReconciliationError(f"{context}: unknown experimental_evidence_grade {grade!r}")
        eligible = _bool_token(row["functional_calibration_eligible"],
                               f"{context}: functional_calibration_eligible")
        bucket = per_family.setdefault(family, {"groups": set(), "genera": set(), "rows": 0})
        if grade not in QUALIFYING_GRADES or not eligible:
            continue
        group = row["independence_group"].strip()
        if not group:
            raise ReadinessReconciliationError(f"{context}: qualifying row without an independence_group")
        organism = row["organism"].strip()
        genus = organism.split()[0] if organism else ""
        if not genus:
            raise ReadinessReconciliationError(f"{context}: qualifying row without an organism/genus")
        bucket["groups"].add(group)  # type: ignore[union-attr]
        bucket["genera"].add(genus)  # type: ignore[union-attr]
        bucket["rows"] = int(bucket["rows"]) + 1  # type: ignore[arg-type]
    return {family: {"groups": sorted(bucket["groups"]),  # type: ignore[arg-type]
                     "genera": sorted(bucket["genera"]),  # type: ignore[arg-type]
                     "rows": int(bucket["rows"])}  # type: ignore[arg-type]
            for family, bucket in per_family.items()}


def read_profile_manifest(path: Path) -> list[dict[str, object]]:
    _, rows = _read_tsv(path, MANIFEST_REQUIRED, "profile manifest")
    profiles: list[dict[str, object]] = []
    seen: set[str] = set()
    for index, row in enumerate(rows):
        context = f"profile manifest row {index + 2}"
        profile_id = row["profile_id"].strip()
        if not profile_id:
            raise ReadinessReconciliationError(f"{context}: empty profile_id")
        if profile_id in seen:
            raise ReadinessReconciliationError(f"profile manifest has a duplicate profile_id: {profile_id}")
        seen.add(profile_id)
        kind = row["profile_kind"].strip()
        if kind not in {"family", "superfamily"}:
            raise ReadinessReconciliationError(f"{context}: unsupported profile_kind {kind!r}")
        family_id = row["phaded_family_id"].strip()
        if kind == "family" and not FAMILY_ID_PATTERN.match(family_id):
            raise ReadinessReconciliationError(
                f"{context}: family profile without a DED_hfam_<n> phaded_family_id: {family_id!r}")
        if kind == "superfamily" and family_id:
            raise ReadinessReconciliationError(
                f"{context}: superfamily profile must not carry a phaded_family_id: {family_id!r}")
        profiles.append({
            "profile_id": profile_id,
            "profile_kind": kind,
            "phaded_superfamily": row["phaded_superfamily"].strip(),
            "phaded_family_id": family_id,
            "training_count": _require_int(row["training_count"], f"{context}: training_count"),
            "model_status": row["model_status"].strip(),
        })
    return profiles


def read_trained_evidence(path: Path) -> dict[str, dict[str, str]]:
    _, rows = _read_tsv(path, TRAINED_EVIDENCE_REQUIRED, "trained-profile evidence")
    evidence: dict[str, dict[str, str]] = {}
    for index, row in enumerate(rows):
        context = f"trained-profile evidence row {index + 2}"
        profile_id = row["profile_id"].strip()
        if not profile_id:
            raise ReadinessReconciliationError(f"{context}: empty profile_id")
        if profile_id in evidence:
            raise ReadinessReconciliationError(f"trained-profile evidence has a duplicate profile_id: {profile_id}")
        _require_int(row["external_bound_negative"], f"{context}: external_bound_negative")
        if not row["external_bound_negative_accessions"].strip():
            raise ReadinessReconciliationError(f"{context}: negative leg without bound accessions")
        if not row["evidence_citation"].strip():
            raise ReadinessReconciliationError(f"{context}: leg values without an evidence_citation")
        decision = row["calibration_decision"].strip()
        if decision not in SUPPORTED_DECISIONS:
            raise ReadinessReconciliationError(f"{context}: unsupported calibration_decision {decision!r}")
        evidence[profile_id] = {key: value.strip() for key, value in row.items() if key is not None}
    return evidence


def read_heldout_validation(path: Path) -> dict[str, dict[str, object]]:
    _, rows = _read_tsv(path, HELDOUT_REQUIRED, "held-out validation")
    heldout: dict[str, dict[str, object]] = {}
    for index, row in enumerate(rows):
        context = f"held-out validation row {index + 2}"
        family = row["family"].strip()
        if not family:
            raise ReadinessReconciliationError(f"{context}: empty family")
        if family in heldout:
            raise ReadinessReconciliationError(f"held-out validation has a duplicate family: {family}")
        verdict = row["pass"].strip().lower()
        heldout[family] = {
            "heldout_positive": 1 if verdict == "pass" else 0,
            "heldout_accession": row["heldout"].strip(),
            "verdict": verdict or "pending",
        }
    return heldout


def read_challenge_validation(path: Path) -> dict[str, dict[str, object]]:
    _, rows = _read_tsv(path, CHALLENGE_REQUIRED, "challenge validation")
    challenge: dict[str, dict[str, object]] = {}
    for index, row in enumerate(rows):
        context = f"challenge validation row {index + 2}"
        family = row["family"].strip()
        if not family:
            raise ReadinessReconciliationError(f"{context}: empty family")
        rejected = _bool_token(row["rejected"], f"{context}: rejected")
        bucket = challenge.setdefault(family, {"count": 0, "rejected": 0})
        bucket["count"] = int(bucket["count"]) + 1  # type: ignore[arg-type]
        if rejected:
            bucket["rejected"] = int(bucket["rejected"]) + 1  # type: ignore[arg-type]
    return {family: {
        "count": int(bucket["count"]),  # type: ignore[arg-type]
        "rejected": int(bucket["rejected"]),  # type: ignore[arg-type]
        "challenge_status": "sufficient" if bucket["count"] == bucket["rejected"] else "not_sufficient",
    } for family, bucket in challenge.items()}


def cross_check_summary(path: Path, qualification: dict[str, dict[str, object]]) -> dict[str, object]:
    """Require this builder's aggregation to reproduce F13's independence summary exactly."""
    if not path.is_file() or path.is_symlink():
        raise ReadinessReconciliationError(f"independence summary is not a regular file: {path}")
    summary = json.loads(path.read_text(encoding="utf-8"))
    checked: list[str] = []
    for key, field in (("family_functional_positive_counts", "groups"),
                       ("family_distinct_genus_counts", "genera")):
        declared = summary.get(key)
        if not isinstance(declared, dict):
            raise ReadinessReconciliationError(f"independence summary has no {key} mapping")
        for family, value in declared.items():
            computed = len(qualification.get(family, {}).get(field, []))  # type: ignore[arg-type]
            if int(value) != computed:
                raise ReadinessReconciliationError(
                    f"independence summary mismatch for {family}.{key}: "
                    f"declared {value}, recomputed {computed}")
            checked.append(family)
            if family not in qualification:
                raise ReadinessReconciliationError(
                    f"independence summary names a family the ledger does not carry: {family}")
        for family in qualification:
            if family not in declared:
                raise ReadinessReconciliationError(
                    f"ledger family missing from independence summary {key}: {family}")
    return {"path": _display_path(path), "status": "matched",
            "checked_family_entries": len(checked)}


def _leg_row(profile: dict[str, object], qualification: dict[str, dict[str, object]],
             *, origin: str, notes: str, qualified_source: str) -> dict[str, str]:
    family_id = str(profile["phaded_family_id"])
    bucket = qualification.get(family_id, {"groups": [], "genera": []}) if family_id else {"groups": [], "genera": []}
    return {
        "qualified_e2_e3_positive_count": str(len(bucket["groups"])),  # type: ignore[arg-type]
        "independent_genus_count": str(len(bucket["genera"])),  # type: ignore[arg-type]
        "qualification_source": qualified_source,
        "readiness_row_origin": origin,
        "readiness_notes": notes,
    }


def build(
    frozen_readiness: str | Path,
    curated_ledger: str | Path,
    profile_manifest: str | Path,
    trained_evidence: str | Path,
    output_dir: str | Path,
    *,
    run_id: str,
    defer_profiles: tuple[str, ...] | list[str] = (),
    defer_reason: str = "",
    heldout_validation: str | Path | None = None,
    challenge_validation: str | Path | None = None,
    independence_summary: str | Path | None = None,
    amendment_report: str | Path | None = None,
) -> dict[str, object]:
    """Build the v2 readiness; any unresolved frame member fails closed."""
    frozen_path = Path(frozen_readiness)
    ledger_path = Path(curated_ledger)
    manifest_path = Path(profile_manifest)
    trained_path = Path(trained_evidence)
    output_dir = Path(output_dir)

    frozen_rows = read_frozen_readiness(frozen_path)
    _, ledger_rows = _read_tsv(ledger_path, LEDGER_REQUIRED, "curated ledger")
    qualification = qualify_ledger(ledger_rows, ledger_path)
    profiles = read_profile_manifest(manifest_path)
    evidence = read_trained_evidence(trained_path)

    ledger_relpath = _display_path(ledger_path)
    qualified_source = "{}::{}".format(ledger_relpath, ",".join(QUALIFICATION_COLUMNS))

    frame_families = [str(profile["phaded_family_id"]) for profile in profiles
                      if profile["profile_kind"] == "family"]
    frame_by_id = {str(profile["profile_id"]): profile for profile in profiles}
    frozen_ids = [row["profile_id"].strip() for row in frozen_rows]

    # ---- fail-closed reconciliation, both directions -------------------------------
    ledger_families = set(qualification)
    missing_profiles = sorted(ledger_families - set(frame_families))
    if missing_profiles:
        raise ReadinessReconciliationError(
            "curated ledger families without a profile in the frame: " + ",".join(missing_profiles))
    without_ledger = sorted(set(frame_families) - ledger_families)
    if without_ledger:
        raise ReadinessReconciliationError(
            "profile families without curated ledger evidence: " + ",".join(without_ledger))
    unknown_frozen = [profile_id for profile_id in frozen_ids if profile_id not in frame_by_id]
    if unknown_frozen:
        raise ReadinessReconciliationError(
            "frozen readiness profiles absent from the frame: " + ",".join(unknown_frozen))
    materialised = [profile for profile in profiles
                    if str(profile["profile_id"]) not in set(frozen_ids)]
    duplicated = sorted(set(evidence) & set(frozen_ids))
    if duplicated:
        raise ReadinessReconciliationError(
            "profiles present both in the frozen readiness and in the trained evidence: "
            + ",".join(duplicated))

    heldout: dict[str, dict[str, object]] = {}
    challenge: dict[str, dict[str, object]] = {}
    unresolved = [str(profile["profile_id"]) for profile in materialised
                  if profile["profile_kind"] == "family" and str(profile["profile_id"]) not in evidence]
    if unresolved:
        raise ReadinessReconciliationError(
            "profiles the frozen readiness does not cover and no cited evidence resolves: "
            + ",".join(sorted(unresolved)))
    if any(profile["profile_kind"] == "family" for profile in materialised):
        if heldout_validation is None or challenge_validation is None:
            raise ReadinessReconciliationError(
                "held-out and challenge validation inputs are required to materialise family rows")
        heldout = read_heldout_validation(Path(heldout_validation))
        challenge = read_challenge_validation(Path(challenge_validation))
        for profile in materialised:
            if profile["profile_kind"] != "family":
                continue
            short = _family_short(str(profile["phaded_family_id"]))
            if short not in heldout:
                raise ReadinessReconciliationError(
                    f"{profile['profile_id']}: no held-out validation row for {short}")
            if short not in challenge:
                raise ReadinessReconciliationError(
                    f"{profile['profile_id']}: no challenge validation rows for {short}")

    # ---- emit rows: frozen verbatim first, then materialised in manifest order ------
    rows: list[dict[str, str]] = []
    for source in frozen_rows:
        row = {field: source[field] for field in FROZEN_FIELDS}
        row.update(_leg_row(frame_by_id[row["profile_id"]], qualification,
                            origin=FROZEN_ROW_ORIGIN, qualified_source=qualified_source,
                            notes="frozen_readiness_row_preserved_verbatim"))
        rows.append(row)

    materialised_notes: dict[str, str] = {}
    for profile in materialised:
        profile_id = str(profile["profile_id"])
        kind = str(profile["profile_kind"])
        family_id = str(profile["phaded_family_id"])
        row = {field: "" for field in FROZEN_FIELDS}
        row["profile_id"] = profile_id
        row["profile_kind"] = kind
        row["phaded_superfamily"] = str(profile["phaded_superfamily"])
        row["phaded_family_id"] = family_id
        row["external_bound_positive"] = "0"
        row["new_family_call_created"] = "false"
        if kind == "family":
            cited = evidence[profile_id]
            short = _family_short(family_id)
            held = heldout[short]
            chal = challenge[short]
            row["existing_experimental_positive"] = str(profile["training_count"])
            row["external_bound_negative"] = cited["external_bound_negative"]
            row["external_bound_challenge"] = str(chal["count"])
            row["heldout_positive"] = str(held["heldout_positive"])
            row["family_resolved_negative_status"] = cited["family_resolved_negative_status"]
            row["challenge_status"] = str(chal["challenge_status"])
            row["calibration_decision"] = cited["calibration_decision"]
            notes = (
                "materialised_from=profile_manifest.tsv(training_count={})"
                "; positive_leg={}"
                "; heldout=heldout_validation.tsv::{}({},verdict={},e_value={})"
                "; challenge=challenge_validation.tsv::{}({} tested,{} rejected,status={})"
                "; negative={}({})"
            ).format(
                profile["training_count"], EXTERNAL_POSITIVE_BASIS, short,
                held["heldout_accession"] or "none", held["verdict"], cited["evidence_citation"],
                short, chal["count"], chal["rejected"], chal["challenge_status"],
                cited["evidence_citation"], cited["external_bound_negative_accessions"],
            )
            origin = MATERIALISED_FAMILY_ORIGIN
        else:
            row["existing_experimental_positive"] = "0"
            row["external_bound_negative"] = "0"
            row["external_bound_challenge"] = "0"
            row["heldout_positive"] = "0"
            row["family_resolved_negative_status"] = POSITIVE_LEG_NOT_ASSESSED
            row["challenge_status"] = POSITIVE_LEG_NOT_ASSESSED
            row["calibration_decision"] = SUPERFAMILY_DECISION
            notes = ("materialised_from=profile_manifest.tsv; " + SUPERFAMILY_LEG_BASIS)
            origin = MATERIALISED_SUPERFAMILY_ORIGIN
        row.update(_leg_row(profile, qualification, origin=origin,
                            qualified_source=qualified_source, notes=notes))
        materialised_notes[profile_id] = notes
        rows.append(row)

    # ---- promotion scope ------------------------------------------------------------
    deferred = [str(profile_id) for profile_id in defer_profiles]
    if deferred and not defer_reason.strip():
        raise ReadinessReconciliationError(
            "a deferred profile requires an explicit --defer-reason; deferral is never implicit")
    emitted_ids = [row["profile_id"] for row in rows]
    unknown_deferrals = sorted(set(deferred) - set(emitted_ids))
    if unknown_deferrals:
        raise ReadinessReconciliationError(
            "deferred profile is not part of the readiness frame: " + ",".join(unknown_deferrals))
    if len(set(deferred)) != len(deferred):
        raise ReadinessReconciliationError("a deferred profile was named twice")
    for row in rows:
        row["promotion_scope"] = SCOPE_DEFERRED if row["profile_id"] in set(deferred) else SCOPE_AUTHORIZED

    # ---- outputs --------------------------------------------------------------------
    output_dir.mkdir(parents=True, exist_ok=True)
    readiness_path = output_dir / OUTPUT_READINESS
    mapping_path = output_dir / OUTPUT_MAPPING
    reconciliation_path = output_dir / OUTPUT_RECONCILIATION
    for path in (readiness_path, mapping_path, reconciliation_path):
        if path.exists():
            raise FileExistsError(path)

    with readiness_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(OUTPUT_FIELDS), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    mapping_source_column = ";".join(QUALIFICATION_COLUMNS)
    mapping_rows = []
    for family_id in sorted(qualification):
        bucket = qualification[family_id]
        profile_id = next((str(profile["profile_id"]) for profile in profiles
                           if str(profile["phaded_family_id"]) == family_id), "")
        superfamily = next((str(profile["phaded_superfamily"]) for profile in profiles
                            if str(profile["phaded_family_id"]) == family_id), "")
        mapping_rows.append({
            "phaded_family_id": family_id,
            "profile_id": profile_id,
            "phaded_superfamily": superfamily,
            "qualified_e2_e3_positive_count": str(len(bucket["groups"])),  # type: ignore[arg-type]
            "independent_genus_count": str(len(bucket["genera"])),  # type: ignore[arg-type]
            "required_qualified_e2_e3_positive_count": str(REQUIRED_QUALIFIED_POSITIVES),
            "required_independent_genus_count": str(REQUIRED_INDEPENDENT_GENERA),
            "qualifying_ledger_row_count": str(bucket["rows"]),
            "contributing_independence_groups": ";".join(bucket["groups"]),  # type: ignore[arg-type]
            "contributing_genera": ";".join(bucket["genera"]),  # type: ignore[arg-type]
            "source_file": ledger_relpath,
            "source_column": mapping_source_column,
            "source_file_sha256": _sha256(ledger_path),
        })
    with mapping_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(MAPPING_FIELDS), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(mapping_rows)

    inputs: dict[str, object] = {
        "frozen_readiness": _describe(frozen_path),
        "curated_ledger": _describe(ledger_path),
        "profile_manifest": _describe(manifest_path),
        "trained_evidence": _describe(trained_path),
    }
    cross_checks: dict[str, object] = {"independence_summary": {"status": "not_supplied"}}
    if heldout_validation is not None:
        inputs["heldout_validation"] = _describe(Path(heldout_validation))
        cross_checks["heldout_validation"] = {"status": "used_for_materialised_rows"}
    if challenge_validation is not None:
        inputs["challenge_validation"] = _describe(Path(challenge_validation))
        cross_checks["challenge_validation"] = {"status": "used_for_materialised_rows"}
    if independence_summary is not None:
        inputs["independence_summary"] = _describe(Path(independence_summary))
        cross_checks["independence_summary"] = cross_check_summary(
            Path(independence_summary), qualification)
    if amendment_report is not None:
        inputs["amendment_report"] = _describe(Path(amendment_report))

    meeting_bar = sorted(family_id for family_id, bucket in qualification.items()
                         if len(bucket["groups"]) >= REQUIRED_QUALIFIED_POSITIVES  # type: ignore[arg-type]
                         and len(bucket["genera"]) >= REQUIRED_INDEPENDENT_GENERA)  # type: ignore[arg-type]
    report = {
        "schema_version": 1,
        "run_id": run_id,
        "task": "build_phaded_calibration_readiness_v2",
        "deterministic": True,
        "timestamp_policy": "no clock value is embedded; the run log and input contract carry timestamps",
        "inputs": inputs,
        "counts": {
            "frame_profiles": len(profiles),
            "frame_families": len(frame_families),
            "frame_superfamilies": len(profiles) - len(frame_families),
            "frozen_readiness_rows": len(frozen_rows),
            "materialised_rows": len(materialised),
            "materialised_family_rows": sum(1 for p in materialised if p["profile_kind"] == "family"),
            "materialised_superfamily_rows": sum(1 for p in materialised if p["profile_kind"] == "superfamily"),
            "emitted_rows": len(rows),
            "ledger_rows": len(ledger_rows),
            "ledger_families": len(qualification),
            "qualifying_ledger_rows": sum(int(bucket["rows"]) for bucket in qualification.values()),  # type: ignore[arg-type]
            "qualifying_independence_groups": sum(len(bucket["groups"]) for bucket in qualification.values()),  # type: ignore[arg-type]
        },
        "coverage": {
            "ledger_families_without_profile": [],
            "profiles_without_ledger_family": [str(profile["profile_id"]) for profile in profiles
                                               if profile["profile_kind"] == "superfamily"],
            "frozen_readiness_profiles": frozen_ids,
            "materialised_profiles": [str(profile["profile_id"]) for profile in materialised],
            "families_meeting_governance_bar": meeting_bar,
        },
        "qualification": {
            "source_file": ledger_relpath,
            "source_column": mapping_source_column,
            "positive_slot_rule": "unique independence_group of E2/E3 rows with functional_calibration_eligible=true",
            "genus_slot_rule": "distinct first token of organism among those rows",
            "required_qualified_e2_e3_positive_count": REQUIRED_QUALIFIED_POSITIVES,
            "required_independent_genus_count": REQUIRED_INDEPENDENT_GENERA,
            "per_family": {family: {
                "qualified_e2_e3_positive_count": len(bucket["groups"]),  # type: ignore[arg-type]
                "independent_genus_count": len(bucket["genera"]),  # type: ignore[arg-type]
                "qualifying_ledger_row_count": int(bucket["rows"]),  # type: ignore[arg-type]
                "contributing_independence_groups": bucket["groups"],
                "contributing_genera": bucket["genera"],
            } for family, bucket in sorted(qualification.items())},
        },
        "cross_checks": cross_checks,
        "materialisation": {
            "family_rule": "frozen profile manifest identity + heldout/challenge validation + cited trained evidence",
            "superfamily_rule": SUPERFAMILY_LEG_BASIS,
            "external_bound_positive_basis": EXTERNAL_POSITIVE_BASIS,
            "notes": materialised_notes,
        },
        "promotion_scope": {
            "default": SCOPE_AUTHORIZED,
            "deferred_value": SCOPE_DEFERRED,
            "authorized_count": sum(1 for row in rows if row["promotion_scope"] == SCOPE_AUTHORIZED),
            "deferred_count": len(deferred),
            "deferred": [{"profile_id": profile_id, "reason": defer_reason.strip()}
                         for profile_id in deferred],
        },
        "append_only": {
            "frozen_fields": list(FROZEN_FIELDS),
            "appended_fields": list(QUALIFICATION_FIELDS + PROVENANCE_FIELDS),
            "frozen_rows_preserved_verbatim": True,
            "frozen_row_order_preserved": True,
            "frozen_values_recomputed": False,
        },
        "outputs": {
            "readiness": {"path": _display_path(readiness_path), "size": readiness_path.stat().st_size,
                          "sha256": _sha256(readiness_path)},
            "mapping": {"path": _display_path(mapping_path), "size": mapping_path.stat().st_size,
                        "sha256": _sha256(mapping_path)},
        },
        "boundaries": [
            "This file is a governance input: it moves gate slots, it does not create evidence.",
            "A gate pass is a candidate statement about reference sequence evidence and is not a PHB/PHA phenotype claim.",
            "No frozen artifact is modified; the frozen readiness keeps its own SHA-256.",
            "The superfamily-resolution legs are not assessed here; superfamily rows keep their frozen planned_not_run decision.",
        ],
    }
    reconciliation_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n",
                                  encoding="utf-8")
    return {"rows": rows, "qualification": qualification, "mapping": mapping_rows,
            "reconciliation": report, "outputs": report["outputs"],
            "promotion_scope": report["promotion_scope"]}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen-readiness", type=Path, required=True)
    parser.add_argument("--curated-ledger", type=Path, required=True)
    parser.add_argument("--profile-manifest", type=Path, required=True)
    parser.add_argument("--trained-evidence", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--heldout-validation", type=Path)
    parser.add_argument("--challenge-validation", type=Path)
    parser.add_argument("--independence-summary", type=Path)
    parser.add_argument("--amendment-report", type=Path)
    parser.add_argument("--defer-profile", action="append", default=[],
                        help="profile excluded from the authorized promotion frame by a design decision")
    parser.add_argument("--defer-reason", default="",
                        help="required when --defer-profile is used; recorded in the reconciliation")
    args = parser.parse_args(argv)
    result = build(
        args.frozen_readiness, args.curated_ledger, args.profile_manifest, args.trained_evidence,
        args.output_dir, run_id=args.run_id, defer_profiles=args.defer_profile,
        defer_reason=args.defer_reason, heldout_validation=args.heldout_validation,
        challenge_validation=args.challenge_validation,
        independence_summary=args.independence_summary,
        amendment_report=args.amendment_report,
    )
    print(json.dumps({
        "counts": result["reconciliation"]["counts"],  # type: ignore[index]
        "coverage": result["reconciliation"]["coverage"],  # type: ignore[index]
        "promotion_scope": result["reconciliation"]["promotion_scope"],  # type: ignore[index]
        "outputs": result["outputs"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
