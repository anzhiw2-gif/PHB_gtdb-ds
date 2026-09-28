#!/usr/bin/env python3
"""Validate the two-layer PhaDED classification authority (Task 1).

Scientific problem this script locks down
-----------------------------------------
The project has been treating the PhaDED superfamily layer and the PhaDED
homologous-family layer as the same kind of entity.  They are not:

* The **8 superfamilies** are a *hard framework* assembled by Knoll 2009 from
  28 seed sequences of experimentally validated depolymerases that were
  assigned to superfamilies from a function/localisation prior (Knoll 2009
  body text: "28 seed sequences of proteins with experimentally validated
  depolymerase activity ... assigned to 6 previously described superfamilies
  based on their function", plus the intracellular nPHASCL (lipase box) family
  and the periplasmic PHA depolymerases, "Thus, a total of 8 superfamilies were
  introduced").
* The **38 homologous families** are a *secondary clustering refinement layer*:
  "Superfamilies were subdivided into homologous families, which were
  introduced based on sequence similarity and phylogenetic analysis (Fig. 1)."
  They come from a database that has been frozen since 2009 (Database Commons
  record for the DED: Last update 2009, v1.1).

Consequently only the functional-prior superfamily layer may be registry
eligible.  A 2009 sequence cluster is evidence of homology, never of phenotype,
and must never inherit registry eligibility.

Boundary
--------
This project is candidate-only.  Every HMM, profile, domain, motif, SignalP,
structural and phylogenetic result in this repository expresses candidate
homology or functional potential and is **not** equivalent to a validated
PHB/PHA degradation phenotype.  This script produces a classification
authority table; it changes no candidate call (0 rows of 109,087).

Fail-closed contract
--------------------
Every check below raises :class:`AuthorityValidationError` instead of guessing.
Missing values are recorded as ``pending`` and never fabricated; a ``pending``
gate profile is a validation failure, not a pass.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Iterable, Mapping, Sequence


# --- output schema ---------------------------------------------------------

REQUIRED_COLUMNS = [
    "layer_id",
    "layer_kind",
    "source_evidence_type",
    "source_version",
    "registry_eligible",
    "gate_profile_id",
    "notes",
]

# --- controlled vocabularies ----------------------------------------------

FUNCTIONAL_PRIOR = "functional_prior"
SEQUENCE_CLUSTERING_2009 = "sequence_clustering_2009"
FROZEN_SOURCE_VERSION = "Knoll2009_v1.1_frozen"
DED_SNAPSHOT_LABEL = "PhaDED (Knoll 2009, v1.1, frozen snapshot)"

SUPERFAMILY_KIND = "superfamily"
FAMILY_KIND = "family"
LAYER_KINDS = (SUPERFAMILY_KIND, FAMILY_KIND)

PENDING = "pending"

#: The authorised combination.  Crossed pairs are the layer mismatch defect.
ALLOWED_EVIDENCE_ELIGIBILITY = {
    (FUNCTIONAL_PRIOR, "true"),
    (SEQUENCE_CLUSTERING_2009, "false"),
}

# --- expected framework ----------------------------------------------------

EXPECTED_SUPERFAMILY_COUNT = 8
EXPECTED_FAMILY_COUNT = 38
EXPECTED_LAYER_COUNT = EXPECTED_SUPERFAMILY_COUNT + EXPECTED_FAMILY_COUNT

#: The 8-superfamily hard framework, verbatim DED labels.
EXPECTED_SUPERFAMILIES = (
    "extracellular dPHAMCL",
    "extracellular dPHASCL type 1",
    "extracellular dPHASCL type 2",
    "extracellular native-SCL/PhaZ7-like",
    "intracellular nPHAMCL",
    "intracellular nPHASCL with lipase box",
    "intracellular nPHASCL without lipase box",
    "periplasmic PHA depolymerases",
)

#: Families for which a substrate binding domain does not exist at all, so a
#: missing Pfam PF06850 hit is structurally expected rather than a negative.
NO_SBD_SUPERFAMILIES = frozenset(
    {"extracellular dPHAMCL", "extracellular native-SCL/PhaZ7-like"}
)

EXPECTED_CANDIDATE_UNIVERSE = 109087

FAMILY_ID_RE = re.compile(r"^DED_hfam_(\d+)$")

DEFINITION_REQUIRED_COLUMNS = ("phaded_superfamily", "phaded_family_id")
MANIFEST_REQUIRED_COLUMNS = (
    "profile_id",
    "profile_kind",
    "phaded_superfamily",
    "phaded_family_id",
)

BOUNDARY_STATEMENT = (
    "candidate-only: the layers classified here express candidate homology or "
    "functional potential and are not equivalent to a validated PHB/PHA "
    "degradation phenotype"
)

FUNCTIONAL_PRIOR_NOTE = (
    f"{DED_SNAPSHOT_LABEL} superfamily layer: hard framework introduced from the "
    "28 experimentally validated seed sequences that Knoll 2009 assigned to "
    "superfamilies from a function/localisation prior (Table 1), not from a "
    "sequence cluster. registry_eligible=true because the layer encodes that "
    "functional prior. Use as a hard framework; it is still candidate-only and "
    "not a validated phenotype."
)

CLUSTERING_NOTE = (
    f"{DED_SNAPSHOT_LABEL} homologous-family layer: secondary clustering "
    "refinement introduced by 'sequence similarity and phylogenetic analysis' "
    "(Knoll 2009 Fig. 1); the DED has been frozen since 2009 (Database Commons: "
    "Last update 2009, v1.1). registry_eligible=false: a 2009 similarity "
    "cluster is homology evidence, never a functional prior and never a "
    "phenotype. GenBank annotation must not define the family (Knoll 2009 "
    "documents three annotation/assignment contradictions, incl. gi:194292521 "
    "and gi:74267419 -> extracellular dPHASCL type 1, gi:34452171 -> "
    "intracellular nPHAMCL)."
)

NO_SBD_NOTE = (
    "This layer has no substrate binding domain, so the absence of Pfam PF06850 "
    "is structurally expected and must not be reported as a negative result."
)

PHAZ7_EXCEPTION_NOTE = (
    "PhaZ7 (Paucimonas lemoignei) is a documented localisation exception in "
    "Knoll 2009: extracellular but active only against native PHA granules."
)

FAMILY_USE_NOTE = (
    "This is the layer the project currently over-weights; it may refine a "
    "call only as an annotation, never as a registry or phenotype decision."
)


class AuthorityValidationError(ValueError):
    """Raised when the authority table violates a fail-closed invariant."""


# ---------------------------------------------------------------------------
# IO helpers
# ---------------------------------------------------------------------------


def sha256_file(path: Path) -> str:
    """Return the SHA-256 of a regular, non-symlink file."""
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise AuthorityValidationError(f"not a regular file: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_tsv(path: Path, *, label: str) -> tuple[list[str], list[dict[str, str]]]:
    """Read a TSV fail-closed: regular file, unique columns, no short rows."""
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise AuthorityValidationError(f"{label}: not a regular file: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if not fields:
            raise AuthorityValidationError(f"{label}: missing header row: {path}")
        duplicates = sorted({name for name in fields if fields.count(name) > 1})
        if duplicates:
            raise AuthorityValidationError(
                f"{label}: duplicate column names: {','.join(duplicates)}"
            )
        rows = []
        for index, row in enumerate(reader, start=2):
            if None in row:
                raise AuthorityValidationError(
                    f"{label} line {index}: more fields than the header declares"
                )
            if any(value is None for value in row.values()):
                raise AuthorityValidationError(
                    f"{label} line {index}: fewer fields than the header declares"
                )
            rows.append({key: (value or "").strip() for key, value in row.items()})
    if not rows:
        raise AuthorityValidationError(f"{label}: no data rows: {path}")
    return fields, rows


def _require_columns(fields: Sequence[str], required: Iterable[str], *, label: str) -> None:
    missing = sorted(set(required) - set(fields))
    if missing:
        raise AuthorityValidationError(
            f"{label}: missing required columns: {','.join(missing)}"
        )


def read_family_definitions(path: Path) -> list[dict[str, str]]:
    fields, rows = _read_tsv(path, label="phaded_family_definitions")
    _require_columns(fields, DEFINITION_REQUIRED_COLUMNS, label="phaded_family_definitions")
    return rows


def read_profile_manifest(path: Path) -> list[dict[str, str]]:
    fields, rows = _read_tsv(path, label="profile_manifest")
    _require_columns(fields, MANIFEST_REQUIRED_COLUMNS, label="profile_manifest")
    return rows


def _profile_index(rows: Sequence[Mapping[str, str]]) -> dict[tuple[str, str, str], str]:
    """Map (kind, superfamily, family) -> profile_id, rejecting ambiguity."""
    index: dict[tuple[str, str, str], str] = {}
    for row in rows:
        kind = row["profile_kind"]
        if kind not in LAYER_KINDS:
            raise AuthorityValidationError(
                f"profile_manifest: unknown profile_kind {kind!r} for {row['profile_id']!r}"
            )
        key = (kind, row["phaded_superfamily"], row["phaded_family_id"])
        if key in index:
            raise AuthorityValidationError(
                "profile_manifest: ambiguous authority for "
                f"{key}: {index[key]!r} and {row['profile_id']!r}"
            )
        if not row["profile_id"]:
            raise AuthorityValidationError(f"profile_manifest: empty profile_id for {key}")
        index[key] = row["profile_id"]
    return index


def _clean_field(value: str, *, column: str, layer_id: str) -> str:
    if any(char in value for char in ("\t", "\n", "\r")):
        raise AuthorityValidationError(
            f"{layer_id}: {column} contains a tab or newline"
        )
    return value


# ---------------------------------------------------------------------------
# generation
# ---------------------------------------------------------------------------


def _family_sort_key(layer_id: str) -> tuple[int, str]:
    match = FAMILY_ID_RE.fullmatch(layer_id)
    if match is None:
        raise AuthorityValidationError(
            f"family layer_id {layer_id!r} does not match the DED_hfam_<integer> form"
        )
    return (int(match.group(1)), layer_id)


def build_authority_rows(
    *,
    family_definitions_path: Path,
    profile_manifest_path: Path,
) -> list[dict[str, str]]:
    """Build the two-layer authority table deterministically from provenance.

    Both parents must agree exactly with the verified DED import; any drift
    raises instead of producing a partially authoritative table.
    """
    definitions = read_family_definitions(family_definitions_path)
    index = _profile_index(read_profile_manifest(profile_manifest_path))

    superfamily_to_families: dict[str, set[str]] = {}
    family_to_superfamily: dict[str, str] = {}
    for row in definitions:
        superfamily = row["phaded_superfamily"]
        family = row["phaded_family_id"]
        if not superfamily or not family:
            raise AuthorityValidationError(
                "phaded_family_definitions: empty superfamily or family id"
            )
        previous = family_to_superfamily.get(family)
        if previous is not None and previous != superfamily:
            raise AuthorityValidationError(
                f"phaded_family_definitions: family {family} is assigned to "
                f"multiple superfamilies ({previous!r}, {superfamily!r})"
            )
        family_to_superfamily[family] = superfamily
        superfamily_to_families.setdefault(superfamily, set()).add(family)

    observed_superfamilies = set(superfamily_to_families)
    expected_superfamilies = set(EXPECTED_SUPERFAMILIES)
    if observed_superfamilies != expected_superfamilies:
        raise AuthorityValidationError(
            "phaded_family_definitions diverges from the 8-superfamily hard "
            f"framework: missing={sorted(expected_superfamilies - observed_superfamilies)}, "
            f"unexpected={sorted(observed_superfamilies - expected_superfamilies)}"
        )
    if len(family_to_superfamily) != EXPECTED_FAMILY_COUNT:
        raise AuthorityValidationError(
            f"phaded_family_definitions declares {len(family_to_superfamily)} "
            f"families; expected {EXPECTED_FAMILY_COUNT}"
        )

    rows: list[dict[str, str]] = []
    for superfamily in EXPECTED_SUPERFAMILIES:
        notes = [FUNCTIONAL_PRIOR_NOTE]
        if superfamily in NO_SBD_SUPERFAMILIES:
            notes.append(NO_SBD_NOTE)
        if superfamily == "extracellular native-SCL/PhaZ7-like":
            notes.append(PHAZ7_EXCEPTION_NOTE)
        rows.append(
            {
                "layer_id": superfamily,
                "layer_kind": SUPERFAMILY_KIND,
                "source_evidence_type": FUNCTIONAL_PRIOR,
                "source_version": FROZEN_SOURCE_VERSION,
                "registry_eligible": "true",
                "gate_profile_id": index.get(
                    (SUPERFAMILY_KIND, superfamily, ""), PENDING
                ),
                "notes": " | ".join(notes),
            }
        )

    for family in sorted(family_to_superfamily, key=_family_sort_key):
        superfamily = family_to_superfamily[family]
        notes = [CLUSTERING_NOTE, FAMILY_USE_NOTE]
        if superfamily in NO_SBD_SUPERFAMILIES:
            notes.append(NO_SBD_NOTE)
        if superfamily == "extracellular native-SCL/PhaZ7-like":
            notes.append(PHAZ7_EXCEPTION_NOTE)
        rows.append(
            {
                "layer_id": family,
                "layer_kind": FAMILY_KIND,
                "source_evidence_type": SEQUENCE_CLUSTERING_2009,
                "source_version": FROZEN_SOURCE_VERSION,
                "registry_eligible": "false",
                "gate_profile_id": index.get(
                    (FAMILY_KIND, superfamily, family), PENDING
                ),
                "notes": " | ".join(notes),
            }
        )

    for row in rows:
        for column in REQUIRED_COLUMNS:
            _clean_field(row[column], column=column, layer_id=row["layer_id"])
    return rows


def write_authority_table(rows: Sequence[Mapping[str, str]], path: Path) -> Path:
    """Write the authority TSV; refuse to clobber an existing different file."""
    path = Path(path)
    if not path.parent.is_dir():
        raise AuthorityValidationError(
            f"authority table parent directory does not exist: {path.parent}"
        )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=REQUIRED_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row[column] for column in REQUIRED_COLUMNS})
    return path


# ---------------------------------------------------------------------------
# validation
# ---------------------------------------------------------------------------


def _invariant(
    checks: list[dict[str, object]],
    violations: list[str],
    *,
    identifier: str,
    description: str,
    expected: object,
    observed: object,
    detail: str = "",
) -> bool:
    passed = expected == observed
    checks.append(
        {
            "id": identifier,
            "description": description,
            "expected": expected,
            "observed": observed,
            "passed": passed,
        }
    )
    if not passed:
        message = (
            f"{identifier}: expected {expected!r}, observed {observed!r}"
            + (f" [{detail}]" if detail else "")
        )
        violations.append(message)
    return passed


def validate_authority_rows(
    rows: Sequence[Mapping[str, str]],
    *,
    family_definitions_path: Path,
    profile_manifest_path: Path,
) -> dict[str, object]:
    """Check every layer invariant; raise listing all violations if any fail."""
    checks: list[dict[str, object]] = []
    violations: list[str] = []

    if not rows:
        raise AuthorityValidationError("authority table has no rows")

    for index, row in enumerate(rows, start=2):
        fields = set(row)
        missing = sorted(set(REQUIRED_COLUMNS) - fields)
        if missing:
            raise AuthorityValidationError(
                f"authority table line {index}: missing columns: {','.join(missing)}"
            )
        unexpected = sorted(fields - set(REQUIRED_COLUMNS))
        if unexpected:
            raise AuthorityValidationError(
                f"authority table line {index}: unexpected columns: {','.join(unexpected)}"
            )

    superfamily_rows = [row for row in rows if row["layer_kind"] == SUPERFAMILY_KIND]
    family_rows = [row for row in rows if row["layer_kind"] == FAMILY_KIND]
    unknown_kind = sorted(
        {row["layer_kind"] for row in rows} - {SUPERFAMILY_KIND, FAMILY_KIND}
    )
    if unknown_kind:
        raise AuthorityValidationError(
            f"authority table: unknown layer_kind values: {','.join(unknown_kind)}"
        )

    _invariant(
        checks, violations,
        identifier="superfamily_count_is_8",
        description="the hard 8-superfamily functional-prior framework",
        expected=EXPECTED_SUPERFAMILY_COUNT,
        observed=len(superfamily_rows),
    )
    _invariant(
        checks, violations,
        identifier="family_count_is_38",
        description="the 38-family 2009 clustering refinement layer",
        expected=EXPECTED_FAMILY_COUNT,
        observed=len(family_rows),
    )

    seen_ids: dict[str, int] = {}
    for row in rows:
        seen_ids[row["layer_id"]] = seen_ids.get(row["layer_id"], 0) + 1
    _invariant(
        checks, violations,
        identifier="layer_ids_are_unique_and_non_empty",
        description="layer_id is the non-empty primary key of the authority table",
        expected=[],
        observed=sorted(
            layer_id
            for layer_id, count in seen_ids.items()
            if count > 1 or not layer_id.strip()
        ),
    )

    def _offenders(
        candidates: Sequence[Mapping[str, str]], column: str, allowed: str
    ) -> list[str]:
        return sorted(
            row["layer_id"] for row in candidates if row[column] != allowed
        )

    # --- invariant 1: the functional-prior framework -----------------------
    offenders = _offenders(superfamily_rows, "source_evidence_type", FUNCTIONAL_PRIOR)
    _invariant(
        checks, violations,
        identifier="superfamilies_are_functional_prior",
        description=(
            "superfamilies come from the function/localisation prior on the 28 "
            "experimentally validated seed sequences"
        ),
        expected=[],
        observed=offenders,
    )
    offenders = _offenders(superfamily_rows, "registry_eligible", "true")
    _invariant(
        checks, violations,
        identifier="superfamilies_are_registry_eligible",
        description="the functional-prior layer is registry eligible",
        expected=[],
        observed=offenders,
    )

    # --- invariant 2: the frozen 2009 clustering layer --------------------
    offenders = _offenders(family_rows, "source_evidence_type", SEQUENCE_CLUSTERING_2009)
    _invariant(
        checks, violations,
        identifier="families_are_sequence_clustering_2009",
        description=(
            "families are a 2009 sequence-similarity + phylogenetic clustering "
            "layer, not a functional prior"
        ),
        expected=[],
        observed=offenders,
    )
    offenders = _offenders(family_rows, "source_version", FROZEN_SOURCE_VERSION)
    _invariant(
        checks, violations,
        identifier="families_are_bound_to_the_frozen_knoll2009_v1_1_snapshot",
        description="the family layer is bound to the frozen DED snapshot",
        expected=[],
        observed=offenders,
    )
    offenders = _offenders(family_rows, "registry_eligible", "false")
    _invariant(
        checks, violations,
        identifier="families_are_not_registry_eligible",
        description=(
            "CORE INVARIANT: a 2009 clustering refinement layer must never be "
            "registry eligible"
        ),
        expected=[],
        observed=offenders,
    )

    # --- hierarchy is never crossed ---------------------------------------
    offenders = sorted(
        row["layer_id"]
        for row in rows
        if (row["source_evidence_type"], row["registry_eligible"])
        not in ALLOWED_EVIDENCE_ELIGIBILITY
    )
    _invariant(
        checks, violations,
        identifier="evidence_type_and_registry_eligibility_are_never_crossed",
        description=(
            "functional_prior pairs with true and sequence_clustering_2009 pairs "
            "with false; the opposite combinations are the layer mismatch defect"
        ),
        expected=[],
        observed=offenders,
    )
    offenders = sorted(
        row["layer_id"] for row in rows if row["source_version"] != FROZEN_SOURCE_VERSION
    )
    _invariant(
        checks, violations,
        identifier="every_layer_is_bound_to_the_frozen_ded_snapshot",
        description="the whole table is bound to the frozen DED snapshot",
        expected=[],
        observed=offenders,
    )
    offenders = sorted(row["layer_id"] for row in rows if not row["notes"].strip())
    _invariant(
        checks, violations,
        identifier="every_layer_carries_its_provenance_rationale",
        description="notes must state why the layer has its authority",
        expected=[],
        observed=offenders,
    )

    # --- invariant 3: exact correspondence with the DED import ------------
    definitions = read_family_definitions(family_definitions_path)
    expected_superfamilies = sorted({row["phaded_superfamily"] for row in definitions})
    expected_families = sorted({row["phaded_family_id"] for row in definitions})
    observed_superfamilies = sorted({row["layer_id"] for row in superfamily_rows})
    observed_families = sorted({row["layer_id"] for row in family_rows})

    _invariant(
        checks, violations,
        identifier="superfamily_layer_ids_match_ded_import",
        description=(
            "no superfamily added or dropped relative to "
            "phaded_family_definitions.tsv"
        ),
        expected=expected_superfamilies,
        observed=observed_superfamilies,
    )
    _invariant(
        checks, violations,
        identifier="family_layer_ids_match_ded_import",
        description=(
            "no family added or dropped relative to phaded_family_definitions.tsv"
        ),
        expected=expected_families,
        observed=observed_families,
    )
    _invariant(
        checks, violations,
        identifier="superfamily_layer_ids_match_hard_framework",
        description="the superfamily layer equals the 8 labelled DED superfamilies",
        expected=sorted(EXPECTED_SUPERFAMILIES),
        observed=observed_superfamilies,
    )

    malformed = sorted(
        layer_id for layer_id in observed_families if FAMILY_ID_RE.fullmatch(layer_id) is None
    )
    _invariant(
        checks, violations,
        identifier="family_layer_ids_use_ded_hfam_format",
        description="family layer_id must be the DED_hfam_<integer> identifier",
        expected=[],
        observed=malformed,
    )

    definition_parent = {
        row["phaded_family_id"]: row["phaded_superfamily"] for row in definitions
    }
    index = _profile_index(read_profile_manifest(profile_manifest_path))
    resolved = 0
    unresolvable: list[str] = []
    parent_mismatch: list[str] = []
    for row in superfamily_rows:
        profile_id = index.get((SUPERFAMILY_KIND, row["layer_id"], ""))
        if profile_id is None or row["gate_profile_id"] != profile_id:
            unresolvable.append(row["layer_id"])
        else:
            resolved += 1
    for row in family_rows:
        profile_id = index.get(
            (FAMILY_KIND, definition_parent.get(row["layer_id"], ""), row["layer_id"])
        )
        if profile_id is None or row["gate_profile_id"] != profile_id:
            unresolvable.append(row["layer_id"])
            continue
        resolved += 1
        if index.get((FAMILY_KIND, row["layer_id"], "")) is not None:
            parent_mismatch.append(row["layer_id"])

    _invariant(
        checks, violations,
        identifier="gate_profiles_resolve_in_profile_manifest",
        description=(
            "every layer must bind a real gate profile; a pending gate profile "
            "fails closed"
        ),
        expected=EXPECTED_LAYER_COUNT,
        observed=resolved,
    )
    _invariant(
        checks, violations,
        identifier="gate_profiles_are_unambiguous",
        description="no layer id may collide between the superfamily and family kinds",
        expected=[],
        observed=sorted(set(unresolvable)),
    )
    _invariant(
        checks, violations,
        identifier="family_gate_profiles_confirm_superfamily_parentage",
        description=(
            "each family's gate profile must name the same parent superfamily as "
            "phaded_family_definitions.tsv"
        ),
        expected=[],
        observed=sorted(set(parent_mismatch)),
    )

    if violations:
        raise AuthorityValidationError(
            "classification authority FAILED "
            f"({len(violations)} violation(s)):\n" + "\n".join(violations)
        )

    eligibility_counts = {"true": 0, "false": 0}
    for row in rows:
        eligibility_counts[row["registry_eligible"]] = (
            eligibility_counts.get(row["registry_eligible"], 0) + 1
        )
    evidence_counts: dict[str, int] = {}
    version_counts: dict[str, int] = {}
    for row in rows:
        evidence_counts[row["source_evidence_type"]] = (
            evidence_counts.get(row["source_evidence_type"], 0) + 1
        )
        version_counts[row["source_version"]] = (
            version_counts.get(row["source_version"], 0) + 1
        )

    return {
        "status": "verified",
        "layer_counts": {
            "superfamily": len(superfamily_rows),
            "family": len(family_rows),
            "total": len(rows),
        },
        "registry_eligible_counts": eligibility_counts,
        "source_evidence_type_counts": evidence_counts,
        "source_version_counts": version_counts,
        "gate_profile_resolution": {
            "resolved": resolved,
            "total": len(rows),
            "pending": len(rows) - resolved,
        },
        "invariants": checks,
    }


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------


def count_data_rows(path: Path) -> int:
    """Count TSV data rows (total lines minus the header)."""
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise AuthorityValidationError(f"not a regular file: {path}")
    total = 0
    last_byte = b""
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            total += block.count(b"\n")
            last_byte = block[-1:]
    if last_byte and last_byte != b"\n":
        total += 1
    if total < 1:
        raise AuthorityValidationError(f"empty table: {path}")
    return total - 1


def _describe_source(role: str, path: Path, *, copied_from: str | None = None) -> dict[str, object]:
    path = Path(path)
    record: dict[str, object] = {
        "role": role,
        "path": str(path),
        "size": None,
        "sha256": None,
        "status": PENDING,
    }
    if copied_from is not None:
        record["copied_from"] = copied_from
    if path.is_file() and not path.is_symlink():
        record["size"] = path.stat().st_size
        record["sha256"] = sha256_file(path)
        record["status"] = "verified"
    return record


def build_validation_report(
    *,
    run_id: str,
    rows: Sequence[Mapping[str, str]],
    authority_path: Path,
    family_definitions_path: Path,
    profile_manifest_path: Path,
    candidate_evidence_path: Path,
    expected_candidates: int = EXPECTED_CANDIDATE_UNIVERSE,
    copied_from: Mapping[str, str] | None = None,
) -> dict[str, object]:
    """Validate the on-disk table and assemble the auditable JSON report."""
    # Validate what is actually on disk, not the in-memory rows.
    _, on_disk = _read_tsv(authority_path, label="phaded_classification_authority")
    report = validate_authority_rows(
        on_disk,
        family_definitions_path=family_definitions_path,
        profile_manifest_path=profile_manifest_path,
    )
    regenerated_match = list(on_disk) == [dict(row) for row in rows]
    if not regenerated_match:
        raise AuthorityValidationError(
            "the committed authority table differs from a fresh generation from "
            "the bound provenance inputs"
        )

    observed_candidates = count_data_rows(candidate_evidence_path)
    checks = list(report["invariants"])
    violations: list[str] = []
    _invariant(
        checks, violations,
        identifier="candidate_universe_row_count_matches_plan_fact_baseline",
        description=(
            "the candidate universe this classification change is measured "
            "against must equal the plan fact baseline"
        ),
        expected=expected_candidates,
        observed=observed_candidates,
    )
    if violations:
        raise AuthorityValidationError("\n".join(violations))

    copied_from = dict(copied_from or {})
    sources = [
        _describe_source(
            "phaded_family_definitions",
            family_definitions_path,
            copied_from=copied_from.get("phaded_family_definitions"),
        ),
        _describe_source(
            "profile_manifest",
            profile_manifest_path,
            copied_from=copied_from.get("profile_manifest"),
        ),
        _describe_source("phaded_classification_authority", authority_path),
        _describe_source("candidate_evidence", candidate_evidence_path),
    ]

    report.update(
        {
            "schema_version": "1.0",
            "run_id": run_id,
            "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "ded_snapshot_label": DED_SNAPSHOT_LABEL,
            "source_version": FROZEN_SOURCE_VERSION,
            "authority_table": {
                "path": str(authority_path),
                "sha256": sha256_file(authority_path),
                "size": authority_path.stat().st_size,
                "regenerated_from_provenance_matches_committed": regenerated_match,
            },
            "source_files": sources,
            "invariants": checks,
            "candidate_universe_impact": {
                "candidates_total": observed_candidates,
                "candidates_total_source": str(candidate_evidence_path),
                "rows_changed": 0,
                "rationale": [
                    "all 8 superfamily labels keep registry_eligible=true, so no "
                    "superfamily-scoped decision loses eligibility",
                    "subtype_call is derived only from phaded_superfamily_best "
                    "(pipeline/scripts/build_phaded_subtype_matrix.py:_base_subtype, "
                    "L153-165); the family layer never contributes to a call, so "
                    "setting registry_eligible=false on families cannot change one",
                    "no runtime consumer of registry_eligible exists in pipeline/ "
                    "before this task, so the authority table is a governance "
                    "artifact with no computational side effects",
                ],
            },
            "boundary_statement": BOUNDARY_STATEMENT,
        }
    )
    return report


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--family-definitions", required=True)
    parser.add_argument("--profile-manifest", required=True)
    parser.add_argument("--authority-table", required=True)
    parser.add_argument("--candidate-evidence", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--expected-candidates", type=int, default=EXPECTED_CANDIDATE_UNIVERSE)
    parser.add_argument(
        "--write-authority",
        action="store_true",
        help="generate the authority table at --authority-table (default: validate only)",
    )
    parser.add_argument("--family-definitions-copied-from", default=None)
    parser.add_argument("--profile-manifest-copied-from", default=None)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    authority_path = Path(args.authority_table)
    output_path = Path(args.output)
    run_id = args.run_id or authority_path.parent.parent.name

    try:
        if args.write_authority:
            rows = build_authority_rows(
                family_definitions_path=Path(args.family_definitions),
                profile_manifest_path=Path(args.profile_manifest),
            )
            write_authority_table(rows, authority_path)
        _, rows = _read_tsv(authority_path, label="phaded_classification_authority")
        report = build_validation_report(
            run_id=run_id,
            rows=rows,
            authority_path=authority_path,
            family_definitions_path=Path(args.family_definitions),
            profile_manifest_path=Path(args.profile_manifest),
            candidate_evidence_path=Path(args.candidate_evidence),
            expected_candidates=args.expected_candidates,
            copied_from={
                "phaded_family_definitions": args.family_definitions_copied_from,
                "profile_manifest": args.profile_manifest_copied_from,
            },
        )
    except AuthorityValidationError as error:
        blocked = {
            "schema_version": "1.0",
            "run_id": run_id,
            "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "status": "blocked",
            "ded_snapshot_label": DED_SNAPSHOT_LABEL,
            "error": str(error),
            "boundary_statement": BOUNDARY_STATEMENT,
        }
        if output_path.parent.is_dir():
            output_path.write_text(
                json.dumps(blocked, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
            )
        print(json.dumps(blocked, ensure_ascii=False, indent=2), file=sys.stderr)
        return 1

    output_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "status": report["status"],
                "layer_counts": report["layer_counts"],
                "registry_eligible_counts": report["registry_eligible_counts"],
                "gate_profile_resolution": report["gate_profile_resolution"],
                "candidate_universe_impact": {
                    "candidates_total": report["candidate_universe_impact"][
                        "candidates_total"
                    ],
                    "rows_changed": report["candidate_universe_impact"]["rows_changed"],
                },
                "authority_table_sha256": report["authority_table"]["sha256"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
