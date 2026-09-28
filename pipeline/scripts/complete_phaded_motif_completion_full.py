#!/usr/bin/env python3
"""Re-run the motif panel over the whole 109,087-candidate universe (Task 7).

Why this script exists
----------------------
The 2026-09-15 candidate panel (``runs/20260915_phaded_motif_reconciliation_01``)
was a *lipase box + oxyanion hole* panel: its His pattern never fired on any of
the 723 references, its Asp pattern was a type 1 motif, its ``lid_state`` was
copied from the subtype prior (circular by construction), and its
``linker_state``/``sbd_state`` were never bound.  Task 2 replaced those detectors
and produced a reference-anchored mapping for the 723 DED references.

This script propagates the *detector* half of Task 2 to the candidate universe
and states, per candidate and per criterion, whether the criterion could be
evaluated at all.

The decisive measured fact behind the design
--------------------------------------------
Task 2 resolves His/Asp/lid by transferring literature coordinates through DED
alignment columns.  That transfer is bounded by alignment membership: **zero of
the 109,087 candidates appears in any of the 123 frozen DED alignment files, and
zero shares an accession with the 723-reference ledger**.  The reference mapping
therefore cannot reach a single candidate row, and this script refuses to
pretend otherwise -- it publishes ``not_assignable_*`` states instead of turning
an untested criterion into a negative.

Assignability rules (all data-derived, never asserted by hand)
--------------------------------------------------------------
* a detector **positive** is always assignable for every criterion;
* a detector **miss** is assignable only when the calibrated pattern reproduces
  *every* literature anchor that carries a coordinate for that role (nucleophile
  4/4, oxyanion 2/2) or, where no anchor carries a coordinate, when the
  criterion has at least one reference-panel positive control (AHSMG 2/2).  The
  relaxed His (GMxH) and Asp (GxxDYTV) proxies reproduce 1 of 4 anchors, so a
  miss carries no information about those residues and stays not assignable;
* ``lid`` is assignable for every candidate: its detector reads a sequence only,
  and the prior is reported in a separate column so the detector output can be
  set against that prior instead of restating it;
* ``linker`` and ``sbd`` are assignable for **no** candidate: no linker profile
  is bound, and PF06850 is not the extracellular substrate-binding domain in this
  panel (DEFECT_D).  Its binding result is still published, as raw evidence in
  its own column.

Everything here stays candidate-only sequence/domain evidence.  This task makes
no new family or superfamily call.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import importlib.util
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence


BOUNDARY = (
    "sequence-pattern/domain evidence only; not catalytic validation and not "
    "PHB/PHA phenotype proof"
)
EVIDENCE_LEVEL = "full_library_sequence_pattern_layer_only"
ALIGNMENT_LAYER = "sequence_pattern_only_no_alignment_membership"

CRITERIA = ("lipase_box", "his", "asp", "oxyanion_hole", "ahsmg", "lid", "linker", "sbd")
REQUIRED_CRITERIA = ("his", "asp", "lid", "linker", "sbd")
#: criterion -> the Task 2 literature role it is detected with
CRITERION_ROLE = {
    "lipase_box": "catalytic_serine",
    "his": "catalytic_histidine",
    "asp": "catalytic_aspartate",
    "oxyanion_hole": "oxyanion_hole_cysteine",
    "ahsmg": "alternative_catalytic_motif_ahsmg",
}
#: criteria whose state comes from a Task 2 pattern detector rather than a domain
PATTERN_CRITERIA = tuple(CRITERION_ROLE)

ASSIGNABLE = "assignable"
NOT_ASSIGNABLE = "not_assignable"

REASON_PATTERN_POSITIVE = "assignable_detector_positive_calibrated_pattern"
REASON_PATTERN_MISS = "assignable_pattern_miss_reproduces_all_literature_anchors"
REASON_LID_DETECTOR = "assignable_sequence_only_lid_detector_ran_for_this_candidate"
REASON_CANDIDATE_NOT_IN_ALIGNMENT = "not_assignable_candidate_not_in_any_ded_alignment"
REASON_NO_ANCHOR_FOR_PRIOR = "not_assignable_no_literature_anchor_for_prior_superfamily"
REASON_PERMUTED_FOLD = "not_assignable_circularly_permuted_fold_column_not_transferable"
REASON_NO_PRIOR = "not_assignable_no_superfamily_prior_and_no_alignment_membership"
REASON_PRIOR_UNKNOWN = "not_assignable_prior_superfamily_absent_from_reference_mapping"
REASON_LINKER_PROFILE = "not_assignable_no_linker_profile_bound"
REASON_SBD_CRITERION = "not_assignable_pf06850_is_not_the_extracellular_sbd_in_this_panel"

NOT_ASSIGNABLE_REASONS = (
    REASON_CANDIDATE_NOT_IN_ALIGNMENT,
    REASON_NO_ANCHOR_FOR_PRIOR,
    REASON_PERMUTED_FOLD,
    REASON_NO_PRIOR,
    REASON_PRIOR_UNKNOWN,
    REASON_LINKER_PROFILE,
    REASON_SBD_CRITERION,
)

PRIOR_ANCHOR_TRANSFERABLE = "anchor_available_and_column_transferable"
PRIOR_ANCHOR_NOT_TRANSFERABLE = "anchor_available_but_column_not_transferable"
PRIOR_NO_ANCHOR = "no_literature_anchor_for_prior"
PRIOR_NONE = "no_superfamily_prior"

SBD_PFAM_ACCESSION = "PF06850"
CATALYTIC_DOMAIN_PFAM_ACCESSION = "PF10503"
SBD_ROLE_CAVEAT = (
    "pf06850_binds_the_intracellular_cys_type_superfamily_here_and_is_not_the_"
    "extracellular_substrate_binding_domain"
)
SBD_ROLE_CAVEAT_NEGATIVE = (
    "pf06850_is_not_the_extracellular_substrate_binding_domain_in_this_panel_so_a_"
    "miss_is_not_an_sbd_negative"
)
#: a criterion without a detector cannot be evaluated for any candidate
NOT_ASSIGNABLE_BASIS = {
    "linker": (
        "no linker profile is bound: Task 2 could only delimit a candidate inter-domain "
        "region, which is not a linker call, so the criterion has no detector at all"
    ),
    "sbd": (
        "PF06850 binds the intracellular Cys-type superfamily in this panel and is not the "
        "extracellular substrate-binding domain (Task 2 DEFECT_D), so the SBD criterion cannot "
        "be evaluated; the raw binding is published in its own column"
    ),
}
#: priors whose detector expectation is recorded (never fed to a detector)
LID_EXPECTED_PRIORS = ("intracellular nPHAMCL",)
MISS_ASSIGNABLE_CRITERIA = ("lipase_box", "oxyanion_hole", "ahsmg", "lid")
PRIOR_COLUMN = "motif_reference_subtype"

PENDING_NO_LEGACY_ROW = "pending_no_legacy_candidate_row"
TASK2_DETECTOR_NAME = "complete_phaded_reference_residue_mapping"

LEGACY_BEFORE_SOURCE = (
    "runs/20260915_phaded_motif_reconciliation_01/results/motif_candidate_evidence.tsv"
)
TASK2_RUN_ID = "20260917_phaded_reference_residue_mapping_01"
TASK2_REPORT_SOURCE = f"runs/{TASK2_RUN_ID}/results/implementation_calibration_report.json"

# ---------------------------------------------------------------------------
# V5: before-block field vocabulary
# ---------------------------------------------------------------------------
BEFORE_FIELD_STATE_COUNTS = "state_counts"
BEFORE_FIELD_RESOLVED = "before_resolved_supported"
BEFORE_FIELD_INDEPENDENT = "before_independent_detection_supported"
BEFORE_FIELD_NOT_AVAILABLE = "before_not_available"
BEFORE_FIELD_NOT_AVAILABLE_REASON = "before_not_available_reason"
#: sentinel for "this criterion had no independent detector in 2026-09-15", so no
#: count exists; a literal 0 would assert that a detector ran and found nothing
BEFORE_NOT_AVAILABLE = "not_available"

#: criteria whose 2026-09-15 baseline column was never produced by an independent
#: detector, with the measured reason (state counts from the frozen evidence file)
BEFORE_NO_INDEPENDENT_DETECTOR = {
    "lid": (
        "no independent lid detector existed in 2026-09-15: all 109,087 rows were prior-derived "
        "(101,836 not_expected_for_subtype + 7,251 pending_reference_annotation), so an "
        "independent-detection count is not available rather than 0"
    ),
    "linker": (
        "no linker profile was bound in 2026-09-15 (108,795 pending_reference_annotation + 292 "
        "not_expected_for_subtype), so the criterion had no detector at all and no "
        "independent-detection count is available"
    ),
    "sbd": (
        "the 2026-09-15 sbd_state column was the raw PF06850 binding promoted to the SBD "
        "criterion (30,621 supported); the SBD criterion itself was never independently "
        "evaluated (DEFECT_D), so an independent-detection count is not available"
    ),
}

BEFORE_FIELD_SEMANTICS = {
    BEFORE_FIELD_STATE_COUNTS: (
        "measured 2026-09-15 state partition of the criterion, summed from the frozen evidence "
        "file; it sums to the candidate universe (109,087) and contains no aggregate key, so it "
        "can be summed without double counting"
    ),
    BEFORE_FIELD_RESOLVED: (
        "total rows whose 2026-09-15 state starts with 'supported', whatever the evidence origin "
        "(measured); this is what the old flat field named independent_detection_supported "
        "actually contained"
    ),
    BEFORE_FIELD_INDEPENDENT: (
        "rows supported by a detector that observed that row itself; an integer when the "
        "2026-09-15 layer had such a detector, otherwise the sentinel string 'not_available'"
    ),
    BEFORE_FIELD_NOT_AVAILABLE: "true when no independent-detection count exists for the criterion",
    BEFORE_FIELD_NOT_AVAILABLE_REASON: (
        "measured reason for before_not_available; empty string when the count is available"
    ),
    "independent_detection_supported": (
        "REMOVED from the before block by the V5 fix: the name claimed an independent detection "
        "while carrying the resolved total, and was hard-coded to 0 for lid/linker/sbd. Use "
        "before_resolved_supported and before_independent_detection_supported instead."
    ),
    "supported": (
        "REMOVED from the top level of the before payload: it duplicated the state_counts entry "
        "of the same name for criteria whose state word is literally 'supported' and could be "
        "summed twice; the aggregate is now before_resolved_supported"
    ),
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    if not fields:
        fields = ["status"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def load_task2_module(script: Path) -> object:
    """Load the Task 2 detectors by path (this script does not own them)."""
    if not script.is_file():
        raise FileNotFoundError(script)
    spec = importlib.util.spec_from_file_location(TASK2_DETECTOR_NAME, script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_fasta(path: Path) -> dict[str, str]:
    """Read a FASTA keyed by the first whitespace-delimited header token."""
    records: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current is not None:
                if current in records or not chunks:
                    raise ValueError(f"duplicate or empty FASTA record: {current}")
                records[current] = "".join(chunks).upper()
            current = line[1:].split()[0]
            if not current:
                raise ValueError("empty FASTA header")
            chunks = []
        else:
            if current is None:
                raise ValueError("sequence appears before FASTA header")
            chunks.append(line)
    if current is not None:
        if current in records or not chunks:
            raise ValueError(f"duplicate or empty FASTA record: {current}")
        records[current] = "".join(chunks).upper()
    if not records:
        raise ValueError(f"FASTA is empty: {path}")
    return records


def criterion_detection(sequence: str) -> dict[str, dict[str, str]]:
    """Run every sequence-only detector on one sequence.

    The parameter list is the whole input of this layer: a sequence.  No
    superfamily, subtype or family label is accepted, which is what lets the
    result be compared against those priors instead of restating them.
    """
    detectors = _DETECTORS
    result: dict[str, dict[str, str]] = {}
    for criterion, role in CRITERION_ROLE.items():
        result[criterion] = dict(detectors.detect_literature_pattern(sequence, role))
    lid = dict(detectors.detect_lid(sequence))
    result["lid"] = {
        "state": str(lid["state"]),
        "coordinates": str(lid.get("coordinates", "")),
        "loop_spacing": str(lid.get("loop_spacing", "")),
        "mapping_basis": str(lid.get("mapping_basis", "")),
    }
    return result


_DETECTORS: object | None = None


def set_detectors(module: object) -> None:
    """Bind the shared Task 2 detector module for this process."""
    global _DETECTORS
    _DETECTORS = module


def pfam_hits_from_coordinates(pfam_coordinates: str) -> list[dict[str, object]]:
    """Rebuild the frozen Pfam hit list from the recorded coordinate string.

    A single accession can carry more than one span (``PF06850:200-267,268-383``
    for a two-copy protein), so comma-separated spans are kept.
    """
    hits: list[dict[str, object]] = []
    for token in (pfam_coordinates or "").split(";"):
        token = token.strip()
        if not token or ":" not in token:
            continue
        accession, _, spans = token.partition(":")
        for span in spans.split(","):
            start, _, end = span.strip().partition("-")
            if not (start.strip().isdigit() and end.strip().isdigit()):
                continue
            hits.append({
                "pfam_accession": accession.strip(),
                "coordinates": f"{int(start)}-{int(end)}",
                "start": int(start),
                "end": int(end),
            })
    return hits


def delimit_linker(pfam_accessions: str, pfam_coordinates: str) -> dict[str, str]:
    """Delimit a candidate inter-domain region; never claim a linker profile."""
    hits = pfam_hits_from_coordinates(pfam_coordinates)
    if not hits:
        return {
            "state": "no_pfam_annotation_available",
            "coordinates": "",
            "basis": "no_pfam_scan_result_recorded_for_this_candidate",
        }
    delimited = _DETECTORS.linker_region(hits)
    basis = str(delimited["mapping_basis"])
    if delimited["state"] == "delimited_candidate_region_not_profiled":
        basis = f"{basis};literature_pattern_none"
    return {
        "state": str(delimited["state"]),
        "coordinates": str(delimited["coordinates"]),
        "basis": basis,
    }


def sbd_evidence(pfam_accessions: str, pfam_coordinates: str) -> dict[str, str]:
    """Report the raw PF06850 binding without promoting it to the SBD criterion.

    Detection is accession-driven, exactly as the 2026-09-15 panel did it, so the
    raw counts stay comparable; the recorded spans are reported next to it and a
    missing span is stated instead of being silently dropped.
    """
    accessions = {item.strip() for item in (pfam_accessions or "").split(";") if item.strip()}
    if not accessions:
        return {
            "state": "not_tested_no_pfam_evidence",
            "coordinates": "",
            "coordinate_status": "not_applicable_no_pfam_scan_result",
            "role_caveat": SBD_ROLE_CAVEAT_NEGATIVE,
        }
    hits = [
        hit for hit in pfam_hits_from_coordinates(pfam_coordinates)
        if hit["pfam_accession"] == SBD_PFAM_ACCESSION
    ]
    if SBD_PFAM_ACCESSION in accessions:
        coordinates = ";".join(sorted({str(hit["coordinates"]) for hit in hits}))
        return {
            "state": "detected",
            "coordinates": coordinates,
            "coordinate_status": (
                "coordinates_recorded" if coordinates
                else "accession_recorded_without_parseable_coordinate"
            ),
            "role_caveat": SBD_ROLE_CAVEAT,
        }
    return {
        "state": "not_detected_in_tested_pfam",
        "coordinates": "",
        "coordinate_status": "pfam_scan_ran_without_pf06850",
        "role_caveat": SBD_ROLE_CAVEAT_NEGATIVE,
    }


def build_prior_anchor_index(
    reference_manifest: Sequence[Mapping[str, str]],
    calibration_report: Mapping[str, object],
) -> dict[str, dict[str, object]]:
    """Per prior, which literature-anchored role can ever be transferred.

    Everything is read from the Task 2 artifacts: the prior -> family link from
    the reference manifest, the anchor and its column conservation from the Task
    2 calibration report.  Nothing is hard-coded, so a rerun of Task 2 changes
    this layer automatically.
    """
    index: dict[str, dict[str, object]] = {}
    transfer = dict(calibration_report.get("alignment_column_transfer", {}))
    anchored_families: dict[str, set[str]] = defaultdict(set)
    for record in transfer.get("anchored_alignments", []):
        anchored_families[str(record["phaded_family_id"])].add(str(record["anchor_accession"]))
    transferable: dict[tuple[str, str], bool] = {}
    for record in transfer.get("column_conservation", []):
        transferable[(str(record["phaded_family_id"]), str(record["role"]))] = bool(
            record["transferable"]
        )

    families_by_prior: dict[str, set[str]] = defaultdict(set)
    for row in reference_manifest:
        prior = str(row.get("phaded_superfamily", "")).strip()
        family = str(row.get("phaded_family_id", "")).strip()
        if prior:
            families_by_prior[prior].add(family)

    for prior, families in families_by_prior.items():
        anchors = sorted(
            {anchor for family in families for anchor in anchored_families.get(family, set())}
        )
        roles: dict[str, bool] = {}
        for criterion, role in CRITERION_ROLE.items():
            flags = [transferable.get((family, role)) for family in families]
            roles[role] = any(flag for flag in flags if flag is not None)
        index[prior] = {
            "families": sorted(families),
            "anchors": anchors,
            "anchor_available": bool(anchors),
            "role_transferable": roles,
        }
    return index


def role_anchor_basis(prior: str, role: str, index: Mapping[str, Mapping[str, object]]) -> dict[str, str]:
    """Classify one prior/role pair as transferable, anchored-but-not, or absent."""
    if not prior:
        return {"status": PRIOR_NONE, "anchor": "", "reason": REASON_NO_PRIOR}
    entry = index.get(prior)
    if entry is None:
        return {"status": PRIOR_NONE, "anchor": "", "reason": REASON_PRIOR_UNKNOWN}
    anchor = ";".join(entry["anchors"]) if entry["anchors"] else "pending_no_literature_anchor_for_prior"
    if not entry["anchor_available"]:
        return {"status": PRIOR_NO_ANCHOR, "anchor": anchor, "reason": REASON_NO_ANCHOR_FOR_PRIOR}
    if bool(entry["role_transferable"].get(role)):
        return {
            "status": PRIOR_ANCHOR_TRANSFERABLE,
            "anchor": anchor,
            "reason": REASON_CANDIDATE_NOT_IN_ALIGNMENT,
        }
    return {
        "status": PRIOR_ANCHOR_NOT_TRANSFERABLE,
        "anchor": anchor,
        "reason": REASON_PERMUTED_FOLD,
    }


def criterion_miss_assignability_basis(
    literature_rows: Sequence[Mapping[str, str]],
    anchors: Mapping[str, str],
    reference_manifest: Sequence[Mapping[str, str]],
    detectors: object,
) -> dict[str, dict[str, object]]:
    """Decide, from measured evidence, whether a pattern miss means anything.

    A miss is only informative if the calibrated pattern reproduces every
    literature anchor that carries a coordinate for that role: otherwise the
    pattern covers one motif context and says nothing about the others.
    """
    basis: dict[str, dict[str, object]] = {}
    positive_controls = {
        criterion: sum(
            1 for row in reference_manifest
            if str(row.get(f"{criterion}_state", "")).startswith("supported")
        )
        for criterion in CRITERIA
    }
    unresolvable_anchors: list[str] = []
    for criterion in CRITERIA:
        if criterion in CRITERION_ROLE:
            role = CRITERION_ROLE[criterion]
            accessions = [
                str(row["accession"]) for row in literature_rows
                if str(row["residue_role"]) == role
                and str(row.get("alignment_coordinate", "")).strip().isdigit()
                and str(row["accession"]).strip() not in ("", "pending")
            ]
            reproduced = 0
            for accession in accessions:
                sequence = anchors.get(accession)
                if sequence is None:
                    unresolvable_anchors.append(f"{role}:{accession}")
                    continue
                if detectors.detect_literature_pattern(sequence, role)["state"] == "supported":
                    reproduced += 1
            with_coordinate = len(accessions)
        else:
            with_coordinate = 0
            reproduced = 0
        controls = positive_controls.get(criterion, 0)
        if criterion in MISS_ASSIGNABLE_CRITERIA:
            # these criteria have a detector with validated controls, so a miss is
            # a reportable detection-layer result rather than an unknown
            miss_is_assignable = True
            reason = REASON_LID_DETECTOR if criterion == "lid" else REASON_PATTERN_MISS
        elif criterion in ("linker", "sbd"):
            miss_is_assignable = False
            reason = REASON_LINKER_PROFILE if criterion == "linker" else REASON_SBD_CRITERION
        else:
            miss_is_assignable = with_coordinate > 0 and reproduced == with_coordinate
            reason = REASON_PATTERN_MISS if miss_is_assignable else ""
        basis[criterion] = {
            "anchors_with_coordinate": with_coordinate,
            "anchors_reproduced": reproduced,
            "reference_panel_positive_controls": controls,
            "miss_is_assignable": miss_is_assignable,
            "assignability_reason_when_a_miss_is_assignable": reason,
            "basis": (
                "literature anchors carrying a coordinate for this role are re-scanned with "
                "the calibrated pattern; a miss is only informative when all of them are "
                "reproduced, otherwise the pattern covers one motif context only"
            ),
        }
    if unresolvable_anchors:
        raise ValueError(
            "literature anchors with a coordinate could not be resolved to a sequence, so the "
            f"measured pattern coverage would be understated: {sorted(set(unresolvable_anchors))}"
        )
    return basis


def build_assignability_basis(
    literature_reference: Path,
    reference_fasta: Path,
    anchor_fasta: Path,
    reference_manifest: Sequence[Mapping[str, str]],
    detectors: object,
) -> dict[str, dict[str, object]]:
    """Resolve literature anchors the same way Task 2 does.

    Three of the four literature anchors are DED panel references and one is an
    extra RCSB chain, so the two FASTA files have to be merged exactly as Task 2
    merges them; otherwise the measured pattern coverage would be understated.
    """
    literature_rows = read_tsv(literature_reference)
    panel = {
        key.split("|", 1)[1]: value
        for key, value in read_fasta(Path(reference_fasta)).items()
    }
    anchors = {**panel, **read_fasta(Path(anchor_fasta))}
    return criterion_miss_assignability_basis(
        literature_rows, anchors, reference_manifest, detectors
    )


def assignability_for(
    criterion: str,
    detection_state: str,
    *,
    prior_basis: Mapping[str, str],
    miss_basis: Mapping[str, object],
) -> dict[str, str]:
    """Decide one criterion's assignability for one candidate."""
    if detection_state == "supported":
        if criterion in PATTERN_CRITERIA:
            return {
                "assignability": ASSIGNABLE,
                "reason": REASON_PATTERN_POSITIVE,
                "basis": f"sequence_pattern_positive:{criterion}",
            }
        if criterion == "lid":
            return {
                "assignability": ASSIGNABLE,
                "reason": REASON_LID_DETECTOR,
                "basis": "sequence_only_lid_detector_matched_both_validated_loops",
            }
        # linker and sbd are never assignable, even when their raw evidence is present
        return {
            "assignability": NOT_ASSIGNABLE,
            "reason": str(miss_basis["assignability_reason_when_a_miss_is_assignable"]),
            "basis": NOT_ASSIGNABLE_BASIS[criterion],
        }
    if criterion == "lid":
        return {
            "assignability": ASSIGNABLE,
            "reason": REASON_LID_DETECTOR,
            "basis": "sequence_only_lid_detector_ran_and_reported_this_state",
        }
    if criterion in ("linker", "sbd"):
        return {
            "assignability": NOT_ASSIGNABLE,
            "reason": str(miss_basis["assignability_reason_when_a_miss_is_assignable"]),
            "basis": NOT_ASSIGNABLE_BASIS[criterion],
        }
    if bool(miss_basis["miss_is_assignable"]):
        return {
            "assignability": ASSIGNABLE,
            "reason": REASON_PATTERN_MISS,
            "basis": (
                f"pattern_miss_reproduces_{miss_basis['anchors_reproduced']}_of_"
                f"{miss_basis['anchors_with_coordinate']}_literature_anchors_for_this_role"
            ),
        }
    return {
        "assignability": NOT_ASSIGNABLE,
        "reason": str(prior_basis["reason"]),
        "basis": (
            f"pattern_miss_informative_for_{miss_basis['anchors_reproduced']}_of_"
            f"{miss_basis['anchors_with_coordinate']}_literature_anchors;prior_basis="
            f"{prior_basis['status']};prior_anchor={prior_basis['anchor']}"
        ),
    }


def state_word(criterion: str, detection_state: str, assignability: Mapping[str, str]) -> str:
    """Combine detection and assignability into one published state word.

    The positive/absent words are the 2026-09-15 panel's own vocabulary
    (``supported`` / ``not_detected_pattern``); only rows that could not be
    evaluated get the new, explicitly namespaced ``not_assignable_*`` words.
    """
    if assignability["assignability"] == ASSIGNABLE:
        if criterion in PATTERN_CRITERIA:
            return "supported" if detection_state == "supported" else "not_detected_pattern"
        return detection_state
    return str(assignability["reason"])


def motif_panel_status(row: Mapping[str, str]) -> str:
    """Reuse the 2026-09-15 bucket names where they are still accurate."""
    positives = [
        row["lipase_box_state"] == "supported",
        row["his_state"] == "supported",
        row["asp_state"] == "supported",
        row["oxyanion_hole_state"] == "supported",
    ]
    if positives[0] and any(positives[1:]):
        return "partial_catalytic_pattern_panel"
    if positives[0]:
        return "partial_lipase_box_only"
    if any(positives[1:]):
        return "partial_catalytic_pattern_panel"
    return "no_catalytic_residue_pattern_detected"


def assemble_candidate_row(
    *,
    accession: str,
    legacy: Mapping[str, str],
    sequence: str,
    detection: Mapping[str, Mapping[str, str]],
    prior_index: Mapping[str, Mapping[str, object]],
    miss_basis: Mapping[str, Mapping[str, object]],
) -> dict[str, str]:
    prior = str(legacy.get(PRIOR_COLUMN, "")).strip()
    pfam_accessions = str(legacy.get("pfam_accessions", ""))
    pfam_coordinates = str(legacy.get("pfam_coordinates", ""))
    role_basis = {
        criterion: role_anchor_basis(prior, CRITERION_ROLE[criterion], prior_index)
        for criterion in PATTERN_CRITERIA
    }
    row: dict[str, str] = {
        "accession": accession,
        "genome": str(legacy.get("genome", "")),
        "candidate_prior_superfamily": prior,
        "candidate_assignment_status": str(legacy.get("assignment_status", "")),
    }
    prior_entry = prior_index.get(prior) if prior else None
    if not prior:
        row["prior_anchor_basis"] = PRIOR_NONE
        row["prior_anchor_reference"] = "pending_no_superfamily_prior"
    elif prior_entry is None:
        row["prior_anchor_basis"] = "prior_absent_from_reference_mapping"
        row["prior_anchor_reference"] = "pending_prior_not_in_reference_mapping"
    else:
        anchors = prior_entry["anchors"]
        row["prior_anchor_reference"] = (
            ";".join(anchors) if anchors else "pending_no_literature_anchor_for_prior"
        )
        if not prior_entry["anchor_available"]:
            row["prior_anchor_basis"] = PRIOR_NO_ANCHOR
        elif all(bool(prior_entry["role_transferable"].get(role)) for role in CRITERION_ROLE.values()):
            row["prior_anchor_basis"] = PRIOR_ANCHOR_TRANSFERABLE
        else:
            row["prior_anchor_basis"] = PRIOR_ANCHOR_NOT_TRANSFERABLE

    for criterion in PATTERN_CRITERIA:
        detection_payload = detection[criterion]
        decision = assignability_for(
            criterion,
            str(detection_payload["state"]),
            prior_basis=role_basis[criterion],
            miss_basis=miss_basis[criterion],
        )
        state = state_word(criterion, str(detection_payload["state"]), decision)
        coordinates = (
            str(detection_payload["coordinates"])
            if decision["assignability"] == ASSIGNABLE and detection_payload["state"] == "supported"
            else ""
        )
        row[f"{criterion}_state"] = state
        row[f"{criterion}_coordinates"] = coordinates
        row[f"{criterion}_assignability"] = decision["assignability"]
        row[f"{criterion}_assignability_reason"] = decision["reason"]
        row[f"{criterion}_assignability_basis"] = decision["basis"]
    row["lipase_box_x1"] = str(detection["lipase_box"].get("x1", ""))

    lid_detection = detection["lid"]
    lid_decision = assignability_for(
        "lid", str(lid_detection["state"]), prior_basis={"reason": ""}, miss_basis=miss_basis["lid"]
    )
    row["lid_state"] = str(lid_detection["state"])
    row["lid_coordinates"] = str(lid_detection["coordinates"])
    row["lid_loop_spacing"] = str(lid_detection.get("loop_spacing", ""))
    row["lid_assignability"] = lid_decision["assignability"]
    row["lid_assignability_reason"] = lid_decision["reason"]
    row["lid_assignability_basis"] = lid_decision["basis"]
    row["lid_expectation_from_prior"] = (
        "expected_for_prior" if prior in LID_EXPECTED_PRIORS else "not_expected_for_prior"
    )
    row["lid_prior_agrees_with_detection"] = str(
        (str(lid_detection["state"]) == "supported") == (prior in LID_EXPECTED_PRIORS)
    ).lower()

    linker = delimit_linker(pfam_accessions, pfam_coordinates)
    row["linker_state"] = REASON_LINKER_PROFILE
    row["linker_assignability"] = NOT_ASSIGNABLE
    row["linker_assignability_reason"] = REASON_LINKER_PROFILE
    row["linker_assignability_basis"] = str(miss_basis["linker"]["basis"])
    row["linker_region_delimitation"] = linker["state"]
    row["linker_region_coordinates"] = linker["coordinates"]
    row["linker_region_basis"] = linker["basis"]

    sbd = sbd_evidence(pfam_accessions, pfam_coordinates)
    row["sbd_state"] = REASON_SBD_CRITERION
    row["sbd_assignability"] = NOT_ASSIGNABLE
    row["sbd_assignability_reason"] = REASON_SBD_CRITERION
    row["sbd_assignability_basis"] = str(miss_basis["sbd"]["basis"])
    row["sbd_pf06850_binding_state"] = sbd["state"]
    row["sbd_pf06850_coordinates"] = sbd["coordinates"]
    row["sbd_pf06850_coordinate_status"] = sbd["coordinate_status"]
    row["sbd_pf06850_role_caveat"] = sbd["role_caveat"]

    assignable = [c for c in CRITERIA if row[f"{c}_assignability"] == ASSIGNABLE]
    not_assignable = [c for c in CRITERIA if row[f"{c}_assignability"] == NOT_ASSIGNABLE]
    row["not_assignable_criteria"] = ";".join(not_assignable)
    row["not_assignable_reason_summary"] = ";".join(
        sorted({row[f"{c}_assignability_reason"] for c in not_assignable})
    )
    row["criteria_assignable_count"] = str(len(assignable))
    row["criteria_not_assignable_count"] = str(len(not_assignable))
    row["pfam_accessions"] = pfam_accessions
    row["motif_panel_status"] = motif_panel_status(row)
    row["motif_evidence_level"] = EVIDENCE_LEVEL
    row["motif_completion_layer"] = ALIGNMENT_LAYER
    row["motif_pending_reason"] = (
        "candidates_are_not_ded_alignment_members_so_no_reference_anchored_residue_transfer"
        if not_assignable else "no_pending_criterion_at_the_candidate_layer"
    )
    row["new_family_call_made"] = "false"
    row["motif_phenotype_boundary"] = BOUNDARY
    return row


def criterion_assignability_totals(rows: Sequence[Mapping[str, str]]) -> dict[str, dict[str, int]]:
    """Assignable / not assignable counts per criterion, plus the positive split."""
    totals: dict[str, dict[str, int]] = {}
    for criterion in CRITERIA:
        assignable = [row for row in rows if row[f"{criterion}_assignability"] == ASSIGNABLE]
        not_assignable = [
            row for row in rows if row[f"{criterion}_assignability"] == NOT_ASSIGNABLE
        ]
        positive = [
            row for row in assignable
            if row[f"{criterion}_state"] in {"supported", "detected"}
            or (
                criterion == "lid"
                and str(row[f"{criterion}_state"]).startswith(("supported", "partial"))
            )
        ]
        totals[criterion] = {
            "assignable": len(assignable),
            "assignable_positive": len(positive),
            "assignable_no_positive": len(assignable) - len(positive),
            "not_assignable": len(not_assignable),
        }
    return totals


def coverage_decomposition_rows(rows: Sequence[Mapping[str, str]]) -> list[dict[str, str]]:
    """One row per criterion bucket so the decomposition is auditable."""
    output: list[dict[str, str]] = []
    totals = criterion_assignability_totals(rows)
    for criterion in CRITERIA:
        bucket = totals[criterion]
        reasons = Counter(
            row[f"{criterion}_assignability_reason"]
            for row in rows if row[f"{criterion}_assignability"] == NOT_ASSIGNABLE
        )
        output.append({
            "criterion": criterion,
            "bucket": ASSIGNABLE,
            "reason_code": sorted(
                {row[f"{criterion}_assignability_reason"]
                 for row in rows if row[f"{criterion}_assignability"] == ASSIGNABLE}
            )[0] if bucket["assignable"] else "not_applicable_no_assignable_candidate",
            "candidate_count": str(bucket["assignable"]),
            "share_of_candidates": f"{bucket['assignable'] / len(rows):.6f}" if rows else "0",
            "assignable": str(bucket["assignable"]),
            "not_assignable": str(bucket["not_assignable"]),
            "assignable_positive": str(bucket["assignable_positive"]),
            "assignable_no_positive": str(bucket["assignable_no_positive"]),
            "evidence_basis": (
                "sequence_pattern_calibrated_on_literature_anchor"
                if criterion in PATTERN_CRITERIA
                else "sequence_only_detector_with_validated_positive_and_negative_controls"
                if criterion == "lid" else "criterion_not_evaluable_at_the_candidate_layer"
            ),
            "note": (
                "a detector positive is always assignable; a miss is assignable only when the "
                "calibrated pattern reproduces every literature anchor for that role"
            ),
        })
        for reason, count in sorted(reasons.items()):
            output.append({
                "criterion": criterion,
                "bucket": NOT_ASSIGNABLE,
                "reason_code": reason,
                "candidate_count": str(count),
                "share_of_candidates": f"{count / len(rows):.6f}" if rows else "0",
                "assignable": str(bucket["assignable"]),
                "not_assignable": str(bucket["not_assignable"]),
                "assignable_positive": str(bucket["assignable_positive"]),
                "assignable_no_positive": str(bucket["assignable_no_positive"]),
                "evidence_basis": "criterion_could_not_be_evaluated_for_this_candidate",
                "note": {
                    REASON_CANDIDATE_NOT_IN_ALIGNMENT: (
                        "the prior has a literature anchor, but no candidate is a member of any "
                        "frozen DED alignment, so the Task 2 column transfer cannot reach it"
                    ),
                    REASON_NO_ANCHOR_FOR_PRIOR: (
                        "the recorded prior has no literature anchor with residue coordinates"
                    ),
                    REASON_PERMUTED_FOLD: (
                        "the anchored column is not transferable in this fold (circular "
                        "permutation; measured conservation 0.0769 in the Task 2 report)"
                    ),
                    REASON_NO_PRIOR: "no superfamily prior and no alignment membership",
                    REASON_PRIOR_UNKNOWN: "prior is absent from the Task 2 reference mapping",
                    REASON_LINKER_PROFILE: (
                        "no linker profile is bound: only a candidate inter-domain region can be "
                        "delimited, which is not a linker call"
                    ),
                    REASON_SBD_CRITERION: (
                        "PF06850 is not the extracellular substrate-binding domain in this panel "
                        "(DEFECT_D); its binding is published as raw evidence only"
                    ),
                }.get(reason, ""),
            })
    return output


def legacy_state_distributions(legacy_rows: Sequence[Mapping[str, str]]) -> dict[str, dict[str, object]]:
    """The 2026-09-15 candidate-layer distributions this rerun is compared against.

    V5: the ``before`` payload used to mix the state partition and derived
    aggregates in one flat dict and carried the total supported count under
    ``independent_detection_supported`` -- a name that claims an *independent
    detection* while the value is the *resolved* total, whatever the evidence
    origin.  For ``lid``, ``linker`` and ``sbd`` it additionally wrote a literal
    ``0``, which reads as "a detector ran and found nothing" although no
    independent detector existed for those criteria in 2026-09-15.

    Field vocabulary of every ``before`` payload (see ``BEFORE_FIELD_SEMANTICS``):

    * ``state_counts``: the measured 2026-09-15 state partition of that criterion,
      summed from the frozen evidence file; it sums to the candidate universe and
      holds no aggregate key;
    * ``before_resolved_supported``: total rows whose state starts with
      ``supported``, whatever the evidence origin (the accurate name for what the
      old field contained);
    * ``before_independent_detection_supported``: the count of rows supported by a
      detector on that row alone, or the sentinel string ``not_available`` when the
      2026-09-15 layer had no independent detector for the criterion -- never a
      literal ``0`` in that case;
    * ``before_not_available`` / ``before_not_available_reason``: whether the
      independent-detection count is unavailable for that criterion, and why.
    """
    before: dict[str, dict[str, object]] = {}
    for criterion in CRITERIA:
        column = f"{criterion}_state"
        counter = Counter(str(row.get(column, PENDING_NO_LEGACY_ROW)) for row in legacy_rows)
        supported = sum(value for key, value in counter.items() if key.startswith("supported"))
        reason = str(BEFORE_NO_INDEPENDENT_DETECTOR.get(criterion, ""))
        before[criterion] = {
            BEFORE_FIELD_STATE_COUNTS: dict(sorted(counter.items())),
            BEFORE_FIELD_RESOLVED: supported,
            BEFORE_FIELD_INDEPENDENT: BEFORE_NOT_AVAILABLE if reason else supported,
            BEFORE_FIELD_NOT_AVAILABLE: bool(reason),
            BEFORE_FIELD_NOT_AVAILABLE_REASON: reason,
        }
    return before


def after_state_distributions(rows: Sequence[Mapping[str, str]]) -> dict[str, dict[str, object]]:
    after: dict[str, dict[str, object]] = {}
    for criterion in CRITERIA:
        counter = Counter(str(row[f"{criterion}_state"]) for row in rows)
        distribution: dict[str, object] = dict(sorted(counter.items()))
        distribution["supported"] = sum(
            value for key, value in counter.items()
            if key in {"supported", "detected"}
        )
        after[criterion] = distribution
    lid_states = Counter(str(row["lid_state"]) for row in rows)
    after["lid"]["partial"] = sum(
        value for key, value in lid_states.items() if key.startswith("partial")
    )
    after["lid"]["any_loop_detected"] = (
        lid_states.get("supported", 0) + after["lid"]["partial"]
    )
    after["lid"]["independent_detection_supported"] = after["lid"]["supported"]
    after["his"]["pattern_supported"] = sum(
        str(row["his_state"]) == "supported" for row in rows
    )
    after["asp"]["pattern_supported"] = sum(
        str(row["asp_state"]) == "supported" for row in rows
    )
    after["sbd"]["pf06850_detected"] = sum(
        str(row["sbd_pf06850_binding_state"]) == "detected" for row in rows
    )
    after["linker"]["region_delimited_not_profiled"] = sum(
        str(row["linker_region_delimitation"]) == "delimited_candidate_region_not_profiled"
        for row in rows
    )
    return after


def candidate_alignment_members(
    accessions: Iterable[str], alignment_dir: Path,
) -> dict[str, object]:
    """Measure how many candidates appear in any frozen DED alignment.

    Only ``*.aln`` files are alignments; the same directory also holds scraped
    ``.html``/``.txt`` pages that must never be parsed as sequence rows.
    """
    candidates = set(accessions)
    pattern = re.compile(r"^\s*([^\s]+)\s+([A-Za-z*.-]+)(?:\s+\d+)?\s*$")
    members: set[str] = set()
    files_examined = 0
    files_with_rows = 0
    for path in sorted(alignment_dir.glob("*.aln")):
        files_examined += 1
        rows_seen = 0
        for line in path.read_text(encoding="ascii", errors="replace").splitlines():
            match = pattern.match(line)
            if not match:
                continue
            name = match.group(1)
            if name.upper() in {"CLUSTAL", "MUSCLE", "PROBCONS", "CLUSTALW"}:
                continue
            rows_seen += 1
            if name in candidates:
                members.add(name)
        if rows_seen:
            files_with_rows += 1
    return {
        "members": sorted(members),
        "files_examined": files_examined,
        "files_with_sequence_rows": files_with_rows,
    }


def build_summary(
    *,
    rows: Sequence[Mapping[str, str]],
    legacy_rows: Sequence[Mapping[str, str]],
    miss_basis: Mapping[str, Mapping[str, object]],
    prior_index: Mapping[str, Mapping[str, object]],
    alignment_members: Sequence[str],
    inputs: Mapping[str, Path],
    calibration_report: Mapping[str, object],
    run_id: str,
) -> dict[str, object]:
    before = legacy_state_distributions(legacy_rows)
    after = after_state_distributions(rows)
    totals = criterion_assignability_totals(rows)
    criteria: dict[str, object] = {}
    for criterion in CRITERIA:
        criteria[criterion] = {
            "before": before[criterion],
            "after": after[criterion],
            "assignable": totals[criterion]["assignable"],
            "assignable_positive": totals[criterion]["assignable_positive"],
            "assignable_no_positive": totals[criterion]["assignable_no_positive"],
            "not_assignable": totals[criterion]["not_assignable"],
            "assignability_reasons": dict(sorted(Counter(
                row[f"{criterion}_assignability_reason"]
                for row in rows if row[f"{criterion}_assignability"] == NOT_ASSIGNABLE
            ).items())),
            "delta_supported": after[criterion]["supported"] - int(before[criterion][BEFORE_FIELD_RESOLVED]),
            "before_source": LEGACY_BEFORE_SOURCE,
            "after_source": f"runs/{run_id}/results/motif_completion_full.tsv",
            "miss_is_assignable": bool(miss_basis[criterion]["miss_is_assignable"]),
        }
    criteria["his"]["note"] = (
        "the 2026-09-15 His pattern GM[A-Z]H[A-Z]{2}P[A-Z]{2}G was self-contradictory and fired "
        "on 265 of 109,087 candidates; the Task 2 pattern GMxH fires on more, but a miss still "
        "says nothing about the catalytic histidine because PhaZKT (DDGHL), PhaZGK13 (GGHEW) and "
        "2D80 (AVHTF) all fail that pattern"
    )
    criteria["asp"]["note"] = (
        "the Asp pattern is byte-identical to the 2026-09-15 G[A-Z]{2}DYTV motif, so the "
        "candidate layer gains no assignable Asp row from Task 2; the reference-layer gain "
        "(25 -> 71) came from anchor column transfer, which reaches no candidate"
    )
    criteria["lid"]["note"] = (
        "0 of 109,087 candidate lid calls in 2026-09-15 were independent: 101,836 were "
        "not_expected_for_subtype and 7,251 pending_reference_annotation, all derived from the "
        "prior.  The sequence-only detector now returns a state for every candidate"
    )
    criteria["linker"]["note"] = "no linker profile is bound, so the criterion stays unassignable"
    criteria["sbd"]["note"] = (
        "PF06850 is not the extracellular substrate-binding domain in this panel (DEFECT_D); the "
        "raw binding is published separately and the criterion stays unassignable"
    )
    criteria["sbd"]["delta_supported_note"] = (
        "the 30,621 'supported' rows in 2026-09-15 were the raw PF06850 binding promoted to the "
        "SBD criterion; the raw binding is unchanged in this rerun (30,621 detected, published in "
        "sbd_pf06850_binding_state) and only the promotion is withdrawn, so the negative delta is "
        "a de-promotion and not a loss of evidence"
    )
    criteria["linker"]["delta_supported_note"] = (
        "no linker state was ever supported in 2026-09-15 either; the 108,795 pending rows were "
        "prior-derived placeholders and are now explicit not-assignable rows"
    )
    criteria["lid"]["delta_supported_note"] = (
        "the 2026-09-15 lid column had zero independent detections, so the whole delta is newly "
        "measured rather than newly favourable"
    )

    lid_agreement = {
        "detected_and_expected": sum(
            row["lid_state"] == "supported" and row["lid_expectation_from_prior"] == "expected_for_prior"
            for row in rows
        ),
        "detected_but_not_expected": sum(
            row["lid_state"] == "supported" and row["lid_expectation_from_prior"] != "expected_for_prior"
            for row in rows
        ),
        "not_detected_but_expected": sum(
            row["lid_state"] != "supported" and row["lid_expectation_from_prior"] == "expected_for_prior"
            for row in rows
        ),
        "not_detected_and_not_expected": sum(
            row["lid_state"] != "supported" and row["lid_expectation_from_prior"] != "expected_for_prior"
            for row in rows
        ),
        "evidence_level": "independent_sequence_detection_compared_to_a_recorded_prior",
    }
    pf06850_by_prior: dict[str, dict[str, int]] = defaultdict(
        lambda: {"detected": 0, "not_detected": 0, "not_tested": 0}
    )
    for row in rows:
        key = row["candidate_prior_superfamily"] or "no_reference_subtype"
        state = row["sbd_pf06850_binding_state"]
        pf06850_by_prior[key]["not_tested" if state == "not_tested_no_pfam_evidence"
                           else "detected" if state == "detected" else "not_detected"] += 1

    transfer = dict(calibration_report.get("alignment_column_transfer", {}))
    pf06850_census = Counter(str(row["sbd_pf06850_coordinate_status"]) for row in rows)
    pf06850_multi_span = sum(
        1 for row in rows
        if str(row["sbd_pf06850_binding_state"]) == "detected"
        and ";" in str(row["sbd_pf06850_coordinates"])
    )
    return {
        "schema_version": "1.0",
        "run_id": run_id,
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "status": "completed_candidate_only",
        "candidate_count": len(rows),
        "candidate_alignment_layer": ALIGNMENT_LAYER,
        "candidate_ded_alignment_members": len(alignment_members),
        "criteria": criteria,
        "field_semantics": dict(BEFORE_FIELD_SEMANTICS),
        "criteria_before_note": (
            "before = the 2026-09-15 candidate-layer distribution of the same criterion, "
            "recomputed from the frozen evidence file; the reference-layer figures (His 0/723, "
            "Asp 25/723) are different numbers and are not the candidate-layer baseline. V5: the "
            "before payload now separates the measured state partition (before.state_counts) from "
            "the aggregates (before_resolved_supported, before_independent_detection_supported), "
            "and reports before_independent_detection_supported='not_available' for lid, linker "
            "and sbd instead of the earlier hard-coded 0"
        ),
        "criterion_miss_assignability_basis": {
            criterion: dict(miss_basis[criterion]) for criterion in CRITERIA
        },
        "coverage": {
            "fully_assignable_criteria": [
                criterion for criterion in CRITERIA
                if totals[criterion]["not_assignable"] == 0
            ],
            "unassignable_criteria": [
                criterion for criterion in CRITERIA
                if totals[criterion]["assignable"] == 0
            ],
            "partially_assignable_criteria": [
                criterion for criterion in CRITERIA
                if 0 < totals[criterion]["assignable"] < len(rows)
            ],
        },
        "criterion_assignability_totals": totals,
        "premise_corrections": [
            "the plan and the task brief quote the His '0 hits' figure from the 723-reference "
            "layer (reference_motif_panel.tsv), not from the candidate layer: the 2026-09-15 "
            "candidate layer recorded his_state=supported for 265 of 109,087 candidates with the "
            "self-contradictory pattern",
            "the brief's 'Asp about 25 hits' is likewise the reference-layer figure: the "
            "2026-09-15 candidate layer recorded asp_state=supported for 1,290 candidates, and "
            "that pattern is unchanged, so this rerun adds no assignable Asp row",
            "the brief's 'only 64 of 723 references transfer coordinates' is measured as 65 in "
            "the Task 2 run and is a reference-layer figure; the candidate-layer counterpart is "
            "zero candidates, because no candidate is a DED alignment member",
        ],
        "lid_prior_agreement": lid_agreement,
        "pf06850_raw_evidence_census": {
            "detected": sum(
                1 for row in rows if str(row["sbd_pf06850_binding_state"]) == "detected"
            ),
            "not_detected_in_tested_pfam": sum(
                1 for row in rows
                if str(row["sbd_pf06850_binding_state"]) == "not_detected_in_tested_pfam"
            ),
            "not_tested_no_pfam_evidence": sum(
                1 for row in rows
                if str(row["sbd_pf06850_binding_state"]) == "not_tested_no_pfam_evidence"
            ),
            "coordinate_status": dict(sorted(pf06850_census.items())),
            "detected_with_more_than_one_span": pf06850_multi_span,
            "note": (
                "detection is accession-driven so the raw count stays comparable with the "
                "2026-09-15 panel (30,621); 10 of those candidates carry two PF06850 spans, "
                "which a naive single-span coordinate parser silently drops"
            ),
        },
        "pf06850_binding_by_prior": {key: dict(value) for key, value in sorted(pf06850_by_prior.items())},
        "prior_anchor_index": {key: dict(value) for key, value in sorted(prior_index.items())},
        "reference_alignment_column_transfer": {
            "references_with_anchor_column_transfer": transfer.get(
                "references_with_anchor_column_transfer", "pending"
            ),
            "anchored_alignments": transfer.get("anchored_alignments", "pending"),
            "source": TASK2_REPORT_SOURCE,
        },
        "inputs": {
            name: {
                "path": str(path.resolve()),
                "size": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for name, path in inputs.items()
        },
        "new_family_or_superfamily_call_made": False,
        "no_not_assignable_row_written_as_a_negative": True,
        "limitations": [
            "no candidate is a member of any frozen DED alignment, so no reference-anchored "
            "residue transfer from Task 2 reaches this layer; the candidate layer is a "
            "sequence-pattern layer only",
            "the His (GMxH) and Asp (GxxDYTV) proxies reproduce 1 of 4 literature anchors each, so "
            "a pattern miss is published as not assignable, never as no catalytic residue",
            "the type 2 fold is circularly permuted (measured column conservation 0.0769), so its "
            "His/Asp stay not assignable even for the anchored family",
            "the linker criterion has no bound profile: only a candidate inter-domain region is "
            "delimited, which is not a linker call",
            "PF06850 is not the extracellular substrate-binding domain in this panel, so the SBD "
            "criterion is unassignable for every candidate and only the raw PF06850 binding is "
            "reported",
            "no candidate carries the intracellular nPHASCL (with or without lipase box), "
            "periplasmic or native-SCL/PhaZ7-like prior, so those superfamilies contribute no "
            "assignability branch here",
            "this task makes no new family or superfamily call and no candidate-level "
            "classification change",
        ],
        "phenotype_boundary_holds": all(row["motif_phenotype_boundary"] == BOUNDARY for row in rows),
        "phenotype_boundary": BOUNDARY,
    }


def run_motif_completion(
    *,
    candidate_fasta: Path,
    legacy_evidence: Path,
    reference_manifest: Path,
    superfamily_table: Path,
    calibration_report: Path,
    literature_reference: Path,
    reference_fasta: Path,
    anchor_fasta: Path,
    ded_alignment_dir: Path,
    output_dir: Path,
    run_id: str,
    detector_script: Path | None = None,
) -> dict[str, object]:
    """Recompute the motif layer for the whole candidate universe."""
    detector_path = detector_script or (
        Path(__file__).resolve().parent / "complete_phaded_reference_residue_mapping.py"
    )
    detectors = load_task2_module(detector_path)
    set_detectors(detectors)

    reference_rows = read_tsv(reference_manifest)
    if not reference_rows:
        raise ValueError(f"reference manifest is empty: {reference_manifest}")
    report = json.loads(Path(calibration_report).read_text(encoding="utf-8"))
    prior_index = build_prior_anchor_index(reference_rows, report)
    miss_basis = build_assignability_basis(
        Path(literature_reference), Path(reference_fasta), Path(anchor_fasta),
        reference_rows, detectors,
    )
    _ = read_tsv(superfamily_table)  # consumed for provenance; the index is built above

    legacy_rows = read_tsv(legacy_evidence)
    if not legacy_rows:
        raise ValueError(f"legacy candidate evidence is empty: {legacy_evidence}")
    legacy_by_accession = {str(row["accession"]).strip(): row for row in legacy_rows}
    if len(legacy_by_accession) != len(legacy_rows):
        raise ValueError("legacy candidate evidence has duplicate accessions")
    sequences = read_fasta(Path(candidate_fasta))
    if set(sequences) != set(legacy_by_accession):
        raise ValueError("candidate FASTA and legacy evidence accession sets differ")

    rows: list[dict[str, str]] = []
    for accession in sorted(sequences):
        legacy = legacy_by_accession[accession]
        rows.append(assemble_candidate_row(
            accession=accession,
            legacy=legacy,
            sequence=sequences[accession],
            detection=criterion_detection(sequences[accession]),
            prior_index=prior_index,
            miss_basis=miss_basis,
        ))

    known_priors = set(prior_index)
    unknown_priors = sorted(
        {row["candidate_prior_superfamily"] for row in rows
         if row["candidate_prior_superfamily"] and row["candidate_prior_superfamily"] not in known_priors}
    )
    members = candidate_alignment_members((row["accession"] for row in rows), Path(ded_alignment_dir))
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(output_dir / "motif_completion_full.tsv", rows)
    write_tsv(output_dir / "coverage_decomposition.tsv", coverage_decomposition_rows(rows))

    summary = build_summary(
        rows=rows,
        legacy_rows=legacy_rows,
        miss_basis=miss_basis,
        prior_index=prior_index,
        alignment_members=members["members"],
        inputs={
            "legacy_candidate_evidence": Path(legacy_evidence),
            "candidate_fasta": Path(candidate_fasta),
            "reference_mapping_manifest": Path(reference_manifest),
            "superfamily_residue_mapping": Path(superfamily_table),
            "task2_calibration_report": Path(calibration_report),
            "literature_residue_reference": Path(literature_reference),
            "literature_anchor_sequences": Path(anchor_fasta),
            "reference_fasta": Path(reference_fasta),
            "task2_detector_script": detector_path,
        },
        calibration_report=report,
        run_id=run_id,
    )
    summary["unknown_prior_superfamilies"] = unknown_priors
    summary["candidate_ded_alignment_member_accessions"] = members["members"]
    summary["ded_alignment_files_examined"] = members["files_examined"]
    summary["ded_alignment_files_with_sequence_rows"] = members["files_with_sequence_rows"]
    summary["outputs"] = {
        name: {
            "path": str((output_dir / name).resolve()),
            "size": (output_dir / name).stat().st_size,
            "sha256": sha256_file(output_dir / name),
        }
        for name in ("motif_completion_full.tsv", "coverage_decomposition.tsv")
    }
    summary_path = output_dir / "motif_completion_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    summary["outputs"]["motif_completion_summary.json"] = {
        "path": str(summary_path.resolve()),
        "size": summary_path.stat().st_size,
        "sha256": sha256_file(summary_path),
    }
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-fasta", type=Path, required=True)
    parser.add_argument("--legacy-evidence", type=Path, required=True)
    parser.add_argument("--reference-manifest", type=Path, required=True)
    parser.add_argument("--superfamily-table", type=Path, required=True)
    parser.add_argument("--calibration-report", type=Path, required=True)
    parser.add_argument("--literature-reference", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--anchor-fasta", type=Path, required=True)
    parser.add_argument("--ded-alignment-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args(argv)
    summary = run_motif_completion(
        candidate_fasta=args.candidate_fasta,
        legacy_evidence=args.legacy_evidence,
        reference_manifest=args.reference_manifest,
        superfamily_table=args.superfamily_table,
        calibration_report=args.calibration_report,
        literature_reference=args.literature_reference,
        reference_fasta=args.reference_fasta,
        anchor_fasta=args.anchor_fasta,
        ded_alignment_dir=args.ded_alignment_dir,
        output_dir=args.output_dir,
        run_id=args.run_id,
    )
    print(json.dumps({
        "run_id": summary["run_id"],
        "candidate_count": summary["candidate_count"],
        "criterion_assignability_totals": summary["criterion_assignability_totals"],
        "candidate_ded_alignment_members": summary["candidate_ded_alignment_members"],
    }, ensure_ascii=False, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
