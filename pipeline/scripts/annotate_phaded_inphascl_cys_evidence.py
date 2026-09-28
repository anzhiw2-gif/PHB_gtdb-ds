#!/usr/bin/env python3
"""Candidate-layer Cys-type annotation evidence for intracellular nPHASCL (no lipase box).

This script publishes a **candidate-layer annotation evidence** table for the
PhaDED superfamily ``intracellular nPHASCL without lipase box`` -- the Cys-His-Asp
family that has no ``Gx1Sx2G`` lipase box. It fits **no** HMM, runs **no**
profile/domain search, opens **no** calibration gate, changes **no**
``subtype_call`` and writes nothing outside its own ``--output-dir``.

Why annotation evidence instead of a fitted profile
--------------------------------------------------
The project rule ``validate_phaded_reference_ledger.py`` (L172-173, enforced by
``test_phaded_reference_ledger``/``test_rejects_annotation_only_record_from_strict_training``)
forbids ``annotation_only``/``pending_review`` ledger records in strict training.
In the frozen 723-row reference ledger the Cys-type superfamily is 275
``annotation_only`` records plus exactly one ``experimental_positive`` record
(``CAJ92291.1``, DED_hfam_65). One trainable sequence cannot support a fitted
profile, and the calibration gate requires >=3 independent positives, so no new
HMM may be trained here. See
``docs/superpowers/plans/2026-09-17-phaded-inphascl-cys-decision-record.md``.

Declared criteria (all measurable on frozen evidence, none trained)
------------------------------------------------------------------
1. ``pf06850_binding`` -- raw Pfam ``PF06850`` binding carried by the frozen
   candidate Pfam architecture evidence. DEFECT_D establishes that in this
   reference set PF06850 marks the intracellular Cys-type superfamily
   (274/276 references) and **not** the extracellular substrate-binding domain
   (0/355 extracellular references). It is therefore used as a *contrast
   calibrated marker*, never as an SBD call.
2. ``lipase_box_absence`` -- ``lipase_box_state == not_detected_pattern``.
   Knoll 2009: "All PHA depolymerases in the DED possess a lipase box ... with
   the exception of the family of intracellular nPHASCL depolymerases (no lipase
   box), which possess a catalytic cysteine instead of the lipase box."
3. ``val_cys_dyad`` -- Knoll 2009 (same paragraph): "all family members of the
   family of intracellular nPHASCL depolymerases (no lipase box) also have a
   hydrophobic residue (almost all valine) at position cysteine-1". The
   descriptor records whether the sequence contains a valine immediately
   followed by a cysteine. It is a supporting descriptor with a measured
   reference-layer contrast, **not** a validated detector: the catalytic cysteine
   coordinate is ``pending_reference_annotation`` for this superfamily.

Evidence boundary: every column below denotes candidate homology or functional
potential only. No column establishes a validated PHB/PHA degradation phenotype,
no column is a family call, and no candidate may be removed on this basis.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
from collections import Counter, OrderedDict
from pathlib import Path


PHENOTYPE_BOUNDARY = (
    "candidate_only; candidate-layer annotation evidence of candidate homology or "
    "functional potential and not a validated PHB/PHA degradation phenotype"
)
CITATION = "Knoll M, et al. BMC Bioinformatics 2009;10:89 (PhaDED; PMC2666664)"

EVIDENCE_LAYER = "candidate_layer_annotation_evidence"
WITHHELD_CLASS = "withheld_falsification_triggered"

PF06850 = "PF06850"
CYS_SUPERFAMILY = "intracellular nPHASCL without lipase box"
EXTRACELLULAR_FALSIFICATION_SUPERFAMILIES = (
    "extracellular dPHASCL type 1",
    "extracellular dPHASCL type 2",
)

# Frozen state vocabularies. Any other token fails closed.
PF06850_BINDING_STATES = (
    "detected",
    "not_detected_in_tested_pfam",
    "not_tested_no_pfam_evidence",
)
LIPASE_BOX_STATES = ("supported", "not_detected_pattern")
VAL_CYS_STATES = ("present", "absent")

CYS_EVIDENCE_CLASSES = tuple(
    f"{pf_state}|{lipase_state}|{val_state}"
    for pf_state in PF06850_BINDING_STATES
    for lipase_state in LIPASE_BOX_STATES
    for val_state in VAL_CYS_STATES
)

# Knoll 2009: the residue at position cysteine-1 is hydrophobic (almost all valine).
HYDROPHOBIC_RESIDUES = frozenset("ACFILMVWY")
VALINE = "V"
CYSTEINE = "C"
GLUTAMINE = "Q"

# Exploratory, *data-derived* descriptor discovered in Task 6 from the frozen reference
# ledger: the tripeptide V-C-Q. It is NOT a literature motif, NOT an HMM and NOT a
# calibrated detector; its measured performance is in-sample only (see
# ``descriptor_discovery_report`` and the ``in_sample`` / ``held_out_set`` flags).
VAL_CYS_GLN_MOTIF = "VCQ"
VAL_CYS_GLN_STATES = ("present", "absent")
EXPLORATORY_ROLE = "exploratory_in_sample_descriptor_discovery"

# Pre-registered falsification thresholds.
PF06850_REFERENCE_SENSITIVITY_FLOOR = 0.90

CANDIDATE_REQUIRED_COLUMNS = (
    "accession",
    "genome",
    "candidate_prior_superfamily",
    "sbd_pf06850_binding_state",
    "sbd_pf06850_coordinates",
    "lipase_box_state",
)

REFERENCE_REQUIRED_COLUMNS = (
    "accession",
    "reference_id",
    "phaded_superfamily",
    "pfam_accessions",
    "lipase_box_state",
)

CANDIDATE_COLUMNS = [
    "accession",
    "genome",
    "candidate_prior_superfamily",
    "sbd_pf06850_binding_state",
    "sbd_pf06850_coordinates",
    "lipase_box_state",
    "val_cys_dyad_state",
    "val_cys_count",
    "val_cys_cysteine_positions",
    "val_cys_first_cysteine_position",
    "val_cys_first_context",
    "hydrophobic_minus1_cys_count",
    "cys_evidence_class",
    "cys_type_candidate_evidence",
    "evidence_layer",
    "family_call_made",
    "hmm_score_present",
    "subtype_call_impact",
    "new_family_call_made",
    "evidence_boundary",
]

REFERENCE_CONTRAST_COLUMNS = [
    "phaded_superfamily",
    "reference_count",
    "pf06850_binding_count",
    "pf06850_binding_rate",
    "lipase_box_absent_count",
    "lipase_box_absent_rate",
    "val_cys_count",
    "val_cys_rate",
    "cys_evidence_classifier_count",
    "cys_evidence_classifier_rate",
    "is_cys_superfamily",
]

DESCRIPTOR_COLUMNS = [
    "accession",
    "candidate_prior_superfamily",
    "val_cys_dyad_state",
    "val_cys_cysteine_positions",
    "val_cys_gln_state",
    "val_cys_gln_count",
    "val_cys_gln_cysteine_positions",
    "val_cys_gln_first_context",
    "evidence_layer",
    "evidence_role",
    "family_call_made",
    "phenotype_boundary",
]

CRITERION_DEFINITIONS = OrderedDict(
    (
        (
            "pf06850_binding",
            {
                "role": "positive marker, contrast calibrated on the frozen reference ledger",
                "observable": "Pfam PF06850 binding carried by the frozen candidate Pfam architecture evidence",
                "literature_basis": "DEFECT_D (Task 2 / Task 7): in this reference set PF06850 marks the intracellular Cys-type superfamily (274/276) and is absent from all 355 extracellular references",
                "load_bearing": True,
            },
        ),
        (
            "lipase_box_absence",
            {
                "role": "required absence",
                "observable": "lipase_box_state == not_detected_pattern under the frozen Gx1Sx2G pattern whose misses are assignable (4/4 literature anchors reproduced)",
                "literature_basis": "Knoll 2009: the intracellular nPHASCL (no lipase box) family possesses a catalytic cysteine instead of the lipase box",
                "load_bearing": True,
            },
        ),
        (
            "val_cys_dyad",
            {
                "role": "supporting descriptor, contrast calibrated on the frozen reference ledger",
                "observable": "the sequence contains a valine immediately followed by a cysteine",
                "literature_basis": "Knoll 2009: all family members have a hydrophobic residue (almost all valine) at position cysteine-1",
                "load_bearing": False,
            },
        ),
    )
)

FALSIFICATION_GATES = (
    "pf06850_reference_specificity",
    "pf06850_reference_sensitivity",
    "extracellular_reference_no_false_positive",
    "lipase_box_miss_assignable",
)


def _require_state(value: object, allowed: tuple[str, ...], what: str) -> str:
    token = "" if value is None else str(value).strip()
    if token not in allowed:
        raise ValueError(f"unknown {what}: {value!r} (allowed: {', '.join(allowed)})")
    return token


def cys_evidence_class(
    pf06850_binding_state: object, lipase_box_state: object, val_cys_state: object
) -> str:
    """Return the declared class label for the three measured primitive states."""
    pf_state = _require_state(pf06850_binding_state, PF06850_BINDING_STATES, "pf06850 binding state")
    lipase_state = _require_state(lipase_box_state, LIPASE_BOX_STATES, "lipase box state")
    val_state = _require_state(val_cys_state, VAL_CYS_STATES, "val-cys dyad state")
    return f"{pf_state}|{lipase_state}|{val_state}"


def is_cys_type_candidate_evidence(pf06850_binding_state: object, lipase_box_state: object) -> bool:
    """Return whether the two load-bearing criteria are both satisfied."""
    pf_state = _require_state(pf06850_binding_state, PF06850_BINDING_STATES, "pf06850 binding state")
    lipase_state = _require_state(lipase_box_state, LIPASE_BOX_STATES, "lipase box state")
    return pf_state == "detected" and lipase_state == "not_detected_pattern"


def _clean_sequence(sequence: object) -> str:
    if sequence is None:
        raise ValueError("sequence is required")
    text = str(sequence).strip().upper()
    if not text:
        raise ValueError("sequence is empty")
    return text


def val_cys_cysteine_positions(sequence: object) -> list[int]:
    """Return 1-based positions of cysteines preceded by a valine (Knoll 2009 cysteine-1)."""
    text = _clean_sequence(sequence)
    positions: list[int] = []
    for offset in range(1, len(text)):
        if text[offset] == CYSTEINE and text[offset - 1] == VALINE:
            positions.append(offset + 1)
    return positions


def hydrophobic_minus1_cysteine_positions(sequence: object) -> list[int]:
    """Return 1-based positions of cysteines preceded by any hydrophobic residue."""
    text = _clean_sequence(sequence)
    positions: list[int] = []
    for offset in range(1, len(text)):
        if text[offset] == CYSTEINE and text[offset - 1] in HYDROPHOBIC_RESIDUES:
            positions.append(offset + 1)
    return positions


def val_cys_gln_cysteine_positions(sequence: object) -> list[int]:
    """Return 1-based positions of the cysteine in every V-C-Q tripeptide.

    Exploratory descriptor only: Knoll 2009 states the cysteine-1 rule but says
    nothing about cysteine+1, so the glutamine is a *measured* coincidence in the
    frozen reference ledger, not a literature motif.
    """
    text = _clean_sequence(sequence)
    positions: list[int] = []
    for offset in range(0, len(text) - 2):
        if text[offset : offset + 3] == VAL_CYS_GLN_MOTIF:
            positions.append(offset + 2)
    return positions


def _context(sequence: str, position: int, flank: int = 4) -> str:
    """Return the sequence window around a 1-based position, unpadded."""
    start = max(0, position - 1 - flank)
    end = min(len(sequence), position - 1 + flank + 1)
    return sequence[start:end]


def _read_tsv(path: Path, required: tuple[str, ...], what: str) -> list[dict[str, str]]:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    rows: list[dict[str, str]] = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        missing = [name for name in required if name not in fields]
        if missing:
            raise ValueError(f"{what} missing required columns: {', '.join(missing)}")
        for row in reader:
            rows.append({name: (row.get(name) or "").strip() for name in required})
    if not rows:
        raise ValueError(f"{what} has no data rows: {path}")
    return rows


def read_candidate_evidence(path: Path) -> list[dict[str, str]]:
    """Read and validate the frozen candidate motif-completion table (read-only)."""
    rows = _read_tsv(Path(path), CANDIDATE_REQUIRED_COLUMNS, "candidate evidence")
    seen: set[str] = set()
    for index, row in enumerate(rows, start=2):
        accession = row["accession"]
        if not accession:
            raise ValueError(f"empty accession at line {index}")
        if accession in seen:
            raise ValueError(f"duplicate accession: {accession}")
        seen.add(accession)
        try:
            _require_state(row["sbd_pf06850_binding_state"], PF06850_BINDING_STATES, "pf06850 binding state")
            _require_state(row["lipase_box_state"], LIPASE_BOX_STATES, "lipase box state")
        except ValueError as error:
            raise ValueError(f"{error} (accession {accession}, line {index})") from error
    return rows


def read_reference_manifest(path: Path) -> list[dict[str, str]]:
    """Read and validate the frozen reference residue-mapping manifest (read-only)."""
    rows = _read_tsv(Path(path), REFERENCE_REQUIRED_COLUMNS, "reference manifest")
    seen: set[str] = set()
    for index, row in enumerate(rows, start=2):
        reference_id = row["reference_id"]
        if not reference_id:
            raise ValueError(f"empty reference_id at line {index}")
        if reference_id in seen:
            raise ValueError(f"duplicate reference_id: {reference_id}")
        seen.add(reference_id)
        try:
            _require_state(row["lipase_box_state"], LIPASE_BOX_STATES, "lipase box state")
        except ValueError as error:
            raise ValueError(f"{error} (reference {reference_id}, line {index})") from error
    return rows


def read_fasta_sequences(path: Path) -> dict[str, str]:
    """Read a FASTA file into ``{first_header_token: sequence}`` (read-only)."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    sequences: dict[str, list[str]] = {}
    name: str | None = None
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                if not name:
                    raise ValueError(f"FASTA header without a name in {path}")
                if name in sequences:
                    raise ValueError(f"duplicate FASTA header: {name}")
                sequences[name] = []
            elif name is None:
                raise ValueError(f"FASTA sequence before the first header in {path}")
            else:
                sequences[name].append(line)
    if not sequences:
        raise ValueError(f"FASTA file has no records: {path}")
    return {key: "".join(value).upper() for key, value in sequences.items()}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _fraction(numerator: int, denominator: int) -> str:
    if denominator <= 0:
        return ""
    return f"{numerator / denominator:.6f}"


def _write_tsv(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    if path.exists():
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def _as_bool(value: object) -> bool | None:
    if isinstance(value, bool):
        return value
    if value is None:
        return None
    token = str(value).strip().lower()
    if token in {"true", "1", "yes"}:
        return True
    if token in {"false", "0", "no"}:
        return False
    return None


def _reference_pf06850(row: dict[str, str]) -> bool:
    accessions = {token.strip() for token in (row.get("pfam_accessions") or "").split(";") if token.strip()}
    return PF06850 in accessions


def _index_reference_sequences(sequences: dict[str, str]) -> dict[str, str]:
    """Index reference FASTA records by ``reference_id``.

    The frozen reference FASTA uses ``<reference_id>|<accession>`` headers, so the
    ``reference_id`` used by the reference manifest is the part before the pipe.
    """
    indexed: dict[str, str] = {}
    for header, sequence in sequences.items():
        reference_id = header.split("|", 1)[0]
        if reference_id in indexed:
            raise ValueError(f"duplicate reference_id in the reference FASTA: {reference_id}")
        indexed[reference_id] = sequence
    return indexed


def build_reference_contrast(
    reference_rows: list[dict[str, str]], reference_sequences: dict[str, str]
) -> list[dict[str, object]]:
    """Measure every declared criterion on the frozen reference ledger, per superfamily."""
    reference_sequences = _index_reference_sequences(reference_sequences)
    grouped: dict[str, list[dict[str, str]]] = OrderedDict()
    for row in reference_rows:
        grouped.setdefault(row["phaded_superfamily"], []).append(row)
    contrast: list[dict[str, object]] = []
    for superfamily, rows in grouped.items():
        total = len(rows)
        pf_count = 0
        absent_count = 0
        val_count = 0
        classifier_count = 0
        for row in rows:
            sequence = reference_sequences.get(row["reference_id"])
            if sequence is None:
                raise ValueError(
                    f"reference {row['reference_id']} has no sequence in the reference FASTA"
                )
            pf_binding = _reference_pf06850(row)
            lipase_absent = row["lipase_box_state"] == "not_detected_pattern"
            val_present = bool(val_cys_cysteine_positions(sequence))
            pf_count += pf_binding
            absent_count += lipase_absent
            val_count += val_present
            classifier_count += pf_binding and lipase_absent and val_present
        contrast.append(
            {
                "phaded_superfamily": superfamily,
                "reference_count": total,
                "pf06850_binding_count": pf_count,
                "pf06850_binding_rate": _fraction(pf_count, total),
                "lipase_box_absent_count": absent_count,
                "lipase_box_absent_rate": _fraction(absent_count, total),
                "val_cys_count": val_count,
                "val_cys_rate": _fraction(val_count, total),
                "cys_evidence_classifier_count": classifier_count,
                "cys_evidence_classifier_rate": _fraction(classifier_count, total),
                "is_cys_superfamily": "true" if superfamily == CYS_SUPERFAMILY else "false",
            }
        )
    return sorted(contrast, key=lambda row: str(row["phaded_superfamily"]))


def evaluate_falsification_gates(
    contrast: list[dict[str, object]], motif_completion_summary: dict[str, object] | None
) -> dict[str, dict[str, object]]:
    """Evaluate the pre-registered gates against the frozen reference ledger."""
    by_superfamily = {str(row["phaded_superfamily"]): row for row in contrast}
    cys_row = by_superfamily.get(CYS_SUPERFAMILY)
    if cys_row is None:
        raise ValueError(
            f"reference ledger does not contain the {CYS_SUPERFAMILY} superfamily"
        )
    observed_sensitivity = float(str(cys_row["pf06850_binding_rate"] or "0"))

    outside = sum(
        int(row["pf06850_binding_count"])
        for name, row in by_superfamily.items()
        if name != CYS_SUPERFAMILY
    )
    specificity = {
        "gate": "Pfam PF06850 must not bind any reference outside the Cys-type superfamily (DEFECT_D must hold)",
        "observed_outside_hits": outside,
        "expected": 0,
        "outcome": "passed" if outside == 0 else "failed",
    }

    sensitivity = {
        "gate": "Pfam PF06850 must bind at least the declared floor of the Cys-type reference superfamily",
        "observed_rate": observed_sensitivity,
        "floor": PF06850_REFERENCE_SENSITIVITY_FLOOR,
        "cys_reference_count": int(cys_row["reference_count"]),
        "cys_pf06850_binding_count": int(cys_row["pf06850_binding_count"]),
        "outcome": "passed" if observed_sensitivity >= PF06850_REFERENCE_SENSITIVITY_FLOOR else "failed",
    }

    extracellular_hits = sum(
        int(by_superfamily[name]["cys_evidence_classifier_count"])
        for name in EXTRACELLULAR_FALSIFICATION_SUPERFAMILIES
        if name in by_superfamily
    )
    extracellular = {
        "gate": (
            "No known extracellular reference (284 dPHASCL type 1 + 71 type 2) may be classified as a "
            "Cys-type candidate"
        ),
        "superfamilies": list(EXTRACELLULAR_FALSIFICATION_SUPERFAMILIES),
        "observed_hits": extracellular_hits,
        "expected": 0,
        "outcome": "passed" if extracellular_hits == 0 else "failed",
        "note": (
            "implied by the specificity gate in the reference layer and therefore not vacuous there; it is "
            "NOT checkable in the candidate layer because no candidate carries an independent family label"
        ),
    }

    if motif_completion_summary is None:
        lipase_assignable = {
            "gate": "the lipase-box criterion miss must be assignable on the frozen evidence",
            "outcome": "not_evaluated_input_pending",
            "note": "no frozen motif-completion summary was supplied to this run",
        }
    else:
        basis = (
            (motif_completion_summary.get("criterion_miss_assignability_basis") or {}).get("lipase_box") or {}
        )
        criterion = ((motif_completion_summary.get("criteria") or {}).get("lipase_box") or {})
        miss_assignable = _as_bool(criterion.get("miss_is_assignable"))
        anchors_with_coordinate = basis.get("anchors_with_coordinate")
        anchors_reproduced = basis.get("anchors_reproduced")
        anchors_ok = (
            anchors_with_coordinate is not None
            and anchors_reproduced is not None
            and int(anchors_with_coordinate) == int(anchors_reproduced)
            and int(anchors_with_coordinate) > 0
        )
        lipase_assignable = {
            "gate": "the lipase-box criterion miss must be assignable on the frozen evidence",
            "miss_is_assignable": miss_assignable,
            "anchors_with_coordinate": anchors_with_coordinate,
            "anchors_reproduced": anchors_reproduced,
            "source": "frozen motif-completion summary (Task 7)",
            "outcome": "passed" if (miss_assignable is True and anchors_ok) else "failed",
        }

    return {
        "pf06850_reference_specificity": specificity,
        "pf06850_reference_sensitivity": sensitivity,
        "extracellular_reference_no_false_positive": extracellular,
        "lipase_box_miss_assignable": lipase_assignable,
    }


def annotate(
    *,
    candidate_evidence: Path,
    reference_manifest: Path,
    candidate_fasta: Path,
    reference_fasta: Path,
    output_dir: Path,
    run_id: str,
    motif_completion_summary: Path | None = None,
) -> dict[str, object]:
    """Write the candidate-layer Cys-type annotation tables and return a manifest."""
    candidate_evidence = Path(candidate_evidence)
    reference_manifest = Path(reference_manifest)
    candidate_fasta = Path(candidate_fasta)
    reference_fasta = Path(reference_fasta)
    output_dir = Path(output_dir)

    evidence_rows = read_candidate_evidence(candidate_evidence)
    reference_rows = read_reference_manifest(reference_manifest)
    candidate_sequences = read_fasta_sequences(candidate_fasta)
    reference_sequences = read_fasta_sequences(reference_fasta)

    contrast = build_reference_contrast(reference_rows, reference_sequences)
    summary_payload: dict[str, object] | None = None
    if motif_completion_summary is not None:
        motif_completion_summary = Path(motif_completion_summary)
        if not motif_completion_summary.is_file():
            raise FileNotFoundError(motif_completion_summary)
        summary_payload = json.loads(motif_completion_summary.read_text(encoding="utf-8-sig"))
    gates = evaluate_falsification_gates(contrast, summary_payload)
    triggered = sorted(name for name, gate in gates.items() if gate["outcome"] == "failed")

    annotated: list[dict[str, object]] = []
    missing_sequences: list[str] = []
    for row in evidence_rows:
        accession = row["accession"]
        sequence = candidate_sequences.get(accession)
        if sequence is None:
            missing_sequences.append(accession)
            continue
        val_positions = val_cys_cysteine_positions(sequence)
        hydrophobic_positions = hydrophobic_minus1_cysteine_positions(sequence)
        val_state = "present" if val_positions else "absent"
        label = cys_evidence_class(
            row["sbd_pf06850_binding_state"], row["lipase_box_state"], val_state
        )
        positive = is_cys_type_candidate_evidence(
            row["sbd_pf06850_binding_state"], row["lipase_box_state"]
        )
        if triggered:
            label = WITHHELD_CLASS
            positive = False
        annotated.append(
            {
                "accession": accession,
                "genome": row["genome"],
                "candidate_prior_superfamily": row["candidate_prior_superfamily"],
                "sbd_pf06850_binding_state": row["sbd_pf06850_binding_state"],
                "sbd_pf06850_coordinates": row["sbd_pf06850_coordinates"],
                "lipase_box_state": row["lipase_box_state"],
                "val_cys_dyad_state": val_state,
                "val_cys_count": len(val_positions),
                "val_cys_cysteine_positions": ";".join(str(p) for p in val_positions),
                "val_cys_first_cysteine_position": val_positions[0] if val_positions else "",
                "val_cys_first_context": _context(sequence, val_positions[0]) if val_positions else "",
                "hydrophobic_minus1_cys_count": len(hydrophobic_positions),
                "cys_evidence_class": label,
                "cys_type_candidate_evidence": "true" if positive else "false",
                "evidence_layer": EVIDENCE_LAYER,
                "family_call_made": "false",
                "hmm_score_present": "false",
                "subtype_call_impact": "none",
                "new_family_call_made": "false",
                "evidence_boundary": PHENOTYPE_BOUNDARY,
            }
        )
    if missing_sequences:
        raise ValueError(
            "candidate evidence rows without a sequence in the candidate FASTA "
            f"({len(missing_sequences)}): {', '.join(sorted(missing_sequences)[:10])}"
        )

    annotated.sort(key=lambda row: str(row["accession"]))
    class_counts = Counter(str(row["cys_evidence_class"]) for row in annotated)
    pf_counts = Counter(str(row["sbd_pf06850_binding_state"]) for row in annotated)
    lipase_counts = Counter(str(row["lipase_box_state"]) for row in annotated)
    prior_counts = Counter(
        str(row["candidate_prior_superfamily"]) or "(no_superfamily_prior)" for row in annotated
    )
    positive_count = sum(1 for row in annotated if row["cys_type_candidate_evidence"] == "true")

    output_dir.mkdir(parents=True, exist_ok=True)
    candidate_table = output_dir / "cys_candidate_evidence.tsv"
    contrast_table = output_dir / "cys_reference_contrast.tsv"
    _write_tsv(candidate_table, CANDIDATE_COLUMNS, annotated)
    _write_tsv(contrast_table, REFERENCE_CONTRAST_COLUMNS, contrast)

    status = "falsification_triggered" if triggered else "completed_candidate_only"
    manifest: dict[str, object] = {
        "schema_version": "1.0",
        "run_id": run_id,
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "status": status,
        "analysis": "inphascl_cys_candidate_annotation_evidence",
        "candidate_count": len(annotated),
        "cys_type_candidate_evidence_count": positive_count,
        "cys_evidence_class_counts": dict(sorted(class_counts.items())),
        "pf06850_binding_state_counts": dict(sorted(pf_counts.items())),
        "lipase_box_state_counts": dict(sorted(lipase_counts.items())),
        "candidate_prior_superfamily_counts": dict(sorted(prior_counts.items())),
        "subtype_call_rows_changed": 0,
        "hmm_fitted": False,
        "hmm_score_present": False,
        "new_family_call_made": False,
        "evidence_layer": EVIDENCE_LAYER,
        "criterion_definitions": CRITERION_DEFINITIONS,
        "falsification_gates": gates,
        "falsification_gates_triggered": triggered,
        "evidence_boundary": PHENOTYPE_BOUNDARY,
        "citation": CITATION,
        "candidate_layer_note": (
            "this table adds evidence columns only; it changes no subtype_call, no family assignment, no "
            "registry entry and no formal scan input, and it contains no HMM score"
        ),
        "gate_note": (
            "the calibration gate (>=3 independent positives + >=1 held-out + >=1 family-resolved negative "
            "+ >=1 challenge + zero unexplained hits) is unchanged and remains unmet for this superfamily"
        ),
        "inputs": {
            "candidate_evidence": {
                "path": str(candidate_evidence.resolve()),
                "size": candidate_evidence.stat().st_size,
                "sha256": _sha256(candidate_evidence),
            },
            "reference_manifest": {
                "path": str(reference_manifest.resolve()),
                "size": reference_manifest.stat().st_size,
                "sha256": _sha256(reference_manifest),
            },
            "candidate_fasta": {
                "path": str(candidate_fasta.resolve()),
                "size": candidate_fasta.stat().st_size,
                "sha256": _sha256(candidate_fasta),
            },
            "reference_fasta": {
                "path": str(reference_fasta.resolve()),
                "size": reference_fasta.stat().st_size,
                "sha256": _sha256(reference_fasta),
            },
        },
        "outputs": {},
    }
    if motif_completion_summary is not None:
        manifest["inputs"]["motif_completion_summary"] = {  # type: ignore[index]
            "path": str(motif_completion_summary.resolve()),
            "size": motif_completion_summary.stat().st_size,
            "sha256": _sha256(motif_completion_summary),
        }

    for path in (candidate_table, contrast_table):
        manifest["outputs"][path.name] = {  # type: ignore[index]
            "path": str(path.resolve()),
            "size": path.stat().st_size,
            "sha256": _sha256(path),
        }

    manifest_path = output_dir / "cys_decision_manifest.json"
    if manifest_path.exists():
        raise FileExistsError(manifest_path)
    with manifest_path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def descriptor_discovery_report(
    reference_rows: list[dict[str, str]], reference_sequences: dict[str, str]
) -> dict[str, object]:
    """Measure the two cysteine-context descriptors on the frozen reference ledger.

    Both descriptors are measured *in sample*: the reference ledger is the only
    labelled data in this workspace and there is no held-out set, so these numbers
    are an in-sample contrast and never a validation.
    """
    reference_sequences = _index_reference_sequences(reference_sequences)
    descriptors = {
        "val_cys_dyad": val_cys_cysteine_positions,
        "val_cys_gln_tripeptide": val_cys_gln_cysteine_positions,
    }
    report: dict[str, object] = {}
    for name, descriptor in descriptors.items():
        buckets: dict[str, dict[str, object]] = {
            "cys_superfamily": {"reference_count": 0, "hits": 0, "hit_accessions": []},
            "other_superfamilies": {"reference_count": 0, "hits": 0, "hit_accessions": []},
        }
        for row in reference_rows:
            sequence = reference_sequences.get(row["reference_id"])
            if sequence is None:
                raise ValueError(
                    f"reference {row['reference_id']} has no sequence in the reference FASTA"
                )
            bucket = (
                "cys_superfamily"
                if row["phaded_superfamily"] == CYS_SUPERFAMILY
                else "other_superfamilies"
            )
            buckets[bucket]["reference_count"] = int(buckets[bucket]["reference_count"]) + 1
            if descriptor(sequence):
                buckets[bucket]["hits"] = int(buckets[bucket]["hits"]) + 1
                buckets[bucket]["hit_accessions"].append(row["accession"])  # type: ignore[union-attr]
        for bucket in buckets.values():
            total = int(bucket["reference_count"])
            bucket["rate"] = _fraction(int(bucket["hits"]), total)
            bucket["hit_accessions"] = sorted(bucket["hit_accessions"])  # type: ignore[arg-type]
        report[name] = {
            **buckets,
            "literature_basis": (
                "Knoll 2009: a hydrophobic residue (almost all valine) at position cysteine-1"
                if name == "val_cys_dyad"
                else "no literature basis: data-derived from the frozen reference ledger in Task 6"
            ),
            "load_bearing": False,
        }
    return {
        "in_sample": True,
        "held_out_set": "none",
        "note": (
            "no held-out cysteine-type data exists in this workspace; the frozen 723-row reference "
            "ledger is both the discovery set and the only measurement set, so these rates are an "
            "in-sample contrast and must not be reported as validation"
        ),
        "descriptors": report,
    }


def annotate_typing_descriptor(
    *,
    candidate_evidence: Path,
    reference_manifest: Path,
    candidate_fasta: Path,
    reference_fasta: Path,
    output_dir: Path,
    run_id: str,
    _allow_without_preregistered_artifact: bool = False,
) -> dict[str, object]:
    """Publish the exploratory cysteine-context descriptor pass.

    This pass writes only new file names. It refuses to run unless the
    pre-registered ``cys_candidate_evidence.tsv`` artifact is already present, and
    it records that artifact's SHA-256 so that a reviewer can confirm the
    pre-registered result was not adjusted after the descriptor was measured.
    """
    candidate_evidence = Path(candidate_evidence)
    reference_manifest = Path(reference_manifest)
    candidate_fasta = Path(candidate_fasta)
    reference_fasta = Path(reference_fasta)
    output_dir = Path(output_dir)

    preregistered = output_dir / "cys_candidate_evidence.tsv"
    checked = preregistered.is_file()
    if not checked and not _allow_without_preregistered_artifact:
        raise FileNotFoundError(
            f"the pre-registered artifact must exist before the exploratory pass: {preregistered}"
        )

    evidence_rows = read_candidate_evidence(candidate_evidence)
    reference_rows = read_reference_manifest(reference_manifest)
    candidate_sequences = read_fasta_sequences(candidate_fasta)
    reference_sequences = read_fasta_sequences(reference_fasta)

    rows: list[dict[str, object]] = []
    missing_sequences: list[str] = []
    for row in evidence_rows:
        accession = row["accession"]
        sequence = candidate_sequences.get(accession)
        if sequence is None:
            missing_sequences.append(accession)
            continue
        dyad = val_cys_cysteine_positions(sequence)
        gln = val_cys_gln_cysteine_positions(sequence)
        rows.append(
            {
                "accession": accession,
                "candidate_prior_superfamily": row["candidate_prior_superfamily"],
                "val_cys_dyad_state": "present" if dyad else "absent",
                "val_cys_cysteine_positions": ";".join(str(p) for p in dyad),
                "val_cys_gln_state": "present" if gln else "absent",
                "val_cys_gln_count": len(gln),
                "val_cys_gln_cysteine_positions": ";".join(str(p) for p in gln),
                "val_cys_gln_first_context": _context(sequence, gln[0]) if gln else "",
                "evidence_layer": EVIDENCE_LAYER,
                "evidence_role": EXPLORATORY_ROLE,
                "family_call_made": "false",
                "phenotype_boundary": PHENOTYPE_BOUNDARY,
            }
        )
    if missing_sequences:
        raise ValueError(
            "candidate evidence rows without a sequence in the candidate FASTA "
            f"({len(missing_sequences)}): {', '.join(sorted(missing_sequences)[:10])}"
        )
    rows.sort(key=lambda row: str(row["accession"]))

    discovery = descriptor_discovery_report(reference_rows, reference_sequences)
    discovery.update(
        {
            "schema_version": "1.0",
            "run_id": run_id,
            "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "status": "completed_candidate_only",
            "analysis": "inphascl_cys_typing_descriptor_discovery",
            "evidence_role": EXPLORATORY_ROLE,
            "preregistered_artifact_checked": checked,
            "preregistered_artifact_sha256": _sha256(preregistered) if checked else None,
            "evidence_boundary": PHENOTYPE_BOUNDARY,
            "promotion_rule": (
                "no descriptor in this report may become a filter for a family call until it has been "
                "measured on a held-out, independent non-Cys negative set"
            ),
        }
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    table = output_dir / "cys_typing_descriptor_evidence.tsv"
    discovery_path = output_dir / "cys_typing_descriptor_discovery.json"
    _write_tsv(table, DESCRIPTOR_COLUMNS, rows)
    if discovery_path.exists():
        raise FileExistsError(discovery_path)
    with discovery_path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(discovery, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")

    gln_bucket = discovery["descriptors"]["val_cys_gln_tripeptide"]  # type: ignore[index]
    manifest = {
        "schema_version": "1.0",
        "run_id": run_id,
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "status": "completed_candidate_only",
        "analysis": "inphascl_cys_typing_descriptor_pass",
        "evidence_role": EXPLORATORY_ROLE,
        "candidate_count": len(rows),
        "val_cys_dyad_present_count": sum(
            1 for row in rows if row["val_cys_dyad_state"] == "present"
        ),
        "val_cys_gln_present_count": sum(
            1 for row in rows if row["val_cys_gln_state"] == "present"
        ),
        "subtype_call_rows_changed": 0,
        "hmm_fitted": False,
        "family_call_made": False,
        "in_sample": True,
        "held_out_set": "none",
        "preregistered_artifact_checked": checked,
        "preregistered_artifact_sha256": _sha256(preregistered) if checked else None,
        "reference_contrast_in_sample": gln_bucket,
        "evidence_boundary": PHENOTYPE_BOUNDARY,
        "inputs": {
            "candidate_evidence": {
                "path": str(candidate_evidence.resolve()),
                "size": candidate_evidence.stat().st_size,
                "sha256": _sha256(candidate_evidence),
            },
            "reference_manifest": {
                "path": str(reference_manifest.resolve()),
                "size": reference_manifest.stat().st_size,
                "sha256": _sha256(reference_manifest),
            },
            "candidate_fasta": {
                "path": str(candidate_fasta.resolve()),
                "size": candidate_fasta.stat().st_size,
                "sha256": _sha256(candidate_fasta),
            },
            "reference_fasta": {
                "path": str(reference_fasta.resolve()),
                "size": reference_fasta.stat().st_size,
                "sha256": _sha256(reference_fasta),
            },
        },
        "outputs": {
            table.name: {
                "path": str(table.resolve()),
                "size": table.stat().st_size,
                "sha256": _sha256(table),
            },
            discovery_path.name: {
                "path": str(discovery_path.resolve()),
                "size": discovery_path.stat().st_size,
                "sha256": _sha256(discovery_path),
            },
        },
    }
    manifest_path = output_dir / "cys_typing_descriptor_manifest.json"
    if manifest_path.exists():
        raise FileExistsError(manifest_path)
    with manifest_path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--candidate-evidence", required=True)
    parser.add_argument("--reference-manifest", required=True)
    parser.add_argument("--candidate-fasta", required=True)
    parser.add_argument("--reference-fasta", required=True)
    parser.add_argument("--motif-completion-summary", default=None)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args(argv)

    try:
        manifest = annotate(
            candidate_evidence=Path(args.candidate_evidence),
            reference_manifest=Path(args.reference_manifest),
            candidate_fasta=Path(args.candidate_fasta),
            reference_fasta=Path(args.reference_fasta),
            output_dir=Path(args.output_dir),
            run_id=args.run_id,
            motif_completion_summary=Path(args.motif_completion_summary)
            if args.motif_completion_summary
            else None,
        )
    except (ValueError, FileNotFoundError, FileExistsError) as error:
        print(json.dumps({"status": "failed", "error": str(error)}, ensure_ascii=False))
        return 2

    print(
        json.dumps(
            {
                "run_id": manifest["run_id"],
                "status": manifest["status"],
                "candidate_count": manifest["candidate_count"],
                "cys_type_candidate_evidence_count": manifest["cys_type_candidate_evidence_count"],
                "subtype_call_rows_changed": manifest["subtype_call_rows_changed"],
                "falsification_gates_triggered": manifest["falsification_gates_triggered"],
                "outputs": sorted(manifest["outputs"]),  # type: ignore[arg-type]
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
    )
    return 3 if manifest["status"] == "falsification_triggered" else 0


if __name__ == "__main__":
    raise SystemExit(main())
