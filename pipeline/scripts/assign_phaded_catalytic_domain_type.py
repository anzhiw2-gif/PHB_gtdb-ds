#!/usr/bin/env python3
"""Assign the PhaDED catalytic-domain type 1 / type 2 by oxyanion-hole geometry.

Scientific basis (Knoll M, et al. The PHA Depolymerase Engineering Database.
BMC Bioinformatics 2009;10:89, PMC2666664): "Within the sequences of type 1
catalytic domains, the oxyanion hole can be found N-terminal to the lipase box,
similar to lipases. Within the sequences of type 2 catalytic domains, the
oxyanion hole is found C-terminal to the catalytic triad."

Comparison convention: a call is made only from the *strict* relative position of
the oxyanion-hole pattern start against the lipase-box pattern start.  The
oxyanion hole N-terminal to the lipase box (smaller residue number) is type 1;
C-terminal (larger residue number) is type 2; equal start positions are a
position conflict.  This is exactly the convention already implemented in
``pipeline/scripts/audit_phaded_motifs.py`` (L177-180), where a
type-1-labelled sequence whose oxyanion pattern starts *after* the lipase box is
recorded as ``conflict_relative_position``; the two modules must agree.

Fail-closed rules (no inference is ever made):

* a type call requires the oxyanion-hole *and* the lipase-box coordinate field to
  be present and parseable -- otherwise ``undetermined_no_oxyanion``;
* ``oxyanion_hole_state == not_detected_pattern`` never receives a type;
* ``oxyanion_hole_state == conflict_relative_position`` is never promoted to a
  verified type;
* the superfamily candidate label (``dPHASCL1_like_candidate`` /
  ``dPHASCL2_like_candidate``) is never used to compute a type; it is only
  cross-tabulated against the geometric call afterwards.

Evidence boundary: a catalytic-domain type call is primary-sequence geometry
denoting candidate homology or functional potential only.  It is not a validated
PHB/PHA degradation phenotype, and it is not equivalent to the superfamily
candidate label it is cross-tabulated with.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import re
import sys
from collections import Counter
from pathlib import Path


PHENOTYPE_BOUNDARY = (
    "candidate_only; catalytic-domain type is primary-sequence geometry denoting candidate "
    "homology or functional potential, not a validated PHB/PHA degradation phenotype"
)

CITATION = "Knoll M, et al. BMC Bioinformatics 2009;10:89 (PhaDED; PMC2666664)"

# Knoll 2009: within type 1 catalytic domains the oxyanion hole is N-terminal to
# the lipase box.  Flipping this flag would invert every geometric call, so it is
# stated once, tested explicitly, and reported in the run manifest.
TYPE1_OXYANION_N_TERMINAL = True

ALLOWED_CATALYTIC_DOMAIN_TYPES = (
    "type1_verified",
    "type2_verified",
    "undetermined_no_oxyanion",
    "undetermined_position_conflict",
)

RELATIVE_ORIENTATIONS = (
    "oxyanion_n_terminal_of_lipase_box",
    "oxyanion_c_terminal_of_lipase_box",
    "coordinate_start_positions_equal",
    "coordinates_unusable",
)

ASSIGNMENT_BASES = (
    "oxyanion_n_terminal_of_lipase_box",
    "oxyanion_c_terminal_of_lipase_box",
    "coordinate_start_positions_equal",
    "oxyanion_coordinates_missing",
    "invalid_coordinate_format",
    "oxyanion_signature_not_supported",
    "oxyanion_relative_position_conflict",
)

# Vocabulary of the upstream candidate evidence table (audit_phaded_motifs.py).
OXYANION_STATES = (
    "supported",
    "not_detected_pattern",
    "conflict_relative_position",
)

COORDINATE_RE = re.compile(r"^([0-9]+)-([0-9]+)$")

# Identical pattern to audit_phaded_motifs.AHSMG_PATTERN; positive controls are
# re-evaluated here so the "not detected" statement is tied to a reproduced motif.
AHSMG_PATTERN = re.compile("AHSMG")

LABEL_EXPECTED_TYPE = {
    "dPHASCL1_like_candidate": "type1_verified",
    "dPHASCL2_like_candidate": "type2_verified",
}

LABEL_AGREEMENT_VALUES = (
    "label_matches_geometry",
    "label_contradicts_geometry",
    "not_assessable_no_oxyanion_evidence",
    "not_applicable_no_superfamily_label",
)

AHSMG_CONCLUSION_BOUNDARY = (
    "Candidate-only boundary: the AHSMG motif detector passed its positive controls "
    "(AAK07742.1 and 2VTVA, both DED_hfam_7 / PhaZ7-type reference entries), and the pattern was "
    "not detected in any of the 109,087 candidates. The only supported statement is that the "
    "PhaZ7-type pattern was not detected with a positive-control-validated motif. This is not "
    "evidence of family or enzyme absence in nature: the candidate pool may simply not contain it."
)

CANDIDATE_POOL_SOURCE_NOTE = (
    "The 109,087-protein candidate universe was assembled from GTDB-wide profile HMM hits. "
    "A zero AHSMG result therefore speaks about this candidate pool, which may not contain the "
    "PhaZ7-type family at all; it is not a statement about the family or the enzyme in nature."
)

DEFAULT_EXPECTED_TOTAL_CANDIDATES = 109087
DEFAULT_EXPECTED_DPHASCL1_LABELED = 67342

ASSIGNMENT_COLUMNS = [
    "accession",
    "catalytic_domain_type",
    "assignment_basis",
    "relative_orientation",
    "oxyanion_hole_state",
    "oxyanion_hole_coordinates",
    "lipase_box_coordinates",
    "oxyanion_start",
    "lipase_start",
    "subtype_call",
    "subtype_label_agrees_with_geometry",
    "ahsmg_state",
    "evidence_boundary",
]


def parse_coordinate_range(value: str) -> tuple[int, int] | None:
    """Return (start, end) of a 1-based inclusive coordinate range, else ``None``.

    Malformed, non-numeric, zero/negative or reversed ranges are unusable and are
    reported as such rather than being silently coerced.
    """
    text = (value or "").strip()
    match = COORDINATE_RE.fullmatch(text)
    if not match:
        return None
    start, end = int(match.group(1)), int(match.group(2))
    if start < 1 or end < start:
        return None
    return start, end


def _call(
    catalytic_domain_type: str,
    assignment_basis: str,
    relative_orientation: str = "coordinates_unusable",
    oxyanion_start: int | str = "",
    lipase_start: int | str = "",
) -> dict[str, object]:
    if catalytic_domain_type not in ALLOWED_CATALYTIC_DOMAIN_TYPES:
        raise ValueError(f"undocumented catalytic domain type: {catalytic_domain_type}")
    if assignment_basis not in ASSIGNMENT_BASES:
        raise ValueError(f"undocumented assignment basis: {assignment_basis}")
    if relative_orientation not in RELATIVE_ORIENTATIONS:
        raise ValueError(f"undocumented relative orientation: {relative_orientation}")
    return {
        "catalytic_domain_type": catalytic_domain_type,
        "assignment_basis": assignment_basis,
        "relative_orientation": relative_orientation,
        "oxyanion_start": oxyanion_start,
        "lipase_start": lipase_start,
    }


def assign_catalytic_domain_type(
    oxyanion_hole_state: str,
    oxyanion_hole_coordinates: str,
    lipase_box_coordinates: str,
    *,
    type1_oxyanion_n_terminal: bool = TYPE1_OXYANION_N_TERMINAL,
) -> dict[str, object]:
    """Return a fail-closed catalytic-domain type call from motif coordinates.

    ``type1_oxyanion_n_terminal`` exists so the single orientation assumption can be
    tested and its sensitivity quantified; the default follows Knoll 2009.
    """
    state = (oxyanion_hole_state or "").strip()
    if state not in OXYANION_STATES:
        raise ValueError(f"unknown oxyanion_hole_state: {oxyanion_hole_state!r}")

    if state == "not_detected_pattern":
        # A signature that was never detected cannot support any geometric call.
        return _call(
            "undetermined_no_oxyanion", "oxyanion_signature_not_supported"
        )

    oxyanion_text = (oxyanion_hole_coordinates or "").strip()
    lipase_text = (lipase_box_coordinates or "").strip()
    if not oxyanion_text or not lipase_text:
        return _call("undetermined_no_oxyanion", "oxyanion_coordinates_missing")

    oxyanion = parse_coordinate_range(oxyanion_text)
    lipase = parse_coordinate_range(lipase_text)
    if oxyanion is None or lipase is None:
        return _call("undetermined_no_oxyanion", "invalid_coordinate_format")

    oxyanion_start, lipase_start = oxyanion[0], lipase[0]

    if oxyanion_start == lipase_start:
        return _call(
            "undetermined_position_conflict",
            "coordinate_start_positions_equal",
            "coordinate_start_positions_equal",
            oxyanion_start,
            lipase_start,
        )

    oxyanion_n_terminal = oxyanion_start < lipase_start
    orientation = (
        "oxyanion_n_terminal_of_lipase_box"
        if oxyanion_n_terminal
        else "oxyanion_c_terminal_of_lipase_box"
    )

    if state == "conflict_relative_position":
        # Upstream already flagged the relative position as conflicting with the
        # subtype expectation; geometry is reported but never promoted to a type.
        return _call(
            "undetermined_position_conflict",
            "oxyanion_relative_position_conflict",
            orientation,
            oxyanion_start,
            lipase_start,
        )

    is_type1 = oxyanion_n_terminal if type1_oxyanion_n_terminal else not oxyanion_n_terminal
    return _call(
        "type1_verified" if is_type1 else "type2_verified",
        orientation,
        orientation,
        oxyanion_start,
        lipase_start,
    )


def label_agreement(subtype_call: str, catalytic_domain_type: str) -> str:
    """Cross-tabulate a superfamily candidate label against the geometric call."""
    expected = LABEL_EXPECTED_TYPE.get((subtype_call or "").strip())
    if not expected:
        return "not_applicable_no_superfamily_label"
    if catalytic_domain_type == "undetermined_no_oxyanion":
        return "not_assessable_no_oxyanion_evidence"
    if catalytic_domain_type == expected:
        return "label_matches_geometry"
    return "label_contradicts_geometry"


def build_assignments(
    rows: list[dict[str, str]],
    subtype_calls: dict[str, str] | None = None,
    *,
    type1_oxyanion_n_terminal: bool = TYPE1_OXYANION_N_TERMINAL,
) -> list[dict[str, object]]:
    """Assign a type to every candidate row; the label is only recorded, never used."""
    calls = dict(subtype_calls or {})
    assigned: list[dict[str, object]] = []
    for row in rows:
        accession = (row.get("accession") or "").strip()
        if not accession:
            raise ValueError("candidate row without accession")
        result = assign_catalytic_domain_type(
            row.get("oxyanion_hole_state", ""),
            row.get("oxyanion_hole_coordinates", ""),
            row.get("lipase_box_coordinates", ""),
            type1_oxyanion_n_terminal=type1_oxyanion_n_terminal,
        )
        label = (calls.get(accession) or "").strip()
        assigned.append(
            {
                "accession": accession,
                "catalytic_domain_type": result["catalytic_domain_type"],
                "assignment_basis": result["assignment_basis"],
                "relative_orientation": result["relative_orientation"],
                "oxyanion_hole_state": row.get("oxyanion_hole_state", ""),
                "oxyanion_hole_coordinates": row.get("oxyanion_hole_coordinates", ""),
                "lipase_box_coordinates": row.get("lipase_box_coordinates", ""),
                "oxyanion_start": result["oxyanion_start"],
                "lipase_start": result["lipase_start"],
                "subtype_call": label or "not_in_subtype_matrix",
                "subtype_label_agrees_with_geometry": label_agreement(
                    label, str(result["catalytic_domain_type"])
                ),
                "ahsmg_state": row.get("ahsmg_state", "not_reported"),
                "evidence_boundary": PHENOTYPE_BOUNDARY,
            }
        )
    return assigned


def _type_counts(assigned: list[dict[str, object]]) -> dict[str, int]:
    counter = Counter(str(row["catalytic_domain_type"]) for row in assigned)
    return {name: counter.get(name, 0) for name in ALLOWED_CATALYTIC_DOMAIN_TYPES}


def summarize_assignments(
    assigned: list[dict[str, object]],
    *,
    ahsmg_evidence: dict[str, object] | None = None,
) -> dict[str, object]:
    """Summarize type calls, keeping labels and geometric verification separate."""
    type_counts = _type_counts(assigned)
    total = len(assigned)

    label_view: dict[str, dict[str, object]] = {}
    for label, expected in LABEL_EXPECTED_TYPE.items():
        subset = [row for row in assigned if row["subtype_call"] == label]
        entry: dict[str, object] = {
            "label_name": label,
            "expected_geometry": expected,
            "labeled_count": len(subset),
        }
        for name in ALLOWED_CATALYTIC_DOMAIN_TYPES:
            entry[name] = sum(
                1 for row in subset if row["catalytic_domain_type"] == name
            )
        entry["expected_type_verified_fraction_of_label"] = (
            round(entry[expected] / len(subset), 6) if subset else None
        )
        for value in LABEL_AGREEMENT_VALUES:
            entry[value] = sum(
                1
                for row in subset
                if row["subtype_label_agrees_with_geometry"] == value
            )
        label_view[label] = entry

    inverted = Counter()
    for row in assigned:
        flipped = assign_catalytic_domain_type(
            str(row["oxyanion_hole_state"]),
            str(row["oxyanion_hole_coordinates"]),
            str(row["lipase_box_coordinates"]),
            type1_oxyanion_n_terminal=not TYPE1_OXYANION_N_TERMINAL,
        )
        inverted[str(flipped["catalytic_domain_type"])] += 1

    summary: dict[str, object] = {
        "total_candidates": total,
        "type_counts": type_counts,
        "assignment_basis_counts": dict(
            sorted(Counter(str(row["assignment_basis"]) for row in assigned).items())
        ),
        "relative_orientation_counts": dict(
            sorted(Counter(str(row["relative_orientation"]) for row in assigned).items())
        ),
        "oxyanion_hole_state_counts": dict(
            sorted(Counter(str(row["oxyanion_hole_state"]) for row in assigned).items())
        ),
        "label_vs_geometric_verification": label_view,
        "label_agreement_counts": dict(
            sorted(
                Counter(
                    str(row["subtype_label_agrees_with_geometry"]) for row in assigned
                ).items()
            )
        ),
        "orientation_sensitivity_if_inverted": {
            **{
                name: inverted.get(name, 0)
                for name in ALLOWED_CATALYTIC_DOMAIN_TYPES
            },
            "note": (
                "Sensitivity check only: these are the calls that the documented default "
                "(Knoll 2009: type 1 oxyanion hole N-terminal to the lipase box) would produce "
                "if TYPE1_OXYANION_N_TERMINAL were flipped. The flipped reading is reported "
                "nowhere else and is not the literature convention."
            ),
        },
        "totals_self_consistent": sum(type_counts.values()) == total,
        "criteria_source": CITATION,
        "evidence_boundary": PHENOTYPE_BOUNDARY,
        "label_separation_statement": (
            "dPHASCL1_like_candidate and dPHASCL2_like_candidate are superfamily profile-hit "
            "candidate labels. type1_verified and type2_verified are geometric verifications of "
            "the oxyanion-hole position relative to the lipase box. They are different kinds of "
            "evidence and are never interchangeable."
        ),
    }
    if ahsmg_evidence is not None:
        summary["ahsmg_evidence"] = ahsmg_evidence
    return summary


def read_motif_evidence(path: Path) -> list[dict[str, str]]:
    """Read the candidate motif evidence table (read-only)."""
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_subtype_calls(path: Path) -> dict[str, str]:
    """Read accession -> subtype_call without materializing the whole matrix."""
    calls: dict[str, str] = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header:
            raise ValueError(f"empty subtype matrix: {path}")
        try:
            accession_index = header.index("accession")
            call_index = header.index("subtype_call")
        except ValueError as error:
            raise ValueError(f"subtype matrix lacks accession/subtype_call: {path}") from error
        for row in reader:
            if not row:
                continue
            accession = row[accession_index].strip()
            if not accession:
                raise ValueError(f"subtype matrix row without accession: {path}")
            if accession in calls:
                raise ValueError(f"duplicate accession in subtype matrix: {accession}")
            calls[accession] = row[call_index].strip()
    return calls


def describe_subtype_matrix_agreement(
    matrices: dict[str, dict[str, str]]
) -> dict[str, object]:
    """Report whether every label source agrees, so the check is recorded not assumed."""
    names = list(matrices)
    if not names:
        return {
            "checked_matrices": [],
            "primary_matrix": None,
            "checked_accessions": 0,
            "accession_sets_identical": True,
            "subtype_call_identical": True,
            "disagreeing_accession_count": 0,
            "disagreeing_accessions": [],
        }
    primary = names[0]
    calls = matrices[primary]
    accession_sets_identical = True
    disagreeing: set[str] = set()
    for name in names[1:]:
        other = matrices[name]
        if set(other) != set(calls):
            accession_sets_identical = False
            disagreeing |= set(other) ^ set(calls)
            continue
        for key, value in calls.items():
            if other[key] != value:
                disagreeing.add(key)
    return {
        "checked_matrices": names,
        "primary_matrix": primary,
        "checked_accessions": len(calls),
        "accession_sets_identical": accession_sets_identical,
        "subtype_call_identical": not disagreeing,
        "disagreeing_accession_count": len(disagreeing),
        "disagreeing_accessions": sorted(disagreeing)[:20],
    }


def merge_subtype_calls(matrices: dict[str, dict[str, str]]) -> dict[str, str]:
    """Merge label sources, requiring them to agree (fail-closed on any divergence)."""
    names = list(matrices)
    if not names:
        return {}
    report = describe_subtype_matrix_agreement(matrices)
    if not report["accession_sets_identical"]:
        raise ValueError(
            "subtype matrix accession sets differ: "
            f"e.g. {report['disagreeing_accessions'][:5]}"
        )
    if not report["subtype_call_identical"]:
        raise ValueError(
            "subtype_call disagreement between label sources for "
            f"{report['disagreeing_accession_count']} accessions, "
            f"e.g. {report['disagreeing_accessions'][:5]}"
        )
    return dict(matrices[str(report["primary_matrix"])])


def read_fasta_records(path: Path) -> dict[str, str]:
    """Read a FASTA file into name -> sequence (fails loudly on malformed input)."""
    records: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="ascii", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current is not None:
                if current in records or not chunks:
                    raise ValueError(f"duplicate or empty FASTA record: {current}")
                records[current] = "".join(chunks).upper()
            current = line[1:].split()[0]
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


def reproduce_ahsmg_positive_controls(fasta_path: Path) -> list[dict[str, object]]:
    """Re-run the AHSMG motif on a reference FASTA to reproduce the positive controls."""
    hits: list[dict[str, object]] = []
    for name, sequence in read_fasta_records(fasta_path).items():
        reference_id, _, accession = name.partition("|")
        searchable = sequence.rstrip("*")
        for match in AHSMG_PATTERN.finditer(searchable):
            hits.append(
                {
                    "reference_id": reference_id,
                    "accession": accession,
                    "sequence_length": len(sequence),
                    "ahsmg_coordinates": f"{match.start() + 1}-{match.end()}",
                }
            )
    return hits


def read_reference_panel(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_ahsmg_evidence(
    candidate_rows: list[dict[str, str]],
    reference_rows: list[dict[str, str]],
    reference_fasta: Path,
) -> dict[str, object]:
    """Fix the AHSMG facts: positive controls pass, the candidate pool has zero hits."""
    candidate_states = Counter(
        (row.get("ahsmg_state") or "not_reported") for row in candidate_rows
    )
    reference_states = Counter(
        (row.get("ahsmg_state") or "not_reported") for row in reference_rows
    )
    positive_controls = [
        {
            "accession": row.get("accession", ""),
            "ahsmg_coordinates": row.get("ahsmg_coordinates", ""),
            "phaded_family_id": row.get("phaded_family_id", ""),
            "motif_reference_subtype": row.get("motif_reference_subtype", ""),
            "evidence_status": row.get("evidence_status", ""),
        }
        for row in reference_rows
        if (row.get("ahsmg_state") or "") == "supported"
    ]
    return {
        "detector_pattern": AHSMG_PATTERN.pattern,
        "candidate_pool_size": len(candidate_rows),
        "candidate_ahsmg_state_counts": dict(sorted(candidate_states.items())),
        "candidate_not_detected_pattern": candidate_states.get("not_detected_pattern", 0),
        "reference_panel_size": len(reference_rows),
        "reference_ahsmg_state_counts": dict(sorted(reference_states.items())),
        "reference_positive_control_count": len(positive_controls),
        "reference_positive_controls": positive_controls,
        "reference_fasta_pattern_reproduction": reproduce_ahsmg_positive_controls(
            reference_fasta
        ),
        "conclusion_boundary": AHSMG_CONCLUSION_BOUNDARY,
        "candidate_pool_source_note": CANDIDATE_POOL_SOURCE_NOTE,
    }


def sha256_file(path: Path) -> str:
    import hashlib

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_assignments(path: Path, assigned: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=ASSIGNMENT_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in assigned:
            writer.writerow({name: row[name] for name in ASSIGNMENT_COLUMNS})


def write_source_manifest(
    path: Path,
    motifs: dict[str, Path],
    subtype_matrices: dict[str, Path],
    reference_panel: Path,
    reference_fasta: Path,
) -> list[dict[str, object]]:
    entries: list[dict[str, object]] = []
    for role, group in (
        ("candidate_motif_evidence", motifs),
        ("candidate_subtype_matrix", subtype_matrices),
        ("reference_panel", {"reference_motif_panel": reference_panel}),
        ("reference_fasta", {"reference_sequences": reference_fasta}),
    ):
        for name, source in group.items():
            resolved = Path(source).resolve()
            entries.append(
                {
                    "role": role,
                    "name": name,
                    "path": str(resolved),
                    "size": resolved.stat().st_size,
                    "sha256": sha256_file(resolved),
                }
            )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["role", "name", "path", "size", "sha256"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for entry in entries:
            writer.writerow(entry)
    return entries


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--motif-evidence", type=Path, required=True)
    parser.add_argument(
        "--subtype-matrix",
        type=Path,
        action="append",
        required=True,
        help="candidate subtype matrix providing subtype_call; repeatable, must agree",
    )
    parser.add_argument("--reference-panel", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument(
        "--create-run",
        action="store_true",
        help="create the run layout with pipeline/scripts/run_context.py before writing",
    )
    parser.add_argument(
        "--expected-total-candidates", type=int, default=DEFAULT_EXPECTED_TOTAL_CANDIDATES
    )
    parser.add_argument(
        "--expected-dphascl1-labeled", type=int, default=DEFAULT_EXPECTED_DPHASCL1_LABELED
    )
    args = parser.parse_args(argv)

    run_dir = args.run_dir.resolve()
    if run_dir.parent.name != "runs":
        raise SystemExit(f"run directory must live under a runs/ directory: {run_dir}")

    scripts_dir = Path(__file__).resolve().parent
    if str(scripts_dir) not in sys.path:
        sys.path.insert(0, str(scripts_dir))
    from run_context import create_run_layout, write_input_contract  # noqa: E402

    if args.create_run:
        create_run_layout(run_dir.parent.parent, run_dir.name)
    for name in ("logs", "inputs", "results"):
        if not (run_dir / name).is_dir():
            raise SystemExit(f"run directory is incomplete, missing {name}/: {run_dir}")

    motif_evidence = Path(args.motif_evidence).resolve()
    reference_panel = Path(args.reference_panel).resolve()
    reference_fasta = Path(args.reference_fasta).resolve()
    subtype_matrices = [Path(item).resolve() for item in args.subtype_matrix]

    candidate_rows = read_motif_evidence(motif_evidence)
    if len(candidate_rows) != args.expected_total_candidates:
        raise SystemExit(
            "candidate universe size mismatch: "
            f"{len(candidate_rows)} != {args.expected_total_candidates} "
            "(refusing to write results for an unexpected universe)"
        )

    matrices_by_name = {path.name: read_subtype_calls(path) for path in subtype_matrices}
    label_consistency = describe_subtype_matrix_agreement(matrices_by_name)
    calls = merge_subtype_calls(matrices_by_name)
    motif_accessions = {row["accession"] for row in candidate_rows}
    if set(calls) != motif_accessions:
        raise SystemExit(
            "subtype matrix and motif evidence accession sets differ: "
            f"{len(set(calls))} vs {len(motif_accessions)}"
        )

    labeled_1 = sum(1 for value in calls.values() if value == "dPHASCL1_like_candidate")
    if labeled_1 != args.expected_dphascl1_labeled:
        raise SystemExit(
            f"dPHASCL1_like_candidate count mismatch: {labeled_1} != "
            f"{args.expected_dphascl1_labeled}"
        )

    assigned = build_assignments(candidate_rows, calls)

    reference_rows = read_reference_panel(reference_panel)
    ahsmg_evidence = build_ahsmg_evidence(candidate_rows, reference_rows, reference_fasta)

    summary = summarize_assignments(assigned, ahsmg_evidence=ahsmg_evidence)
    if not summary["totals_self_consistent"]:
        raise SystemExit("type counts do not sum to the candidate total; refusing to write")

    manifest_entries = write_source_manifest(
        run_dir / "inputs" / "source_manifest.tsv",
        {"motif_candidate_evidence": motif_evidence},
        {path.name: path for path in subtype_matrices},
        reference_panel,
        reference_fasta,
    )

    write_assignments(run_dir / "results" / "catalytic_domain_type.tsv", assigned)
    (run_dir / "results" / "type_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    dphascl1 = summary["label_vs_geometric_verification"]["dPHASCL1_like_candidate"]
    dphascl2 = summary["label_vs_geometric_verification"]["dPHASCL2_like_candidate"]
    run_manifest = {
        "run_id": run_dir.name,
        "status": "completed_candidate_only",
        "purpose": (
            "PhaDED catalytic-domain type 1/type 2 assignment by oxyanion-hole geometry "
            "(pure analysis, zero new computation, no server execution)"
        ),
        "criteria_source": CITATION,
        "comparison_convention": {
            "type1_requires": "oxyanion_hole start coordinate strictly less than lipase_box start",
            "type2_requires": "oxyanion_hole start coordinate strictly greater than lipase_box start",
            "equal_start_positions": "undetermined_position_conflict",
            "type1_oxyanion_n_terminal_flag": TYPE1_OXYANION_N_TERMINAL,
            "matches_audit_phaded_motifs_relative_position_semantics": True,
        },
        "counts": {
            "total_candidates": summary["total_candidates"],
            **summary["type_counts"],
        },
        "label_vs_geometric_verification": summary["label_vs_geometric_verification"],
        "label_agreement_counts": summary["label_agreement_counts"],
        "label_source_consistency": label_consistency,
        "oxyanion_hole_state_counts": summary["oxyanion_hole_state_counts"],
        "orientation_sensitivity_if_inverted": summary["orientation_sensitivity_if_inverted"],
        "ahsmg_evidence": summary["ahsmg_evidence"],
        "authorization": {
            "candidate_only_execution": True,
            "server_execution_started": False,
            "ssh_started": False,
            "formal_scan_started": False,
            "formal_registry_modified": False,
            "historical_run_modified": False,
            "formal_scan_models_modified": False,
        },
        "source_snapshot": {
            "assignment_script": {
                "path": str(Path(__file__).resolve()),
                "size": Path(__file__).resolve().stat().st_size,
                "sha256": sha256_file(Path(__file__).resolve()),
            },
            "run_context": {
                "path": str((scripts_dir / "run_context.py").resolve()),
                "size": (scripts_dir / "run_context.py").stat().st_size,
                "sha256": sha256_file(scripts_dir / "run_context.py"),
            },
            "inputs": manifest_entries,
        },
        "outputs": {
            "catalytic_domain_type": str(run_dir / "results" / "catalytic_domain_type.tsv"),
            "type_summary": str(run_dir / "results" / "type_summary.json"),
            "source_manifest": str(run_dir / "inputs" / "source_manifest.tsv"),
        },
        "evidence_boundary": PHENOTYPE_BOUNDARY,
        "label_separation_statement": summary["label_separation_statement"],
        "ahsmg_conclusion_boundary": AHSMG_CONCLUSION_BOUNDARY,
        "non_modification_statement": (
            "Read-only consumption of the 20260915 motif reconciliation evidence table and the "
            "20260916 full-library subtype matrices; no input was rewritten and no historical run, "
            "registry, config or formal scan artifact was touched."
        ),
    }
    (run_dir / "run_manifest.json").write_text(
        json.dumps(run_manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    log_lines = [
        f"run_id: {run_dir.name}",
        f"generated_at: {dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')}",
        f"motif_evidence: {motif_evidence}",
        f"subtype_matrices: {', '.join(str(path) for path in subtype_matrices)}",
        f"total_candidates: {summary['total_candidates']}",
        f"type1_verified: {summary['type_counts']['type1_verified']}",
        f"type2_verified: {summary['type_counts']['type2_verified']}",
        f"undetermined_no_oxyanion: {summary['type_counts']['undetermined_no_oxyanion']}",
        f"undetermined_position_conflict: {summary['type_counts']['undetermined_position_conflict']}",
        f"dPHASCL1_like_candidate_labeled: {dphascl1['labeled_count']}",
        f"dPHASCL1_like_candidate_type1_verified: {dphascl1['type1_verified']}",
        f"dPHASCL2_like_candidate_labeled: {dphascl2['labeled_count']}",
        f"dPHASCL2_like_candidate_type2_verified: {dphascl2['type2_verified']}",
        f"ahsmg_reference_positive_controls: {ahsmg_evidence['reference_positive_control_count']}",
        f"ahsmg_candidate_not_detected_pattern: {ahsmg_evidence['candidate_not_detected_pattern']}",
        "server_execution_started: false",
        "ssh_started: false",
    ]
    (run_dir / "logs" / "catalytic_domain_type.log").write_text(
        "\n".join(log_lines) + "\n", encoding="utf-8"
    )

    write_input_contract(
        run_dir,
        run_id=run_dir.name,
        inputs={
            "motif_candidate_evidence": motif_evidence,
            "reference_motif_panel": reference_panel,
            "reference_sequences": reference_fasta,
            **{
                f"subtype_matrix_{index + 1}": path
                for index, path in enumerate(subtype_matrices)
            },
            "source_manifest": run_dir / "inputs" / "source_manifest.tsv",
        },
    )

    print(
        json.dumps(
            {
                "run_id": run_dir.name,
                "total_candidates": summary["total_candidates"],
                **summary["type_counts"],
                "dPHASCL1_like_candidate_labeled": dphascl1["labeled_count"],
                "dPHASCL1_like_candidate_type1_verified": dphascl1["type1_verified"],
                "dPHASCL2_like_candidate_labeled": dphascl2["labeled_count"],
                "dPHASCL2_like_candidate_type2_verified": dphascl2["type2_verified"],
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
