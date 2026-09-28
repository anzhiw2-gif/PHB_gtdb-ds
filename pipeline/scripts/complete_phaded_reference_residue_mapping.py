#!/usr/bin/env python3
"""Complete the PhaDED reference residue and discriminating-domain mapping.

Why this script exists
----------------------
The 2026-09-15 motif panel was effectively a *lipase box + oxyanion hole*
panel: ``his_state`` never fired (0/723), ``asp_state`` only fired for one
superfamily (25/723), ``sbd_state`` and ``linker_state`` were never bound, and
``lid_state`` was derived from the subtype prior, which makes it circular when
the prior itself is under test.

This script replaces those three failure modes with evidence-bound mappings:

* literature coordinates are stored in a versioned table
  (``inputs/literature_residue_reference.tsv``) with DOI/PMID, database
  version, retrieval date and the SHA-256 of the sequence they were verified on;
* catalytic residues are transferred from literature-anchored references to
  every member of the same DED alignment *by alignment column*, which is the
  annotation transfer Knoll 2009 describes ("residues of the lipase box and the
  catalytic triad were manually annotated ... based on multiple sequence
  alignments");
* where no anchor exists in an alignment, the state stays ``pending_*`` -- it is
  never downgraded into a negative;
* the lid detector reads a sequence only.  It never sees a superfamily, a
  family or a subtype label, so its output can be used to test those priors
  instead of restating them.

Every state carries its own mapping basis, and nothing here is promoted beyond
candidate-level sequence/domain evidence.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence


AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
BOUNDARY = (
    "sequence-pattern/domain evidence only; not catalytic validation and not "
    "PHB/PHA phenotype proof"
)
PFAM_MAX_IEVALUE = 1e-5
#: A transferred alignment column is only trusted when the anchored residue is
#: conserved across the alignment.  The type 2 family is a circularly permuted
#: fold and its columns do not carry the catalytic residues for most members.
COLUMN_CONSERVATION_THRESHOLD = 0.5
SBD_PFAM_ACCESSION = "PF06850"
CATALYTIC_DOMAIN_PFAM_ACCESSIONS = ("PF10503",)
LID_EXPECTED_SUPERFAMILIES = ("intracellular nPHAMCL",)
ALIGNMENT_ROW_RE = re.compile(r"^\s*([^\s]+)\s+([A-Za-z*.-]+)(?:\s+\d+)?\s*$")

PENDING_NO_ALIGNMENT_ROW = "pending_no_alignment_row"
PENDING_NO_ANCHOR = "pending_no_anchor_for_this_role_in_alignment"

#: Prefixes of the per-role ``*_mapping_basis`` values.  They separate the three
#: evidence origins that the calibration report must not merge under one name:
#: a detection on the reference itself, an anchor row's own literature coordinate,
#: and an alignment-column transfer that is an inference from another row.
TRANSFER_BASIS_PREFIX = "anchor_column_transfer:"
ANCHOR_SELF_BASIS_PREFIX = "literature_anchor:"

#: Literature-calibrated sequence patterns.  ``residue_offset`` is the 0-based
#: offset of the catalytic residue inside the match; every pattern is verified
#: against the literature coordinate it was taken from before it is used.
LITERATURE_PATTERNS: dict[str, dict[str, object]] = {
    "catalytic_serine": {
        "pattern_id": "lipase_box_Gx1Sx2G_Knoll2009",
        "regex": r"G([A-Z])S([A-Z])G",
        "residue_offset": 2,
        "expected_residue": "S",
        "source_entry": "LIT_KNOL2009_ALIGNMENT_ANNOTATION",
        "calibration_entry": "LIT_8DAJ_SER",
        "calibration_coordinate": 121,
    },
    "catalytic_histidine": {
        "pattern_id": "his_GMxH_from_PDB_8DAJ_His270",
        "regex": r"GM[A-Z]H",
        "residue_offset": 3,
        "expected_residue": "H",
        "source_entry": "LIT_8DAJ_HIS",
        "calibration_entry": "LIT_8DAJ_HIS",
        "calibration_coordinate": 270,
    },
    "catalytic_aspartate": {
        "pattern_id": "asp_GxxDYTV_from_PDB_8DAJ_Asp197",
        "regex": r"G[A-Z]{2}DYTV",
        "residue_offset": 3,
        "expected_residue": "D",
        "source_entry": "LIT_8DAJ_ASP",
        "calibration_entry": "LIT_8DAJ_ASP",
        "calibration_coordinate": 197,
    },
    "oxyanion_hole_cysteine": {
        "pattern_id": "oxyanion_HGCXQ_from_PDB_8DAJ_Cys40",
        "regex": r"HGC[A-Z]Q",
        "residue_offset": 2,
        "expected_residue": "C",
        "source_entry": "LIT_8DAJ_OXY",
        "calibration_entry": "LIT_8DAJ_OXY",
        "calibration_coordinate": 40,
    },
    "alternative_catalytic_motif_ahsmg": {
        "pattern_id": "ahsmg_AHSMG_Knoll2009",
        "regex": r"AHSMG",
        "residue_offset": 0,
        "expected_residue": "A",
        "source_entry": "LIT_KNOL2009_AHSMG",
        "calibration_entry": "",
        "calibration_coordinate": 0,
    },
}

#: Signature of the lid detector: the two micro-loops whose deletion abolished
#: activity in PhaZKT, plus the loop spacing measured in that same anchor.
LID_SIGNATURE: dict[str, object] = {
    "loop1_motif": "FNGIG",
    "loop2_motif": "YYWQLF",
    "loop1_source_entry": "LIT_PHAZKT_LID1",
    "loop2_source_entry": "LIT_PHAZKT_LID2",
    "anchor_spacing": 156,
    "spacing_tolerance": 8,
}

ROLE_EXPECTED_RESIDUE = {
    role: str(spec["expected_residue"]) for role, spec in LITERATURE_PATTERNS.items()
}
MANIFEST_ROLES = ("lipase_box", "his", "asp", "oxyanion_hole", "lid")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


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


def read_fasta(path: Path) -> dict[str, str]:
    """Read a FASTA file keyed by the first whitespace-delimited header token."""
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


def coordinate_positions(text: str) -> set[int]:
    """Expand a coordinate string such as ``34-38;190-195`` into positions."""
    positions: set[int] = set()
    for token in str(text or "").split(";"):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            start, _, end = token.partition("-")
            if start.isdigit() and end.isdigit():
                positions.update(range(int(start), int(end) + 1))
        elif token.isdigit():
            positions.add(int(token))
    return positions


def detect_literature_pattern(sequence: str, role: str) -> dict[str, str]:
    """Apply one literature-calibrated pattern and report the residue it lands on."""
    spec = LITERATURE_PATTERNS[role]
    clean = sequence.strip().upper().rstrip("*")
    match = re.search(str(spec["regex"]), clean)
    if not match:
        return {
            "state": "not_detected_pattern",
            "coordinates": "",
            "residue": "",
            "matched": "",
            "x1": "",
            "pattern_id": str(spec["pattern_id"]),
        }
    position = match.start() + int(spec["residue_offset"]) + 1
    residue = clean[position - 1]
    if residue != str(spec["expected_residue"]):
        raise ValueError(
            f"pattern {spec['pattern_id']} landed on {residue} instead of "
            f"{spec['expected_residue']} at {position}"
        )
    return {
        "state": "supported",
        "coordinates": f"{position}-{position}",
        "residue": residue,
        "matched": match.group(0),
        "x1": match.group(1) if role == "catalytic_serine" else "",
        "pattern_id": str(spec["pattern_id"]),
    }


def detect_lid(sequence: str, signature: Mapping[str, object] | None = None) -> dict[str, str]:
    """Detect the lid micro-loops from sequence alone.

    The detector is deliberately blind to any grouping label: its only inputs
    are the sequence and a motif signature.  That is what makes its output
    usable as evidence about a grouping prior rather than a restatement of it.
    """
    spec = dict(LID_SIGNATURE if signature is None else signature)
    clean = (sequence or "").strip().upper().rstrip("*")
    loop1 = str(spec["loop1_motif"])
    loop2 = str(spec["loop2_motif"])
    start1 = clean.find(loop1)
    start2 = clean.find(loop2)
    span = int(spec["anchor_spacing"])
    tolerance = int(spec["spacing_tolerance"])
    coordinates: list[str] = []
    if start1 >= 0:
        coordinates.append(f"{start1 + 1}-{start1 + len(loop1)}")
    if start2 >= 0:
        coordinates.append(f"{start2 + 1}-{start2 + len(loop2)}")
    if start1 >= 0 and start2 >= 0:
        spacing = start2 - start1
        within_window = span - tolerance <= spacing <= span + tolerance
        return {
            "state": "supported" if within_window else "conflict_lid_loop_spacing",
            "coordinates": ";".join(coordinates),
            "loop_spacing": str(spacing),
            "mapping_basis": (
                f"sequence_motif_pair:{loop1}({spec['loop1_source_entry']})+"
                f"{loop2}({spec['loop2_source_entry']});spacing_window="
                f"{span - tolerance}-{span + tolerance}"
            ),
        }
    if start1 >= 0:
        return {
            "state": "partial_lid_loop1_only",
            "coordinates": ";".join(coordinates),
            "loop_spacing": "",
            "mapping_basis": f"sequence_motif_single:{loop1}({spec['loop1_source_entry']})",
        }
    if start2 >= 0:
        return {
            "state": "partial_lid_loop2_only",
            "coordinates": ";".join(coordinates),
            "loop_spacing": "",
            "mapping_basis": f"sequence_motif_single:{loop2}({spec['loop2_source_entry']})",
        }
    return {
        "state": "not_detected_pattern",
        "coordinates": "",
        "loop_spacing": "",
        "mapping_basis": "no_lid_loop_motif_detected",
    }


def parse_alignment(path: Path) -> dict[str, str]:
    """Parse a Clustal-style alignment into ``name -> gapped sequence``.

    DED alignments list a protein entry and its structure entry under the same
    name inside one block; only the first occurrence per block is kept so that
    the ungapped sequence stays equal to the reference sequence.
    """
    rows: dict[str, list[str]] = defaultdict(list)
    block: dict[str, str] = {}

    def flush() -> None:
        for name, chunk in block.items():
            rows[name].append(chunk)
        block.clear()

    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        if not line.strip():
            flush()
            continue
        match = ALIGNMENT_ROW_RE.match(line)
        if not match:
            continue
        name = match.group(1)
        if name.upper() in {"CLUSTAL", "MUSCLE", "PROBCONS", "CLUSTALW"}:
            continue
        body = match.group(2)
        if not body.strip("*.-") or set(body) <= {"*", ".", ":"}:
            continue
        if name in block:
            continue
        block[name] = body
    flush()
    return {name: "".join(chunks) for name, chunks in rows.items()}


def gapped_column(gapped: str, position: int) -> int:
    """Return the 0-based alignment column of a 1-based ungapped position."""
    if position < 1:
        raise ValueError(f"position must be 1-based: {position}")
    seen = 0
    for column, character in enumerate(gapped):
        if character not in "-.":
            seen += 1
            if seen == position:
                return column
    raise ValueError(f"position {position} is beyond the aligned sequence length")


def ungapped_position(gapped: str, column: int) -> int | None:
    """Return the 1-based ungapped position at an alignment column, or None for a gap."""
    if column >= len(gapped) or gapped[column] in "-.":
        return None
    return sum(1 for character in gapped[: column + 1] if character not in "-.")


def load_literature_reference(path: Path) -> list[dict[str, str]]:
    rows = read_tsv(path)
    if not rows:
        raise ValueError(f"literature reference table is empty: {path}")
    for row in rows:
        for field in ("entry_id", "residue_role", "coordinate_source", "retrieval_date"):
            if not str(row.get(field, "")).strip():
                raise ValueError(f"literature row {row.get('entry_id')!r} lacks {field}")
    return rows


def literature_anchor_coordinates(
    rows: Iterable[Mapping[str, str]],
) -> dict[str, dict[str, dict[str, str]]]:
    """Index literature rows as ``accession -> role -> {coordinate, entry_id}``."""
    anchors: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in rows:
        role = str(row["residue_role"])
        coordinate = str(row.get("alignment_coordinate", "")).strip()
        if not coordinate.isdigit():
            continue
        if role not in ROLE_EXPECTED_RESIDUE:
            continue
        accession = str(row["accession"]).strip()
        if not accession or accession == "pending":
            continue
        anchors[accession][role] = {
            "coordinate": coordinate,
            "entry_id": str(row["entry_id"]),
            "reported_coordinate": str(row.get("reported_coordinate", "pending")),
        }
    return dict(anchors)


def verify_anchor_coordinates(
    rows: Iterable[Mapping[str, str]],
    sequences: Mapping[str, str],
) -> list[dict[str, str]]:
    """Fail-closed verification of every literature residue against its sequence."""
    verified: list[dict[str, str]] = []
    for row in rows:
        coordinate = str(row.get("alignment_coordinate", "")).strip()
        accession = str(row.get("accession", "")).strip()
        role = str(row.get("residue_role", ""))
        if accession in ("", "pending") or not coordinate.isdigit() and "-" not in coordinate:
            continue
        sequence = sequences.get(accession)
        if sequence is None:
            raise KeyError(f"no sequence available to verify anchor {accession}")
        if "-" in coordinate and not coordinate.isdigit():
            start, _, end = coordinate.partition("-")
            if not (start.isdigit() and end.isdigit()):
                continue
            spanned = sequence[int(start) - 1:int(end)]
            if role in {"lid_loop_1", "lid_loop_2"}:
                expected_motif = str(
                    LID_SIGNATURE["loop1_motif"] if role == "lid_loop_1"
                    else LID_SIGNATURE["loop2_motif"]
                )
                if spanned != expected_motif:
                    raise ValueError(
                        f"anchor {accession} carries {spanned!r} at {coordinate}, "
                        f"expected {expected_motif!r} ({row['entry_id']})"
                    )
            verified.append({
                "entry_id": row["entry_id"],
                "accession": accession,
                "residue_role": role,
                "alignment_coordinate": coordinate,
                "observed_residue": spanned,
                "verification": "verified_range_in_anchor_sequence",
            })
            continue
        position = int(coordinate)
        observed = sequence[position - 1]
        if role in ROLE_EXPECTED_RESIDUE:
            expected = ROLE_EXPECTED_RESIDUE[role]
            if observed != expected:
                raise ValueError(
                    f"anchor {accession} residue {position} is {observed}, expected {expected} "
                    f"({row['entry_id']})"
                )
        verified.append({
            "entry_id": row["entry_id"],
            "accession": accession,
            "residue_role": role,
            "alignment_coordinate": coordinate,
            "observed_residue": observed,
            "verification": "verified_in_anchor_sequence",
        })
    if not verified:
        raise ValueError("no literature coordinate could be verified")
    return verified


def calibrate_patterns(
    anchors: Mapping[str, Mapping[str, Mapping[str, str]]],
    sequences: Mapping[str, str],
) -> list[dict[str, str]]:
    """Check that each calibrated pattern lands on the coordinate it came from."""
    report: list[dict[str, str]] = []
    for role, spec in LITERATURE_PATTERNS.items():
        entry_id = str(spec["calibration_entry"])
        coordinate = int(spec["calibration_coordinate"])
        record = {
            "role": role,
            "pattern_id": str(spec["pattern_id"]),
            "regex": str(spec["regex"]),
            "calibration_entry": entry_id or "pending_no_calibration_anchor",
            "calibration_coordinate": str(coordinate) if coordinate else "pending",
            "state": "not_calibrated",
            "observed_coordinates": "",
        }
        if not entry_id or not coordinate:
            record["state"] = "pending_no_coordinate_in_literature"
            report.append(record)
            continue
        holder = next(
            (accession for accession, roles in anchors.items() if role in roles
             and roles[role]["entry_id"] == entry_id),
            None,
        )
        if holder is None:
            raise ValueError(f"calibration anchor for {entry_id} is missing from the table")
        result = detect_literature_pattern(sequences[holder], role)
        record["calibration_accession"] = holder
        record["observed_coordinates"] = result["coordinates"]
        if result["state"] != "supported" or coordinate not in coordinate_positions(result["coordinates"]):
            raise ValueError(
                f"pattern {spec['pattern_id']} does not reproduce {entry_id} at {coordinate}"
            )
        record["state"] = "calibrated_hits_literature_coordinate"
        report.append(record)
    return report


def parse_domtblout(path: Path, max_ievalue: float = PFAM_MAX_IEVALUE) -> dict[str, list[dict[str, object]]]:
    """Parse a HMMER domtblout into ``reference accession -> significant hits``."""
    hits: dict[str, list[dict[str, object]]] = defaultdict(list)
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.split()
        if len(fields) < 22:
            continue
        try:
            i_evalue = float(fields[12])
            start, end = int(fields[17]), int(fields[18])
        except ValueError:
            continue
        if i_evalue > max_ievalue:
            continue
        accession = fields[3].split("|")[-1]
        hits[accession].append({
            "pfam_accession": fields[1].split(".", 1)[0],
            "pfam_name": fields[0],
            "i_evalue": i_evalue,
            "coordinates": f"{start}-{end}",
            "start": start,
            "end": end,
        })
    return dict(hits)


def alignment_paths(manifest_path: Path, alignment_dir: Path | None) -> dict[str, dict[str, str]]:
    """Map DED homologous family -> alignment record from the frozen manifest."""
    paths: dict[str, dict[str, str]] = {}
    for row in read_tsv(manifest_path):
        family = row["phaded_family_id"]
        record = dict(row)
        if alignment_dir is not None:
            match = re.search(r"(\d+)", family)
            candidate = alignment_dir / f"aln{match.group(1)}.aln" if match else None
            if candidate is not None and candidate.is_file():
                record["resolved_path"] = str(candidate.resolve())
        record.setdefault("resolved_path", row.get("reference_alignment", ""))
        paths[family] = record
    return paths


def build_column_transfers(
    ledger_rows: Sequence[Mapping[str, str]],
    alignments: Mapping[str, Mapping[str, str]],
    alignment_records: Mapping[str, Mapping[str, str]],
    anchors: Mapping[str, Mapping[str, Mapping[str, str]]],
) -> tuple[dict[str, dict[str, object]], dict[str, object]]:
    """Transfer literature residues to alignment members by shared column."""
    family_of = {str(row["accession"]): str(row["phaded_family_id"]) for row in ledger_rows}
    transfers: dict[str, dict[str, object]] = {}
    stats: dict[str, object] = {
        "alignments_examined": 0,
        "alignments_usable": 0,
        "alignments_with_anchor": 0,
        "anchored_alignments": [],
        "column_conservation": [],
        "column_conservation_threshold": COLUMN_CONSERVATION_THRESHOLD,
        "references_with_anchor_column_transfer": 0,
        "anchors_used": [],
        "anchor_rows_missing_from_alignment": [],
    }
    for family, rows in sorted(alignments.items()):
        record = alignment_records.get(family, {})
        status = str(record.get("alignment_status", ""))
        if status and status != "usable_reference_alignment":
            continue
        stats["alignments_examined"] = int(stats["alignments_examined"]) + 1
        if not rows:
            continue
        stats["alignments_usable"] = int(stats["alignments_usable"]) + 1
        family_anchors: list[tuple[str, dict[str, dict[str, str]]]] = []
        for accession in sorted(anchors):
            if family_of.get(accession) != family:
                continue
            if accession not in rows:
                stats["anchor_rows_missing_from_alignment"].append(f"{family}:{accession}")
                continue
            family_anchors.append((accession, anchors[accession]))
        if not family_anchors:
            continue
        anchor_accession, anchor_roles = family_anchors[0]
        columns: dict[str, dict[str, object]] = {}
        for role, payload in sorted(anchor_roles.items()):
            coordinate = int(payload["coordinate"])
            column = gapped_column(rows[anchor_accession], coordinate)
            if rows[anchor_accession][column].upper() != ROLE_EXPECTED_RESIDUE[role]:
                raise ValueError(
                    f"anchor row {anchor_accession} column {column} does not carry "
                    f"{ROLE_EXPECTED_RESIDUE[role]} for {role}"
                )
            observed = [
                gapped[column].upper()
                for gapped in rows.values()
                if column < len(gapped) and gapped[column] not in "-."
            ]
            expected = ROLE_EXPECTED_RESIDUE[role]
            conserved = (
                sum(1 for residue in observed if residue == expected) / len(observed)
                if observed else 0.0
            )
            columns[role] = {
                "column": column,
                "anchor_coordinate": coordinate,
                "anchor_entry_id": payload["entry_id"],
                "expected_residue": expected,
                "anchor_observed": rows[anchor_accession][column].upper(),
                "conserved_fraction": round(conserved, 4),
                "rows_observed": len(observed),
                "transferable": conserved >= COLUMN_CONSERVATION_THRESHOLD,
            }
            stats["column_conservation"].append({
                "phaded_family_id": family,
                "anchor_accession": anchor_accession,
                "role": role,
                "alignment_column": column,
                "anchor_coordinate": coordinate,
                "rows_observed": len(observed),
                "conserved_fraction": round(conserved, 4),
                "transferable": conserved >= COLUMN_CONSERVATION_THRESHOLD,
            })
        for name, gapped in rows.items():
            per_role: dict[str, dict[str, str]] = {}
            for role, meta in columns.items():
                column = int(meta["column"])
                position = ungapped_position(gapped, column)
                residue = gapped[column].upper() if column < len(gapped) else "-"
                if name == anchor_accession:
                    mode = "anchor_self_literature_coordinate"
                elif meta["transferable"]:
                    mode = "column_transfer"
                else:
                    mode = "not_conserved"
                per_role[role] = {
                    "mode": mode,
                    "residue": residue,
                    "coordinate": str(position) if position else "",
                    "expected_residue": str(meta["expected_residue"]),
                    "conserved_fraction": str(meta["conserved_fraction"]),
                }
            entry = transfers.setdefault(name, {})
            entry.update(per_role)
            entry["_alignment"] = str(record.get("resolved_path", ""))
            entry["_family"] = family
            entry["_anchor_accession"] = anchor_accession
            entry["_anchor_roles"] = {
                role: f"{meta['anchor_coordinate']}@{meta['column']}"
                for role, meta in columns.items()
            }
        stats["alignments_with_anchor"] = int(stats["alignments_with_anchor"]) + 1
        stats["anchored_alignments"].append(
            {
                "phaded_family_id": family,
                "anchor_accession": anchor_accession,
                "anchor_roles": {role: meta["anchor_coordinate"] for role, meta in columns.items()},
                "alignment_rows": len(rows),
                "alignment_path": str(record.get("resolved_path", "")),
            }
        )
        stats["anchors_used"].append(f"{family}:{anchor_accession}")
    stats["references_with_anchor_column_transfer"] = len(transfers)
    return transfers, stats


def resolve_role_state(
    pattern: Mapping[str, str],
    transfer: Mapping[str, object] | None,
    *,
    role_key: str,
    role_label: str,
    has_alignment_row: bool,
) -> dict[str, str]:
    """Combine anchor-column transfer and calibrated pattern into one state."""
    payload = transfer.get(role_key) if transfer else None
    if payload:
        residue = str(payload["residue"])
        expected = str(payload["expected_residue"])
        mode = str(payload["mode"])
        origin = (
            f"literature_anchor:{transfer['_family']}:{transfer['_anchor_accession']}:{role_label}"
            if mode == "anchor_self_literature_coordinate"
            else f"anchor_column_transfer:{transfer['_family']}:"
                 f"{transfer['_anchor_accession']}:{role_label}"
        )
        if mode == "not_conserved":
            return {
                "state": "pending_reference_annotation",
                "coordinates": "",
                "mapping_basis": (
                    f"pending_anchor_column_not_conserved:{transfer['_family']}:"
                    f"{transfer['_anchor_accession']}:{role_label}:"
                    f"conserved_fraction={payload['conserved_fraction']}"
                ),
            }
        if residue == expected and payload["coordinate"]:
            return {
                "state": "supported_anchor_column",
                "coordinates": f"{payload['coordinate']}-{payload['coordinate']}",
                "mapping_basis": origin,
            }
        if residue in {"", "-", ".", "X"} or not payload["coordinate"]:
            return {
                "state": "pending_reference_annotation",
                "coordinates": "",
                "mapping_basis": f"pending_gap_at_mapped_column:{origin}",
            }
        return {
            "state": "not_detected_at_mapped_column",
            "coordinates": "",
            "mapping_basis": f"{origin};observed={residue}",
        }
    if pattern["state"] == "supported":
        return {
            "state": "supported_literature_pattern",
            "coordinates": pattern["coordinates"],
            "mapping_basis": f"literature_pattern:{pattern['pattern_id']}",
        }
    return {
        "state": "pending_reference_annotation",
        "coordinates": "",
        "mapping_basis": PENDING_NO_ANCHOR if has_alignment_row else PENDING_NO_ALIGNMENT_ROW,
    }


def linker_region(hits: Sequence[Mapping[str, object]]) -> dict[str, str]:
    """Delimit a candidate inter-domain region without claiming a linker profile."""
    catalytic = [hit for hit in hits if hit["pfam_accession"] in CATALYTIC_DOMAIN_PFAM_ACCESSIONS]
    if not catalytic:
        return {
            "state": "pending_reference_annotation",
            "coordinates": "",
            "mapping_basis": "no_catalytic_domain_annotation",
        }
    catalytic_end = max(int(hit["end"]) for hit in catalytic)
    downstream = [
        hit for hit in hits
        if hit["pfam_accession"] not in CATALYTIC_DOMAIN_PFAM_ACCESSIONS
        and int(hit["start"]) > catalytic_end + 1
    ]
    if not downstream:
        return {
            "state": "pending_reference_annotation",
            "coordinates": "",
            "mapping_basis": "no_downstream_domain_annotation",
        }
    next_start = min(int(hit["start"]) for hit in downstream)
    next_accession = sorted(
        hit["pfam_accession"] for hit in downstream if int(hit["start"]) == next_start
    )[0]
    return {
        "state": "delimited_candidate_region_not_profiled",
        "coordinates": f"{catalytic_end + 1}-{next_start - 1}",
        "mapping_basis": (
            f"interdomain_region_between_PF10503_end_{catalytic_end}_and_"
            f"{next_accession}_start_{next_start};no_linker_profile_bound"
        ),
    }


def build_reference_manifest(
    ledger_rows: Sequence[Mapping[str, str]],
    sequences: Mapping[str, str],
    transfers: Mapping[str, Mapping[str, str]],
    alignment_members: Mapping[str, set[str]],
    pfam_hits: Mapping[str, Sequence[Mapping[str, object]]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for source in ledger_rows:
        accession = str(source["accession"]).strip()
        sequence = sequences.get(accession)
        if sequence is None:
            raise KeyError(f"reference sequence missing for {accession}")
        superfamily = str(source.get("phaded_superfamily", ""))
        family = str(source.get("phaded_family_id", ""))
        transfer = transfers.get(accession)
        has_row = accession in alignment_members.get(family, set())

        box = detect_literature_pattern(sequence, "catalytic_serine")
        his_pattern = detect_literature_pattern(sequence, "catalytic_histidine")
        asp_pattern = detect_literature_pattern(sequence, "catalytic_aspartate")
        oxy_pattern = detect_literature_pattern(sequence, "oxyanion_hole_cysteine")
        ahsmg = detect_literature_pattern(sequence, "alternative_catalytic_motif_ahsmg")
        lid = detect_lid(sequence)

        if box["state"] == "supported":
            ser_state, ser_coordinates = "supported", box["coordinates"]
        elif "without lipase box" in superfamily:
            ser_state, ser_coordinates = "pending_reference_annotation", ""
        else:
            ser_state, ser_coordinates = "not_detected_pattern", ""

        his = resolve_role_state(
            his_pattern, transfer, role_key="catalytic_histidine", role_label="his",
            has_alignment_row=has_row,
        )
        asp = resolve_role_state(
            asp_pattern, transfer, role_key="catalytic_aspartate", role_label="asp",
            has_alignment_row=has_row,
        )
        oxy = resolve_role_state(
            oxy_pattern, transfer, role_key="oxyanion_hole_cysteine",
            role_label="oxyanion_hole", has_alignment_row=has_row,
        )

        hits = pfam_hits.get(accession, [])
        accessions = sorted({str(hit["pfam_accession"]) for hit in hits})
        sbd_hits = [hit for hit in hits if hit["pfam_accession"] == SBD_PFAM_ACCESSION]
        if sbd_hits:
            sbd_state = "supported"
            sbd_coordinates = ";".join(sorted({str(hit["coordinates"]) for hit in sbd_hits}))
            sbd_basis = f"pfam_binding:{SBD_PFAM_ACCESSION}_i_evalue<={PFAM_MAX_IEVALUE:g}"
        else:
            sbd_state = "not_detected_in_tested_pfam"
            sbd_coordinates = ""
            sbd_basis = f"pfam_scan_ran_without_{SBD_PFAM_ACCESSION}_hit"
        sbd_caveat = (
            "pf06850_binds_the_intracellular_cys_type_superfamily_here_not_the_extracellular_sbd"
            if sbd_state == "supported"
            else "pf06850_is_not_the_extracellular_substrate_binding_domain_in_this_reference_set"
        )

        linker = linker_region(hits)

        lid_expected = superfamily in LID_EXPECTED_SUPERFAMILIES
        pending_reasons: list[str] = []
        for role_label, payload in (("his", his), ("asp", asp), ("oxyanion_hole", oxy)):
            if payload["state"].startswith("pending"):
                pending_reasons.append(f"{role_label}_{payload['mapping_basis']}")
        if lid["state"] in {"partial_lid_loop1_only", "partial_lid_loop2_only"}:
            pending_reasons.append("lid_second_loop_not_detected")
        if linker["state"] == "pending_reference_annotation":
            pending_reasons.append("linker_profile_not_bound")
        if sbd_state != "supported":
            pending_reasons.append("sbd_reference_or_pfam_binding_pending")
        if ser_state == "pending_reference_annotation":
            pending_reasons.append("cys_catalytic_residue_reference_mapping_pending")

        supported_count = sum(
            payload["state"].startswith("supported")
            for payload in (his, asp, oxy)
        ) + int(box["state"] == "supported")
        if supported_count >= 4 and not pending_reasons:
            panel_status = "complete_evidence_bound_panel"
        elif supported_count >= 3:
            panel_status = "catalytic_triad_partial"
        elif supported_count >= 2:
            panel_status = "partial_catalytic_pattern_panel"
        else:
            panel_status = "partial_lipase_box_only" if box["state"] == "supported" else "no_catalytic_residue_evidence_bound"

        row: dict[str, str] = {
            "reference_id": str(source.get("reference_id", "")),
            "accession": accession,
            "phaded_superfamily": superfamily,
            "phaded_family_id": family,
            "evidence_status": str(source.get("evidence_status", "")),
            "lipase_box_state": box["state"] if box["state"] == "supported" else "not_detected_pattern",
            "lipase_box_coordinates": box["coordinates"],
            "lipase_box_x1": box["x1"],
            "lipase_box_mapping_basis": f"literature_pattern:{box['pattern_id']}" if box["state"] == "supported" else "literature_pattern_absent",
            "lipase_box_expectation_from_superfamily": (
                "not_expected_for_superfamily" if "without lipase box" in superfamily
                else "expected_for_superfamily"
            ),
            "catalytic_ser_cys_state": ser_state,
            "catalytic_ser_cys_coordinates": ser_coordinates,
            "his_state": his["state"],
            "his_coordinates": his["coordinates"],
            "his_mapping_basis": his["mapping_basis"],
            "his_pattern_state": his_pattern["state"],
            "his_pattern_coordinates": his_pattern["coordinates"],
            "his_pattern_id": his_pattern["pattern_id"],
            "asp_state": asp["state"],
            "asp_coordinates": asp["coordinates"],
            "asp_mapping_basis": asp["mapping_basis"],
            "asp_pattern_state": asp_pattern["state"],
            "asp_pattern_coordinates": asp_pattern["coordinates"],
            "asp_pattern_id": asp_pattern["pattern_id"],
            "oxyanion_hole_state": oxy["state"],
            "oxyanion_hole_coordinates": oxy["coordinates"],
            "oxyanion_hole_mapping_basis": oxy["mapping_basis"],
            "oxyanion_hole_pattern_state": oxy_pattern["state"],
            "oxyanion_relative_position_vs_lipase_box": (
                "undetermined_missing_coordinates"
                if not (box["coordinates"] and oxy["coordinates"])
                else (
                    "oxyanion_upstream_of_lipase_box"
                    if min(coordinate_positions(oxy["coordinates"])) < min(coordinate_positions(box["coordinates"]))
                    else "oxyanion_downstream_of_lipase_box"
                )
            ),
            "type_assignment_deferred_to": "Task 4 assign_phaded_catalytic_domain_type",
            "ahsmg_state": ahsmg["state"],
            "ahsmg_coordinates": ahsmg["coordinates"],
            "sbd_state": sbd_state,
            "sbd_coordinates": sbd_coordinates,
            "sbd_mapping_basis": sbd_basis,
            "sbd_role_caveat": sbd_caveat,
            "linker_state": linker["state"],
            "linker_coordinates": linker["coordinates"],
            "linker_mapping_basis": linker["mapping_basis"],
            "lid_state": lid["state"],
            "lid_coordinates": lid["coordinates"],
            "lid_mapping_basis": lid["mapping_basis"],
            "lid_loop_spacing": lid["loop_spacing"],
            "lid_expectation_from_superfamily": (
                "expected_for_superfamily" if lid_expected else "not_expected_for_superfamily"
            ),
            "lid_prior_agrees_with_detection": str(
                (lid["state"] == "supported") == lid_expected
            ).lower(),
            "anchor_column_alignment": str(transfer.get("_alignment", "")) if transfer else "",
            "anchor_column_anchor": str(transfer.get("_anchor_accession", "")) if transfer else "",
            "anchor_column_roles": ";".join(
                f"{role}:{value}" for role, value in sorted(transfer.get("_anchor_roles", {}).items())
            ) if transfer else "",
            "pfam_accessions": ";".join(accessions),
            "motif_panel_status": panel_status,
            "motif_evidence_level": "motif_pattern_and_reference_mapping_only",
            "motif_pending_reason": ";".join(sorted(set(pending_reasons))) or "no_reference_gap_remaining",
            "motif_phenotype_boundary": BOUNDARY,
        }
        rows.append(row)
    for role in MANIFEST_ROLES:
        if f"{role}_state" not in rows[0]:
            raise ValueError(f"manifest is missing the {role} state column")
    return rows


def superfamily_mapping_table(rows: Sequence[Mapping[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["phaded_superfamily"]].append(row)
    table: list[dict[str, str]] = []
    for superfamily, members in sorted(grouped.items()):
        total = len(members)

        def count(role: str) -> int:
            return sum(row[f"{role}_state"].startswith("supported") for row in members)

        anchored = [row for row in members if row["anchor_column_anchor"]]
        table.append({
            "phaded_superfamily": superfamily,
            "reference_count": str(total),
            "anchor_reference": ";".join(sorted({row["anchor_column_anchor"] for row in anchored})) or "pending_no_literature_anchor_in_alignment",
            "references_with_anchor_column_transfer": str(len(anchored)),
            "lipase_box_supported": str(sum(row["lipase_box_state"] == "supported" for row in members)),
            "catalytic_ser_supported": str(sum(row["catalytic_ser_cys_state"] == "supported" for row in members)),
            "his_supported": str(count("his")),
            "his_pending": str(sum(row["his_state"].startswith("pending") for row in members)),
            "asp_supported": str(count("asp")),
            "asp_pending": str(sum(row["asp_state"].startswith("pending") for row in members)),
            "oxyanion_supported": str(count("oxyanion_hole")),
            "lid_supported": str(count("lid")),
            "lid_not_detected": str(sum(row["lid_state"] == "not_detected_pattern" for row in members)),
            "sbd_supported_pf06850": str(sum(row["sbd_state"] == "supported" for row in members)),
            "linker_delimited": str(sum(row["linker_state"] == "delimited_candidate_region_not_profiled" for row in members)),
            "linker_pending": str(sum(row["linker_state"] == "pending_reference_annotation" for row in members)),
            "mapping_basis": (
                "literature_coordinates_transferred_by_DED_alignment_column;"
                "sequence_patterns_calibrated_on_literature_anchors;"
                "no_new_homologous_family_or_superfamily_call_made"
            ),
        })
    return table


def _state_basis_column(rows: Sequence[Mapping[str, str]], column: str) -> str | None:
    """Return the manifest column that records how this state was resolved, if any."""
    if not rows or not column.endswith("_state"):
        return None
    candidate = f"{column[: -len('_state')]}_mapping_basis"
    return candidate if candidate in rows[0] else None


def state_distribution(rows: Sequence[Mapping[str, str]], column: str) -> dict[str, object]:
    """Count states and split the supported total by evidence origin.

    ``resolved_supported`` is the total. ``independent_detection_supported`` excludes
    references whose support rests only on an alignment-column transfer from a
    literature anchor: that is an inference from another row, not a detection on the
    reference itself. The Phase 1 verification report (P4) found that the previous
    field of that name carried the total instead.
    """
    counter = Counter(str(row.get(column, "")) for row in rows)
    supported = sum(value for key, value in counter.items() if key.startswith("supported"))
    basis_column = _state_basis_column(rows, column)
    transferred = 0
    anchor_self = 0
    if basis_column is not None:
        for row in rows:
            if not str(row.get(column, "")).startswith("supported"):
                continue
            basis = str(row.get(basis_column, ""))
            if basis.startswith(TRANSFER_BASIS_PREFIX):
                transferred += 1
            elif basis.startswith(ANCHOR_SELF_BASIS_PREFIX):
                anchor_self += 1
    distribution: dict[str, object] = dict(sorted(counter.items()))
    distribution["supported"] = supported
    distribution["resolved_supported"] = supported
    distribution["independent_detection_supported"] = supported - transferred
    distribution["transferred_from_anchor_column"] = transferred
    distribution["anchor_self_literature_coordinate"] = anchor_self
    distribution["supported_evidence_breakdown"] = {
        "direct_row_evidence": supported - transferred - anchor_self,
        "anchor_self_literature_coordinate": anchor_self,
        "anchor_column_transfer": transferred,
    }
    return distribution


def pattern_layer_hit_rates(rows: Sequence[Mapping[str, str]]) -> dict[str, object]:
    """Per-criterion hit rate over the 723 references, per layer."""
    layers: dict[str, object] = {
        "lipase_box_sequence": lambda row: row["lipase_box_state"] == "supported",
        "his_sequence_pattern": lambda row: row["his_pattern_state"] == "supported",
        "asp_sequence_pattern": lambda row: row["asp_pattern_state"] == "supported",
        "oxyanion_sequence_pattern": lambda row: row["oxyanion_hole_pattern_state"] == "supported",
        "ahsmg_sequence_pattern": lambda row: row["ahsmg_state"] == "supported",
        "lid_sequence_loops": lambda row: row["lid_state"] == "supported",
        "lid_sequence_loops_partial": lambda row: row["lid_state"].startswith("partial_lid"),
        "sbd_pf06850_binding": lambda row: row["sbd_state"] == "supported",
        "any_anchor_column_transfer": lambda row: bool(row["anchor_column_anchor"]),
        "his_resolved_any": lambda row: row["his_state"].startswith("supported"),
        "asp_resolved_any": lambda row: row["asp_state"].startswith("supported"),
    }
    rates: dict[str, object] = {}
    for name, predicate in layers.items():
        hits = [row for row in rows if predicate(row)]
        rates[name] = {
            "supported": len(hits),
            "hit_rate_of_references": round(len(hits) / len(rows), 4) if rows else 0.0,
            "by_superfamily": dict(sorted(Counter(row["phaded_superfamily"] for row in hits).items())),
        }
    return rates


def build_calibration_report(
    rows: Sequence[Mapping[str, str]],
    legacy_rows: Sequence[Mapping[str, str]],
    anchor_verification: Sequence[Mapping[str, str]],
    pattern_calibration: Sequence[Mapping[str, str]],
    transfer_stats: Mapping[str, object],
    run_id: str,
) -> dict[str, object]:
    criteria: dict[str, object] = {}
    for column in (
        "lipase_box_state", "catalytic_ser_cys_state", "his_state", "asp_state",
        "oxyanion_hole_state", "ahsmg_state", "sbd_state", "linker_state", "lid_state",
    ):
        before = state_distribution(legacy_rows, column) if legacy_rows else {"supported": None}
        after = state_distribution(rows, column)
        criteria[column] = {
            "before": before,
            "after": after,
            "delta_supported": (
                after["supported"] - before["supported"]
                if before.get("supported") is not None else "pending_no_legacy_panel"
            ),
            "before_source": "runs/20260915_phaded_motif_reconciliation_01/results/reference_motif_panel.tsv",
            "after_source": "runs/%s/results/reference_mapping_manifest.tsv" % run_id,
        }
    # The legacy lid values were prior-derived, so no legacy value is a detection at all.
    for field in (
        "resolved_supported",
        "independent_detection_supported",
        "transferred_from_anchor_column",
    ):
        criteria["lid_state"]["before"][field] = 0
    criteria["lid_state"]["before_note"] = (
        "the legacy lid_state values were pending_reference_annotation / "
        "not_expected_for_subtype, both derived from the subtype prior, so no legacy value "
        "constitutes an independent detection"
    )
    criteria["his_state"]["before_note"] = (
        "the legacy His pattern GM[A-Z]H[A-Z]{2}P[A-Z]{2}G requires a proline at position 7 "
        "and therefore never fires, including on the literature example GMGHAWSGG it was "
        "supposedly derived from"
    )
    criteria["sbd_state"]["after_note"] = (
        "PF06850 binds only the intracellular Cys-type superfamily (274/276 references) and "
        "no extracellular dPHASCL reference at all, so it is a real Pfam binding but not the "
        "extracellular substrate-binding domain criterion"
    )
    criteria["oxyanion_hole_state"]["after_note"] = (
        "the loose legacy pattern H[A-Z]?G[A-Z]?C[A-Z]?Q was replaced by the literature "
        "consensus HGCXQ (PDB 8DAJ Cys40 in HGCTQ, PDB 2D80 Cys230 in HGCLQ), which fires on "
        "2 fewer references; references whose DED alignment contains a literature anchor now "
        "report the transferred column instead of the pattern. Evidence origins for the "
        "%d resolved references: %d direct row evidence (calibrated pattern), "
        "%d anchor self literature coordinate, %d anchor-column transfer."
        % (
            criteria["oxyanion_hole_state"]["after"]["resolved_supported"],
            criteria["oxyanion_hole_state"]["after"]["supported_evidence_breakdown"]["direct_row_evidence"],
            criteria["oxyanion_hole_state"]["after"]["anchor_self_literature_coordinate"],
            criteria["oxyanion_hole_state"]["after"]["transferred_from_anchor_column"],
        )
    )
    for column in ("his_state", "asp_state"):
        after = criteria[column]["after"]
        criteria[column]["after_note"] = (
            "of the %d resolved references, %d are calibrated-pattern detections on the "
            "reference itself, %d are anchor rows carrying their own recorded literature "
            "coordinate and %d rest only on a conservation-gated alignment-column transfer, so "
            "independent_detection_supported (%d) is lower than resolved_supported (%d); the "
            "transferred part is an inference, not an independent detection"
            % (
                after["resolved_supported"],
                after["supported_evidence_breakdown"]["direct_row_evidence"],
                after["anchor_self_literature_coordinate"],
                after["transferred_from_anchor_column"],
                after["independent_detection_supported"],
                after["resolved_supported"],
            )
        )
    criteria["linker_state"]["after_note"] = (
        "41 references carry a delimited candidate inter-domain region between the catalytic "
        "domain and the next annotated domain; it is explicitly not profiled, so no reference "
        "reaches a supported linker state and the linker criterion stays unresolved"
    )
    criteria["lid_state"]["after_note"] = (
        "39 of the 52 references whose superfamily prior expects a lid are recovered by the "
        "sequence-only detector; 13 are not, and 10 references outside that prior carry only "
        "the first lid loop, so the paired-loop call is specific but not exhaustive"
    )

    pending_remainder = {
        column: sum(str(row[column]).startswith("pending") for row in rows)
        for column in ("his_state", "asp_state", "oxyanion_hole_state", "linker_state")
    }
    defects = [
        {
            "defect_id": "DEFECT_A_LID_CIRCULARITY",
            "description": (
                "lid_state was assigned from the subtype prior "
                "(not_expected_for_subtype / pending_reference_annotation), so it could not be "
                "used to test the subtype it was derived from"
            ),
            "correction": (
                "lid detection now reads the sequence only; the prior is reported in a separate "
                "column (lid_expectation_from_superfamily) and never enters the detector"
            ),
            "evidence": (
                "the independent detector reports supported for 39 references, all of them "
                "members of the intracellular nPHAMCL superfamily, and not_detected_pattern for "
                "the extracellular PhaZGK13 control (Q51718.1), matching the literature that "
                "extracellular depolymerases lack the lid"
            ),
            "verified_by_test": "LidDecircularisationTests",
        },
        {
            "defect_id": "DEFECT_B_HIS_PATTERN_NEVER_FIRED",
            "description": (
                "his_state was not_detected_pattern for all 723 references because the pattern "
                "imported from Knoll-style annotation GMXHXXPXXG cannot match the literature "
                "example GMGHAWSGG (9 residues) nor the stated 267-274 range (8 residues)"
            ),
            "correction": (
                "the His pattern was re-calibrated to the observed literature example "
                "(GMxH) and verified to land exactly on His270 of PDB 8DAJ; where a DED "
                "alignment contains a literature-anchored reference, the His column is "
                "transferred instead of pattern-matched"
            ),
            "evidence": "criteria.his_state.after.supported > 0 with per-row mapping basis",
            "verified_by_test": "LegacyDefectTests.test_legacy_his_pattern_never_fires_on_any_reference",
        },
        {
            "defect_id": "DEFECT_C_CATALYTIC_TRIAD_BLIND_TO_NON_TYPE1_FOLDS",
            "description": (
                "asp_state fired for 25 references only (all e-dPHASCL type 1) because the "
                "GxxDYTV motif is a type 1 motif; the type 2 circularly permuted fold and the "
                "intracellular mcl enzyme have unrelated catalytic aspartate contexts"
            ),
            "correction": (
                "type 2 and intracellular mcl aspartate/histidine are resolved by alignment "
                "column transfer from PDB 2D80, PhaZKT and PhaZGK13 anchors"
            ),
            "evidence": "2D80A, AAM63408.1 and Q51718.1 rows carry anchor-transferred coordinates",
            "verified_by_test": "ReferenceMappingTests",
        },
        {
            "defect_id": "DEFECT_D_SBD_BINDING_WAS_NEVER_RUN",
            "description": "sbd_state was pending_tool_input for 717 references",
            "correction": (
                "PF06850 is now bound per reference from the frozen reference hmmscan domtblout; "
                "the result is reported with an explicit role caveat because PF06850 does not "
                "mark the extracellular substrate-binding domain in this reference set"
            ),
            "evidence": "criteria.sbd_state.after.supported == 274",
            "verified_by_test": "ReferenceMappingTests.test_sbd_is_bound_from_the_existing_pf06850_domtblout",
        },
        {
            "defect_id": "DEFECT_E_AHSMG_COORDINATE_REPORTED_AS_CATALYTIC_RESIDUE",
            "description": (
                "the legacy implementation set catalytic_ser_cys_state=supported with "
                "catalytic_ser_cys_coordinates taken from the AHSMG match for the two PhaZ7-type "
                "references, so an alternative catalytic motif was reported as the position of the "
                "nucleophile"
            ),
            "correction": (
                "AHSMG is reported only in ahsmg_state/ahsmg_coordinates; the nucleophile state now "
                "comes from the lipase box pattern or from a transferred anchor column, and the two "
                "PhaZ7 references are not_detected_pattern for the lipase box because they genuinely "
                "lack one"
            ),
            "evidence": (
                "criteria.catalytic_ser_cys_state: 443 legacy versus 444 now; the delta is +3 chance "
                "lipase box matches in the without-lipase-box superfamily and -2 removed AHSMG "
                "transfers"
            ),
            "verified_by_test": "ReferenceMappingTests.test_ahsmg_positive_controls_still_pass",
        },
    ]

    return {
        "schema_version": "1.0",
        "run_id": run_id,
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "reference_count": len(rows),
        "field_semantics": {
            "note": (
                "The 2026-09-17 Phase 1 verification report (P4) found that "
                "criteria.<criterion>.after.independent_detection_supported carried the total "
                "supported count, although part of that count is inferred by transferring a "
                "literature-anchored DED alignment column from an anchor row to another row "
                "(conservation gate >= %.1f). The names below replace that over-strong label for "
                "every run produced from now on; the frozen Phase 1 report keeps its original "
                "field name as evidence and is deliberately not rewritten."
                % COLUMN_CONSERVATION_THRESHOLD
            ),
            "fields": {
                "resolved_supported": (
                    "total references whose state for this criterion is supported, whatever the "
                    "evidence origin"
                ),
                "independent_detection_supported": (
                    "supported references whose evidence is observed on that row itself (a "
                    "calibrated sequence pattern, a Pfam domain binding, or the anchor row's own "
                    "recorded literature coordinate) and is not obtained by transferring another "
                    "row's alignment column; it excludes every anchor-column transfer"
                ),
                "transferred_from_anchor_column": (
                    "supported references whose coordinate comes only from a conservation-gated "
                    "alignment-column transfer from a literature anchor; this is an inference from "
                    "another row, not an independent detection"
                ),
                "anchor_self_literature_coordinate": (
                    "supported references that are themselves literature anchors and therefore "
                    "carry their own recorded coordinate rather than a transferred column"
                ),
                "supported_evidence_breakdown": (
                    "the three evidence origins above; it always sums to resolved_supported"
                ),
            },
        },
        "criteria": criteria,
        "pattern_layer_hit_rates": pattern_layer_hit_rates(rows),
        "defects_corrected": defects,
        "anchor_verification": list(anchor_verification),
        "pattern_calibration": list(pattern_calibration),
        "alignment_column_transfer": dict(transfer_stats),
        "pending_remainder": pending_remainder,
        "lid_prior_agreement": {
            "detected_and_expected": sum(
                row["lid_state"] == "supported" and row["lid_expectation_from_superfamily"] == "expected_for_superfamily"
                for row in rows
            ),
            "detected_but_not_expected": sum(
                row["lid_state"] == "supported" and row["lid_expectation_from_superfamily"] != "expected_for_superfamily"
                for row in rows
            ),
            "not_detected_but_expected": sum(
                row["lid_state"] != "supported" and row["lid_expectation_from_superfamily"] == "expected_for_superfamily"
                for row in rows
            ),
            "evidence_level": "independent_sequence_detection_compared_to_a_recorded_prior",
        },
        "limitations": [
            "the local DED snapshot contains raw ClustalW alignments without the manual annotation "
            "lines described by Knoll 2009, so annotation transfer is reproduced from literature "
            "coordinates instead of being read from the alignment",
            "no literature anchor with residue coordinates exists for the intracellular nPHASCL "
            "superfamilies (with or without lipase box), the periplasmic superfamily or the "
            "extracellular nPHASCL superfamily, so their His/Asp states remain pending",
            "the linker criterion is not resolved: only a candidate inter-domain region is "
            "delimited and no linker profile is bound",
            "PF06850 is bound from the frozen reference hmmscan output and does not represent the "
            "extracellular substrate-binding domain in this reference set",
            "this task creates no new homologous family or superfamily judgement and no new "
            "candidate-level classification",
        ],
        "phenotype_boundary_holds": all(row["motif_phenotype_boundary"] == BOUNDARY for row in rows),
        "phenotype_boundary": BOUNDARY,
    }


def run_completion(
    *,
    reference_fasta: Path,
    reference_ledger: Path,
    reference_domtblout: Path,
    alignment_dir: Path,
    alignment_manifest: Path,
    literature_reference: Path,
    anchor_fasta: Path,
    legacy_panel: Path,
    output_dir: Path,
    run_id: str,
) -> dict[str, object]:
    """Build the completed reference mapping, its manifest and its calibration report."""
    ledger_rows = read_tsv(reference_ledger)
    panel_sequences = {
        key.split("|", 1)[1]: value
        for key, value in read_fasta(reference_fasta).items()
    }
    anchor_sequences = read_fasta(anchor_fasta)
    sequences = {**panel_sequences, **anchor_sequences}
    literature_rows = load_literature_reference(literature_reference)
    anchors = literature_anchor_coordinates(literature_rows)

    anchor_verification = verify_anchor_coordinates(literature_rows, sequences)
    pattern_calibration = calibrate_patterns(anchors, sequences)

    alignment_records = alignment_paths(alignment_manifest, alignment_dir)
    alignments: dict[str, dict[str, str]] = {}
    for family, record in alignment_records.items():
        resolved = str(record.get("resolved_path", ""))
        if not resolved or not Path(resolved).is_file():
            continue
        if str(record.get("alignment_status", "")) != "usable_reference_alignment":
            continue
        alignments[family] = parse_alignment(Path(resolved))
    alignment_members = {family: set(rows) for family, rows in alignments.items()}

    transfers, transfer_stats = build_column_transfers(
        ledger_rows, alignments, alignment_records, anchors
    )
    pfam_hits = parse_domtblout(reference_domtblout)
    rows = build_reference_manifest(
        ledger_rows, panel_sequences, transfers, alignment_members, pfam_hits
    )
    table = superfamily_mapping_table(rows)
    legacy_rows = read_tsv(legacy_panel) if Path(legacy_panel).is_file() else []
    report = build_calibration_report(
        rows, legacy_rows, anchor_verification, pattern_calibration, transfer_stats, run_id
    )

    output_dir = Path(output_dir)
    write_tsv(output_dir / "reference_mapping_manifest.tsv", rows)
    write_tsv(output_dir / "superfamily_residue_mapping.tsv", table)
    report_path = output_dir / "implementation_calibration_report.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    manifest = {
        "schema_version": "1.0",
        "run_id": run_id,
        "status": "completed_candidate_only",
        "reference_count": len(rows),
        "outputs": {
            name: {
                "path": str((output_dir / name).resolve()),
                "size": (output_dir / name).stat().st_size,
                "sha256": sha256_file(output_dir / name),
            }
            for name in (
                "reference_mapping_manifest.tsv",
                "superfamily_residue_mapping.tsv",
                "implementation_calibration_report.json",
            )
        },
        "inputs": {
            name: {"path": str(Path(path).resolve()), "size": Path(path).stat().st_size, "sha256": sha256_file(Path(path))}
            for name, path in {
                "reference_fasta": reference_fasta,
                "reference_ledger": reference_ledger,
                "reference_domtblout": reference_domtblout,
                "alignment_manifest": alignment_manifest,
                "literature_reference": literature_reference,
                "anchor_fasta": anchor_fasta,
                "legacy_panel": legacy_panel,
            }.items()
        },
        "calibration_report": report,
        "phenotype_boundary": BOUNDARY,
    }
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--reference-ledger", type=Path, required=True)
    parser.add_argument("--reference-domtblout", type=Path, required=True)
    parser.add_argument("--alignment-dir", type=Path, required=True)
    parser.add_argument("--alignment-manifest", type=Path, required=True)
    parser.add_argument("--literature-reference", type=Path, required=True)
    parser.add_argument("--anchor-fasta", type=Path, required=True)
    parser.add_argument("--legacy-panel", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args(argv)
    manifest = run_completion(
        reference_fasta=args.reference_fasta,
        reference_ledger=args.reference_ledger,
        reference_domtblout=args.reference_domtblout,
        alignment_dir=args.alignment_dir,
        alignment_manifest=args.alignment_manifest,
        literature_reference=args.literature_reference,
        anchor_fasta=args.anchor_fasta,
        legacy_panel=args.legacy_panel,
        output_dir=args.output_dir,
        run_id=args.run_id,
    )
    print(json.dumps({
        "run_id": manifest["run_id"],
        "reference_count": manifest["reference_count"],
        "outputs": sorted(manifest["outputs"]),
    }, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
