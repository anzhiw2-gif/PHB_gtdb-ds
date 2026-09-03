#!/usr/bin/env python3
"""Split the ePhaZ seed library into curated and discovery layers.

The splitter is deliberately conservative: a sequence enters the curated
layer only when the evidence bundle explicitly marks experimental support,
typical architecture and complete length.  All other valid sequences remain
in ``ePhaZ_broad_discovery`` and are recorded with a reason.  In particular,
short sequences are reviewed, never silently discarded.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


CORE_LAYER = "ePhaZ_curated_core"
BROAD_LAYER = "ePhaZ_broad_discovery"
OUTPUT_NAMES = (
    f"{CORE_LAYER}.faa",
    f"{BROAD_LAYER}.faa",
    "ePhaZ_layer_manifest.tsv",
    "ePhaZ_short_sequence_review.tsv",
)
PROVENANCE_NAME = "ePhaZ_split_manifest.json"

# Include common ambiguity symbols because these can occur in real protein
# records; punctuation, digits and stop characters are rejected.
VALID_AA = frozenset("ABCDEFGHIKLMNPQRSTVWXYZJUOB")
ACCESSION_RE = re.compile(r"^[^\s|]+$")


class SplitError(ValueError):
    """Raised when the seed/evidence contract cannot be proven."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _delimiter(path: Path) -> str:
    """Detect TSV/CSV input while remaining deterministic for one-column data."""
    first = path.read_text(encoding="utf-8-sig").splitlines()[:1]
    if not first:
        raise SplitError(f"empty tabular input: {path}")
    if "\t" in first[0]:
        return "\t"
    if "," in first[0]:
        return ","
    raise SplitError(f"manifest must be TSV or CSV: {path}")


def _read_table(path: Path, required: Iterable[str], label: str) -> List[Dict[str, str]]:
    if not path.is_file() or path.is_symlink():
        raise SplitError(f"{label} is missing or not a regular file: {path}")
    delimiter = _delimiter(path)
    try:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            fields = [field.strip() for field in (reader.fieldnames or []) if field]
            if not set(required).issubset(fields):
                missing = sorted(set(required) - set(fields))
                raise SplitError(f"{label} missing required columns: {','.join(missing)}")
            rows = []
            for row in reader:
                rows.append({(key or "").strip(): (value or "").strip() for key, value in row.items()})
    except UnicodeError as error:
        raise SplitError(f"cannot decode {label}: {error}") from error
    return rows


def _read_fasta(path: Path) -> List[Tuple[str, str, str]]:
    if not path.is_file() or path.is_symlink():
        raise SplitError(f"seed FASTA is missing or not a regular file: {path}")
    records: List[Tuple[str, str, str]] = []
    header: str | None = None
    sequence: List[str] = []

    def finish() -> None:
        nonlocal header, sequence
        if header is None:
            return
        accession = header[1:].split("|", 1)[0].strip()
        if not accession or not ACCESSION_RE.fullmatch(accession):
            raise SplitError(f"invalid FASTA accession in header: {header}")
        seq = "".join(sequence).upper()
        if not seq:
            raise SplitError(f"FASTA record has an empty sequence: {header}")
        invalid = sorted(set(seq) - VALID_AA)
        if invalid:
            raise SplitError(f"invalid amino-acid symbol(s) for {accession}: {','.join(invalid)}")
        records.append((accession, header, seq))

    try:
        lines = path.read_text(encoding="utf-8-sig").splitlines()
    except UnicodeError as error:
        raise SplitError(f"cannot decode seed FASTA: {error}") from error
    for raw in lines:
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            header, sequence = line, []
        else:
            if header is None:
                raise SplitError(f"sequence occurs before first FASTA header: {path}")
            sequence.append(line)
    finish()
    if not records:
        raise SplitError(f"FASTA contains no records: {path}")
    accessions = [record[0] for record in records]
    duplicates = sorted({acc for acc in accessions if accessions.count(acc) > 1})
    if duplicates:
        raise SplitError("duplicate FASTA accession(s): " + ",".join(duplicates))
    return records


def _canonical(value: str | None) -> str:
    return (value or "").strip().lower().replace("-", "_").replace(" ", "_")


def _is_experimental(value: str) -> bool:
    # Exact level is preferred; these explicit labels are accepted for
    # compatibility with prior curated evidence tables.
    return _canonical(value) in {
        "experimental",
        "experimentally_supported",
        "biochemical",
        "functional_validation",
        "literature_supported_experimental",
    }


def _is_typical(value: str) -> bool:
    return _canonical(value) in {"typical", "canonical", "pass", "complete", "confirmed"}


def _is_complete(value: str) -> bool:
    return _canonical(value) in {"complete", "full", "pass", "confirmed"}


def _manifest_rows(path: Path) -> Dict[str, Dict[str, str]]:
    rows = _read_table(path, {"accession"}, "seed manifest")
    result: Dict[str, Dict[str, str]] = {}
    for row in rows:
        accession = row.get("accession", "").strip()
        if not accession:
            raise SplitError("seed manifest contains an empty accession")
        if accession in result:
            raise SplitError(f"duplicate accession in seed manifest: {accession}")
        family = row.get("family") or row.get("query_group", "")
        if not family:
            raise SplitError(f"seed manifest has no family/query_group: {accession}")
        result[accession] = row
    if not result:
        raise SplitError("seed manifest contains no records")
    return result


def _evidence_rows(path: Path) -> Dict[str, Dict[str, str]]:
    rows = _read_table(path, {"accession"}, "curated evidence")
    columns = set(rows[0]) if rows else set()
    aliases = (
        {"evidence_level", "evidence"},
        {"architecture_status", "architecture"},
        {"completeness_status", "completeness"},
    )
    missing = ["/".join(sorted(options)) for options in aliases if not columns.intersection(options)]
    if missing:
        raise SplitError("curated evidence missing required field(s): " + ",".join(missing))
    result: Dict[str, Dict[str, str]] = {}
    for row in rows:
        accession = row.get("accession", "").strip()
        if not accession:
            raise SplitError("curated evidence contains an empty accession")
        if accession in result:
            raise SplitError(f"duplicate accession in curated evidence: {accession}")
        result[accession] = row
    if not result:
        raise SplitError("curated evidence contains no records")
    return result


def _write_fasta(path: Path, records: Sequence[Tuple[str, str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for _, header, sequence in records:
            handle.write(header + "\n" + sequence + "\n")


def _safe_output_path(outdir: Path) -> None:
    if outdir.exists():
        raise SplitError(f"output directory already exists: {outdir}")
    if outdir.name in {"", ".", ".."}:
        raise SplitError(f"unsafe output directory: {outdir}")
    # The caller may choose a temporary/synthetic directory in tests.  Reject
    # symlink components so a run cannot unexpectedly write elsewhere.
    current = outdir
    while True:
        if current.is_symlink():
            raise SplitError(f"output path contains a symlink: {current}")
        parent = current.parent
        if parent == current:
            break
        current = parent


def split_ephaz_seeds(
    seed_fasta: Path | str,
    manifest: Path | str,
    curated_evidence: Path | str,
    outdir: Path | str,
    short_review_tsv: Path | str | None = None,
) -> Mapping[str, object]:
    """Create the two ePhaZ layers and auditable review/checksum files."""
    seed_fasta, manifest, curated_evidence, outdir = map(Path, (seed_fasta, manifest, curated_evidence, outdir))
    _safe_output_path(outdir)
    if short_review_tsv is not None:
        requested_review = Path(short_review_tsv)
        if requested_review.resolve().parent != outdir.resolve():
            raise SplitError("short_review_tsv must be inside outdir")
        if requested_review.name != OUTPUT_NAMES[3]:
            raise SplitError(f"short_review_tsv must be named {OUTPUT_NAMES[3]}")

    records = _read_fasta(seed_fasta)
    manifest_by_acc = _manifest_rows(manifest)
    evidence_by_acc = _evidence_rows(curated_evidence)
    fasta_accessions = {record[0] for record in records}
    missing_manifest = sorted(fasta_accessions - set(manifest_by_acc))
    if missing_manifest:
        raise SplitError("FASTA accessions missing from seed manifest: " + ",".join(missing_manifest))
    extra_manifest = sorted(set(manifest_by_acc) - fasta_accessions)
    if extra_manifest:
        raise SplitError("manifest accessions missing from seed FASTA: " + ",".join(extra_manifest))

    layered: Dict[str, List[Tuple[str, str, str]]] = {CORE_LAYER: [], BROAD_LAYER: []}
    layer_rows: List[Dict[str, str]] = []
    short_rows: List[Dict[str, str]] = []
    for accession, header, sequence in records:
        meta = manifest_by_acc[accession]
        family_value = (meta.get("family") or meta.get("query_group") or "").lower()
        family_token = family_value.replace("-", "").replace("_", "")
        if "ephaz" not in family_token:
            raise SplitError(f"manifest family is not ePhaZ for {accession}: {family_value}")
        evidence = evidence_by_acc.get(accession, {})
        evidence_level = evidence.get("evidence_level") or evidence.get("evidence") or meta.get("evidence", "")
        architecture = evidence.get("architecture_status") or evidence.get("architecture") or meta.get("architecture", "")
        completeness = evidence.get("completeness_status") or evidence.get("completeness") or meta.get("completeness", "")
        length = len(sequence)
        reasons: List[str] = []
        if length < 200:
            reasons.append("short_sequence")
        if not _is_experimental(evidence_level):
            reasons.append("evidence_missing")
        if not _is_typical(architecture):
            reasons.append("architecture_pending")
        if not _is_complete(completeness):
            reasons.append("completeness_pending")
        is_core = not reasons
        layer = CORE_LAYER if is_core else BROAD_LAYER
        layered[layer].append((accession, header, sequence))
        review_status = ""
        if length < 200:
            review_status = "architecture_pending" if not _is_typical(architecture) else "review_required"
        row = {
            "accession": accession,
            "layer": layer,
            "length": str(length),
            "organism": meta.get("organism", ""),
            "protein_name": meta.get("protein_name", ""),
            "reviewed": meta.get("reviewed", ""),
            "evidence_level": evidence_level,
            "evidence_reference": evidence.get("reference") or evidence.get("references") or meta.get("evidence", ""),
            "architecture_status": architecture,
            "completeness_status": completeness,
            "review_status": review_status,
            "reason": ";".join(reasons) if reasons else "curated_core_criteria_pass",
            "input_header": header[1:],
        }
        layer_rows.append(row)
        if length < 200:
            short_rows.append(row.copy())

    if not layered[CORE_LAYER] or not layered[BROAD_LAYER]:
        raise SplitError("both ePhaZ layers must contain at least one sequence")

    # Stage all files in a sibling temporary directory, then rename only after
    # every validation and checksum has succeeded.  Failed runs leave no output.
    outdir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{outdir.name}.", dir=str(outdir.parent)))
    try:
        _write_fasta(staging / OUTPUT_NAMES[0], layered[CORE_LAYER])
        _write_fasta(staging / OUTPUT_NAMES[1], layered[BROAD_LAYER])
        fields = [
            "accession", "layer", "length", "organism", "protein_name", "reviewed",
            "evidence_level", "evidence_reference", "architecture_status", "completeness_status",
            "review_status", "reason", "input_header",
        ]
        with (staging / OUTPUT_NAMES[2]).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(layer_rows)
        with (staging / OUTPUT_NAMES[3]).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(short_rows)
        # Keep a compact provenance record next to the layer outputs.  This
        # binds the result to all three inputs and states the exact gate used
        # for curated membership without claiming phenotype validation.
        provenance = {
            "schema": "ephaz-layered-seeds-v1",
            "inputs": {
                "seed_fasta": {"path": str(seed_fasta), "bytes": seed_fasta.stat().st_size, "sha256": sha256_file(seed_fasta)},
                "manifest": {"path": str(manifest), "bytes": manifest.stat().st_size, "sha256": sha256_file(manifest)},
                "curated_evidence": {"path": str(curated_evidence), "bytes": curated_evidence.stat().st_size, "sha256": sha256_file(curated_evidence)},
            },
            "criteria": {
                "minimum_length_for_core": 200,
                "layer_core": ["evidence_level=experimental", "architecture_status=typical", "completeness_status=complete", "length>=200"],
                "short_sequence_policy": "retain_in_broad_and_review",
            },
            "counts": {CORE_LAYER: len(layered[CORE_LAYER]), BROAD_LAYER: len(layered[BROAD_LAYER]), "short": len(short_rows)},
            "outputs": list(OUTPUT_NAMES) + [PROVENANCE_NAME, "sha256.tsv"],
        }
        (staging / PROVENANCE_NAME).write_text(json.dumps(provenance, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")

        checksums = []
        for name in OUTPUT_NAMES:
            target = staging / name
            checksums.append({"path": name, "bytes": str(target.stat().st_size), "sha256": sha256_file(target)})
        checksums.append({"path": PROVENANCE_NAME, "bytes": str((staging / PROVENANCE_NAME).stat().st_size), "sha256": sha256_file(staging / PROVENANCE_NAME)})
        with (staging / "sha256.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["path", "bytes", "sha256"], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(checksums)
        staging.replace(outdir)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise

    return {
        "outdir": str(outdir),
        "outputs": {name: str(outdir / name) for name in OUTPUT_NAMES + (PROVENANCE_NAME, "sha256.tsv")},
        "counts": {CORE_LAYER: len(layered[CORE_LAYER]), BROAD_LAYER: len(layered[BROAD_LAYER]), "short": len(short_rows)},
        "sha256": {row["path"]: row["sha256"] for row in checksums},
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-fasta", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--curated-evidence", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--short-review-tsv")
    args = parser.parse_args(argv)
    result = split_ephaz_seeds(**vars(args))
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
