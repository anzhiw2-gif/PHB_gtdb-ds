#!/usr/bin/env python3
"""Build provenance-bound PhaDED family/superfamily profiles.

Discovery construction and functional calibration are separate decisions
(2026-09-28 evidence-model redesign, Task 5):

* this builder emits exactly two model layers.  ``reference_query_only`` covers one or
  two unique eligible sequences (or a family that cannot form a stable alignment);
  ``discovery_hmm_uncalibrated`` covers an HMM actually built here from at least three
  unique quality-controlled eligible sequences.  ``sequence_family_hmm_validated``
  (independent held-out/confounder validation) and ``calibrated_candidate_model``
  (explicitly authorized finalize action) are **never** emitted by this builder;
* the discovery layer only recalls and scores.  It never filters, deletes or demotes a
  candidate and never makes a family call, so nothing in its manifest may be read as
  functional discrimination.  In particular an ``intracellular nPHASCL with lipase box``
  (with-lipase) family still builds, but only as an uncalibrated discovery layer that
  stays in deferred structural review downstream;
* a direct contradictory experimental negative blocks *functional calibration*, not
  discovery construction.  It never silently removes the sequence from the reference
  audit: every ledger row of every profile is listed in ``profile_training_audit.tsv``
  with an explicit training/exclusion reason and, when applicable, the calibration
  block reason;
* MAFFT output is not bit-reproducible (MAFFT 7.525 yields different alignments for
  identical input), so ``bit_reproducible`` is always ``false`` and the concrete
  alignment SHA-256 is authoritative.  ``training_set_sha256`` is the real SHA-256 of
  the exact training FASTA that was written; profiles without a training set record
  ``pending``.

Discovery training eligibility: a row trains the discovery layer when its
``discovery_training_eligible`` column says ``true``.  The frozen 723-row ledger
predates that column, so when the column is **absent** eligibility falls back to
``evidence_status in {"experimental_positive", "annotation_only"}`` (quality-controlled
experimental and annotation-only references).  A row may never be eligible without a
sequence SHA-256.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Iterable


ALLOWED_EVIDENCE = {
    "experimental_positive", "experimental_negative", "challenge_control",
    "annotation_only", "pending_review",
}
ALLOWED_SUPERFAMILIES = {
    "intracellular nPHASCL without lipase box",
    "intracellular nPHASCL with lipase box",
    "periplasmic PHA depolymerases",
    "intracellular nPHAMCL",
    "extracellular dPHASCL type 1",
    "extracellular dPHASCL type 2",
    "extracellular native-SCL/PhaZ7-like",
    "extracellular dPHAMCL",
}
REQUIRED_COLUMNS = {
    "reference_id", "accession", "sequence_sha256", "phaded_superfamily",
    "phaded_family_id", "evidence_status", "experimental_assay",
    "experimental_result",
}

#: Optional ledger column added by the 2026-09-28 redesign.  Absent in the frozen ledger.
ELIGIBILITY_COLUMN = "discovery_training_eligible"
#: Documented fallback used when ``ELIGIBILITY_COLUMN`` is absent from the ledger.
DISCOVERY_TRAINING_ELIGIBLE_STATUSES = {"experimental_positive", "annotation_only"}
_TRUE_TOKENS = {"true", "1", "yes"}
_FALSE_TOKENS = {"false", "0", "no"}

#: A discovery HMM needs at least this many *unique* eligible sequences.
MINIMUM_DISCOVERY_TRAINING_SEQUENCES = 3
#: Layer names owned by this builder.  The validated/calibrated layers are never emitted here.
MODEL_LAYER_REFERENCE_QUERY_ONLY = "reference_query_only"
MODEL_LAYER_DISCOVERY = "discovery_hmm_uncalibrated"
BUILDER_MODEL_LAYERS = (MODEL_LAYER_REFERENCE_QUERY_ONLY, MODEL_LAYER_DISCOVERY)
NON_BUILDER_MODEL_LAYERS = ("sequence_family_hmm_validated", "calibrated_candidate_model")
#: Functional calibration statuses this builder may emit; it never calibrates a function.
FUNCTIONAL_CALIBRATION_NOT_CALIBRATED = "not_function_calibrated"
FUNCTIONAL_CALIBRATION_BLOCKED = "blocked_contradictory_experimental_negative"
#: with-lipase families stay a deferred structural-review layer (AGENTS.md).
DEFERRED_STRUCTURE_REVIEW_SUPERFAMILY = "intracellular nPHASCL with lipase box"

MANIFEST_COLUMNS = [
    "profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id",
    "training_accessions", "training_count", "model_status", "model_reason",
    "alignment_sha256", "hmm_sha256", "mafft_version", "hmmer_version",
    "alignment_path", "hmm_path",
    "model_layer", "functional_calibration_status", "training_sequence_count",
    "experimental_anchor_count", "annotation_only_count", "training_set_sha256",
    "bit_reproducible",
]
TRAINING_AUDIT_COLUMNS = [
    # ``model_layer``/``functional_calibration_status`` label every audit output row, so a
    # discovery-layer row is identifiable on its own (AGENTS.md: every discovery-layer
    # output carries its layer name).
    "profile_id", "profile_kind", "model_layer", "functional_calibration_status",
    "accession", "evidence_status", "discovery_training_eligible", "training_selected",
    "exclusion_reason", "blocks_functional_calibration", "block_reason",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"ledger is not a regular file: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        missing = sorted(REQUIRED_COLUMNS - set(fields))
        if missing:
            raise ValueError("ledger missing required columns: " + ",".join(missing))
        return fields, list(reader)


def _read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"reference FASTA is not a regular file: {path}")
    records: dict[str, tuple[str, str]] = {}
    current: str | None = None
    accession: str | None = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current is not None:
                if current in records or not chunks:
                    raise ValueError(f"invalid FASTA record: {current}")
                records[current] = (accession or "", "".join(chunks))
            header = line[1:].strip().split("|")
            if len(header) != 2 or not header[0] or not header[1]:
                raise ValueError("FASTA header must be >reference_id|accession")
            current, accession, chunks = header[0], header[1], []
        else:
            if current is None:
                raise ValueError("FASTA sequence appears before header")
            chunks.append(line)
    if current is not None:
        if current in records or not chunks:
            raise ValueError(f"invalid FASTA record: {current}")
        records[current] = (accession or "", "".join(chunks))
    if not records:
        raise ValueError(f"empty FASTA: {path}")
    return records


def _slug(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")
    return value or "unnamed"


def _profile_id(kind: str, value: str) -> str:
    """Keep paths readable while making punctuation-colliding IDs unique."""
    digest = hashlib.sha256(value.encode("utf-8")).hexdigest()[:12]
    return f"{kind}_{_slug(value)}_{digest}"


def _run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=True, capture_output=True, text=True)


def tool_version(executable: str, *, kind: str) -> str:
    flag = "--version" if kind == "mafft" else "-h"
    result = subprocess.run([executable, flag], check=False, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"cannot determine {kind} version for {executable}: exit {result.returncode}")
    return next((line.strip() for line in (result.stdout + "\n" + result.stderr).splitlines() if line.strip()), "pending")


def write_fasta(path: Path, records: dict[str, tuple[str, str]], accessions: Iterable[str]) -> None:
    by_accession = {accession: sequence for accession, sequence in records.values()}
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for accession in accessions:
            handle.write(f">{accession}\n{by_accession[accession]}\n")


def discovery_training_eligible(row: dict[str, str]) -> bool:
    """Decide whether one ledger row may train a *discovery-layer* HMM.

    The explicit ``discovery_training_eligible`` column is authoritative when present.
    The frozen ledger predates it, so an absent (or empty) value falls back to
    ``evidence_status in {"experimental_positive", "annotation_only"}``: quality-controlled
    direct experimental positives and annotation-only references both qualify for the
    discovery layer, which only recalls and never makes a family call.  An unrecognized
    declared value fails closed with an error instead of being read as either answer.
    """
    accession = str(row.get("accession") or "").strip() or "(unknown)"
    declared = str(row.get(ELIGIBILITY_COLUMN) or "").strip().lower()
    if declared:
        if declared in _TRUE_TOKENS:
            return True
        if declared in _FALSE_TOKENS:
            return False
        raise ValueError(
            f"{accession}: {ELIGIBILITY_COLUMN}={declared!r} is not a recognized flag; "
            "expected true or false"
        )
    return str(row.get("evidence_status") or "").strip() in DISCOVERY_TRAINING_ELIGIBLE_STATUSES


def _require_eligible_sequence_hash(row: dict[str, str]) -> None:
    if not str(row.get("sequence_sha256") or "").strip():
        accession = str(row.get("accession") or "").strip() or "(unknown)"
        raise ValueError(
            f"{accession}: discovery training eligibility requires a sequence_sha256"
        )


def select_discovery_training_rows(records: list[dict[str, str]]) -> list[dict[str, str]]:
    """Return the discovery-training rows in a deterministic, order-independent order.

    Selection is by eligibility only; the sort key ``(sequence_sha256, accession)`` makes
    the result independent of ledger row order.  A claimed-eligible row without a
    sequence SHA-256 is refused rather than silently dropped.
    """
    selected = [row for row in records if discovery_training_eligible(row)]
    for row in selected:
        _require_eligible_sequence_hash(row)
    return sorted(selected, key=lambda row: (row["sequence_sha256"], row["accession"]))


def unique_discovery_training_rows(records: list[dict[str, str]]) -> list[dict[str, str]]:
    """Collapse eligible rows sharing one sequence SHA-256 into one training sequence.

    Two ledger rows with the same hash are the *same* training sequence and must never be
    double-counted towards the three-sequence threshold.  The representative is the first
    row of :func:`select_discovery_training_rows` (lowest accession within that hash), so
    the outcome does not depend on ledger row order.
    """
    unique: dict[str, dict[str, str]] = {}
    for row in select_discovery_training_rows(records):
        unique.setdefault(row["sequence_sha256"], row)
    return [unique[digest] for digest in sorted(unique)]


def build_hmm(
    training_fasta: Path,
    alignment: Path,
    hmm: Path,
    mafft_bin: str,
    hmmbuild_bin: str,
    *,
    mafft_threads: int | None = None,
) -> None:
    command = [mafft_bin, "--auto"]
    if mafft_threads is not None:
        if mafft_threads < 1:
            raise ValueError("mafft_threads must be positive")
        command.extend(["--thread", str(mafft_threads)])
    command.append(str(training_fasta))
    result = _run(command)
    alignment.write_text(result.stdout, encoding="ascii", newline="\n")
    _run([hmmbuild_bin, str(hmm), str(alignment)])
    if not alignment.is_file() or not alignment.stat().st_size:
        raise RuntimeError(f"MAFFT did not create {alignment}")
    if not hmm.is_file() or not hmm.stat().st_size:
        raise RuntimeError(f"hmmbuild did not create {hmm}")


def _load_rows(ledger_path: Path, fasta_path: Path) -> tuple[list[dict[str, str]], dict[str, tuple[str, str]]]:
    _, rows = _read_tsv(ledger_path)
    fasta = _read_fasta(fasta_path)
    seen_accessions: set[str] = set()
    seen_refs: set[str] = set()
    seen_hashes: dict[str, str] = {}
    family_superfamilies: dict[str, set[str]] = defaultdict(set)
    for index, row in enumerate(rows, start=2):
        accession = (row.get("accession") or "").strip()
        reference_id = (row.get("reference_id") or "").strip()
        superfamily = (row.get("phaded_superfamily") or "").strip()
        family = (row.get("phaded_family_id") or "").strip()
        status = (row.get("evidence_status") or "").strip()
        if not accession or not reference_id or not family:
            raise ValueError(f"row {index}: accession, reference_id, and family are required")
        if accession in seen_accessions or reference_id in seen_refs:
            raise ValueError(f"duplicate accession or reference_id: {accession}")
        if superfamily not in ALLOWED_SUPERFAMILIES:
            raise ValueError(f"{accession}: invalid PhaDED superfamily")
        if status not in ALLOWED_EVIDENCE:
            raise ValueError(f"{accession}: invalid evidence_status")
        if status in {"experimental_positive", "experimental_negative"} and (
            not row.get("experimental_assay", "").strip() or not row.get("experimental_result", "").strip()
        ):
            raise ValueError(f"{accession}: experimental call lacks assay/result")
        if reference_id not in fasta:
            raise ValueError(f"{accession}: missing FASTA reference_id {reference_id}")
        fasta_accession, sequence = fasta[reference_id]
        if fasta_accession != accession:
            raise ValueError(f"{accession}: FASTA accession mismatch")
        digest = (row.get("sequence_sha256") or "").strip().lower()
        # Eligibility is a claim about a concrete sequence: it must carry its hash.
        if discovery_training_eligible(row):
            _require_eligible_sequence_hash(row)
        observed = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        if digest != observed:
            raise ValueError(f"{accession}: FASTA sequence SHA-256 mismatch")
        previous = seen_hashes.get(digest)
        if previous is not None:
            raise ValueError(f"duplicate sequence hash across profile partitions: {previous}, {accession}")
        seen_hashes[digest] = accession
        family_superfamilies[family].add(superfamily)
        seen_accessions.add(accession)
        seen_refs.add(reference_id)
    for family, superfamilies in family_superfamilies.items():
        if len(superfamilies) > 1:
            raise ValueError(f"mixed superfamily for family {family}: {sorted(superfamilies)}")
    if set(seen_refs) != set(fasta):
        missing = sorted(set(seen_refs) - set(fasta))
        extra = sorted(set(fasta) - set(seen_refs))
        raise ValueError(f"FASTA-to-ledger mismatch: missing={missing}, extra={extra}")
    return rows, fasta


def _training_audit_rows(
    profile_id: str,
    kind: str,
    model_layer: str,
    functional_calibration_status: str,
    records: list[dict[str, str]],
    selected_accessions: set[str],
) -> list[dict[str, str]]:
    """One audit row per ledger row per profile, with an explicit reason either way.

    A contradictory direct experimental negative is recorded here as blocking functional
    calibration instead of being dropped from the reference audit.
    """
    audit: list[dict[str, str]] = []
    for row in sorted(records, key=lambda item: item["accession"]):
        accession = row["accession"]
        eligible = discovery_training_eligible(row)
        if accession in selected_accessions:
            exclusion_reason = ""
        elif eligible:
            exclusion_reason = "duplicate_sequence_sha256_counts_as_one_training_sequence"
        elif str(row.get(ELIGIBILITY_COLUMN) or "").strip():
            exclusion_reason = "declared_discovery_training_eligible_false"
        else:
            exclusion_reason = "evidence_status_not_discovery_training_eligible"
        blocking = row["evidence_status"] == "experimental_negative"
        audit.append({
            "profile_id": profile_id,
            "profile_kind": kind,
            "model_layer": model_layer,
            "functional_calibration_status": functional_calibration_status,
            "accession": accession,
            "evidence_status": row["evidence_status"],
            "discovery_training_eligible": "true" if eligible else "false",
            "training_selected": "true" if accession in selected_accessions else "false",
            "exclusion_reason": exclusion_reason,
            "blocks_functional_calibration": "true" if blocking else "false",
            "block_reason": "contradictory_experimental_negative" if blocking else "",
        })
    return audit


def _profile_record(
    profile_id: str,
    kind: str,
    superfamily: str,
    family: str,
    records: list[dict[str, str]],
    fasta: dict[str, tuple[str, str]],
    profile_dir: Path,
    mafft_bin: str,
    hmmbuild_bin: str,
    software_versions: dict[str, str],
    mafft_threads: int | None = None,
) -> tuple[dict[str, object], list[dict[str, str]]]:
    """Build one profile record plus its per-row training audit.

    The discovery HMM is built when at least three *unique* eligible sequences exist.
    ``training_accessions`` is written in accession order, which is exactly the order of
    the training FASTA bytes bound by ``training_set_sha256``.
    """
    eligible = unique_discovery_training_rows(records)
    training_rows = sorted(eligible, key=lambda row: row["accession"])
    accessions = [row["accession"] for row in training_rows]
    experimental_anchor_count = sum(
        row["evidence_status"] == "experimental_positive" for row in training_rows
    )
    annotation_only_count = sum(
        row["evidence_status"] == "annotation_only" for row in training_rows
    )
    blocking_negatives = sorted(
        row["accession"] for row in records if row["evidence_status"] == "experimental_negative"
    )
    functional_calibration_status = (
        FUNCTIONAL_CALIBRATION_BLOCKED if blocking_negatives
        else FUNCTIONAL_CALIBRATION_NOT_CALIBRATED
    )
    if len(accessions) >= MINIMUM_DISCOVERY_TRAINING_SEQUENCES:
        status = "trained"
        model_layer = MODEL_LAYER_DISCOVERY
        reason = "at_least_three_unique_discovery_training_eligible_sequences"
    else:
        status = "reference_only"
        model_layer = MODEL_LAYER_REFERENCE_QUERY_ONLY
        reason = "fewer_than_three_unique_discovery_training_eligible_sequences"
    if blocking_negatives:
        reason = f"{reason};contradictory_experimental_negative_blocks_functional_calibration"
    alignment_path = profile_dir / f"{profile_id}.alignment.faa"
    hmm_path = profile_dir / f"{profile_id}.hmm"
    alignment_sha = hmm_sha = ""
    training_set_sha = "pending"
    mafft_version = software_versions.get("mafft", "pending")
    hmmer_version = software_versions.get("hmmer", "pending")
    if status == "trained":
        training_path = profile_dir / f"{profile_id}.training.faa"
        write_fasta(training_path, fasta, accessions)
        training_set_sha = sha256(training_path)
        build_hmm(
            training_path, alignment_path, hmm_path, mafft_bin, hmmbuild_bin,
            mafft_threads=mafft_threads,
        )
        alignment_sha, hmm_sha = sha256(alignment_path), sha256(hmm_path)
        if mafft_version == "pending":
            mafft_version = tool_version(mafft_bin, kind="mafft")
        if hmmer_version == "pending":
            hmmer_version = tool_version(hmmbuild_bin, kind="hmmer")
    else:
        alignment_path = hmm_path = Path("")
    record: dict[str, object] = {
        "profile_id": profile_id, "profile_kind": kind,
        "phaded_superfamily": superfamily, "phaded_family_id": family,
        "training_accessions": ";".join(accessions), "training_count": len(accessions),
        "model_status": status, "model_reason": reason,
        "model_layer": model_layer,
        "functional_calibration_status": functional_calibration_status,
        "training_sequence_count": len(accessions),
        "experimental_anchor_count": experimental_anchor_count,
        "annotation_only_count": annotation_only_count,
        "training_set_sha256": training_set_sha,
        "alignment_sha256": alignment_sha, "hmm_sha256": hmm_sha,
        "mafft_version": mafft_version, "hmmer_version": hmmer_version,
        # MAFFT 7.525 is not bit-reproducible: never claim it, the alignment hash is authoritative.
        "bit_reproducible": "false",
        "alignment_path": str(alignment_path) if status == "trained" else "",
        "hmm_path": str(hmm_path) if status == "trained" else "",
    }
    audit = _training_audit_rows(
        profile_id, kind, model_layer, functional_calibration_status, records, set(accessions)
    )
    return record, audit


def build_profiles(
    ledger_path: str | Path,
    fasta_path: str | Path,
    output_dir: str | Path,
    *,
    mafft_bin: str = "mafft",
    hmmbuild_bin: str = "hmmbuild",
    software_versions: dict[str, str] | None = None,
    expected_family_count: int | None = None,
    mafft_threads: int | None = None,
) -> dict[str, object]:
    ledger_path, fasta_path, output_dir = map(Path, (ledger_path, fasta_path, output_dir))
    rows, fasta = _load_rows(ledger_path, fasta_path)
    family_keys = {(row["phaded_family_id"], row["phaded_superfamily"]) for row in rows}
    if expected_family_count is not None and len({family for family, _ in family_keys}) != expected_family_count:
        raise ValueError(
            f"ledger contains {len({family for family, _ in family_keys})} families; "
            f"expected {expected_family_count}"
        )
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    profile_dir = output_dir / "profiles"
    profile_dir.mkdir()
    versions = software_versions or {}
    grouped_family: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    grouped_superfamily: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped_family[(row["phaded_family_id"], row["phaded_superfamily"])].append(row)
        grouped_superfamily[row["phaded_superfamily"]].append(row)
    profiles: list[dict[str, object]] = []
    training_audit: list[dict[str, str]] = []
    for (family, superfamily), group in sorted(grouped_family.items()):
        record, audit = _profile_record(
            _profile_id("family", family), "family", superfamily, family, group,
            fasta, profile_dir, mafft_bin, hmmbuild_bin, versions,
            mafft_threads,
        )
        profiles.append(record)
        training_audit.extend(audit)
    for superfamily, group in sorted(grouped_superfamily.items()):
        record, audit = _profile_record(
            _profile_id("superfamily", superfamily), "superfamily", superfamily, "", group,
            fasta, profile_dir, mafft_bin, hmmbuild_bin, versions,
            mafft_threads,
        )
        profiles.append(record)
        training_audit.extend(audit)
    manifest_path = output_dir / "profile_manifest.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(profiles)
    audit_path = output_dir / "profile_training_audit.tsv"
    with audit_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TRAINING_AUDIT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(training_audit)
    family_profiles = [profile for profile in profiles if profile["profile_kind"] == "family"]
    result: dict[str, object] = {
        "profile_manifest": str(manifest_path), "profiles": profiles,
        "profile_count": len(profiles), "ledger_sha256": sha256(ledger_path),
        "reference_fasta_sha256": sha256(fasta_path),
        "training_audit_path": str(audit_path), "training_audit": training_audit,
    }
    if len(family_profiles) == 1:
        result["family_status"] = family_profiles[0]["model_status"]
        result["hmm_path"] = family_profiles[0]["hmm_path"]
        result["family_model_layer"] = family_profiles[0]["model_layer"]
        result["family_functional_calibration_status"] = family_profiles[0]["functional_calibration_status"]
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--mafft", default="mafft")
    parser.add_argument("--hmmbuild", default="hmmbuild")
    parser.add_argument("--expected-family-count", type=int, default=38)
    parser.add_argument(
        "--mafft-threads", type=int, default=None,
        help="maximum MAFFT threads per alignment (set explicitly on shared servers)",
    )
    args = parser.parse_args(argv)
    result = build_profiles(
        args.ledger, args.fasta, args.output_dir,
        mafft_bin=args.mafft, hmmbuild_bin=args.hmmbuild,
        expected_family_count=args.expected_family_count,
        mafft_threads=args.mafft_threads,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
