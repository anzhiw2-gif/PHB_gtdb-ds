#!/usr/bin/env python3
"""Fail-closed validation for an experimental ePhaZ bridge evidence ledger."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Sequence


REQUIRED_COLUMNS = {
    "candidate_id", "accession", "protein_name", "organism", "sequence_source",
    "sequence_url", "sequence_length", "primary_reference", "pmid", "doi",
    "experimental_evidence_type", "substrate", "extracellular_evidence",
    "architecture_evidence", "completeness_status", "core_overlap_status",
    "iPhaZ_challenge_status", "decision", "decision_reason", "review_date", "notes",
}
ACCEPT = "accept_bridge_candidate"
EXPERIMENTAL_TOKENS = ("purified", "recombinant", "cloned_gene", "biochemical", "characterization")


class LedgerError(ValueError):
    """Raised when an accepted bridge row lacks its required evidence."""


def _required(row: dict[str, str], name: str) -> str:
    value = (row.get(name) or "").strip()
    if not value:
        raise LedgerError(f"{row.get('candidate_id', '<unknown>')}: missing {name}")
    return value


def validate_ledger(path: str | Path) -> list[str]:
    """Validate accepted rows and return their accessions in ledger order."""
    path = Path(path)
    if not path.is_file() or path.is_symlink():
        raise LedgerError(f"ledger is not a regular file: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        missing = sorted(REQUIRED_COLUMNS - fields)
        if missing:
            raise LedgerError("ledger missing required columns: " + ",".join(missing))
        rows = list(reader)
    seen: set[str] = set()
    accepted: list[str] = []
    for row in rows:
        accession = _required(row, "accession")
        if accession in seen:
            raise LedgerError(f"duplicate accession: {accession}")
        seen.add(accession)
        if (row.get("decision") or "").strip() != ACCEPT:
            continue
        for name in (
            "candidate_id", "protein_name", "organism", "sequence_source", "sequence_url",
            "primary_reference", "pmid", "doi", "experimental_evidence_type", "substrate",
            "extracellular_evidence", "architecture_evidence", "completeness_status",
            "core_overlap_status", "iPhaZ_challenge_status", "decision_reason", "review_date",
        ):
            _required(row, name)
        try:
            length = int(_required(row, "sequence_length"))
        except ValueError as exc:
            raise LedgerError(f"{accession}: sequence_length must be an integer") from exc
        if length < 200:
            raise LedgerError(f"{accession}: accepted bridge sequence is shorter than 200 aa")
        evidence = row["experimental_evidence_type"].lower()
        if not any(token in evidence for token in EXPERIMENTAL_TOKENS):
            raise LedgerError(f"{accession}: accepted row lacks direct experimental enzyme evidence")
        if row["completeness_status"].strip().lower() != "complete":
            raise LedgerError(f"{accession}: accepted row is not complete")
        if row["iPhaZ_challenge_status"].strip().lower() != "not_challenge":
            raise LedgerError(f"{accession}: challenge sequence cannot enter bridge")
        accepted.append(accession)
    if not accepted:
        raise LedgerError("ledger contains no accepted bridge candidates")
    return accepted


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ledger")
    args = parser.parse_args(argv)
    accepted = validate_ledger(args.ledger)
    print(f"validated {len(accepted)} accepted ePhaZ bridge candidates")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
