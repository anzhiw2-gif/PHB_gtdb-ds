#!/usr/bin/env python3
"""Build independent ePhaZ validation and structural-neighbor control panels.

The builder consumes browser-captured UniProt/NCBI FASTA responses and a human
curation manifest. It never infers experimental function from an annotation;
the manifest must carry the evidence decision for every accession.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

MANIFEST_FIELDS = (
    "accession",
    "panel",
    "substrate",
    "organism",
    "length",
    "reviewed",
    "pmid",
    "doi",
    "response_kind",
    "response_file",
    "evidence_summary",
)

ALLOWED_PANELS = {
    "independent_experimental_positive",
    "mcl_pha_experimental_positive",
    "mcl_pha_non_phb_negative",
    "intracellular_non_ephaz_negative",
    "annotation_only_near_neighbor_negative",
    "fragment_or_incomplete_negative",
}

PROTECTED_ACCESSIONS = {
    # bridge training proteins
    "Q51871",
    "Q5SLU4",
    # existing ePhaZ core/seed validation proteins
    "B2NHN2",
    "O05527",
    "P12625",
    "Q51718",
}

PANEL_FILES = {
    "independent_experimental_positive": "independent_experimental_positive.faa",
    "mcl_pha_experimental_positive": "mcl_pha_experimental_positive.faa",
    "mcl_pha_non_phb_negative": "mcl_pha_non_phb_negative.faa",
    "intracellular_non_ephaz_negative": "intracellular_non_ephaz_negative.faa",
    "annotation_only_near_neighbor_negative": "annotation_only_near_neighbor_negative.faa",
    "fragment_or_incomplete_negative": "fragment_or_incomplete_negative.faa",
}


class PanelBuildError(ValueError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_response(path: Path) -> str:
    if not path.is_file():
        raise PanelBuildError(f"missing raw response: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
        value = payload["value"]
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise PanelBuildError(f"invalid raw response JSON: {path}") from exc
    if not isinstance(value, str) or not value.strip():
        raise PanelBuildError(f"empty raw response: {path}")
    return value


def _bound_input_path(inputs: Path, name: str) -> Path:
    candidate = Path(name)
    if candidate.is_absolute() or ".." in candidate.parts or candidate.name != name:
        raise PanelBuildError(f"response path must be a basename inside inputs: {name!r}")
    resolved = (inputs / candidate).resolve()
    if resolved.parent != inputs.resolve():
        raise PanelBuildError(f"response path escapes inputs: {name!r}")
    return resolved


def _find_repo_root(run_dir: Path) -> Path:
    """Locate the checkout root from a run directory, including test fixtures."""
    for candidate in (run_dir, *run_dir.parents):
        if (candidate / "pipeline" / "seeds" / "controls" / "positive.faa").is_file():
            return candidate
    if run_dir.parent.name == "runs":
        return run_dir.parent.parent
    return run_dir.parent


def _header_accession(header: str) -> str:
    fields = header[1:].strip().split("|")
    if fields and fields[0].lower() in {"sp", "tr", "gi", "ref"} and len(fields) >= 2:
        return fields[1].split()[0]
    return fields[0].split()[0]


def _parse_fasta_text(text: str) -> Tuple[str, str]:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not any(line.startswith(">") for line in lines):
        # NCBI GenBank pages expose the translated CDS in a wrapped qualifier.
        match = re.search(r'/protein_id="([^"]+)".*?/translation="([A-Za-z\s]+)"', text, re.DOTALL)
        if match:
            accession = match.group(1)
            sequence = re.sub(r"\s+", "", match.group(2)).upper()
            if sequence and re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYBXZJUO]+", sequence):
                return accession, sequence
    try:
        header_index = next(i for i, line in enumerate(lines) if line.startswith(">"))
    except StopIteration as exc:
        raise PanelBuildError("raw response contains no FASTA record") from exc
    header = lines[header_index]
    if sum(1 for line in lines[header_index:] if line.startswith(">")) != 1:
        raise PanelBuildError("raw response contains multiple FASTA records")
    sequence: List[str] = []
    for line in lines[header_index + 1 :]:
        if line.startswith(">"):
            break
        if re.fullmatch(r"[A-Za-z*.-]+", line):
            sequence.append(line.replace("-", "").replace(".", "").upper())
    seq = "".join(sequence).rstrip("*")
    if not seq or not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYBXZJUO]+", seq):
        raise PanelBuildError(f"invalid FASTA sequence for {header}")
    accession = _header_accession(header)
    return accession, seq


def _parse_manifest(path: Path) -> List[dict]:
    if not path.is_file():
        raise PanelBuildError(f"missing candidate manifest: {path}")
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    if not lines:
        raise PanelBuildError("candidate manifest is empty")
    header = lines[0].split("\t")
    if tuple(header) != MANIFEST_FIELDS:
        raise PanelBuildError("candidate manifest header does not match required schema")
    rows: List[dict] = []
    seen = set()
    for line_no, line in enumerate(lines[1:], 2):
        if not line.strip():
            continue
        values = line.split("\t")
        if len(values) != len(MANIFEST_FIELDS):
            raise PanelBuildError(f"manifest line {line_no} has wrong field count")
        row = dict(zip(MANIFEST_FIELDS, values))
        acc = row["accession"].strip()
        panel = row["panel"].strip()
        if not acc or acc in seen:
            raise PanelBuildError(f"missing or duplicate accession at line {line_no}: {acc}")
        if acc in PROTECTED_ACCESSIONS:
            raise PanelBuildError(f"protected accession cannot enter independent panels: {acc}")
        if panel not in ALLOWED_PANELS:
            raise PanelBuildError(f"unsupported panel {panel!r} at line {line_no}")
        if not row["response_file"].strip() or not row["response_kind"].strip():
            raise PanelBuildError(f"missing response binding at line {line_no}")
        if panel in {"independent_experimental_positive", "mcl_pha_experimental_positive"} and not row["pmid"].strip():
            raise PanelBuildError(f"experimental positive lacks PMID at line {line_no}")
        seen.add(acc)
        rows.append(row)
    if not rows:
        raise PanelBuildError("candidate manifest has no data rows")
    return rows


def _write_fasta(path: Path, records: Iterable[Tuple[dict, str]]) -> int:
    count = 0
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row, sequence in records:
            handle.write(f">{row['accession']}|{row['panel']}|{row['organism']}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset : offset + 80] + "\n")
            count += 1
    return count


def _sequence_hashes(path: Path) -> Dict[str, str]:
    """Return accession -> normalized sequence hash for a local FASTA."""
    if not path.is_file():
        return {}
    records: Dict[str, str] = {}
    accession = None
    chunks: List[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if accession is not None:
                if accession in records:
                    raise PanelBuildError(f"duplicate accession in protected FASTA: {accession}")
                records[accession] = hashlib.sha256("".join(chunks).encode()).hexdigest()
            accession = _header_accession(line)
            chunks = []
        elif line.strip():
            chunks.append(line.strip().replace("-", "").replace(".", "").upper())
    if accession is not None:
        if accession in records:
            raise PanelBuildError(f"duplicate accession in protected FASTA: {accession}")
        records[accession] = hashlib.sha256("".join(chunks).encode()).hexdigest()
    return records


def build_panels(run_dir: Path | str, manifest_path: Path | str) -> dict:
    run_dir = Path(run_dir)
    manifest_path = Path(manifest_path)
    inputs = run_dir / "inputs"
    results = run_dir / "results"
    logs = run_dir / "logs"
    for directory in (inputs, results, logs):
        directory.mkdir(parents=True, exist_ok=True)
    rows = _parse_manifest(manifest_path)
    parsed: List[Tuple[dict, str]] = []
    seen = set()
    repo_root = _find_repo_root(run_dir)
    protected_paths = [
        repo_root / "runs" / "20260830_ephaz_bridge_curation_02" / "results" / "ephaz_bridge_candidate.faa",
        repo_root / "runs" / "20260830_ephaz_bridge_curation_02" / "inputs" / "ephaz_curated_core.faa",
        repo_root / "pipeline" / "seeds" / "controls" / "positive.faa",
    ]
    protected_hashes: Dict[str, str] = {}
    protected_sources = []
    for protected_path in protected_paths:
        if not protected_path.is_file() or protected_path.stat().st_size == 0:
            raise PanelBuildError(f"missing protected independence source: {protected_path}")
        hashes = _sequence_hashes(protected_path)
        if not hashes:
            raise PanelBuildError(f"protected independence source is not valid FASTA: {protected_path}")
        protected_sources.append({
            "path": str(protected_path),
            "bytes": protected_path.stat().st_size if protected_path.is_file() else None,
            "sha256": sha256_file(protected_path) if protected_path.is_file() else "pending",
            "accession_count": len(hashes),
        })
        protected_hashes.update(hashes)
    for row in rows:
        response_path = _bound_input_path(inputs, row["response_file"])
        response_text = _read_response(response_path)
        found_acc, sequence = _parse_fasta_text(response_text)
        expected = row["accession"].strip()
        if found_acc != expected:
            raise PanelBuildError(f"response accession mismatch: expected {expected}, found {found_acc}")
        declared_length = row["length"].strip()
        if declared_length and declared_length.isdigit() and int(declared_length) != len(sequence):
            raise PanelBuildError(f"length mismatch for {expected}: manifest {declared_length}, FASTA {len(sequence)}")
        if expected in seen:
            raise PanelBuildError(f"duplicate parsed accession: {expected}")
        sequence_hash = hashlib.sha256(sequence.encode()).hexdigest()
        if sequence_hash in protected_hashes.values():
            matching = sorted(acc for acc, digest in protected_hashes.items() if digest == sequence_hash)
            raise PanelBuildError(f"sequence-identical to protected training/core accession(s): {','.join(matching)}")
        seen.add(expected)
        parsed.append((row, sequence))

    # Keep excluded but informative sequences in a non-denominator challenge panel.
    exclusion_path = inputs / "excluded_candidates.tsv"
    challenge_records: List[Tuple[dict, str]] = []
    excluded_accessions = set()
    if exclusion_path.is_file():
        exclusion_lines = exclusion_path.read_text(encoding="utf-8").splitlines()
        for line in exclusion_lines[1:]:
            if not line.strip():
                continue
            values = line.split("\t")
            if len(values) != 4:
                raise PanelBuildError("excluded candidate manifest has wrong field count")
            accession, exclusion_class, reason, reference = values
            if not accession or not exclusion_class or not reason or not reference:
                raise PanelBuildError("excluded candidate manifest contains an empty field")
            if accession in seen:
                raise PanelBuildError(f"excluded accession overlaps candidate manifest: {accession}")
            response_file = reference.split(";", 1)[0].strip()
            response_path = _bound_input_path(inputs, response_file)
            found_acc, sequence = _parse_fasta_text(_read_response(response_path))
            if found_acc != accession:
                raise PanelBuildError(f"excluded response accession mismatch: {accession} vs {found_acc}")
            if accession in excluded_accessions:
                raise PanelBuildError(f"duplicate excluded accession: {accession}")
            excluded_accessions.add(accession)
            challenge_records.append((
                {
                    "accession": accession,
                    "panel": f"excluded_challenge_{exclusion_class}",
                    "organism": "excluded candidate",
                },
                sequence,
            ))

    counts = Counter(row["panel"] for row, _ in parsed)
    for panel, filename in PANEL_FILES.items():
        _write_fasta(results / filename, ((row, seq) for row, seq in parsed if row["panel"] == panel))
    challenge_path = results / "ephaz_excluded_challenge.faa"
    _write_fasta(challenge_path, challenge_records)

    ledger_path = results / "ephaz_panel_evidence.tsv"
    ledger_fields = list(MANIFEST_FIELDS) + ["sequence_length", "normalized_sequence_sha256", "source_sha256", "independence_status", "decision"]
    with ledger_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(ledger_fields) + "\n")
        for row, sequence in parsed:
            response_path = _bound_input_path(inputs, row["response_file"])
            values = [row[field] for field in MANIFEST_FIELDS]
            values += [
                str(len(sequence)),
                hashlib.sha256(sequence.encode()).hexdigest(),
                sha256_file(response_path),
                "independent_of_bridge_and_core",
                "accept_panel",
            ]
            handle.write("\t".join(values) + "\n")

    raw_index_path = results / "raw_response_index.tsv"
    with raw_index_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("accession\tresponse_file\tresponse_kind\tbytes\tsha256\treason\n")
        indexed = set()
        for row, _ in parsed:
            response_path = _bound_input_path(inputs, row["response_file"])
            indexed.add(response_path.name)
            handle.write(
                f"{row['accession']}\t{row['response_file']}\t{row['response_kind']}\t"
                f"{response_path.stat().st_size}\t{sha256_file(response_path)}\tmanifest_panel\n"
            )
        for record, _ in challenge_records:
            response_file = next(
                value.split(";", 1)[0].strip()
                for value in [line.split("\t")[3] for line in exclusion_lines[1:] if line.strip() and line.split("\t")[0] == record["accession"]]
            )
            response_path = _bound_input_path(inputs, response_file)
            handle.write(
                f"{record['accession']}\t{response_file}\texcluded_challenge\t"
                f"{response_path.stat().st_size}\t{sha256_file(response_path)}\texcluded_challenge\n"
            )
        for response_path in sorted(inputs.glob("*_eval.json")):
            if response_path.name in indexed:
                continue
            handle.write(
                f"\t{response_path.name}\tliterature_or_database_evidence\t"
                f"{response_path.stat().st_size}\t{sha256_file(response_path)}\n"
            )

    outputs = {}
    for path in sorted(results.iterdir()):
        if path.is_file():
            outputs[path.name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    count_payload = {panel: counts.get(panel, 0) for panel in PANEL_FILES}
    count_payload["expanded_ephaz_near_neighbor_negative"] = sum(
        count_payload[panel]
        for panel in (
            "intracellular_non_ephaz_negative",
            "mcl_pha_non_phb_negative",
            "annotation_only_near_neighbor_negative",
        )
    )
    exclusions = {
        "path": str(exclusion_path),
        "bytes": exclusion_path.stat().st_size if exclusion_path.is_file() else None,
        "sha256": sha256_file(exclusion_path) if exclusion_path.is_file() else "pending",
    }
    contract = {
        "run_id": run_dir.name,
        "purpose": "independent experimental ePhaZ validation and expanded structural-neighbor controls",
        "frozen_hmm_unchanged": True,
        "full_scan_authorized": False,
        "protected_accessions": sorted(PROTECTED_ACCESSIONS),
        "protected_sequence_sources": protected_sources,
        "candidate_manifest": {
            "path": str(manifest_path),
            "bytes": manifest_path.stat().st_size,
            "sha256": sha256_file(manifest_path),
        },
        "exclusion_manifest": exclusions,
        "counts": count_payload,
        "formal_negative_denominator": {
            "intracellular_non_ephaz_negative": count_payload["intracellular_non_ephaz_negative"],
            "mcl_pha_non_phb_negative": count_payload["mcl_pha_non_phb_negative"],
            "annotation_only_near_neighbor_negative": count_payload["annotation_only_near_neighbor_negative"],
            "total": sum(count_payload[key] for key in ("intracellular_non_ephaz_negative", "mcl_pha_non_phb_negative", "annotation_only_near_neighbor_negative")),
        },
        "challenge_only_count": count_payload["fragment_or_incomplete_negative"] + len(challenge_records),
        "outputs": outputs,
    }
    (run_dir / "input_contract.json").write_text(json.dumps(contract, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return {"counts": contract["counts"], "outputs": outputs, "protected_accessions": contract["protected_accessions"]}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    args = parser.parse_args()
    result = build_panels(args.run_dir, args.manifest)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
