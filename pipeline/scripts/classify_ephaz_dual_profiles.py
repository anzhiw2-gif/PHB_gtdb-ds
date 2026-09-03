#!/usr/bin/env python3
"""Classify sequences with independent PHB and MCL HMM profiles."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path


def parse_fasta_accession(header: str) -> str:
    """Return the stable accession from a FASTA header or HMMER target name."""
    token = header.lstrip("> ").split()[0]
    fields = token.split("|")
    if len(fields) >= 2 and fields[0].lower() in {"sp", "tr", "ref", "gi"}:
        return fields[1]
    return fields[0]


def parse_tblout(text: str) -> dict[str, dict[str, float]]:
    hits: dict[str, dict[str, float]] = {}
    for raw in text.splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        fields = raw.split()
        if len(fields) < 6:
            continue
        try:
            accession = parse_fasta_accession(fields[0])
            evalue, bitscore = float(fields[4]), float(fields[5])
        except (ValueError, IndexError):
            continue
        if accession not in hits or bitscore > hits[accession]["bitscore"]:
            hits[accession] = {"evalue": evalue, "bitscore": bitscore}
    return hits


def parse_domtblout(text: str) -> dict[str, dict[str, float]]:
    """Parse HMMER domain output and merge non-overlapping HMM intervals."""
    grouped: dict[str, dict[str, object]] = {}
    for raw in text.splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        fields = raw.split()
        if len(fields) < 22:
            continue
        try:
            accession = parse_fasta_accession(fields[0])
            qlen = int(fields[5])
            evalue = float(fields[12])
            bitscore = float(fields[13])
            start, end = int(fields[15]), int(fields[16])
        except (ValueError, IndexError):
            continue
        if qlen <= 0 or end < start:
            continue
        item = grouped.setdefault(accession, {"evalue": evalue, "bitscore": bitscore, "qlen": qlen, "intervals": []})
        item["evalue"] = min(float(item["evalue"]), evalue)
        item["bitscore"] = max(float(item["bitscore"]), bitscore)
        item["intervals"].append((start, end))
    hits: dict[str, dict[str, float]] = {}
    for accession, item in grouped.items():
        intervals = sorted(item["intervals"])
        covered = 0
        current_start = current_end = None
        for start, end in intervals:
            if current_start is None:
                current_start, current_end = start, end
            elif start > current_end + 1:
                covered += current_end - current_start + 1
                current_start, current_end = start, end
            else:
                current_end = max(current_end, end)
        if current_start is not None:
            covered += current_end - current_start + 1
        hits[accession] = {
            "evalue": float(item["evalue"]),
            "bitscore": float(item["bitscore"]),
            "coverage": min(1.0, covered / int(item["qlen"])),
        }
    return hits


def _accepted(hit: dict[str, float] | None, threshold: float, min_domain_coverage: float | None = None) -> bool:
    if hit is None or hit.get("evalue", float("inf")) > threshold:
        return False
    return min_domain_coverage is None or hit.get("coverage", 0.0) >= min_domain_coverage


def classify_accession(phb_hit: dict[str, float] | None, mcl_hit: dict[str, float] | None, threshold: float = 1e-5, margin: float = 10.0, min_domain_coverage: float | None = None) -> str:
    phb_ok = _accepted(phb_hit, threshold, min_domain_coverage)
    mcl_ok = _accepted(mcl_hit, threshold, min_domain_coverage)
    if phb_ok and not mcl_ok:
        return "PHB_like"
    if mcl_ok and not phb_ok:
        return "MCL_like"
    if not phb_ok and not mcl_ok:
        return "no_hit"
    phb_score = phb_hit["bitscore"] if phb_hit else float("-inf")
    mcl_score = mcl_hit["bitscore"] if mcl_hit else float("-inf")
    if abs(phb_score - mcl_score) < margin:
        return "ambiguous"
    return "PHB_like" if phb_score > mcl_score else "MCL_like"


def _is_hit_map(value: object) -> bool:
    return isinstance(value, dict) and all(
        isinstance(item, dict) and {"evalue", "bitscore"}.issubset(item)
        for item in value.values()
    )


def _profile_label(name: str) -> str:
    if name.startswith("MCL_"):
        return name
    if name.startswith("mcl_"):
        return f"MCL_{name[4:]}"
    return f"MCL_{name}"


def _best_mcl_profile(
    profiles: dict[str, dict[str, dict[str, float]]], accession: str,
    threshold: float = 1e-5, min_domain_coverage: float | None = None,
) -> tuple[str | None, dict[str, float] | None]:
    candidates = [
        (name, hits[accession])
        for name, hits in profiles.items()
        if accession in hits
    ]
    accepted = [item for item in candidates if _accepted(item[1], threshold, min_domain_coverage)]
    if accepted:
        candidates = accepted
    if not candidates:
        return None, None
    return min(candidates, key=lambda item: (item[1]["evalue"], -item[1]["bitscore"], item[0]))


def classify_hits(
    phb_hits: dict[str, dict[str, float]],
    mcl_hits: dict[str, dict[str, float]] | dict[str, dict[str, dict[str, float]]],
    threshold: float = 1e-5,
    margin: float = 10.0,
    accessions: list[str] | None = None,
    min_domain_coverage: float | None = None,
) -> list[dict[str, object]]:
    """Classify hits from one PHB profile and one or more named MCL profiles."""
    profiles = mcl_hits if not _is_hit_map(mcl_hits) else {"mcl": mcl_hits}
    rows = []
    universe = set(accessions) if accessions is not None else set(phb_hits) | {item for hits in profiles.values() for item in hits}
    for accession in sorted(universe):
        phb = phb_hits.get(accession)
        mcl_name, mcl = _best_mcl_profile(profiles, accession, threshold, min_domain_coverage)
        base = classify_accession(phb, mcl, threshold, margin, min_domain_coverage)
        classification = _profile_label(mcl_name) if mcl_name and mcl_name != "mcl" and base == "MCL_like" else base
        row: dict[str, object] = {
            "accession": accession,
            "classification": classification,
            "phb_evalue": phb.get("evalue", "") if phb else "",
            "phb_bitscore": phb.get("bitscore", "") if phb else "",
            "mcl_evalue": mcl.get("evalue", "") if mcl else "",
            "mcl_bitscore": mcl.get("bitscore", "") if mcl else "",
            "mcl_coverage": mcl.get("coverage", "") if mcl else "",
            "best_mcl_subfamily": mcl_name or "",
        }
        for name, hits in profiles.items():
            hit = hits.get(accession)
            slug = name.lower().replace("-", "_")
            row[f"mcl_{slug}_evalue"] = hit.get("evalue", "") if hit else ""
            row[f"mcl_{slug}_bitscore"] = hit.get("bitscore", "") if hit else ""
        rows.append(row)
    return rows


def _run(hmm: Path, fasta: Path, work: Path, executable: str, cpu: int, label: str | None = None) -> dict[str, dict[str, float]]:
    tbl = work / f"{label or hmm.stem}.tbl"
    dom = work / f"{label or hmm.stem}.domtbl"
    subprocess.run([executable, "--tblout", str(tbl), "--domtblout", str(dom), "-E", "1e-2", "--cpu", str(cpu), str(hmm), str(fasta)], check=True, capture_output=True, text=True)
    hits = parse_tblout(tbl.read_text(encoding="utf-8"))
    for accession, domain in parse_domtblout(dom.read_text(encoding="utf-8")).items():
        if accession in hits:
            hits[accession]["coverage"] = domain["coverage"]
    return hits


def classify_fasta(
    phb_hmm: Path | str,
    mcl_hmm: Path | str | None,
    fasta: Path | str,
    out: Path | str,
    threshold: float = 1e-5,
    margin: float = 10.0,
    cpu: int = 4,
    hmmsearch_bin: str | None = None,
    mcl_profiles: dict[str, Path | str] | None = None,
    min_domain_coverage: float | None = None,
) -> list[dict[str, object]]:
    executable = hmmsearch_bin or shutil.which("hmmsearch")
    if not executable:
        raise ValueError("hmmsearch executable not found")
    if mcl_hmm is None and not mcl_profiles:
        raise ValueError("at least one MCL HMM profile is required")
    if cpu < 1 or not 0 < threshold or (min_domain_coverage is not None and not 0 <= min_domain_coverage <= 1):
        raise ValueError("cpu must be positive and threshold must be > 0")
    output = Path(out)
    output.parent.mkdir(parents=True, exist_ok=True)
    work = output.parent / f".{output.stem}_work"
    work.mkdir(exist_ok=False)
    try:
        phb_hits = _run(Path(phb_hmm), Path(fasta), work, executable, cpu)
        mcl_hits = {}
        if mcl_hmm is not None:
            mcl_hits["mcl"] = _run(Path(mcl_hmm), Path(fasta), work, executable, cpu, "mcl")
        for name, profile in (mcl_profiles or {}).items():
            mcl_hits[name] = _run(Path(profile), Path(fasta), work, executable, cpu, name)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    accessions = []
    for line in Path(fasta).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            accessions.append(parse_fasta_accession(line[1:]))
    rows = classify_hits(phb_hits, mcl_hits, threshold, margin, accessions, min_domain_coverage)
    with output.open("w", encoding="utf-8", newline="") as handle:
        fields = ["accession", "classification", "phb_evalue", "phb_bitscore", "mcl_evalue", "mcl_bitscore", "mcl_coverage", "best_mcl_subfamily"]
        fields.extend(sorted({key for row in rows for key in row if key.startswith("mcl_") and key not in fields}))
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    return rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phb-hmm", required=True, type=Path)
    parser.add_argument("--mcl-hmm", type=Path)
    parser.add_argument("--mcl-profile", action="append", default=[], metavar="NAME=PATH")
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--threshold", type=float, default=1e-5)
    parser.add_argument("--margin", type=float, default=10.0)
    parser.add_argument("--cpu", type=int, default=4)
    parser.add_argument("--hmmsearch-bin")
    parser.add_argument("--min-domain-coverage", type=float)
    args = parser.parse_args(argv)
    profiles = {}
    for item in args.mcl_profile:
        if "=" not in item:
            parser.error("--mcl-profile must use NAME=PATH")
        name, path = item.split("=", 1)
        if not name or not path or name in profiles:
            parser.error("--mcl-profile names must be non-empty and unique")
        profiles[name] = Path(path)
    classify_fasta(args.phb_hmm, args.mcl_hmm, args.fasta, args.out, args.threshold, args.margin, args.cpu, args.hmmsearch_bin, profiles, args.min_domain_coverage)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
