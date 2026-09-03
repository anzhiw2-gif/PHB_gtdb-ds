#!/usr/bin/env python3
"""Stratify and deterministically sample ePhaZ/iPhaZ ambiguous candidates."""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


REQUIRED_COMPETITION = {"accession", "length", "ephaz_bitscore", "iphaz_bitscore", "assignment"}


def read_competition(path: Path | str) -> list[dict[str, str]]:
    path = Path(path)
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = REQUIRED_COMPETITION - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"competition TSV missing columns: {sorted(missing)}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    if not rows:
        raise ValueError("competition TSV is empty")
    if any(row["assignment"] != "ambiguous" for row in rows):
        raise ValueError("input competition TSV must contain only ambiguous rows")
    if len({row["accession"] for row in rows}) != len(rows):
        raise ValueError("competition TSV contains duplicate accessions")
    return rows


def parse_domtblout(paths: list[Path | str]) -> dict[str, float]:
    """Return merged HMM-coordinate coverage per full target identifier."""
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    lengths: dict[str, float] = {}
    for path in paths:
        with Path(path).open(encoding="utf-8") as handle:
            for raw in handle:
                if not raw.strip() or raw.lstrip().startswith("#"):
                    continue
                fields = raw.split()
                if len(fields) < 21:
                    continue
                try:
                    target = fields[0]
                    qlen = float(fields[5])
                    # HMMER domtblout: HMM coordinates are columns 16-17
                    # (zero-based 15-16); envelope coordinates (19-20) are
                    # broader and would overstate model coverage.
                    start, end = int(fields[15]), int(fields[16])
                    if start < 1 or end < start or qlen <= 0:
                        continue
                except (IndexError, ValueError, ZeroDivisionError):
                    continue
                intervals[target].append((start, end))
                lengths[target] = qlen
    coverage: dict[str, float] = {}
    for target, spans in intervals.items():
        merged: list[list[int]] = []
        for start, end in sorted(spans):
            if merged and start <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        covered = sum(end - start + 1 for start, end in merged)
        coverage[target] = min(1.0, max(0.0, covered / lengths[target]))
    return coverage


def read_neighborhood(path: Path | str | None) -> dict[str, str]:
    """Map full candidate IDs to marker-present or record-no-marker."""
    if path is None:
        return {}
    result: dict[str, str] = {}
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        if not {"genome", "hit_locus"}.issubset(fields):
            raise ValueError("neighborhood TSV requires genome and hit_locus columns")
        marker_field = "nearby_markers" if "nearby_markers" in fields else "marker_family" if "marker_family" in fields else None
        for row in reader:
            genome = (row.get("genome") or "").strip()
            locus = (row.get("hit_locus") or "").strip()
            if not genome or not locus:
                continue
            key = locus if "|" in locus else f"{genome}|{locus}"
            marker = (row.get(marker_field) or "").strip() if marker_field else ""
            current = result.get(key)
            value = "marker_present" if marker else "record_no_marker"
            if current != "marker_present":
                result[key] = value
    return result


def _delta_bin(value: float) -> str:
    for lower, upper in ((0, 2), (2, 5), (5, 10), (10, 15), (15, 20)):
        if lower <= value < upper:
            return f"{lower}-<{upper}"
    return ">=20"


def _length_bin(value: int) -> str:
    if value < 250:
        return "<250"
    if value < 400:
        return "250-399"
    if value < 600:
        return "400-599"
    return ">=600"


def _coverage_bin(ephaz: float | None, iphaz: float | None) -> str:
    values = [value for value in (ephaz, iphaz) if value is not None]
    if not values:
        return "pending"
    value = max(values)
    if value < 0.5:
        return "<0.5"
    if value < 0.8:
        return "0.5-<0.8"
    return ">=0.8"


def enrich_row(row: dict[str, str], coverage: dict[str, dict[str, float]], neighborhood: dict[str, str]) -> dict[str, str]:
    accession = row["accession"]
    ephaz = coverage.get(accession, {}).get("ephaz")
    iphaz = coverage.get(accession, {}).get("iphaz")
    length = int(row["length"])
    delta = abs(float(row["ephaz_bitscore"]) - float(row["iphaz_bitscore"]))
    out = dict(row)
    out.update({
        "delta_abs": f"{delta:.6f}",
        "delta_abs_bin": _delta_bin(delta),
        "length_bin": _length_bin(length),
        "ephaz_domain_coverage": "" if ephaz is None else f"{ephaz:.6f}",
        "iphaz_domain_coverage": "" if iphaz is None else f"{iphaz:.6f}",
        "domain_coverage_bin": _coverage_bin(ephaz, iphaz),
        "neighborhood_bin": neighborhood.get(accession, "no_record"),
    })
    out["stratum"] = "|".join((out["delta_abs_bin"], out["length_bin"], out["domain_coverage_bin"], out["neighborhood_bin"]))
    return out


def sample_rows(rows: list[dict[str, str]], per_stratum: int, seed: str) -> list[dict[str, str]]:
    if per_stratum < 1:
        raise ValueError("per_stratum must be positive")
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[row["stratum"]].append(row)
    selected = []
    for stratum in sorted(groups):
        ranked = sorted(groups[stratum], key=lambda row: hashlib.sha256(f"{seed}|{stratum}|{row['accession']}".encode()).hexdigest())
        selected.extend(ranked[:per_stratum])
    return sorted(selected, key=lambda row: (row["stratum"], row["accession"]))


def read_fasta(path: Path | str) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header = None
    sequence: list[str] = []
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith(">"):
            if header is not None:
                key = header.split(None, 1)[0]
                if key in records:
                    raise ValueError(f"duplicate FASTA accession: {key}")
                records[key] = (header, "".join(sequence))
            header, sequence = line[1:], []
        elif line:
            if header is None:
                raise ValueError("FASTA sequence precedes header")
            sequence.append(line)
    if header is not None:
        key = header.split(None, 1)[0]
        if key in records:
            raise ValueError(f"duplicate FASTA accession: {key}")
        records[key] = (header, "".join(sequence))
    if not records:
        raise ValueError(f"empty FASTA: {path}")
    return records


def fasta_text(records: dict[str, tuple[str, str]], accessions: list[str]) -> str:
    missing = set(accessions) - set(records)
    if missing:
        raise ValueError(f"sample accessions missing from FASTA: {sorted(missing)[:3]}")
    return "".join(f">{records[accession][0]}\n{records[accession][1]}\n" for accession in accessions)


def _record(path: Path) -> dict[str, object]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": digest}


def _expand_paths(values: list[str]) -> list[Path]:
    paths = [Path(item) for value in values for item in glob.glob(value)]
    if not paths:
        raise ValueError("no domtblout files matched")
    return sorted(set(paths))


def stratify(competition_tsv: Path | str, ephaz_domtblout: list[str], iphaz_domtblout: list[str], neighborhood_tsv: Path | str | None, outdir: Path | str, seed: str = "20260829", per_stratum: int = 10, ambiguous_faa: Path | str | None = None) -> dict[str, object]:
    rows = read_competition(competition_tsv)
    ephaz_paths = _expand_paths(ephaz_domtblout)
    iphaz_paths = _expand_paths(iphaz_domtblout)
    ephaz_cov = parse_domtblout(ephaz_paths)
    iphaz_cov = parse_domtblout(iphaz_paths)
    coverage = {acc: {"ephaz": ephaz_cov[acc], "iphaz": iphaz_cov[acc]} for acc in set(ephaz_cov) | set(iphaz_cov)}
    neighborhood = read_neighborhood(neighborhood_tsv)
    enriched = [enrich_row(row, coverage, neighborhood) for row in rows]
    sampled = sample_rows(enriched, per_stratum, seed)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=False)
    columns = list(enriched[0])
    for name, values in (("ambiguous_stratified.tsv", enriched), ("ambiguous_sample.tsv", sampled)):
        with (outdir / name).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(values)
    summary = []
    counts = Counter(row["stratum"] for row in enriched)
    selected_counts = Counter(row["stratum"] for row in sampled)
    for stratum in sorted(counts):
        summary.append({"stratum": stratum, "population": counts[stratum], "sampled": selected_counts[stratum]})
    with (outdir / "strata_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["stratum", "population", "sampled"], delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(summary)
    outputs = ["ambiguous_stratified.tsv", "ambiguous_sample.tsv", "strata_summary.tsv"]
    inputs = {"competition_tsv": _record(Path(competition_tsv)), "ephaz_domtblout": [_record(path) for path in ephaz_paths], "iphaz_domtblout": [_record(path) for path in iphaz_paths], "neighborhood_tsv": None if neighborhood_tsv is None else _record(Path(neighborhood_tsv))}
    if ambiguous_faa is not None:
        records = read_fasta(ambiguous_faa)
        sample_faa = outdir / "ambiguous_sample.faa"
        sample_faa.write_text(fasta_text(records, [row["accession"] for row in sampled]), encoding="utf-8")
        inputs["ambiguous_faa"] = _record(Path(ambiguous_faa))
        outputs.append("ambiguous_sample.faa")
    metadata = {"schema_version": 1, "seed": seed, "per_stratum": per_stratum, "population": len(enriched), "sampled": len(sampled), "strata": len(summary), "inputs": inputs, "outputs": {name: _record(outdir / name) for name in outputs}}
    (outdir / "sampling_metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return metadata


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--competition-tsv", required=True)
    parser.add_argument("--ephaz-domtblout", action="append", required=True)
    parser.add_argument("--iphaz-domtblout", action="append", required=True)
    parser.add_argument("--neighborhood-tsv")
    parser.add_argument("--ambiguous-faa")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--seed", default="20260829")
    parser.add_argument("--per-stratum", type=int, default=10)
    args = parser.parse_args(argv)
    try:
        stratify(**vars(args))
    except (OSError, ValueError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
