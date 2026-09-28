"""Audit predicted dPHASCL2 structures against domain-coordinate evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def _values(pdb: Path) -> dict[int, float]:
    values: dict[int, float] = {}
    with pdb.open(encoding="ascii", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            try:
                residue = int(line[22:26])
                score = float(line[60:66])
            except ValueError:
                continue
            values.setdefault(residue, score)
    return values


def _stats(values: list[float]) -> tuple[int, float, float, float, float]:
    if not values:
        return 0, 0.0, 0.0, 0.0, 0.0
    low = sum(value < 70.0 for value in values) / len(values)
    very_high = sum(value >= 90.0 for value in values) / len(values)
    longest = current = 0
    for value in values:
        if value < 70.0:
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return len(values), sum(values) / len(values), low, very_high, float(longest)


def _low_runs(values: dict[int, float]) -> str:
    runs: list[str] = []
    start = None
    previous = None
    for residue in sorted(values):
        if values[residue] < 70.0:
            if start is None or previous is None or residue != previous + 1:
                if start is not None:
                    runs.append(f"{start}-{previous}")
                start = residue
            previous = residue
        elif start is not None:
            runs.append(f"{start}-{previous}")
            start = previous = None
    if start is not None:
        runs.append(f"{start}-{previous}")
    return ";".join(runs)


def _mean_pairs(matrix: list[list[float]], left: list[int], right: list[int], exclude_diagonal: bool = False) -> float:
    values = [float(matrix[i - 1][j - 1]) for i in left for j in right if not (exclude_diagonal and i == j)]
    return round(sum(values) / len(values), 3) if values else 0.0


def summarize_pae(matrix: list[list[float]], sequence_length: int, domain_start: int, domain_end: int, nterm_end: int = 30) -> dict[str, object]:
    domain = list(range(domain_start, domain_end + 1))
    nterm = list(range(1, min(nterm_end, sequence_length) + 1))
    outside = [i for i in range(1, sequence_length + 1) if i not in domain]
    return {
        "pae_domain_internal_mean": _mean_pairs(matrix, domain, domain, exclude_diagonal=True),
        "pae_nterm_domain_mean": _mean_pairs(matrix, nterm, domain),
        "pae_domain_outside_mean": _mean_pairs(matrix, domain, outside),
    }


def summarize_pdb(pdb: Path, sequence_length: int, domain_start: int, domain_end: int, nterm_end: int = 30) -> dict[str, object]:
    values = _values(pdb)
    domain = [values[i] for i in range(domain_start, domain_end + 1) if i in values]
    nterm = [values[i] for i in range(1, min(nterm_end, sequence_length) + 1) if i in values]
    outside = [values[i] for i in range(1, sequence_length + 1) if i not in range(domain_start, domain_end + 1) and i in values]
    d_count, d_mean, d_low, d_high, d_longest = _stats(domain)
    n_count, n_mean, n_low, n_high, n_longest = _stats(nterm)
    o_count, o_mean, o_low, o_high, o_longest = _stats(outside)
    return {
        "pdb": str(pdb),
        "pdb_sha256": hashlib.sha256(pdb.read_bytes()).hexdigest(),
        "sequence_length": sequence_length,
        "pdb_residue_count": len(values),
        "pdb_coverage": len(values) / sequence_length if sequence_length else 0.0,
        "low_confidence_runs": _low_runs(values),
        "domain_start": domain_start,
        "domain_end": domain_end,
        "domain_residue_count": d_count,
        "domain_plddt_mean": round(d_mean, 3),
        "domain_low_plddt_fraction": round(d_low, 4),
        "domain_high_plddt_fraction": round(d_high, 4),
        "domain_longest_low_run": int(d_longest),
        "nterm_end": min(nterm_end, sequence_length),
        "nterm_residue_count": n_count,
        "nterm_plddt_mean": round(n_mean, 3),
        "nterm_low_plddt_fraction": round(n_low, 4),
        "nterm_high_plddt_fraction": round(n_high, 4),
        "nterm_longest_low_run": int(n_longest),
        "outside_domain_residue_count": o_count,
        "outside_domain_plddt_mean": round(o_mean, 3),
        "outside_domain_low_plddt_fraction": round(o_low, 4),
        "outside_domain_high_plddt_fraction": round(o_high, 4),
        "outside_domain_longest_low_run": int(o_longest),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    records = manifest.get("candidates", manifest if isinstance(manifest, list) else [])
    rows = []
    for item in records:
        domain_start = int(item.get("domain_start", item["pfam_start"]))
        domain_end = int(item.get("domain_end", item["pfam_end"]))
        nterm_end = int(item.get("nterm_end", 30))
        row = summarize_pdb(Path(item["pdb"]), int(item["sequence_length"]), domain_start, domain_end, nterm_end)
        if item.get("score_json"):
            score_data = json.loads(Path(item["score_json"]).read_text(encoding="utf-8"))
            row.update(summarize_pae(score_data["pae"], int(item["sequence_length"]), domain_start, domain_end, nterm_end))
        row.update({k: v for k, v in item.items() if k not in row})
        rows.append(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    args.json.write_text(json.dumps({"row_count": len(rows), "rows": rows}, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
