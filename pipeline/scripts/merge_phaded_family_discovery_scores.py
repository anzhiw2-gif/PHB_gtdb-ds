#!/usr/bin/env python3
"""Merge the 33 family discovery-layer HMM pool scores into a superfamily-primary
classification with an explicit ambiguity tier.

The discovery-layer HMMs are broad and, because the 38 DED families are
sequence-clustering subdivisions *within* superfamilies, family-level assignment
is unreliable (most candidates hit 2-22 families). The classification is
therefore anchored at the SUPERFAMILY level (the functional prior), with:

* ``superfamily_confidence`` in {unique, ambiguous_2, ambiguous_3} — how many
  distinct superfamilies the candidate hits;
* family hits reported as a transparent list (``families_hit`` and
  ``families_in_best_superfamily``), never forced into a single "best family".

Experimental-positive provenance is carried at both superfamily and family level.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

MODEL_LAYER = "discovery_hmm_uncalibrated"


def parse_tblout(path: Path) -> dict[str, dict[str, object]]:
    out: dict[str, dict[str, object]] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 6:
                continue
            target = fields[0]
            try:
                e = float(fields[4])
                score = float(fields[5])
            except ValueError:
                continue
            row = out.get(target)
            if row is None:
                out[target] = {"best_e": e, "best_score": score, "n_rows": 1}
            else:
                row["n_rows"] = int(row["n_rows"]) + 1
                if e < float(row["best_e"]):
                    row["best_e"] = e
                    row["best_score"] = score
    return out


def read_provenance(path: Path) -> dict[str, dict[str, str]]:
    prov: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            prov[row["phaded_family_id"]] = row
    return prov


def read_family_superfamily(path: Path) -> dict[str, str]:
    fam2sf: dict[str, str] = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            fam2sf[row["phaded_family_id"]] = row["phaded_superfamily"]
    return fam2sf


def read_candidate_universe(path: Path) -> set[str]:
    ids: set[str] = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    return ids


def merge(
    pool_scores_dir: Path,
    provenance_manifest: Path,
    family_definitions: Path,
    candidate_faa: Path,
    out_dir: Path,
) -> dict[str, object]:
    families = sorted(p.name[:-4] for p in pool_scores_dir.glob("DED_hfam_*.tbl"))
    prov = read_provenance(provenance_manifest)
    fam2sf = read_family_superfamily(family_definitions)
    universe = read_candidate_universe(candidate_faa)

    # candidate -> {family: best_e}
    per_candidate: dict[str, dict[str, float]] = defaultdict(dict)
    family_hits: dict[str, int] = {}
    for fam in families:
        hits = parse_tblout(pool_scores_dir / f"{fam}.tbl")
        family_hits[fam] = len(hits)
        for target, rec in hits.items():
            per_candidate[target][fam] = float(rec["best_e"])

    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) Long-form claim (many-to-many), now with superfamily per family.
    long_path = out_dir / "family_pool_claims_long.tsv"
    with long_path.open("w", encoding="utf-8", newline="\n") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow([
            "protein_id", "phaded_family_id", "phaded_superfamily",
            "best_evalue", "n_families_hit",
            "family_n_experimental_positive", "model_layer",
        ])
        for cand in sorted(per_candidate):
            hits = per_candidate[cand]
            for fam in sorted(hits, key=lambda f: hits[f]):
                p = prov[fam]
                w.writerow([
                    cand, fam, fam2sf[fam], "%.3g" % hits[fam], str(len(hits)),
                    p["n_experimental_positive"], MODEL_LAYER,
                ])

    # 2) Superfamily-primary classification per candidate.
    cls_path = out_dir / "superfamily_classification.tsv"
    n_unique = n_amb2 = n_amb3 = 0
    with cls_path.open("w", encoding="utf-8", newline="\n") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow([
            "protein_id", "superfamily_claim", "superfamily_confidence",
            "n_superfamilies_hit", "best_superfamily", "best_evalue",
            "n_families_hit", "families_hit", "families_in_best_superfamily",
            "best_superfamily_experimental_positive_count", "model_layer",
        ])
        for cand in sorted(per_candidate):
            hits = per_candidate[cand]
            best_fam = min(hits, key=hits.get)
            best_sf = fam2sf[best_fam]
            sfs_hit = sorted({fam2sf[f] for f in hits})
            n_sf = len(sfs_hit)
            confidence = "unique" if n_sf == 1 else (f"ambiguous_{n_sf}")
            if n_sf == 1:
                n_unique += 1
            elif n_sf == 2:
                n_amb2 += 1
            else:
                n_amb3 += 1
            claim = best_sf if n_sf == 1 else "ambiguous:" + "|".join(sfs_hit)
            fams_in_best_sf = sorted(
                (f for f in hits if fam2sf[f] == best_sf), key=lambda f: hits[f]
            )
            sf_ep = sum(
                int(prov[f]["n_experimental_positive"])
                for f in fams_in_best_sf
            )
            w.writerow([
                cand, claim, confidence, str(n_sf), best_sf, "%.3g" % hits[best_fam],
                str(len(hits)), ";".join(sorted(hits, key=lambda f: hits[f])),
                ";".join(fams_in_best_sf), str(sf_ep), MODEL_LAYER,
            ])

    # 3) Superfamily summary.
    sf_summary_path = out_dir / "superfamily_classification_summary.tsv"
    sf_claimed: dict[str, int] = defaultdict(int)
    for cand, hits in per_candidate.items():
        best_fam = min(hits, key=hits.get)
        sf_claimed[fam2sf[best_fam]] += 1
    with sf_summary_path.open("w", encoding="utf-8", newline="\n") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow([
            "phaded_superfamily", "claimed_candidates",
            "experimental_positive_in_reference_only_families",
            "model_layer",
        ])
        for sf in sorted(sf_claimed, key=lambda s: -sf_claimed[s]):
            ep = sum(
                int(prov[f]["n_experimental_positive"])
                for f in prov if fam2sf.get(f) == sf
            )
            w.writerow([sf, str(sf_claimed[sf]), str(ep), MODEL_LAYER])

    n_claimed = len(per_candidate)
    n_unclaimed = len(universe - set(per_candidate))

    return {
        "families": len(families),
        "candidate_universe": len(universe),
        "claimed_candidates": n_claimed,
        "unclaimed_candidates": n_unclaimed,
        "superfamily_confidence_counts": {
            "unique": n_unique, "ambiguous_2": n_amb2, "ambiguous_3": n_amb3,
        },
        "outputs": {
            "long": str(long_path),
            "classification": str(cls_path),
            "superfamily_summary": str(sf_summary_path),
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pool-scores", type=Path, required=True)
    parser.add_argument("--provenance-manifest", type=Path, required=True)
    parser.add_argument("--family-definitions", type=Path, required=True)
    parser.add_argument("--candidate-faa", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    print(json.dumps(
        merge(args.pool_scores, args.provenance_manifest, args.family_definitions,
              args.candidate_faa, args.out_dir),
        indent=2,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
