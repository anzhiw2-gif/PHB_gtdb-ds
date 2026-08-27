#!/usr/bin/env python3
"""Build compact, provenance-tracked source tables for Scheme A figures.

The script reads recorded project outputs and writes figure-specific TSV files.
It does not modify any screening or cluster result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


CORE_FAMILIES = ("ePhaZ", "iPhaZ", "OH", "ArchPhaZ_hydrolase")
ALL_FIGURE_FAMILIES = (*CORE_FAMILIES, "ArchPhaZ_patatin")
GTDB_REPRESENTATIVES = 199_923
FUNNEL_STAGES = (
    ("HMM hit rows", 6_769_772, "hit rows", "All family-model hit records"),
    ("Filtered HMM hit rows", 6_767_839, "hit rows", "Coverage-filtered hit records"),
    ("Unique protein candidates", 6_642_556, "proteins", "Unique protein identifiers"),
    ("Validated candidates", 534_314, "records", "Validated records across five screening groups"),
    ("Strict tier1 core candidates", 74_339, "sequences", "Four core candidate families; patatin excluded"),
)
AVAILABLE_MARKERS = ("PhaC", "PhaE", "PhaJ", "BdhA", "phasin", "PHA_gran_rgn")


class FigureDataError(ValueError):
    """Raised when an input cannot support the approved figure contract."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_table(path: Path) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        raise FigureDataError(f"missing or empty table: {path}")
    frame = pd.read_csv(path, sep=None, engine="python", dtype=str, keep_default_na=False)
    if frame.empty or not len(frame.columns):
        raise FigureDataError(f"empty table: {path}")
    return frame


def require_columns(frame: pd.DataFrame, columns: set[str], label: str) -> None:
    missing = columns.difference(frame.columns)
    if missing:
        raise FigureDataError(f"{label} missing required columns: {', '.join(sorted(missing))}")


def count_fasta_headers(path: Path) -> int:
    if not path.is_file() or path.stat().st_size == 0:
        raise FigureDataError(f"missing or empty FASTA: {path}")
    with path.open(encoding="utf-8") as handle:
        count = sum(line.startswith(">") for line in handle)
    if not count:
        raise FigureDataError(f"FASTA contains no records: {path}")
    return count


def core_union(tier_df: pd.DataFrame) -> int:
    require_columns(tier_df, {"genome", "family"}, "tier table")
    return int(tier_df.loc[tier_df["family"].isin(CORE_FAMILIES), "genome"].nunique())


def figure1_genome_coverage(tier_df: pd.DataFrame, total_genomes: int = GTDB_REPRESENTATIVES) -> pd.DataFrame:
    core_genomes = core_union(tier_df)
    return pd.DataFrame([
        {
            "stage": "GTDB representative genomes",
            "value": total_genomes,
            "unit": "genomes",
            "note": "GTDB R232 input census",
        },
        {
            "stage": "Genomes with >=1 tier1 core candidate",
            "value": core_genomes,
            "unit": "genomes",
            "note": "Unique union across four core candidate families",
        },
    ])


def core_phylum_totals(tier_df: pd.DataFrame) -> pd.DataFrame:
    require_columns(tier_df, {"genome", "family", "phylum"}, "tier table")
    core = tier_df.loc[tier_df["family"].isin(CORE_FAMILIES), ["genome", "phylum"]].drop_duplicates()
    return (
        core.groupby("phylum", sort=True)["genome"].nunique()
        .rename("genomes").sort_values(ascending=False, kind="stable").reset_index()
    )


def seed_family(query_group: str) -> str:
    group = query_group.lower()
    if group.startswith("e-phaz"):
        return "ePhaZ"
    if group.startswith("i-phaz"):
        return "iPhaZ"
    if "oligomer" in group:
        return "OH"
    if "3hb" in group:
        return "BdhA"
    if "phasin" in group:
        return "phasin"
    return query_group


def evidence_class(value: str) -> str:
    tokens = [token.strip().lower() for token in value.split(";") if token.strip()]
    has_pmid = any(token.startswith("pmid:") or token.isdigit() for token in tokens)
    has_doi = any(token.startswith("doi:") or token.startswith("10.") for token in tokens)
    if has_pmid and has_doi:
        return "PMID and DOI"
    if has_pmid:
        return "PMID only"
    if has_doi:
        return "DOI only"
    return "No linked ID"


def extract_seed_tables(seed_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    require_columns(seed_df, {"query_group", "reviewed", "split", "evidence"}, "seed manifest")
    seeds = seed_df.copy()
    seeds["family"] = seeds["query_group"].map(seed_family)
    seeds["review_status"] = seeds["reviewed"].str.lower().eq("true").map({True: "reviewed", False: "unreviewed"})
    family = (
        seeds.groupby(["family", "review_status"], sort=True)
        .size().rename("seeds").reset_index()
    )
    split = seeds.groupby(["family", "split"], sort=True).size().rename("seeds").reset_index()
    seeds["evidence_source"] = seeds["evidence"].map(evidence_class)
    evidence = (
        seeds.groupby(["family", "evidence_source"], sort=True)
        .size().rename("seeds").reset_index()
    )
    return family, evidence, split


def normalize_context(context_df: pd.DataFrame) -> pd.DataFrame:
    base = {"genome", "hit_locus", "hit_family"}
    require_columns(context_df, base, "cluster context")
    if "marker_family" in context_df.columns:
        markers = context_df[["genome", "hit_locus", "hit_family", "marker_family"]].copy()
    elif "nearby_markers" in context_df.columns:
        markers = context_df[["genome", "hit_locus", "hit_family", "nearby_markers"]].copy()
        markers["marker_family"] = markers["nearby_markers"].str.split(",")
        markers = markers.explode("marker_family")
        markers["marker_family"] = markers["marker_family"].str.strip()
        markers = markers.drop(columns="nearby_markers")
    else:
        raise FigureDataError("cluster context requires marker_family or nearby_markers")
    markers = markers.loc[markers["marker_family"].isin(AVAILABLE_MARKERS)]
    return markers.drop_duplicates().reset_index(drop=True)


def neighborhood_rates(context_df: pd.DataFrame, audit_df: pd.DataFrame) -> pd.DataFrame:
    require_columns(audit_df, {"genome", "locus", "family", "status"}, "cluster locus audit")
    if audit_df.empty or audit_df["status"].eq("").any():
        raise FigureDataError("cluster locus audit is incomplete")
    analyzed = audit_df.loc[audit_df["status"] == "analyzed", ["genome", "locus", "family"]].drop_duplicates()
    if analyzed.empty:
        raise FigureDataError("cluster locus audit contains no analyzed candidate loci")
    denominators = analyzed.groupby("family").size().rename("candidate_loci")
    context = normalize_context(context_df)
    supported = (
        context.groupby(["hit_family", "marker_family"], sort=True)
        .size().rename("supported_loci").reset_index()
    )
    supported["candidate_loci"] = supported["hit_family"].map(denominators).fillna(0).astype(int)
    supported = supported.loc[supported["candidate_loci"] > 0].copy()
    supported["support_rate"] = supported["supported_loci"] / supported["candidate_loci"]
    return supported.sort_values(["hit_family", "marker_family"], kind="stable").reset_index(drop=True)


def top_signal_phyla(signal_phylum_df: pd.DataFrame, limit: int = 15) -> pd.DataFrame:
    require_columns(signal_phylum_df, {"phylum", "genomes"}, "SignalP phylum table")
    ranked = signal_phylum_df.copy()
    ranked["genomes"] = pd.to_numeric(ranked["genomes"], errors="raise")
    return ranked.sort_values("genomes", ascending=False, kind="stable").head(limit).reset_index(drop=True)


def table_record(path: Path, frame: pd.DataFrame | None = None) -> dict[str, object]:
    return {
        "path": str(path),
        "sha256": sha256(path),
        "rows": int(len(frame)) if frame is not None else count_fasta_headers(path),
    }


def write_table(frame: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, encoding="utf-8", lineterminator="\n")
    return path


def build_figure_data(args: argparse.Namespace) -> dict[str, Path]:
    seed_path = Path(args.seed_manifest)
    seed_stats_path = Path(args.seed_stats)
    tier_path = Path(args.tier_table)
    phylum_path = Path(args.tier_phylum)
    signal_path = Path(args.signal_summary)
    signal_phylum_path = Path(args.signal_phylum)
    context_path = Path(args.cluster_context)
    audit_path = Path(args.cluster_audit)
    summary_path = Path(args.cluster_summary)
    fasta_dir = Path(args.tier_fasta_dir)
    output_dir = Path(args.output_dir)

    seed_df = read_table(seed_path)
    tier_df = read_table(tier_path)
    phylum_df = read_table(phylum_path)
    signal_df = read_table(signal_path)
    signal_phylum_df = read_table(signal_phylum_path)
    context_df = read_table(context_path)
    audit_df = read_table(audit_path)
    summary_df = read_table(summary_path)
    require_columns(phylum_df, {"phylum", "family", "genomes"}, "tier phylum table")
    require_columns(signal_df, {"metric", "value"}, "SignalP summary")
    require_columns(signal_phylum_df, {"phylum", "genomes"}, "SignalP phylum table")
    require_columns(summary_df, {"hit_family", "marker_family"}, "cluster summary")

    provenance_inputs = {
        "seed_manifest": table_record(seed_path, seed_df),
        "seed_stats": {"path": str(seed_stats_path), "sha256": sha256(seed_stats_path)},
        "tier_table": table_record(tier_path, tier_df),
        "tier_phylum": table_record(phylum_path, phylum_df),
        "signal_summary": table_record(signal_path, signal_df),
        "signal_phylum": table_record(signal_phylum_path, signal_phylum_df),
        "cluster_context": table_record(context_path, context_df),
        "cluster_audit": table_record(audit_path, audit_df),
        "cluster_summary": table_record(summary_path, summary_df),
    }

    seed_family_df, seed_evidence_df, seed_split_df = extract_seed_tables(seed_df)
    outputs = {
        "figure0_seed_family": write_table(seed_family_df, output_dir / "figure0_seed_family.tsv"),
        "figure0_seed_evidence": write_table(seed_evidence_df, output_dir / "figure0_seed_evidence.tsv"),
        "figure0_seed_split": write_table(seed_split_df, output_dir / "figure0_seed_split.tsv"),
        "figure1_funnel": write_table(pd.DataFrame(FUNNEL_STAGES, columns=["stage", "value", "unit", "note"]), output_dir / "figure1_funnel.tsv"),
        "figure1_genome_coverage": write_table(figure1_genome_coverage(tier_df), output_dir / "figure1_genome_coverage.tsv"),
    }

    core = tier_df.loc[tier_df["family"].isin(CORE_FAMILIES), ["genome", "family", "copies"]].copy()
    core["copies"] = pd.to_numeric(core["copies"], errors="raise")
    fasta_rows = []
    for family in CORE_FAMILIES:
        fasta_path = fasta_dir / f"{family}_tier1.faa"
        sequence_count = count_fasta_headers(fasta_path)
        provenance_inputs[f"tier1_fasta_{family}"] = table_record(fasta_path)
        fasta_rows.append({
            "family": family,
            "sequences": sequence_count,
            "genomes": int(core.loc[core["family"] == family, "genome"].nunique()),
        })
    outputs["figure2_core_scale"] = write_table(pd.DataFrame(fasta_rows), output_dir / "figure2_core_scale.tsv")
    outputs["figure2_union"] = write_table(pd.DataFrame([{
        "core_family_union_genomes": core_union(tier_df),
        "gtdb_representatives": GTDB_REPRESENTATIVES,
        "coverage_fraction": core_union(tier_df) / GTDB_REPRESENTATIVES,
    }]), output_dir / "figure2_union.tsv")

    core_phylum = phylum_df.loc[phylum_df["family"].isin(CORE_FAMILIES)].copy()
    core_phylum["genomes"] = pd.to_numeric(core_phylum["genomes"], errors="raise")
    totals = core_phylum_totals(tier_df)
    leading_phyla = list(totals.head(15)["phylum"])
    matrix = core_phylum.loc[core_phylum["phylum"].isin(leading_phyla)].copy()
    matrix["family"] = pd.Categorical(matrix["family"], categories=CORE_FAMILIES, ordered=True)
    matrix = matrix.sort_values(["phylum", "family"], kind="stable")
    outputs["figure3_phylum_matrix"] = write_table(matrix, output_dir / "figure3_phylum_matrix.tsv")
    outputs["figure3_phylum_totals"] = write_table(totals.head(15), output_dir / "figure3_phylum_totals.tsv")

    signal_classes = signal_df.loc[signal_df["metric"].isin(["pred_SP", "pred_LIPO", "pred_TAT", "pred_TATLIPO", "pred_OTHER"]), ["metric", "value"]].copy()
    signal_classes["signal_class"] = signal_classes["metric"].str.removeprefix("pred_")
    signal_classes["sequences"] = pd.to_numeric(signal_classes["value"], errors="raise")
    signal_classes = signal_classes[["signal_class", "sequences"]]
    outputs["figure4_signal_classes"] = write_table(signal_classes, output_dir / "figure4_signal_classes.tsv")
    outputs["figure4_signal_phyla"] = write_table(top_signal_phyla(signal_phylum_df), output_dir / "figure4_signal_phyla.tsv")

    rates = neighborhood_rates(context_df, audit_df)
    outputs["figure5_neighborhood_rate"] = write_table(rates, output_dir / "figure5_neighborhood_rate.tsv")
    patatin = rates.loc[rates["hit_family"] == "ArchPhaZ_patatin"].copy()
    audit_analyzed = audit_df.loc[(audit_df["family"] == "ArchPhaZ_patatin") & (audit_df["status"] == "analyzed")]
    patatin["candidate_genomes"] = int(audit_analyzed["genome"].nunique())
    patatin["supported_genomes"] = [
        int(normalize_context(context_df).loc[
            (normalize_context(context_df)["hit_family"] == "ArchPhaZ_patatin") &
            (normalize_context(context_df)["marker_family"] == marker), "genome"
        ].nunique())
        for marker in patatin["marker_family"]
    ]
    outputs["figure5_patatin_context"] = write_table(patatin, output_dir / "figure5_patatin_context.tsv")

    provenance = {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "cluster_snapshot_id": args.cluster_snapshot_id,
        "core_families": list(CORE_FAMILIES),
        "available_neighborhood_markers": list(AVAILABLE_MARKERS),
        "inputs": provenance_inputs,
        "outputs": {name: {"file": path.name, "sha256": sha256(path), "rows": int(len(read_table(path)))} for name, path in outputs.items()},
    }
    provenance_path = output_dir / "provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    outputs["provenance"] = provenance_path
    return outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-manifest", required=True)
    parser.add_argument("--seed-stats", required=True)
    parser.add_argument("--tier-table", required=True)
    parser.add_argument("--tier-phylum", required=True)
    parser.add_argument("--tier-fasta-dir", required=True)
    parser.add_argument("--signal-summary", required=True)
    parser.add_argument("--signal-phylum", required=True)
    parser.add_argument("--cluster-context", required=True)
    parser.add_argument("--cluster-audit", required=True)
    parser.add_argument("--cluster-summary", required=True)
    parser.add_argument("--cluster-snapshot-id", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def main() -> None:
    outputs = build_figure_data(parse_args())
    print(f"Wrote {len(outputs)} Scheme A figure data artifacts.")


if __name__ == "__main__":
    main()
