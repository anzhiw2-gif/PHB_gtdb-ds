import csv
import json
from pathlib import Path

from pipeline.scripts.stratified_oh_arch_audit import audit


def test_audit_stratifies_family_phylum_and_coverage(tmp_path: Path):
    tier = tmp_path / "tier.tsv"
    hits = tmp_path / "hits.tsv"
    faa = tmp_path / "family.faa"
    out = tmp_path / "out"
    tier.write_text(
        "genome\tfamily\tcopies\tgtdb_taxonomy\tphylum\tclass\n"
        "GCF_A\tOH\t1\td__Bacteria;p__Firmicutes;c__X\tFirmicutes\tX\n"
        "GCF_B\tOH\t1\td__Archaea;p__Euryarchaeota;c__Y\tEuryarchaeota\tY\n"
        "GCF_C\tArchPhaZ_hydrolase\t1\td__Bacteria;p__Proteobacteria;c__Z\tProteobacteria\tZ\n",
        encoding="utf-8",
    )
    hits.write_text(
        "family\tprotein\tE-value\tscore\tcov\n"
        "OH\tGCF_A|p1\t1e-20\t100\t0.7\n"
        "OH\tGCF_B|p2\t1e-30\t110\t0.9\n"
        "ArchPhaZ_hydrolase\tGCF_C|p3\t1e-25\t90\t0.3\n",
        encoding="utf-8",
    )
    faa.write_text(
        ">GCF_A|p1\nM" + "A" * 20 + "\n>GCF_B|p2\nM" + "B" * 20 + "\n>GCF_C|p3\nM" + "C" * 20 + "\n",
        encoding="utf-8",
    )
    faa2 = tmp_path / "family2.faa"
    faa2.write_text(">unused|p9\nAAAA\n", encoding="utf-8")
    audit(tier, hits, [faa, faa2], out, per_stratum=1, seed=7)
    rows = list(csv.DictReader((out / "sample.tsv").open(), delimiter="\t"))
    assert {r["family"] for r in rows} == {"OH", "ArchPhaZ_hydrolase"}
    assert {r["coverage_tier"] for r in rows} == {"oh_review", "oh_high", "arch_exploratory"}
    assert all(r["sequence_status"] == "verified" for r in rows)
    summary = json.loads((out / "summary.json").read_text())
    assert summary["status"] == "completed"
