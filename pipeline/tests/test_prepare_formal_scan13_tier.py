import csv

from pipeline.scripts.prepare_formal_scan13_tier import prepare


def test_prepare_keeps_broad_separate_and_applies_oh_coverage(tmp_path):
    hits = tmp_path / "hits.tsv"
    hits.write_text(
        "family\tprotein\tE-value\tscore\tcov\n"
        "ePhaZ_curated_core\tGCA_1|p1\t1e-20\t90\t1.0\n"
        "ePhaZ_broad_discovery\tGCA_1|p2\t1e-20\t80\t1.0\n"
        "iPhaZ\tGCA_1|p2\t1e-30\t100\t1.0\n"
        "OH\tGCF_2|p3\t1e-20\t70\t0.5\n"
        "OH\tGCF_2|p4\t1e-20\t70\t0.7\n"
        "ArchPhaZ_hydrolase\tGCF_3|p5\t1e-20\t70\t1.0\n",
        encoding="utf-8",
    )
    registry = tmp_path / "registry.tsv"
    registry.write_text(
        "model\thmm_source\tthreshold\tmin_cov\treport_group\n"
        "ePhaZ_curated_core\tx\te-5\t0\tePhaZ\n"
        "ePhaZ_broad_discovery\tx\te-5\t0\tePhaZ\n"
        "iPhaZ\tx\te-5\t0\tiPhaZ\n"
        "OH\tx\te-5\t0.6\tOH\n"
        "ArchPhaZ_hydrolase\tx\te-5\t0\tArchPhaZ_hydrolase\n",
        encoding="utf-8",
    )
    out = tmp_path / "out"
    prepare(hits, registry, out)

    core = list(csv.DictReader((out / "hits_filtered.tsv").open(), delimiter="\t"))
    assert {(row["family"], row["protein"]) for row in core} == {
        ("ePhaZ", "GCA_1|p1"), ("iPhaZ", "GCA_1|p2"),
        ("OH", "GCF_2|p4"), ("ArchPhaZ_hydrolase", "GCF_3|p5"),
    }
    broad = list(csv.DictReader((out / "broad_discovery.tsv").open(), delimiter="\t"))
    assert [(row["protein"], row["source_model"]) for row in broad] == [
        ("GCA_1|p2", "ePhaZ_broad_discovery")
    ]
