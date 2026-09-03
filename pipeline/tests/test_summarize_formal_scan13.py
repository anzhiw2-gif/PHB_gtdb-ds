import csv
from pathlib import Path

from pipeline.scripts.summarize_formal_scan13 import main


def test_scan13_summary_filters_registry_and_unions(tmp_path, monkeypatch):
    hits = tmp_path / "hits.tsv"
    hits.write_text(
        "family\tprotein\tE-value\tscore\tcov\n"
        "ePhaZ_curated_core\tGCA_1|p1\t1e-20\t100\t0.9\n"
        "ePhaZ_broad_discovery\tGCA_1|p2\t1e-10\t90\t0.1\n"
        "OH\tGCF_2|p3\t1e-10\t80\t0.5\n"
        "OH\tGCF_2|p4\t1e-10\t80\t0.7\n",
        encoding="utf-8",
    )
    registry = tmp_path / "registry.tsv"
    registry.write_text(
        "model\thmm_source\tthreshold\tmin_cov\treport_group\n"
        "ePhaZ_curated_core\tx\te-5\t0.0\tePhaZ\n"
        "ePhaZ_broad_discovery\tx\te-5\t0.0\tePhaZ\n"
        "OH\tx\te-5\t0.6\tOH\n",
        encoding="utf-8",
    )
    out = tmp_path / "out"
    monkeypatch.setattr("sys.argv", ["summarize", "--hits", str(hits), "--registry", str(registry), "--outdir", str(out)])
    main()
    rows = list(csv.DictReader((out / "model_summary.tsv").open(), delimiter="\t"))
    assert {r["model"]: int(r["accepted_hits"]) for r in rows} == {
        "ePhaZ_curated_core": 1, "ePhaZ_broad_discovery": 1, "OH": 1
    }
    union = (out / "genome_union_summary.tsv").read_text(encoding="utf-8")
    assert "core_union\t2" in union
