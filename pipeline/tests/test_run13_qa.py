import csv
import json
from pathlib import Path

from pipeline.scripts.run13_qa import audit


def test_run13_qa_counts_registry_and_oh_filter(tmp_path: Path):
    parent = tmp_path / "parent"
    tier = tmp_path / "tier"
    out = tmp_path / "out"
    (parent / "results").mkdir(parents=True)
    (parent / "logs").mkdir()
    (tier / "results" / "tables").mkdir(parents=True)
    (tier / "logs").mkdir()
    (parent / "results" / "hits_all.tsv").write_text(
        "family\tprotein\tE-value\tscore\tcov\n"
        "OH\tGCF_1|p1\t1e-10\t80\t0.7\n"
        "OH\tGCF_1|p2\t1e-10\t70\t0.5\n"
        "ePhaZ_curated_core\tGCF_2|p3\t1e-20\t90\t0.9\n"
        "OH\tGCF_1|p1\t1e-10\t80\t0.7\n",
        encoding="utf-8",
    )
    (parent / "results" / "scan_manifest.json").write_text('{"status":"completed"}\n')
    (parent / "logs" / "parallel.joblog").write_text(
        "Seq Host Starttime Runtime Send Receive Exitval Signal Command\n"
        "1 host start 1 1 1 0 0 cmd\n"
    )
    (tier / "results" / "tables" / "tier1_genome_family.tsv").write_text(
        "genome\tfamily\tmodels\nGCF_1\tOH\tOH\nGCF_2\tePhaZ\tePhaZ_curated_core\n",
        encoding="utf-8",
    )
    registry = tmp_path / "registry.tsv"
    registry.write_text(
        "model\thmm_source\tthreshold\tmin_cov\treport_group\n"
        "OH\tx\te-5\t0.6\tOH\n"
        "ePhaZ_curated_core\tx\te-5\t0.0\tePhaZ\n",
        encoding="utf-8",
    )

    audit(parent, tier, registry, out)

    rows = {r["model"]: r for r in csv.DictReader((out / "model_counts.tsv").open(), delimiter="\t")}
    assert rows["OH"]["raw_rows"] == "3"
    assert rows["OH"]["accepted_rows"] == "2"
    assert rows["OH"]["coverage_filtered_rows"] == "1"
    assert json.loads((out / "qa_summary.json").read_text())["status"] == "completed"
    report = (out / "QA_REPORT.md").read_text(encoding="utf-8")
    assert "candidate-only" in report


def test_run13_qa_uses_scan_manifest_for_task_completion(tmp_path: Path):
    parent = tmp_path / "parent"; tier = tmp_path / "tier"; out = tmp_path / "out"
    (parent / "results").mkdir(parents=True); (parent / "logs").mkdir(); (tier / "results" / "tables").mkdir(parents=True)
    (parent / "results" / "hits_all.tsv").write_text("family\tprotein\tE-value\tscore\tcov\n", encoding="utf-8")
    (parent / "results" / "scan_manifest.json").write_text('{"status":"completed","task_total":1000,"task_completed":1000,"reused_tasks":297}\n')
    (tier / "results" / "tables" / "tier1_genome_family.tsv").write_text("genome\tfamily\tmodels\n", encoding="utf-8")
    registry = tmp_path / "registry.tsv"; registry.write_text("model\thmm_source\tthreshold\tmin_cov\treport_group\nOH\tx\te-5\t0.6\tOH\n", encoding="utf-8")
    audit(parent, tier, registry, out)
    summary = json.loads((out / "qa_summary.json").read_text())
    assert summary["jobs"]["manifest_task_total"] == 1000
    assert summary["jobs"]["manifest_reused_tasks"] == 297
