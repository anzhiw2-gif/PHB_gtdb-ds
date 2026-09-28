import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "pipeline" / "scripts" / "annotate_archaea_domains.py"


def load_module():
    spec = importlib.util.spec_from_file_location("annotate_archaea_domains", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_validate_candidate_table_requires_high_and_review_layers(tmp_path):
    module = load_module()
    candidate = tmp_path / "candidate.tsv"
    candidate.write_text(
        "genome\tlocus\tlayer\tphylum\tnearby_markers\n"
        "g1\tp1\tPhaZh1_like_high\tHalobacteriota\tBdhA\n"
        "g2\tp2\tPhaZh1_like_review\tThermoproteota\tPhaJ\n",
        encoding="utf-8",
    )
    rows = module.read_candidate_table(candidate)
    assert [row["layer"] for row in rows] == ["PhaZh1_like_high", "PhaZh1_like_review"]


def test_validate_candidate_table_accepts_utf8_bom_header(tmp_path):
    module = load_module()
    candidate = tmp_path / "candidate.tsv"
    candidate.write_text(
        "genome\tlocus\tlayer\tphylum\tnearby_markers\n"
        "g1\tp1\tPhaZh1_like_high\tHalobacteriota\tBdhA\n",
        encoding="utf-8-sig",
    )
    assert module.read_candidate_table(candidate)[0]["genome"] == "g1"


def test_validate_candidate_table_rejects_duplicate_locus_and_unknown_layer(tmp_path):
    module = load_module()
    candidate = tmp_path / "candidate.tsv"
    candidate.write_text(
        "genome\tlocus\tlayer\tphylum\tnearby_markers\n"
        "g1\tp1\tPhaZh1_like_high\tHalobacteriota\tBdhA\n"
        "g1\tp1\tarchaeal_patatin_exploratory\tHalobacteriota\t\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="duplicate locus"):
        module.read_candidate_table(candidate)


def test_model_manifest_requires_bound_hash_and_role(tmp_path):
    module = load_module()
    model = tmp_path / "ArchPhaZ.hmm"
    model.write_text("HMMER3/f [3.4 | 2023]\n", encoding="ascii")
    manifest = tmp_path / "models.tsv"
    manifest.write_text(
        "family\tpath\tversion\tsha256\trole\n"
        f"ArchPhaZ_patatin\t{model}\tHMMER3.4\t{module.sha256_file(model)}\tpositive_domain\n",
        encoding="utf-8",
    )
    rows = module.read_model_manifest(manifest)
    assert rows[0]["family"] == "ArchPhaZ_patatin"
    assert rows[0]["role"] == "positive_domain"


def test_model_manifest_rejects_unbound_hash(tmp_path):
    module = load_module()
    model = tmp_path / "ArchPhaZ.hmm"
    model.write_text("HMMER3/f [3.4 | 2023]\n", encoding="ascii")
    manifest = tmp_path / "models.tsv"
    manifest.write_text(
        "family\tpath\tversion\tsha256\trole\n"
        f"ArchPhaZ_patatin\t{model}\tHMMER3.4\tpending\tpositive_domain\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="sha256"):
        module.read_model_manifest(manifest)


def test_classify_domain_status_keeps_conflicts_and_pending_separate():
    module = load_module()
    assert module.classify_domain_status({"ArchPhaZ_patatin": 1}, set()) == "domain_supported"
    assert module.classify_domain_status(
        {"ArchPhaZ_hydrolase": 1}, set(), {"ArchPhaZ_patatin", "ArchPhaZ_hydrolase"}
    ) == "domain_supported"
    assert module.classify_domain_status({"ArchPhaZ_patatin": 1, "PhaC": 1}, {"PhaC"}) == "domain_conflict"
    assert module.classify_domain_status({}, set()) == "no_registered_domain_hit"


def test_parse_domtblout_uses_query_accession_and_domain_family():
    module = load_module()
    line = (
        "g1|p1 - 321 ArchPhaZ_patatin - 317 1e-40 120 0 "
        "1 1 1e-40 1e-40 120 0 1 300 21 320 15 330 0.99 -\n"
    )
    parsed = module.parse_domtblout([line])
    assert parsed == [{"accession": "g1|p1", "family": "ArchPhaZ_patatin", "i_evalue": 1e-40, "bitscore": 120.0}]


def test_parse_domtblout_can_bind_manifest_family_over_hmm_name():
    module = load_module()
    line = (
        "g1|p1 - 321 ArchPhaZ_patatin_aln - 317 1e-40 120 0 "
        "1 1 1e-40 1e-40 120 0 1 300 21 320 15 330 0.99 -\n"
    )
    parsed = module.parse_domtblout([line], family_override="ArchPhaZ_patatin")
    assert parsed[0]["family"] == "ArchPhaZ_patatin"


def test_parse_domtblout_excludes_hits_above_registered_i_evalue_cutoff():
    module = load_module()
    line = (
        "g1|p1 - 321 ArchPhaZ_patatin_aln - 317 1e-4 120 0 "
        "1 1 1e-4 1e-4 120 0 1 300 21 320 15 330 0.99 -\n"
    )
    assert module.parse_domtblout([line], max_i_evalue=1e-5) == []
