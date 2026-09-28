import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
PREPARE = ROOT / "pipeline" / "scripts" / "prepare_archaea_phazh1_candidate.py"
CLASSIFY = ROOT / "pipeline" / "scripts" / "classify_archaea_phazh1_candidates.py"
BUILD = ROOT / "pipeline" / "scripts" / "build_archaea_phazh1_candidate.py"
CALIBRATE = ROOT / "pipeline" / "scripts" / "calibrate_archaea_phazh1_candidate.py"


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_build_commands_bind_seed_alignment_and_hmm():
    module = load(BUILD, "build_archaea_phazh1_candidate")
    mafft, hmmbuild = module.build_commands(Path("seed.faa"), Path("candidate.sto"), Path("candidate.hmm"), "mafft-x", "hmmbuild-x")
    assert mafft == ["mafft-x", "--auto", "seed.faa"]
    assert hmmbuild == ["hmmbuild-x", "--amino", "candidate.hmm", "candidate.sto"]


def test_calibration_command_has_explicit_cpu_and_no_alignment():
    module = load(CALIBRATE, "calibrate_archaea_phazh1_candidate")
    command = module.search_command(Path("model.hmm"), Path("panel.faa"), Path("out.tbl"), cpu=60)
    assert command == ["hmmsearch", "--noali", "--cpu", "60", "--tblout", "out.tbl", "model.hmm", "panel.faa"]


def test_parse_domtblout_calculates_target_coverage():
    module = load(CLASSIFY, "classify_archaea_phazh1_candidates_domtbl")
    domtbl = "target|locus - 400 model - 305 1e-20 100 0 1 1 1e-20 1e-20 100 0 1 300 21 320 15 330 0.99 -\n"
    parsed = module.parse_domtblout(domtbl.splitlines())
    assert parsed[0]["target"] == "target|locus"
    assert parsed[0]["coverage"] == 0.75


def test_classify_requires_archaea_and_uses_neighborhood_layers():
    module = load(CLASSIFY, "classify_archaea_phazh1_candidates")
    high = module.classify_record({
        "domain": "Archaea", "coverage": "0.91", "evalue": "1e-30",
        "complete": "true", "nearby_markers": "BdhA;PhaJ",
        "patatin_conflict": "false",
    })
    assert high == "PhaZh1_like_high"
    review = module.classify_record({
        "domain": "Archaea", "coverage": "0.67", "evalue": "1e-20",
        "complete": "true", "nearby_markers": "PhaC",
        "patatin_conflict": "false",
    })
    assert review == "PhaZh1_like_review"
    bacteria = module.classify_record({
        "domain": "Bacteria", "coverage": "0.95", "evalue": "1e-50",
        "complete": "true", "nearby_markers": "BdhA",
        "patatin_conflict": "false",
    })
    assert bacteria == "non_PhaZh1_patatin"


def test_classify_fails_closed_for_missing_neighborhood_and_conflict():
    module = load(CLASSIFY, "classify_archaea_phazh1_candidates_missing")
    missing = module.classify_record({
        "domain": "Archaea", "coverage": "0.91", "evalue": "1e-30",
        "complete": "true", "nearby_markers": "", "patatin_conflict": "false",
    })
    assert missing == "archaeal_patatin_exploratory"
    conflict = module.classify_record({
        "domain": "Archaea", "coverage": "0.95", "evalue": "1e-40",
        "complete": "true", "nearby_markers": "BdhA", "patatin_conflict": "true",
    })
    assert conflict == "non_PhaZh1_patatin"


def test_prepare_run_records_required_accessions_and_pending_gtdb(tmp_path):
    module = load(PREPARE, "prepare_archaea_phazh1_candidate")
    records = [
        {"accession": "I3RBH0", "sequence": "M" * 321, "organism": "Haloferax mediterranei", "evidence": "Liu 2015"},
        {"accession": "M1XPT2", "sequence": "M" * 320, "organism": "Natronomonas moolapensis", "evidence": "UniProt reviewed"},
        {"accession": "SUPPORT1", "sequence": "M" * 300, "organism": "Methanobacterium formicicum", "evidence": "supporting homologue"},
    ]
    negative = [{"accession": "NEG1", "sequence": "M" * 280, "label": "non_PhaZh1_patatin"}]
    run = module.prepare_run(tmp_path, "20260904_test_archaea", records, negative)
    contract = json.loads((run / "input_contract.json").read_text(encoding="utf-8"))
    assert contract["status"] == "pending"
    assert (run / "inputs" / "seed_panel.faa").is_file()
    assert (run / "inputs" / "negative_panel.faa").is_file()
    manifest = (run / "inputs" / "reference_manifest.tsv").read_text(encoding="utf-8")
    assert "I3RBH0" in manifest and "M1XPT2" in manifest


def test_prepare_rejects_missing_required_accession(tmp_path):
    module = load(PREPARE, "prepare_archaea_phazh1_candidate_missing")
    with pytest.raises(ValueError, match="I3RBH0"):
        module.prepare_run(tmp_path, "20260904_test_archaea_missing", [
            {"accession": "M1XPT2", "sequence": "M" * 320, "organism": "x", "evidence": "x"},
        ], [])
