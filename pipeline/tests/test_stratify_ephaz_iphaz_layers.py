import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "stratify_ephaz_iphaz_layers.py"


def load():
    spec = importlib.util.spec_from_file_location("stratify_ephaz_iphaz_layers", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_ephaz_layers_keep_core_and_signal_context_distinct():
    module = load()
    assert module.layer_for({"family": "ePhaZ", "source_layer": "ePhaZ_curated_core", "signal_type": "SP"}) == "ePhaZ_curated_core_secreted"
    assert module.layer_for({"family": "ePhaZ", "source_layer": "ePhaZ_curated_core", "signal_type": "OTHER"}) == "ePhaZ_curated_core_nonsecreted"
    assert module.layer_for({"family": "ePhaZ", "source_layer": "ePhaZ_broad_discovery", "signal_type": "SP"}) == "ePhaZ_broad_discovery"


def test_iphaz_tiers_are_explicit_and_missing_metadata_is_review():
    module = load()
    assert module.layer_for({"family": "iPhaZ", "source_layer": "tier1", "signal_type": ""}) == "iPhaZ_tier1"
    assert module.layer_for({"family": "iPhaZ", "source_layer": "tier2", "signal_type": ""}) == "iPhaZ_tier2_review"
    assert module.layer_for({"family": "iPhaZ", "source_layer": "", "signal_type": ""}) == "iPhaZ_review_pending"


def test_tier2_ephaz_is_review_layer_and_competition_overrides():
    module = load()
    assert module.layer_for({"family": "ePhaZ", "source_layer": "tier2", "signal_type": "SP"}) == "ePhaZ_tier2_review"
    assert module.layer_for({"family": "ePhaZ", "source_layer": "ePhaZ_curated_core", "assignment": "ambiguous", "signal_type": "OTHER"}) == "ePhaZ_competition_review"
    assert module.layer_for({"family": "iPhaZ", "source_layer": "tier1", "assignment": "ambiguous", "signal_type": ""}) == "iPhaZ_competition_review"


def test_unknown_family_fails_closed():
    module = load()
    with pytest.raises(ValueError, match="unsupported family"):
        module.layer_for({"family": "OH", "source_layer": "tier1", "signal_type": ""})


def test_stratify_writes_genome_counts_without_duplicate_scans(tmp_path):
    module = load()
    ephaz = tmp_path / "e.faa"
    iphaz = tmp_path / "i.faa"
    ephaz.write_text(">G1|p1\nAAAA\n>G1|p2\nCCCC\n", encoding="utf-8")
    iphaz.write_text(">G1|p3\nGGGG\n", encoding="utf-8")
    output = tmp_path / "out"
    module.stratify(ephaz, iphaz, output)
    rows = (output / "genome_layers.tsv").read_text(encoding="utf-8").splitlines()
    assert rows == [
        "genome\tfamily\tlayers\tprotein_count",
        "G1\tePhaZ\tePhaZ_curated_core_nonsecreted\t2",
        "G1\tiPhaZ\tiPhaZ_tier1\t1",
    ]


def test_stratify_deduplicates_accessions_with_highest_evidence_layer(tmp_path):
    module = load()
    ephaz = tmp_path / "e.faa"
    iphaz = tmp_path / "i.faa"
    ephaz.write_text(">G1|p1\nAAAA\n", encoding="utf-8")
    iphaz.write_text(">G1|p1\nGGGG\n", encoding="utf-8")
    ephaz_tier2 = tmp_path / "e2.faa"
    iphaz_tier2 = tmp_path / "i2.faa"
    ephaz_tier2.write_text(">G1|p1\nAAAA\n>G1|p2\nCCCC\n", encoding="utf-8")
    iphaz_tier2.write_text(">G1|p1\nGGGG\n>G1|p3\nTTTT\n", encoding="utf-8")
    output = tmp_path / "out"
    module.stratify(ephaz, iphaz, output, ephaz_tier2=ephaz_tier2, iphaz_tier2=iphaz_tier2)
    rows = (output / "protein_layers.tsv").read_text(encoding="utf-8").splitlines()
    assert sum("G1|p1" in row for row in rows) == 2
    assert sum("G1|p2" in row for row in rows) == 1
    assert sum("G1|p3" in row for row in rows) == 1
    by_key = {(parts[2], parts[0]): parts[-1] for parts in (row.split("\t") for row in rows[1:])}
    assert by_key[("ePhaZ", "G1|p1")] == "ePhaZ_curated_core_nonsecreted"
    assert by_key[("iPhaZ", "G1|p1")] == "iPhaZ_tier1"
