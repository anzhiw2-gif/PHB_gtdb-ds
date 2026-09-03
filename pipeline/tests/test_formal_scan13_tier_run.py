from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_tier_run_creates_a_new_auditable_run_and_keeps_broad_separate():
    script = (ROOT / "scripts" / "formal_scan13_tier_processing.sh").read_text(encoding="utf-8")
    assert "create_run_layout" in script
    assert "write_input_contract" in script
    assert "prepare_formal_scan13_tier.py" in script
    assert "parallel_extract_sequences.py" in script
    assert "parallel -j 20" in script
    assert "broad_discovery.tsv" in script
    assert "07b_extract_seqs.py" in script
    assert "08_validate.py" in script
    assert "08c_tier_rescore.py" in script
    assert "--hmm-cpu must be 1..60" in script
    assert 'mv "$BUILD" "$RUN_ROOT/results/tier_processing"' in script
