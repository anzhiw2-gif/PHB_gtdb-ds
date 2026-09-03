import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "08c_tier_rescore.py"
SPEC = importlib.util.spec_from_file_location("tier_rescore", SCRIPT)
tier = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tier)


def test_hmmsearch_uses_requested_cpu(monkeypatch, tmp_path):
    captured = {}

    def fake_run(command, **kwargs):
        captured["command"] = command

    monkeypatch.setattr(tier.subprocess, "run", fake_run)
    tier.hmmsearch("model.hmm", "input.faa", str(tmp_path / "out.tbl"), "1e-20", cpu=17)
    assert captured["command"][captured["command"].index("--cpu") + 1] == "17"
