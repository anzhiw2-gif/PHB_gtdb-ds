"""Regression guard: the reframed figure funnel must not mislabel the prescreen step as DIAMOND.

Frozen scan 13 (`20260901_formal_frozen_scan_13`) contains no DIAMOND step
(``scan_manifest.json`` records a pure hmmsearch of 100 shards x 10 models).
The 6,769,772 number is the raw HMM hit-row count of an earlier 900-tbl
aggregate (``results/logs/rerun_candidates.log``), NOT a DIAMOND blastp output.
Any figure-funnel definition that still names "DIAMOND" is therefore wrong and
must be caught here.
"""

import unittest
from pathlib import Path

PIPELINE = Path(__file__).resolve().parents[1]

FIGURE_SCRIPTS = [
    PIPELINE / "scripts" / "plot_progress_report_0907.py",
    PIPELINE / "scripts" / "plot_reframed_nature_figure_suite.py",
]

FORBIDDEN_TOKENS = ["DIAMOND", "diamond", "blastp"]


class FigureFunnelLabelTests(unittest.TestCase):
    def test_no_figure_script_mentions_diamond_or_blastp(self):
        offenders = []
        for script in FIGURE_SCRIPTS:
            if not script.is_file():
                self.skipTest(f"script missing: {script}")
            text = script.read_text(encoding="utf-8")
            for token in FORBIDDEN_TOKENS:
                if token in text:
                    offenders.append(f"{script.name}: {token}")
        self.assertEqual(
            offenders,
            [],
            "figure funnel mislabels the prescreen step; "
            "scan 13 has no DIAMOND/blastp step",
        )


if __name__ == "__main__":
    unittest.main()
