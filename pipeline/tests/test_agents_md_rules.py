"""Pin the operating rules that ``AGENTS.md`` must state.

``AGENTS.md`` is read by humans and agents, not by the pipeline, so the rules have no
mechanical enforcement anywhere else in the repository (grep confirms the only reference to
it in ``pipeline/`` is ``finalize_phaded_subtype_reconciliation.py``, which merely records the
file as a provenance input).  These tests give the amendment a regression guard: a rule that
was agreed in the
``20260917_phaded_agents_md_amendment_01`` run cannot silently disappear again.

Evidence behind each rule is recorded in
``runs/20260917_phaded_agents_md_amendment_01/results/amendment_rationale.md``.
"""

from __future__ import annotations

import unittest
from pathlib import Path

AGENTS_MD = Path(__file__).resolve().parents[2] / "AGENTS.md"

# Rules that existed before the amendment and must survive it verbatim in substance.
PRE_EXISTING_RULE_ANCHORS = [
    "运行必须使用新的 `runs/<run_id>/`",
    "`input_contract.json`",
    "不得伪造哈希",
    "`deploy/<run_id>/`",
    "不等同于已验证 PHB 降解表型",
    "保留历史运行残留和失败证据",
    "先写失败测试，再实现",
    "`git diff --check`",
    "发布前核对本地、GitHub、服务器 deploy",
]

# The authorization sentence that the operator removed from the worktree copy of AGENTS.md
# (still present in HEAD:AGENTS.md).  Its blanket form must not come back, because the
# amended rules carve out local commits while keeping remote pushes authorized.
REMOVED_BLANKET_CLAUSE = "提交或推送都必须得到明确授权"


class AgentsMdRulesTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.text = AGENTS_MD.read_text(encoding="utf-8")
        cls.lines = cls.text.splitlines()

    def test_is_present_and_stays_a_rule_file(self):
        self.assertTrue(self.text.strip(), "AGENTS.md must not be empty")
        self.assertLessEqual(
            len(self.lines),
            30,
            "AGENTS.md is an operating-rule file, not a report; keep it under 30 lines",
        )
        rules = [line for line in self.lines if line.startswith("- ")]
        self.assertGreaterEqual(len(rules), 12, "the amended rule set must keep every rule")

    def test_preserves_every_pre_existing_rule(self):
        for anchor in PRE_EXISTING_RULE_ANCHORS:
            with self.subTest(anchor=anchor):
                self.assertIn(anchor, self.text)

    def test_states_the_discovery_layer_rule(self):
        self.assertIn("`discovery_hmm_uncalibrated`", self.text)
        self.assertIn("`annotation_only`", self.text)
        self.assertIn("发现层只用于**召回**", self.text)
        self.assertIn("不得**用于筛选、删除或降级任何候选", self.text)

    def test_states_that_discovery_layer_never_reaches_the_registry(self):
        self.assertIn("`pipeline/config/formal_scan_models.tsv`", self.text)
        self.assertIn("不得产生 family 判定", self.text)

    def test_discovery_layer_may_only_score_existing_intermediates(self):
        # Cross-check against the "no new GTDB full-library scan" constraint: the two rules
        # coexist only because the discovery layer is confined to existing intermediates.
        self.assertIn("既有中间产物", self.text)
        self.assertIn("全库重扫仍需**单独授权**", self.text)

    def test_states_the_model_reproducibility_rule(self):
        self.assertIn("MAFFT 7.525", self.text)
        self.assertIn("必须把**具体比对哈希**写进 manifest", self.text)
        self.assertIn("不得声称逐位可复现", self.text)

    def test_states_the_failure_evidence_archival_rule(self):
        self.assertIn("1,016 GB", self.text)
        self.assertIn("归档目录", self.text)
        self.assertIn("SHA-256 清单", self.text)
        self.assertIn("**禁止删除**", self.text)

    def test_authorization_rule_is_explicit_about_local_commit_and_remote_push(self):
        self.assertIn("本地 `git commit` 允许", self.text)
        self.assertIn("每次必须在交接文档记录", self.text)
        self.assertIn("**push 到 GitHub 仍需操作者明确授权**", self.text)

    def test_authorization_rule_confirmed_by_operator_not_left_pending(self):
        # The operator confirmed option A on 2026-09-17; the "awaiting
        # confirmation" marker must no longer be present.
        self.assertNotIn("待操作者一行确认", self.text)

    def test_blanket_clause_does_not_return_and_contradict_the_carve_out(self):
        self.assertNotIn(
            REMOVED_BLANKET_CLAUSE,
            self.text,
            "restoring the blanket clause would contradict the local-commit carve-out",
        )
        # The remaining actions are still gated on explicit authorization.
        self.assertIn("重算、安装、配置变更、删除与推送仍须明确授权", self.text)
        self.assertIn("只读审计应保持只读", self.text)


if __name__ == "__main__":
    unittest.main()
