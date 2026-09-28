# 古菌 PhaZh1 候选重建实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在独立的 candidate-only run 中，依据 *Haloferax mediterranei* PhaZh1 实验证据，构建一个带邻域联合判定的古菌 patatin/PhaZh1-like 候选层，并保留不可区分的广谱 patatin 结果。

**Architecture:** 先固定 UniProt accession、现有古菌 tier 输入和负对照的 SHA-256，再用小型候选 HMM 做同源召回；最终分类器同时读取 domain、coverage、序列完整性、GTDB domain、邻域 marker 和冲突注释，输出 high/review/exploratory/non-PhaZh1 四层。所有结果只写入 `runs/20260904_archaea_phazh1_candidate_01/`，不修改正式 registry 或正式扫描模型。

**Tech Stack:** Python 3 标准库、HMMER 3.4 (`hmmalign`/`hmmbuild`/`hmmsearch`)、现有 TSV/FASTA 解析器、SHA-256 provenance。

## Global Constraints

- 运行目录必须包含 `logs/`、`inputs/`、`results/` 和 `input_contract.json`。
- 新 run 只允许 candidate-only；不得覆盖 `results/`、历史 run、`pipeline/config/formal_scan_models.tsv`。
- HMM/domain/邻域结果只表示候选同源或功能潜力，不等同于实验验证的 PHB 降解表型。
- 缺失输入、无法核验的 accession 或未绑定的 taxonomy 必须 fail-closed，并记录 `pending`。
- 代码改动遵循 TDD；完成后运行相关测试、`compileall` 和 `git diff --check`。

### Task 1: Run scaffold and input contract

**Files:**
- Create: `runs/20260904_archaea_phazh1_candidate_01/logs/`
- Create: `runs/20260904_archaea_phazh1_candidate_01/inputs/`
- Create: `runs/20260904_archaea_phazh1_candidate_01/results/`
- Create: `runs/20260904_archaea_phazh1_candidate_01/input_contract.json`
- Test: `pipeline/tests/test_archaea_phazh1_candidate.py`

- [x] Write failing tests for required run directories, candidate-only status, and SHA-256/pending input fields.
- [x] Run the focused test and verify failure.
- [x] Implement the minimal manifest/scaffold helpers.
- [x] Re-run focused tests and verify pass.

### Task 2: Reference and control panel

**Files:**
- Create: `pipeline/scripts/prepare_archaea_phazh1_candidate.py`
- Create: `runs/20260904_archaea_phazh1_candidate_01/inputs/reference_manifest.tsv`
- Create: `runs/20260904_archaea_phazh1_candidate_01/inputs/negative_panel.faa`
- Create: `runs/20260904_archaea_phazh1_candidate_01/inputs/seed_panel.faa`
- Modify: `pipeline/tests/test_archaea_phazh1_candidate.py`

- [x] Add tests requiring `I3RBH0`, `M1XPT2`, source URL/date, complete sequence, and explicit negative-panel labels.
- [x] Fetch accession records through the UniProt REST endpoint, save raw responses only in the new run, and record response hashes.
- [x] Assemble a small positive/supporting panel and non-PhaZh1 patatin/phospholipase controls without silently claiming experimental validation.
- [x] Fail closed if a required sequence is absent or malformed.

### Task 3: Candidate HMM and calibration

**Files:**
- Create: `pipeline/scripts/build_archaea_phazh1_candidate.py`
- Create: `runs/20260904_archaea_phazh1_candidate_01/results/archaea_PhaZh1_like_candidate.sto`
- Create: `runs/20260904_archaea_phazh1_candidate_01/results/archaea_PhaZh1_like_candidate.hmm`
- Create: `runs/20260904_archaea_phazh1_candidate_01/results/calibration.tsv`

- [x] Test deterministic command construction and rejection of empty/invalid panels.
- [x] Build alignment/HMM once from the fixed panel; record HMMER versions and hashes.
- [x] Calibrate E-value and coverage on the positive/supporting and negative panels; retain only thresholds with zero negative hits and improved positive recall.

### Task 4: Neighborhood joint classifier

**Files:**
- Create: `pipeline/scripts/classify_archaea_phazh1_candidates.py`
- Create: `runs/20260904_archaea_phazh1_candidate_01/results/candidate_layers.tsv`
- Create: `runs/20260904_archaea_phazh1_candidate_01/results/classifier_summary.tsv`
- Modify: `pipeline/tests/test_archaea_phazh1_candidate.py`

- [x] Test domain filtering, coverage tiers, missing-neighborhood fail-closed behavior, positive marker support, and phospholipase conflict handling.
- [x] Implement four labels: `PhaZh1_like_high`, `PhaZh1_like_review`, `archaeal_patatin_exploratory`, `non_PhaZh1_patatin`.
- [ ] Join existing run-13/tier evidence only through explicit column names and bound input hashes.

### Task 5: Audit report and verification

**Files:**
- Create: `runs/20260904_archaea_phazh1_candidate_01/results/provenance.json`
- Create: `docs/T141_20260904_archaea_phazh1_candidate_01_status.md`

- [x] Generate counts by layer, archaeal phylum, marker support, and conflict reason.
- [x] State that the model is candidate-only and that PhaZh1 experimental evidence is limited to the cited accession/paper.
- [x] Run focused tests, all pipeline tests, `python -m compileall pipeline/scripts`, and `git diff --check`.
- [x] Confirm no formal model, registry, historical run, or GTDB scan was modified or started.
