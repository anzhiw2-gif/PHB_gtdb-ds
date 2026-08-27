# Research Change Split Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans (recommended) to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Split the current unpublished research workspace into auditable script, test, result, and documentation commits without mixing large machine-local artifacts into Git.

**Architecture:** Preserve the existing governance baseline and classify changes by ownership boundary. Commit source/configuration first, then tests, then lightweight derived evidence, and finally reports/plans. Large generated figures, staging ledgers, raw SignalP output, and oversized tables remain local until a Release/Zenodo packaging decision.

**Tech Stack:** Python, Bash, unittest, TSV/JSON/Newick evidence files, Markdown, Git.

## Global Constraints

- Do not modify or delete unrelated user work.
- Do not commit raw GTDB data, bulk generated FASTA, raw SignalP output, TIFFs, or staging ledgers over 25 MB.
- Every commit must pass relevant tests, `compileall`, and `git diff --check`.
- Research evidence remains candidate/function-potential evidence, not experimentally verified PHB phenotype.

---

### Task 1: Script implementation commit

**Files:**
- Modify: `pipeline/config/params.yaml`, `pipeline/dev/resume_0611.sh`, and modified `pipeline/scripts/*.py|*.sh`
- Create: new pipeline scripts under `pipeline/scripts/`
- Create: `pipeline/seeds/controls/*`

- [x] Stage only source/configuration/control-set files.
- [x] Run `python -m compileall -q pipeline/scripts pipeline/tests`.
- [x] Run the existing governance tests `test_run_context` and `test_run_isolation`.
- [x] Run `git diff --cached --check` and commit as `科研：提交脚本与配置改造` (`84dccc5`).

### Task 2: Test commit

**Files:**
- Create: `pipeline/tests/test_export_archphaz_hydrolase_archaea.py`
- Create: `pipeline/tests/test_scheme_a_figure_data.py`
- Create: `pipeline/tests/test_scheme_a_figure_exports.py`
- Create: `pipeline/tests/test_scheme_a_figure_qa.py`
- Create: `pipeline/tests/test_smoke.py`
- Create: `pipeline/tests/test_smoke.sh`

- [x] Stage only test files.
- [x] Run all available local tests with the current research workspace dependencies (`53/53` passed).
- [x] Run `compileall` and `git diff --cached --check`.
- [x] Commit as `科研：加入结果与管线回归测试` (`b3ae977`).

### Task 3: Lightweight result commit

**Files:**
- Modify: `results/tables/tier1_phylum_distribution.tsv`, `results/trees_tier1/iPhaZ.treefile`
- Create: small TSV/JSON outputs under `results/tables/`, `results/trees_tier1/`, `results/archphaz_hydrolase_archaea/`, and `results/figures/scheme_a/source_data/`

- [x] Exclude `results/tables/tier1_genome_family.tsv` (about 25 MB), TIFF/large PNG outputs, staging TSVs, raw SignalP output, and other large generated files.
- [x] Validate TSV/JSON/Newick files are nonempty and inspect their provenance fields.
- [x] Run `git diff --cached --check` and commit as `科研：提交轻量结果与图源证据` (`8c4f90d`).

### Task 4: Documentation commit

**Files:**
- Modify: `CHANGELOG.md`, `docs/final_results_report.md`, `docs/审查报告.md`, `docs/项目审核报告.md`
- Create: `docs/T141_RECOMPUTE_HANDOFF.md`, current audit report, and `docs/superpowers/{plans,specs}/*`

- [x] Stage report, handoff, plan, and specification files only.
- [x] Check that claims distinguish completed, preliminary, pending, and unpublished work.
- [x] Run `git diff --cached --check` and commit as `文档：归档科研改动与复算交接`.

### Task 5: Final verification and remote decision

- [x] Run the full local test suite, `compileall`, and whitespace checks (`53/53` passed).
- [x] Review four commit diffs and confirm no large/local-only artifact was included.
- [x] Report the final commit list and remaining local artifacts; push only with explicit authorization.
