# Repository Governance Cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove browser automation residue, establish reproducible GitHub CI checks, and publish a clear repository license boundary.

**Architecture:** Keep runtime and research files unchanged. Delete the tracked `.playwright-cli/` artifact directory, use one GitHub Actions workflow to run the repository's Python and shell checks on Ubuntu, and add a root MIT license with a note that third-party datasets and dependencies retain upstream terms.

**Tech Stack:** GitHub Actions, Python 3.12, `unittest`, Bash, MIT license.

## Global Constraints

- Do not modify unrelated uncommitted research files.
- CI must not require GTDB inputs, HMMER, or server access.
- Missing/third-party data licensing must not be represented as covered by the repository license.

---

### Task 1: Remove browser automation residue

**Files:**
- Delete: `.playwright-cli/` tracked page snapshots
- Archive: local `.playwright-cli` logs outside the repository
- Modify: `.gitignore` only if needed to retain the existing ignore rule

- [x] Verify the exact tracked target with `git ls-files .playwright-cli`.
- [x] Remove the tracked snapshots with `git rm -r .playwright-cli`.
- [x] Archive local logs at `D:\PHB_gtdb-ds-archive\playwright-cli_20260825` and verify no working-tree `.playwright-cli` files remain.

### Task 2: Harden GitHub Actions checks

**Files:**
- Create or modify: `.github/workflows/ci.yml`

- [x] Run Python compilation for `pipeline/scripts` and `pipeline/tests`.
- [x] Run the committed governance baseline tests (`test_run_context` and `test_run_isolation`). The audit/figure suites remain tied to the separate uncommitted research snapshot.
- [x] Configure shell syntax checks without external bioinformatics inputs; the untracked smoke test is intentionally excluded from this public baseline, and Bash syntax is delegated to Ubuntu CI because Bash is unavailable on this Windows host.
- [x] Run `git diff --check HEAD^ HEAD` in CI against the submitted commit.

### Task 3: Add license boundary

**Files:**
- Create: `LICENSE`
- Modify: `README.md` with a concise license section

- [x] Add the standard MIT license text for repository-authored code and documentation.
- [x] State that GTDB, database, software, and other third-party materials remain under their own terms.
- [x] Verify the license text and README link are present.

### Task 4: Verify and commit

- [x] Run the local test, compile, and whitespace checks; Bash syntax is delegated to CI.
- [x] Review the staged diff to ensure only this governance change is included.
- [x] Create one focused governance commit (amended after review; verify the final SHA with `git log -1`).
