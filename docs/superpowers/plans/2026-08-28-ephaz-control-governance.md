# ePhaZ Control Governance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reclassify the six cross-family-confounded controls, recalibrate ePhaZ controls by evidence class, and audit the existing no-signal-peptide ePhaZ candidates against ePhaZ and iPhaZ HMMs without launching the 268G screening.

**Architecture:** A deterministic control-governance script will produce a new control-class manifest and separate FASTA/TSV panels while preserving all six sequences. A second audit script will parse SignalP and HMMER outputs, classify no-signal-peptide candidates by model competition, and write only to a new dated run. Existing historical runs and `results/` remain read-only.

**Tech Stack:** Python 3.10 standard library, HMMER 3.4 on T141, existing TSV/FASTA contracts, `unittest`.

## Global Constraints

- Use a new `runs/<run_id>/`; never overwrite historical `results/` or server runs.
- Keep the six FASTA sequences; do not promote them to `ePhaZ_curated_core` or ePhaZ bridge/reference.
- Record paths, sizes, versions, and SHA-256; use `pending` for unavailable values.
- Server execution must use dated `deploy/<run_id>/` source bindings.
- HMM/domain/SignalP evidence indicates candidate homology or functional potential, not validated PHB-degradation phenotype.
- Run focused tests, `compileall`, and `git diff --check`; no commit or push in this turn.

### Task 1: Control classification contract

**Files:**
- Create: `pipeline/scripts/reclassify_ephaz_controls.py`
- Create: `pipeline/tests/test_reclassify_ephaz_controls.py`
- Create: `pipeline/seeds/ephaz_control_classes.tsv`

- [ ] Write failing tests for six fixed accessions mapping to `cross_family_confounded`, preservation in a challenge FASTA, exclusion from ePhaZ positives, and fail-closed behavior for missing accessions.
- [ ] Run the focused test and confirm failure because the module/output contract is absent.
- [ ] Implement deterministic TSV/FASTA parsing, exact accession validation, class assignment, SHA-256 manifest, and outputs `ephaz_positive_controls.faa`, `iphaZ_like_challenge.faa`, `ephaz_control_classes.tsv`, and `control_governance.json`.
- [ ] Run focused tests and confirm pass.

### Task 2: Evidence-class calibration

**Files:**
- Create: `pipeline/scripts/calibrate_ephaz_control_classes.py`
- Create: `pipeline/tests/test_calibrate_ephaz_control_classes.py`

- [ ] Write failing parser/classification tests proving separate `ePhaZ_curated_core`, `ePhaZ_architecture_remote`, and `iPhaZ_like_challenge` metrics.
- [ ] Run the focused test and confirm failure.
- [ ] Implement HMMER tblout/domtblout parsing, per-class TP/FN/FP/TN summaries, model/evidence hashes, and fail-closed control-manifest binding.
- [ ] Run focused tests and confirm pass.

### Task 3: No-signal-peptide competition audit

**Files:**
- Create: `pipeline/scripts/audit_ephaz_iphaZ_competition.py`
- Create: `pipeline/tests/test_audit_ephaz_iphaZ_competition.py`

- [ ] Write failing tests for SignalP `OTHER` selection, HMM score parsing, tie/ambiguous assignment, and missing-record rejection.
- [ ] Run the focused test and confirm failure.
- [ ] Implement deterministic extraction of no-signal-peptide candidates, ePhaZ/iPhaZ score comparison, architecture flags, and TSV/JSON provenance outputs.
- [ ] Run focused tests and confirm pass.

### Task 4: Dated T141 smoke run and documentation

**Files:**
- Create: `docs/T141_20260828_ephaz_control_governance_status.md`
- Create: `runs/20260828_ephaz_control_governance_01/input_contract.json`
- Modify: `README.md` only after result acceptance; no full-screening counts in this turn.

- [ ] Stage dated deploy files and inputs without touching historical runs.
- [ ] Run only control-class calibration and no-signal-peptide competition audit on T141.
- [ ] Verify output completeness, SHA-256, tests, `compileall`, and `git diff --check`.
- [ ] Document the six controls as cross-family challenges and record whether full screening remains blocked.
