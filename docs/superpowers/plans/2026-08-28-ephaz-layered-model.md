# ePhaZ Layered Model Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task.

**Goal:** Build, validate, and calibrate separate curated-core and broad-discovery ePhaZ HMM layers without deleting short sequences.

**Architecture:** A fail-closed Python splitter creates auditable layer FASTA/TSV outputs. The shell HMM builder accepts an explicit family list and emits independent alignment/HMM artifacts. A parameterized calibration command evaluates each layer with evidence-aware controls.

**Tech Stack:** Python 3 standard library, Bash, HMMER, MAFFT, CD-HIT, unittest.

## Global Constraints

- Runs use a new safe `runs/<run_id>/`; historical results remain untouched.
- Every input/output artifact records path, size, and SHA-256; unavailable values are `pending`.
- All `<200 aa` ePhaZ records remain auditable and are never deleted by length alone.
- HMM/domain evidence is not a phenotype claim.
- Write failing tests before production code; no commit or push without explicit authorization.

### Task 1: Layering contract and tests

**Files:** Create `pipeline/tests/test_split_ephaz.py`.

- [ ] Add synthetic FASTA/manifest fixtures covering core, broad, missing evidence, architecture pending, and short records.
- [ ] Assert disjoint/covering layers, all short records retained, short records excluded from core, evidence gating, fail-closed malformed inputs, and SHA-256 outputs.
- [ ] Add builder/calibration interface contract tests using shell text and synthetic calibration output.
- [ ] Run `python -m unittest pipeline.tests.test_split_ephaz -v`; expected initial failures because the splitter/interfaces do not exist.

### Task 2: Splitter implementation

**Files:** Create `pipeline/scripts/split_ephaz_seeds.py`.

- [ ] Implement accession-aware FASTA parsing, manifest binding, explicit evidence/architecture/length rules, short review TSV, layer manifest, and checksums.
- [ ] Fail closed on duplicates, missing sequence/manifest/evidence fields, invalid amino acids, and unsafe output paths.
- [ ] Run Task 1 tests and `python -m compileall -q pipeline/scripts/split_ephaz_seeds.py`.

### Task 3: HMM builder interface

**Files:** Modify `pipeline/scripts/04b_build_hmms_v2.sh`; add focused shell contract tests if needed.

- [ ] Add `--families` and `--family-list` support while preserving the existing nine-family default.
- [ ] Permit `ePhaZ_curated_core` and `ePhaZ_broad_discovery` FASTA names and emit matching alignment/HMM names.
- [ ] Refuse duplicate family names, missing inputs, and unsafe output locations.

### Task 4: Layer calibration

**Files:** Create `pipeline/scripts/calibrate_ephaz_layers.py`; add `pipeline/tests/test_calibrate_ephaz_layers.py`.

- [ ] Parameterize HMM/control/output paths and layer-specific positive-control selections.
- [ ] Emit layer-specific summary/hit tables with metrics, coverage, evidence level, and SHA-256/tool metadata.
- [ ] Ensure no fallback to legacy unsplit `ePhaZ.hmm`.

### Task 5: Dated T141 run

**Files:** New runtime directories only: `runs/20260828_ephaz_split_01/` and `deploy/20260828_ephaz_split_01/`.

- [ ] Bind current source and copy `_07` seed inputs read-only into the new run.
- [ ] Run splitter, two alignments/HMMs, and two calibrations on T141 through dated deploy scripts.
- [ ] Record manifests, checksums, logs, tool versions, and an explicit no-full-screening status.

### Task 6: Verification and handoff

- [ ] Run `python -m unittest discover -s pipeline/tests -v`, `python -m compileall -q pipeline/scripts pipeline/tests`, and `git diff --check`.
- [ ] Verify remote HMMs with `hmmstat`, layer disjointness/coverage, and all 61 short records.
- [ ] Report whether calibration supports a subsequent full screening rerun; do not launch it in this task.
