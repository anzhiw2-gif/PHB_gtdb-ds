# PhaDED Evidence Model and Candidate Catalog Redesign Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Preserve the frozen 2026-09-21 results while building a versioned, high-recall PhaDED sequence catalog that separates discovery homology, sequence-family validation, experimental evidence, and functional calibration.

**Architecture:** Introduce a shared evidence ontology, extend the reference ledger, allow quality-controlled annotation-only records to train discovery HMMs, independently validate sequence-family models, and generate a new candidate catalog with one primary disposition per protein. Existing 40,690 outputs remain immutable comparison inputs; every recalculation writes a new auditable run and is executed only after the relevant authorization gate.

**Tech Stack:** Python 3.12, standard-library `csv`/`json`/`hashlib`, HMMER 3.4, MAFFT 7.525, SignalP 6, pandas/numpy/scipy for downstream statistics, `unittest`, Bash entrypoints, GitHub Actions.

## Global Constraints

- Never overwrite historical `results/`, `runs/`, or server directories; every execution uses a new auditable `runs/<run_id>/` with `logs/`, `inputs/`, `results/`, and `input_contract.json`.
- Missing paths, versions, sizes, or SHA-256 values are recorded as `pending`; hashes are never invented.
- Server work executes only from dated `deploy/<run_id>/` bound to the exact source snapshot.
- Server task threads are calculated immediately before launch as `min(40, total_cores - one_minute_load - 10)`; abort when the result is below 1.
- HMM, domain, SignalP, motif, structure, phylogeny, classification, and growth outputs remain candidate-only and never prove phenotype.
- The fixed discovery layer name is `discovery_hmm_uncalibrated`; this layer only recalls and scores, never filters, deletes, demotes, or makes a family call.
- Discovery profiles never enter `pipeline/config/formal_scan_models.tsv`; any GTDB full-library rescan requires separate authorization.
- `DED_hfam_2` pool-external hits remain in `runs/20260919_phaded_with_lipase_deferred_archive_01` as a preserved deferred layer and never enter candidate counts without structure evidence plus downstream gates.
- Every profile binds its exact training-set, alignment, and HMM SHA-256. MAFFT output is treated as non-bit-reproducible and its concrete alignment hash is authoritative.
- Code changes follow failing-test-first development. Every task ends with focused tests, applicable suite, `compileall`, and `git diff --check`.
- Local commits are allowed only with a handoff entry. Push, recomputation, installation, configuration changes, deletion, and server execution require explicit authorization.
- Existing user changes and untracked research files are preserved. No bulk staging, deletion, reset, or cleanup is allowed.

## System Boundaries and File Map

**New shared components**

- `pipeline/scripts/phaded_evidence_schema.py` — controlled vocabularies and row validators.
- `pipeline/scripts/evaluate_phaded_sequence_models.py` — sequence-only held-out/confounder evaluation and promotion to `sequence_family_hmm_validated`.
- `pipeline/scripts/build_phaded_candidate_catalog_v2.py` — produces one-row-per-protein v2 catalog and primary dispositions.
- `pipeline/scripts/reconcile_phaded_candidate_flow.py` — closes the 616, 52, 4,301/6,214, and other count transitions.
- `pipeline/scripts/compare_phaded_catalog_versions.py` — v1/v2 impact and sensitivity reports.
- `pipeline/scripts/server_resources.py` — dynamic thread calculation and auditable resource record.
- `pipeline/config/phaded_sequence_model_gates.tsv` — preregistered sequence-model gates, separate from functional calibration.

**Existing components to modify**

- `pipeline/scripts/validate_phaded_reference_ledger.py`
- `pipeline/scripts/build_phaded_profiles.py`
- `pipeline/scripts/finalize_phaded_reference_panel_calibration.py`
- `pipeline/scripts/filter_phaded_high_confidence.py` — compatibility path only; old output semantics remain frozen.
- `pipeline/scripts/06_screen.sh`
- `pipeline/scripts/run_pipeline.sh`
- `pipeline/config/params.yaml`
- `deploy/20260920_phaded_grodon_growth_01/scripts/build_grodon_manifest.py` — logic is ported into a new dated deploy during authorized execution; the historical deploy remains untouched.
- `.github/workflows/ci.yml`
- `README.md`, `CHANGELOG.md`, `AGENTS.md`, and the new status/handoff document.

**New tests**

- `pipeline/tests/test_phaded_evidence_schema.py`
- `pipeline/tests/test_evaluate_phaded_sequence_models.py`
- `pipeline/tests/test_build_phaded_candidate_catalog_v2.py`
- `pipeline/tests/test_reconcile_phaded_candidate_flow.py`
- `pipeline/tests/test_compare_phaded_catalog_versions.py`
- `pipeline/tests/test_server_resources.py`
- `pipeline/tests/test_hmmer_database_size.py`

## Authorization Checkpoints

The implementation may proceed through code, tests, documentation, and local synthetic fixtures. Stop and request explicit authorization before each of the following:

1. Curating and committing the existing untracked 2026-09-21 research snapshot.
2. Creating any real analysis run or dated deploy.
3. Rebuilding MAFFT alignments or HMMs from the 723-reference collection.
4. Running SignalP, Foldseek, structure prediction, HMMER against GTDB shards, or gRodon.
5. Recomputing the candidate catalog from real project data.
6. Installing/changing environments, modifying server configuration, deleting/archiving data, or pushing to GitHub.

---

### Task 1: Freeze the authority map before editing scientific logic

**Files:**
- Create: `docs/T141_20260928_phaded_evidence_model_redesign_status.md`
- Modify: `docs/T141_20260917_project_handoff.md`
- Inspect only: `docs/T141_20260921_phaded_comprehensive_review.md`

**Interfaces:**
- Consumes: current `main` commit, `git status`, the frozen v1 run paths and their manifests.
- Produces: a reviewed list of files that constitute the 2026-09-21 research snapshot and a source-of-truth table for v1 versus v2.

- [ ] **Step 1: Record the immutable v1 authorities**

Write a status table containing these exact authorities:

```markdown
| Topic | Frozen authority | Status |
|---|---|---|
| candidate universe | runs/20260905_ephaz_iphaz_layering_01/results/layers_retry_02/protein_layers.tsv | frozen-v1 |
| pool-in high confidence | runs/20260919_phaded_nucleophile_unify_01/results/high_confidence_candidates.tsv | frozen-v1 |
| final v1 merge | runs/20260920_phaded_high_confidence_pool_merge_01/results/merged_high_confidence_candidates.tsv | frozen-v1 |
| with-lipase deferred | runs/20260919_phaded_with_lipase_deferred_archive_01 | frozen-deferred |
| v2 candidate catalog | not created | pending-authorization |
```

- [ ] **Step 2: Inventory tracked and untracked research files without staging**

Run:

```powershell
git status --short
git ls-files pipeline/scripts pipeline/tests docs
git ls-files --others --exclude-standard pipeline/scripts pipeline/tests docs
```

Expected: a human-reviewable list; no index change.

- [ ] **Step 3: Compute and record hashes for the three v1 authority files**

Run:

```powershell
Get-FileHash -Algorithm SHA256 'runs\20260905_ephaz_iphaz_layering_01\results\layers_retry_02\protein_layers.tsv'
Get-FileHash -Algorithm SHA256 'runs\20260919_phaded_nucleophile_unify_01\results\high_confidence_candidates.tsv'
Get-FileHash -Algorithm SHA256 'runs\20260920_phaded_high_confidence_pool_merge_01\results\merged_high_confidence_candidates.tsv'
```

Expected: three non-empty SHA-256 values copied verbatim into the status document.

- [ ] **Step 4: Stop for explicit snapshot-commit authorization**

Present the exact proposed staging list. Do not run `git add -A`. If authorized, stage only reviewed files and add a handoff entry containing the commit hash.

- [ ] **Step 5: Verify documentation-only changes**

Run:

```powershell
git diff --check
```

Expected: exit 0.

---

### Task 2: Enforce dynamic server resources and complete run contracts

**Files:**
- Create: `pipeline/scripts/server_resources.py`
- Create: `pipeline/tests/test_server_resources.py`
- Modify: `pipeline/scripts/run_pipeline.sh`
- Modify: `pipeline/scripts/05_predict_proteins.sh`
- Modify: `pipeline/scripts/06_screen.sh`
- Modify: `pipeline/config/params.yaml`
- Modify: `pipeline/tests/test_run_isolation.py`

**Interfaces:**
- Consumes: integer `total_cores`, float `load_1m`, reserve 10, hard cap 40.
- Produces: `compute_thread_limit(total_cores: int, load_1m: float, reserve: int = 10, hard_cap: int = 40) -> int` and a JSON resource record.

- [ ] **Step 1: Write failing resource-limit tests**

```python
def test_compute_thread_limit_reserves_ten_and_caps_at_forty():
    assert module.compute_thread_limit(80, 12.2) == 40
    assert module.compute_thread_limit(32, 10.1) == 11

def test_compute_thread_limit_aborts_when_no_capacity_remains():
    with self.assertRaisesRegex(RuntimeError, "insufficient free capacity"):
        module.compute_thread_limit(16, 6.2)
```

Run:

```powershell
python -m unittest pipeline.tests.test_server_resources -v
```

Expected: FAIL because `server_resources.py` does not exist.

- [ ] **Step 2: Implement the resource calculation**

```python
def compute_thread_limit(total_cores: int, load_1m: float, reserve: int = 10, hard_cap: int = 40) -> int:
    if total_cores < 1 or load_1m < 0 or reserve < 0 or hard_cap < 1:
        raise ValueError("invalid resource measurement")
    available = int(total_cores - load_1m - reserve)
    if available < 1:
        raise RuntimeError("insufficient free capacity after reserve")
    return min(hard_cap, available)
```

The CLI reads `nproc` and `/proc/loadavg`, prints JSON with measurements, and exits non-zero when no capacity remains.

- [ ] **Step 3: Replace every 70-thread default**

`run_pipeline.sh` must call the helper immediately before prediction and screening and pass the returned value. Standalone stage scripts default to 40 but reject values above the freshly measured limit when `--server` is active. `params.yaml` changes both 70 values to 40 and documents that 40 is a ceiling, not a fixed launch value.

- [ ] **Step 4: Extend isolation tests**

Add assertions that the driver contains `server_resources.py`, `nproc`, `/proc/loadavg`, and no `THREADS_PREDICT=70`, `THREADS_SCREEN=70`, or `THREADS=70` default.

- [ ] **Step 5: Run verification and commit**

```powershell
python -m unittest pipeline.tests.test_server_resources pipeline.tests.test_run_isolation -v
python -m compileall -q pipeline/scripts pipeline/tests
git diff --check
git add pipeline/scripts/server_resources.py pipeline/tests/test_server_resources.py pipeline/scripts/run_pipeline.sh pipeline/scripts/05_predict_proteins.sh pipeline/scripts/06_screen.sh pipeline/config/params.yaml pipeline/tests/test_run_isolation.py docs/T141_20260917_project_handoff.md
git commit -m "fix: enforce dynamic server thread limits"
```

Expected: tests pass; commit hash is appended to the handoff document.

---

### Task 3: Add the shared evidence ontology

**Files:**
- Create: `pipeline/scripts/phaded_evidence_schema.py`
- Create: `pipeline/tests/test_phaded_evidence_schema.py`
- Create: `pipeline/config/phaded_sequence_model_gates.tsv`

**Interfaces:**
- Consumes: strings from reference ledgers, profile manifests, and candidate evidence tables.
- Produces: `validate_reference_evidence(row)`, `validate_candidate_evidence(row)`, and controlled vocabulary constants used by later tasks.

- [ ] **Step 1: Write failing schema tests**

```python
def test_ahsmg_is_a_ser_motif_not_a_nucleophile_identity():
    row = module.normalize_nucleophile("extracellular native-SCL/PhaZ7-like", "", "supported")
    self.assertEqual(row, {"nucleophile_identity": "ser", "motif_class": "AHSMG"})

def test_signalp_does_not_become_localization_truth():
    row = module.normalize_localization("OTHER", "historical_family_label_only")
    self.assertEqual(row["transport_signal_prediction"], "OTHER")
    self.assertEqual(row["localization_evidence"], "historical_family_label_only")

def test_model_layers_are_disjoint_from_functional_status():
    module.validate_model_state("sequence_family_hmm_validated", "not_function_calibrated")
```

Run:

```powershell
python -m unittest pipeline.tests.test_phaded_evidence_schema -v
```

Expected: FAIL because the module is missing.

- [ ] **Step 2: Implement controlled vocabularies**

```python
MODEL_LAYERS = {
    "reference_query_only",
    "discovery_hmm_uncalibrated",
    "sequence_family_hmm_validated",
    "calibrated_candidate_model",
}
NUCLEOPHILE_IDENTITIES = {"ser", "cys", "unresolved"}
MOTIF_CLASSES = {"GxSxG", "AHSMG", "Cys-associated", "other", "unresolved"}
EVIDENCE_GRADES = {"E3", "E2", "E1", "A"}
PRIMARY_DISPOSITIONS = {
    "core_sequence_homolog",
    "probable_sequence_homolog",
    "remote_homolog_candidate",
    "function_unresolved",
    "deferred_structure_review",
    "excluded_input_quality",
}
```

Validators reject `nucleophile_identity=ahsmg`, localization inferred solely from SignalP, and discovery-layer rows with a family call or exclusion action.

- [ ] **Step 3: Add the sequence-model gate configuration**

Create a TSV with columns:

```text
gate_version	min_training_sequences	min_heldout_sequences	min_heldout_recall	max_unexplained_confounder_hits	require_unique_best_model	require_alignment_hash
sequence-gate-v1	3	1	0.8	0	true	true
```

For held-out panels smaller than five, all held-out sequences must be recovered; this rule is encoded in the evaluator, not inferred from `0.8` rounding.

- [ ] **Step 4: Verify and commit**

```powershell
python -m unittest pipeline.tests.test_phaded_evidence_schema -v
python -m compileall -q pipeline/scripts pipeline/tests
git diff --check
git add pipeline/scripts/phaded_evidence_schema.py pipeline/tests/test_phaded_evidence_schema.py pipeline/config/phaded_sequence_model_gates.tsv docs/T141_20260917_project_handoff.md
git commit -m "feat: define PhaDED evidence ontology"
```

---

### Task 4: Extend the 723-reference ledger with evidence quality and independence

**Files:**
- Modify: `pipeline/scripts/validate_phaded_reference_ledger.py`
- Modify: `pipeline/tests/test_validate_phaded_reference_ledger.py`
- Create: `pipeline/config/phaded_reference_evidence_columns.tsv`

**Interfaces:**
- Consumes: the existing reference ledger plus new evidence columns.
- Produces: a validated ledger in which training eligibility, experimental quality, and independence are explicit and machine-checkable.

- [ ] **Step 1: Write failing evidence-grade tests**

```python
def test_direct_experimental_record_requires_complete_provenance(self):
    row = valid_row(evidence_status="experimental_positive", experimental_evidence_grade="E3")
    row["study_id"] = ""
    with self.assertRaisesRegex(ValueError, "E3.*study_id"):
        module.validate_reference_row(row)

def test_annotation_only_record_can_be_discovery_training_eligible(self):
    row = valid_row(evidence_status="annotation_only", experimental_evidence_grade="A")
    row["discovery_training_eligible"] = "true"
    module.validate_reference_row(row)
```

- [ ] **Step 2: Require the new columns**

The config and validator require:

```text
accession_version
gene_id
study_id
experiment_unit_id
independence_group
experimental_evidence_grade
assay_directness
experimental_substrate_class
experimental_system
negative_control_description
catalytic_mutant_evidence
localization_evidence
sequence_integrity
discovery_training_eligible
functional_calibration_eligible
```

E3/E2 require study, experiment unit, substrate, system, direct result, and evidence source. E1 cannot be counted as a direct functional positive. A may train discovery models only when sequence integrity and family consistency pass.

- [ ] **Step 3: Add duplicate-independence checks**

The validator reports counts by accession, sequence SHA, gene, experiment unit, study, genus, and independence group. Functional positive counts use unique `independence_group`, never row count.

- [ ] **Step 4: Run synthetic tests**

```powershell
python -m unittest pipeline.tests.test_validate_phaded_reference_ledger -v
python -m compileall -q pipeline/scripts pipeline/tests
git diff --check
```

Expected: all tests pass without modifying the real 723-row ledger.

- [ ] **Step 5: Commit code and schema**

```powershell
git add pipeline/scripts/validate_phaded_reference_ledger.py pipeline/tests/test_validate_phaded_reference_ledger.py pipeline/config/phaded_reference_evidence_columns.tsv docs/T141_20260917_project_handoff.md
git commit -m "feat: grade PhaDED experimental evidence"
```

- [ ] **Step 6: Stop for authorization before curating the real ledger**

Real accession/literature curation is a separate evidence-acquisition run. No missing field is inferred; unavailable values are `pending`.

---

### Task 5: Separate HMM construction from functional calibration

**Files:**
- Modify: `pipeline/scripts/build_phaded_profiles.py`
- Modify: `pipeline/tests/test_build_phaded_profiles.py`
- Modify: `pipeline/scripts/finalize_phaded_reference_panel_calibration.py`
- Modify: `pipeline/tests/test_finalize_phaded_reference_panel_calibration.py`

**Interfaces:**
- Consumes: validated reference ledger and FASTA.
- Produces: profile manifest with `model_layer`, training composition, training-set/alignment/HMM hashes, and independent `functional_calibration_status`.

- [ ] **Step 1: Replace the “three experimental positives to build” test**

```python
def test_annotation_only_family_builds_discovery_hmm_when_three_eligible_sequences_exist(self):
    rows = [annotation_row("A1"), annotation_row("A2"), annotation_row("A3")]
    result = build_with_fake_tools(rows)
    family = family_profile(result)
    self.assertEqual(family["model_layer"], "discovery_hmm_uncalibrated")
    self.assertEqual(family["experimental_anchor_count"], 0)
    self.assertEqual(family["annotation_only_count"], 3)

def test_two_sequence_family_remains_reference_query_only(self):
    result = build_with_fake_tools([annotation_row("A1"), annotation_row("A2")])
    self.assertEqual(family_profile(result)["model_layer"], "reference_query_only")
```

Expected initial result: FAIL because current `_profile_record` trains only with three experimental positives.

- [ ] **Step 2: Implement deterministic training selection**

Add:

```python
def select_discovery_training_rows(records: list[dict[str, str]]) -> list[dict[str, str]]:
    selected = [row for row in records if row["discovery_training_eligible"] == "true"]
    return sorted(selected, key=lambda row: (row["sequence_sha256"], row["accession"]))
```

Build an HMM when at least three unique eligible sequences exist. One or two sequences remain `reference_query_only`. A contradictory direct experimental negative blocks functional calibration but does not silently remove the sequence from the reference audit.

- [ ] **Step 3: Expand profile manifest provenance**

Add exact fields:

```text
model_layer
functional_calibration_status
training_sequence_count
experimental_anchor_count
annotation_only_count
training_accessions
training_set_sha256
alignment_sha256
hmm_sha256
mafft_version
hmmer_version
bit_reproducible
```

`bit_reproducible` is always `false` for MAFFT-built alignments; the actual alignment hash is mandatory.

- [ ] **Step 4: Preserve the functional gate but count qualified independent evidence**

Modify the readiness input to include `qualified_e2_e3_positive_count` and `independent_genus_count`. The family gate becomes:

```python
return (
    kind == "family"
    and qualified_e2_e3_positive_count >= 3
    and independent_genus_count >= 3
    and heldout >= 1
    and negative >= 1
    and challenge >= 1
    and family_resolved_negative_status == "sufficient"
    and challenge_status == "sufficient"
)
```

The output status for a passed but unpromoted profile is `candidate_gate_passed_not_promoted`; `calibrated_candidate_model` is emitted only by an explicitly authorized finalize action.

- [ ] **Step 5: Run focused verification**

```powershell
python -m unittest pipeline.tests.test_build_phaded_profiles pipeline.tests.test_finalize_phaded_reference_panel_calibration -v
python -m compileall -q pipeline/scripts pipeline/tests
git diff --check
```

- [ ] **Step 6: Commit without building real models**

```powershell
git add pipeline/scripts/build_phaded_profiles.py pipeline/tests/test_build_phaded_profiles.py pipeline/scripts/finalize_phaded_reference_panel_calibration.py pipeline/tests/test_finalize_phaded_reference_panel_calibration.py docs/T141_20260917_project_handoff.md
git commit -m "feat: separate discovery HMMs from functional calibration"
```

---

### Task 6: Independently validate sequence-family models

**Files:**
- Create: `pipeline/scripts/evaluate_phaded_sequence_models.py`
- Create: `pipeline/tests/test_evaluate_phaded_sequence_models.py`
- Consume: `pipeline/config/phaded_sequence_model_gates.tsv`

**Interfaces:**
- Consumes: profile manifest, held-out assignments, nearest-family/competing-fold challenge scores, and preregistered gates.
- Produces: `sequence_model_evaluation.tsv`, `sequence_model_gate_detail.tsv`, and a manifest. It never changes functional calibration status.

- [ ] **Step 1: Write failing gate tests**

```python
def test_sequence_gate_does_not_require_experimental_positive(self):
    result = module.evaluate_model(record(
        training_sequence_count=8,
        heldout_count=5,
        heldout_recovered=5,
        unexplained_confounder_hits=0,
        unique_best_model=True,
        alignment_sha256="a" * 64,
        experimental_anchor_count=0,
    ), gates())
    self.assertEqual(result["model_layer"], "sequence_family_hmm_validated")
    self.assertEqual(result["functional_calibration_status"], "not_function_calibrated")

def test_effn_below_one_is_diagnostic_not_an_automatic_failure():
    item = record(eff_nseq=0.7)
    self.assertNotIn("eff_nseq_below_one", module.blocking_reasons(item, gates()))
```

- [ ] **Step 2: Implement evaluation metrics**

Compute held-out recall, exact unexplained confounder hit count, unique-best-model status, score margin, coverage, family assignment stability, and `EFFN` as a diagnostic field. For fewer than five held-out records, require all to be recovered; otherwise use the configured minimum recall.

- [ ] **Step 3: Prevent functional overclaiming**

Every output contains:

```python
{
    "sequence_claim": "sequence-defined family homology",
    "phenotype_boundary": "does not establish PHB/PHA degradation",
    "functional_calibration_status": input_row["functional_calibration_status"],
}
```

- [ ] **Step 4: Run and commit**

```powershell
python -m unittest pipeline.tests.test_evaluate_phaded_sequence_models -v
python -m compileall -q pipeline/scripts pipeline/tests
git diff --check
git add pipeline/scripts/evaluate_phaded_sequence_models.py pipeline/tests/test_evaluate_phaded_sequence_models.py pipeline/config/phaded_sequence_model_gates.tsv docs/T141_20260917_project_handoff.md
git commit -m "feat: validate PhaDED sequence family models"
```

---

### Task 7: Build a v2 candidate catalog with corrected field semantics

**Files:**
- Create: `pipeline/scripts/build_phaded_candidate_catalog_v2.py`
- Create: `pipeline/tests/test_build_phaded_candidate_catalog_v2.py`
- Modify: `pipeline/scripts/filter_phaded_high_confidence.py` only to add a deprecation/boundary note; do not change frozen-v1 behavior.

**Interfaces:**
- Consumes: candidate universe, profile scores/evaluations, domain geometry, motifs, SignalP predictions, localization evidence, InterPro/Pfam provenance, competition flags, and explicit demotion/deferred sets.
- Produces: `phaded_candidate_catalog_v2.tsv`, `catalog_summary.json`, and `catalog_manifest.json`.

- [ ] **Step 1: Write failing semantic tests**

```python
def test_ahsmg_candidate_is_ser_with_ahsmg_motif():
    row = build_row(superfamily=PHAZ7, ahsmg_state="supported")
    out = module.classify_row(row)
    self.assertEqual(out["nucleophile_identity"], "ser")
    self.assertEqual(out["motif_class"], "AHSMG")

def test_sbd_does_not_define_type_one_or_type_two():
    row = build_row(catalytic_domain_type="type2_verified", sbd_state="supported")
    out = module.classify_row(row)
    self.assertEqual(out["catalytic_domain_type"], "type2_verified")
    self.assertEqual(out["accessory_domain_architecture"], "SBD:supported")

def test_signalp_other_is_not_equivalent_to_experimental_intracellular():
    out = module.classify_row(build_row(signalp_class="OTHER", localization_evidence="unknown"))
    self.assertEqual(out["transport_signal_prediction"], "OTHER")
    self.assertEqual(out["localization_evidence"], "unknown")
```

- [ ] **Step 2: Encode special-family dispositions**

```python
def primary_disposition(row: dict[str, str]) -> str:
    if row["superfamily"] == "intracellular nPHASCL with lipase box":
        return "deferred_structure_review"
    if row["superfamily"] == "intracellular nPHAMCL":
        return "function_unresolved"
    if row["model_layer"] == "sequence_family_hmm_validated" and row["assignment_unique"] == "true":
        return "core_sequence_homolog"
    if row["model_layer"] in {"sequence_family_hmm_validated", "discovery_hmm_uncalibrated"}:
        return "probable_sequence_homolog"
    return "remote_homolog_candidate"
```

The complete implementation also handles input-quality exclusion and never routes scientific uncertainty to `excluded_input_quality`.

- [ ] **Step 3: Make localization a compatibility annotation**

SignalP discordance adds `localization_discordant` to `evidence_flags`; it does not by itself delete a sequence homolog. A separate functional view may require localization support, but that view cannot change the master catalog disposition.

- [ ] **Step 4: Wire the 52 Cys-export records deterministically**

The explicit demotion input adds `localization_discordant` and prevents these rows from entering a functional-support view. Tests assert all 52 input IDs are accounted for exactly once when the real run is authorized.

- [ ] **Step 5: Run synthetic tests and commit**

```powershell
python -m unittest pipeline.tests.test_build_phaded_candidate_catalog_v2 pipeline.tests.test_phaded_evidence_schema -v
python -m compileall -q pipeline/scripts pipeline/tests
git diff --check
git add pipeline/scripts/build_phaded_candidate_catalog_v2.py pipeline/tests/test_build_phaded_candidate_catalog_v2.py pipeline/scripts/filter_phaded_high_confidence.py docs/T141_20260917_project_handoff.md
git commit -m "feat: add evidence-separated PhaDED candidate catalog"
```

---

### Task 8: Close the candidate flow ledger and all count discrepancies

**Files:**
- Create: `pipeline/scripts/reconcile_phaded_candidate_flow.py`
- Create: `pipeline/tests/test_reconcile_phaded_candidate_flow.py`

**Interfaces:**
- Consumes: candidate universe and every disposition/hold/demotion/deferred input with accession keys.
- Produces: `candidate_flow.tsv`, `candidate_flow_summary.json`, and `overlap_matrix.tsv`.

- [ ] **Step 1: Write failing partition tests**

```python
def test_every_candidate_has_one_primary_disposition():
    rows, summary = module.reconcile(universe={"A", "B", "C"}, records=[
        record("A", "core_sequence_homolog"),
        record("B", "function_unresolved"),
        record("C", "remote_homolog_candidate"),
    ])
    self.assertEqual(summary["unaccounted"], 0)
    self.assertEqual(summary["multiply_disposed"], 0)

def test_overlapping_reason_sets_are_not_added_as_unique_proteins():
    summary = module.summarize_reason_sets({"PhaC": {"A", "B"}, "iPhaZ": {"B", "C"}})
    self.assertEqual(summary["union_count"], 3)
    self.assertEqual(summary["sum_of_reason_counts"], 4)
```

- [ ] **Step 2: Implement exact accounting**

The summary must separately report:

```text
universe_count
primary_disposition_counts
unaccounted
multiply_disposed
demoted_from_v1
hold_reason_counts_nonexclusive
hold_reason_union_count
```

The implementation exits non-zero when `unaccounted` or `multiply_disposed` is non-zero.

- [ ] **Step 3: Add named reconciliation checks**

Machine-readable checks explain:

- the 616 records in `109087 - 36611 - 71860`;
- whether the 52 demotions join hold or form a separate disposition;
- why PhaC 1,970 plus iPhaZ 4,244 relates to a 4,301-protein union;
- the 66-record difference between 4,507 eligible positives and 4,441 manifest positives;
- motif run coverage versus motif resolvability;
- exact counts of superfamily profiles, family definitions, trained HMM files, and reference-only entries.

- [ ] **Step 4: Run and commit**

```powershell
python -m unittest pipeline.tests.test_reconcile_phaded_candidate_flow -v
python -m compileall -q pipeline/scripts pipeline/tests
git diff --check
git add pipeline/scripts/reconcile_phaded_candidate_flow.py pipeline/tests/test_reconcile_phaded_candidate_flow.py docs/T141_20260917_project_handoff.md
git commit -m "feat: close PhaDED candidate flow accounting"
```

---

### Task 9: Audit and normalize HMMER shard database size

**Files:**
- Create: `pipeline/tests/test_hmmer_database_size.py`
- Modify: `pipeline/scripts/06_screen.sh`
- Modify: `pipeline/scripts/06_validate_screen_manifest.py`
- Modify: `pipeline/tests/test_formal_scan13_tier_run.py`
- Reuse: `pipeline/scripts/parse_phaded_cys_targeted_recall.py::rescale_evalue`

**Interfaces:**
- Consumes: verified total target sequence count and per-shard target counts.
- Produces: shard commands with a common `-Z`, plus manifest fields proving the scale.

- [ ] **Step 1: Audit scan-13 commands and logs read-only**

Search the frozen run/deploy evidence for `hmmsearch`, `-Z`, `--domZ`, `tblout`, `domtblout`, shard sequence counts, and raw bitscores. Record one of two outcomes in the status document:

```text
normalizable_from_frozen_scores
requires_authorized_rescan
```

Do not infer the outcome from current source code; use the frozen executed command evidence.

- [ ] **Step 2: Write a failing command test**

```python
def test_every_shard_uses_the_same_full_database_z():
    command = module.build_hmmsearch_command("model.hmm", "shard.faa", total_targets=292000000)
    self.assertIn("-Z", command)
    self.assertEqual(command[command.index("-Z") + 1], "292000000")
```

- [ ] **Step 3: Modify the screening entrypoint**

`06_screen.sh` accepts `--database-size-z`, validates it as a positive integer, and passes the same value to every HMMER call:

```bash
hmmsearch --tblout "$out" --domtblout "${out%.tbl}.dom" \
  -Z "$DATABASE_SIZE_Z" -E "$EVAL" --cpu 1 "$hmm" "$shard"
```

Do not add `--domZ` unless a separate verified decision shows that downstream logic uses a domain E-value requiring that scale.

- [ ] **Step 4: Bind the scale in the manifest**

The screen manifest records `database_size_Z`, `database_size_basis`, every shard sequence count, the total, and the exact command template. Validation fails if shard runs use different Z values.

- [ ] **Step 5: Verify and commit**

```powershell
python -m unittest pipeline.tests.test_hmmer_database_size pipeline.tests.test_formal_scan13_tier_run -v
python -m compileall -q pipeline/scripts pipeline/tests
git diff --check
git add pipeline/tests/test_hmmer_database_size.py pipeline/scripts/06_screen.sh pipeline/scripts/06_validate_screen_manifest.py pipeline/tests/test_formal_scan13_tier_run.py docs/T141_20260928_phaded_evidence_model_redesign_status.md docs/T141_20260917_project_handoff.md
git commit -m "fix: bind HMMER shard searches to full database size"
```

- [ ] **Step 6: Stop at the execution decision**

If frozen scores are normalizable, propose a new reconciliation run with bound inputs and no HMM search. If not, request explicit authorization for a new dated deploy and GTDB rescan; do not silently mix old and new E-values.

---

### Task 10: Create v1/v2 impact and rule-sensitivity reports

**Files:**
- Create: `pipeline/scripts/compare_phaded_catalog_versions.py`
- Create: `pipeline/tests/test_compare_phaded_catalog_versions.py`

**Interfaces:**
- Consumes: frozen v1 40,690 table, v2 catalog, and scenario definitions.
- Produces: accession-level changes, per-family count deltas, reason summaries, and localization-sensitivity scenarios.

- [ ] **Step 1: Write failing comparison tests**

```python
def test_comparison_reports_added_removed_and_relabelled():
    result = module.compare(
        v1={"A": "high_confidence", "B": "high_confidence"},
        v2={"B": "function_unresolved", "C": "core_sequence_homolog"},
    )
    self.assertEqual(result["removed"], {"A"})
    self.assertEqual(result["added"], {"C"})
    self.assertEqual(result["relabelled"], {"B"})
```

- [ ] **Step 2: Implement three preregistered scenarios**

```text
current_v1_rules
localization_uncertain_not_excluded
sequence_validated_models_only
```

For each scenario report total proteins, genomes, superfamily/family composition, transport-signal composition, and accession-level differences. Report 76:24 only as candidate composition under the named rule.

- [ ] **Step 3: Verify and commit**

```powershell
python -m unittest pipeline.tests.test_compare_phaded_catalog_versions -v
python -m compileall -q pipeline/scripts pipeline/tests
git diff --check
git add pipeline/scripts/compare_phaded_catalog_versions.py pipeline/tests/test_compare_phaded_catalog_versions.py docs/T141_20260917_project_handoff.md
git commit -m "feat: compare PhaDED catalog rule scenarios"
```

---

### Task 11: Prepare targeted evidence runs in scientific priority order

**Files:**
- Modify after authorization: `pipeline/scripts/run_phaded_full_signalp.py`
- Modify after authorization: `pipeline/scripts/complete_phaded_reference_residue_mapping.py`
- Modify after authorization: `pipeline/scripts/prepare_phaded_structure_review.py`
- Create after authorization: `deploy/20260928_phaded_reference_signalp_v2_01/`
- Create after authorization: `runs/20260928_phaded_reference_signalp_v2_01/{logs,inputs,results}` and `input_contract.json`
- Create after authorization: `deploy/20260928_phaded_cys_mapping_v2_01/`
- Create after authorization: `runs/20260928_phaded_cys_mapping_v2_01/{logs,inputs,results}` and `input_contract.json`
- Create after authorization: `deploy/20260928_phaded_structure_review_v2_01/`
- Create after authorization: `runs/20260928_phaded_structure_review_v2_01/{logs,inputs,results}` and `input_contract.json`

**Interfaces:**
- Consumes: validated 723-reference ledger, v2 catalog, exact target accession lists.
- Produces: bounded evidence tables that annotate v2 without deleting candidates.

- [ ] **Step 1: Prepare authorization packets without executing tools**

Prepare exact target lists and anticipated outputs for:

1. 723-reference SignalP calibration, stratified by experimental localization, functional-only evidence, and inherited historical label;
2. 11,610 SP-less type-1-like candidates, using catalytic geometry and independent accessory-domain evidence;
3. 29,974 Cys candidates, using reference-anchored residue mapping with explicit substitution/truncation/uncertain states;
4. 987 nPHAMCL-like candidates, with a competition panel rather than a single 8YNV reference;
5. stratified representatives from the with-lipase deferred pool by sequence cluster, taxonomy, length, completeness, score, and architecture.

- [ ] **Step 2: Define success criteria before requesting execution**

Each run must state its decision question, target population, sampling or full-set basis, tool/version, false-positive competitors, output schema, and how results alter only evidence flags/dispositions.

- [ ] **Step 3: Request separate authorization for each compute class**

SignalP, HMM building, Foldseek/structure prediction, and GTDB searching are separate approvals. Approval of one does not authorize the others.

- [ ] **Step 4: After authorization, create the dated run and deploy**

Use the current date plus a descriptive suffix; validate the exact path before creation. Bind source, reference FASTA/ledger, HMMs, environment, commands, threads, sizes, and SHA-256 values. Missing items are `pending`.

- [ ] **Step 5: Verify each bounded run**

Run its focused tests plus:

```powershell
python -m compileall -q pipeline/scripts pipeline/tests
git diff --check
```

On the server, record `nproc`, `/proc/loadavg`, calculated threads, and dated deploy identity before starting.

---

### Task 12: Correct gRodon labels, units, overlap, and dependence

**Files:**
- Create during authorized execution: `deploy/20260928_phaded_grodon_reanalysis_v2_01/scripts/build_grodon_manifest.py`
- Create during authorized execution: `deploy/20260928_phaded_grodon_reanalysis_v2_01/scripts/grodon_group_stats.py`
- Modify: `pipeline/tests/test_build_grodon_manifest.py`
- Modify: `pipeline/tests/test_grodon_group_stats.py`
- Modify: `pipeline/scripts/plot_phaded_grodon_growth_nature.py`
- Modify: `pipeline/tests/test_plot_phaded_grodon_growth_nature.py`

**Interfaces:**
- Consumes: v2 genome-level carrier table, complete final-positive exclusion set, compatible gRodon prediction records.
- Produces: candidate-carrier versus candidate-not-detected analysis with genus-level primary estimates and sensitivity analyses.

- [ ] **Step 1: Write failing overlap and label tests**

```python
def test_final_positive_and_control_sets_are_disjoint():
    with self.assertRaisesRegex(ValueError, "positive/control overlap"):
        module.build_manifest(
            carriers={"G1": carrier_row()},
            exclusion=set(),
            tax_by_gid={"G1": tax(), "G2": tax()},
            seed=42,
            max_per_genus=0,
            file_exists=lambda _: True,
        )

def test_manifest_uses_candidate_carrier_language():
    rows, _ = valid_manifest()
    self.assertEqual({row["candidate_detection_status"] for row in rows}, {
        "candidate_gene_carrier", "candidate_not_detected_under_defined_search"
    })
```

- [ ] **Step 2: Verify the growth-rate formula**

The dated deploy must expose one tested conversion function. If gRodon returns minimum doubling time in hours:

```python
def growth_rate_per_h(doubling_time_h: float) -> float:
    if doubling_time_h <= 0:
        raise ValueError("doubling time must be positive")
    return math.log(2.0) / doubling_time_h
```

If the actual R output is already a growth rate, record that field and do not apply a second conversion. The run manifest states the verified interpretation.

- [ ] **Step 3: Make genus-level mean difference the primary estimand**

Bootstrap and permutation resample genera as the dependence unit. Median, sign test, per-genome analyses, and intracellular/extracellular comparisons are secondary. Replace “no difference” with “no difference detected under the current design” unless an equivalence margin is preregistered.

- [ ] **Step 4: Audit the 713 reused predictions**

Require equal tool version, model/settings, temperature handling, CDS preparation, ribosomal marker procedure, and output formula. Incompatible records are excluded from the primary result and retained in a sensitivity table.

- [ ] **Step 5: Quantify selection coverage**

Compare included, no-control, missing-FASTA, and failed-prediction genomes by taxonomy, candidate family, group, and available quality metrics. Explain the 4,507-to-4,441 difference with exact dispositions.

- [ ] **Step 6: Run tests before requesting gRodon execution**

```powershell
python -m unittest pipeline.tests.test_build_grodon_manifest pipeline.tests.test_grodon_group_stats pipeline.tests.test_plot_phaded_grodon_growth_nature -v
python -m compileall -q pipeline deploy/20260928_phaded_grodon_reanalysis_v2_01
git diff --check
```

Do not create the deploy or execute gRodon until separately authorized.

---

### Task 13: Expand CI and publish a versioned scientific release candidate

**Files:**
- Modify: `.github/workflows/ci.yml`
- Modify: `README.md`
- Modify: `CHANGELOG.md`
- Modify: `AGENTS.md`
- Modify: `docs/T141_20260921_phaded_comprehensive_review.md` only by adding a supersession banner; preserve its historical content.
- Create: `docs/T141_20260928_phaded_evidence_model_redesign_status.md`
- Modify: `docs/T141_20260917_project_handoff.md`

**Interfaces:**
- Consumes: completed tasks, focused test commands, v2 manifests and comparison report.
- Produces: a release-blocking CI gate and one current project authority statement.

- [ ] **Step 1: Expand the offline CI suite**

Run all tracked tests that use synthetic/local fixtures:

```yaml
- name: Run PhaDED governance and evidence tests
  run: >-
    python -m unittest
    pipeline.tests.test_run_context
    pipeline.tests.test_run_isolation
    pipeline.tests.test_server_resources
    pipeline.tests.test_phaded_evidence_schema
    pipeline.tests.test_validate_phaded_reference_ledger
    pipeline.tests.test_build_phaded_profiles
    pipeline.tests.test_evaluate_phaded_sequence_models
    pipeline.tests.test_build_phaded_candidate_catalog_v2
    pipeline.tests.test_reconcile_phaded_candidate_flow
    pipeline.tests.test_hmmer_database_size
    pipeline.tests.test_compare_phaded_catalog_versions
    -v
```

- [ ] **Step 2: Correct current-status language**

README and CHANGELOG must state:

```text
PhaDED outputs are a versioned sequence-homology candidate resource.
Discovery and sequence-family HMMs do not prove PHB/PHA degradation.
hfam_70 passed the project candidate gate but is not finalized, registered, or released as a function-calibrated production model.
```

Replace “million false positives” with “million-scale unresolved hits excluded from main counts pending discrimination.” Replace “degrader/non-degrader” with candidate carrier/not detected under defined search.

- [ ] **Step 3: Mark T141 as historical rather than rewriting history**

Add a top banner pointing to the 2026-09-28 status and design. Do not silently alter the old 40,690 narrative or historical error record.

- [ ] **Step 4: Run the complete local gate**

```powershell
python -m unittest discover -s pipeline/tests -p 'test_*.py' -v
python -m compileall -q pipeline/scripts pipeline/tests
git diff --check
```

Expected: all applicable tests pass. Data/tool-dependent tests are explicitly skipped with a reason, never removed or weakened.

- [ ] **Step 5: Review branch diff**

```powershell
git status --short
git diff --stat
git diff --check
```

Expected: only reviewed source, tests, configuration, and documentation appear; historical runs remain untouched.

- [ ] **Step 6: Commit and record the handoff**

```powershell
git add .github/workflows/ci.yml README.md CHANGELOG.md AGENTS.md docs/T141_20260921_phaded_comprehensive_review.md docs/T141_20260928_phaded_evidence_model_redesign_status.md docs/T141_20260917_project_handoff.md
git commit -m "docs: publish PhaDED evidence model redesign status"
```

Record every implementation commit in the handoff. Do not push until the operator explicitly authorizes it.

## Final Acceptance Review

- [ ] Frozen v1 paths and hashes are recorded and unchanged.
- [ ] The complete candidate universe has exactly one v2 primary disposition per accession.
- [ ] Experimental evidence grade and independence are machine-checkable.
- [ ] Annotation-only records can train discovery HMMs without implying functional calibration.
- [ ] `sequence_family_hmm_validated` is based on held-out/confounder performance, not experimental count.
- [ ] AHSMG is represented as a Ser motif class; SignalP is not localization truth; SBD does not define catalytic type.
- [ ] nPHAMCL and with-lipase remain unresolved/deferred rather than being deleted or called validated functions.
- [ ] The 52 demotions and all count discrepancies are reproducible from accession-level ledgers.
- [ ] HMMER shard E-values use a recorded common database size or are explicitly blocked pending rescan.
- [ ] gRodon labels, overlap assertion, units, reuse compatibility, and genus-level dependence are verified.
- [ ] Every new real run has a dated deploy, input contract, source/environment binding, logs, results, and canonical manifest.
- [ ] Relevant tests, `compileall`, and `git diff --check` pass.
- [ ] README, CHANGELOG, T141 banner, status, handoff, local commit, GitHub, deploy, and run manifests agree before release.
