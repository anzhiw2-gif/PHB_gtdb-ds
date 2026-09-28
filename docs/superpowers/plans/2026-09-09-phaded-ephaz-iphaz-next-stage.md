# PhaDED-guided ePhaZ/iPhaZ Next-stage Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `subagent-driven-development` (recommended) or an execution-plan workflow to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** On the accepted GTDB candidate set, produce a PhaDED-grounded, family-level candidate annotation that is explicit about uncertainty; independently repair and calibrate the four-seed `ePhaZ_curated_core` before it is ever used as a general classifier or as the basis of another GTDB-wide scan.

## Execution Status (2026-09-13)

The checklist below is updated from the actual dated runs and manifests. A
completed calibration task may conclude `reference_only_insufficient_panel`;
that is a completed governance result, not a successful model promotion.

| Task | Status | Authoritative evidence |
|---|---|---|
| Task 0 | completed; historical contract discrepancies reconciled in a new amendment without editing old contracts | `runs/20260909_phaded_reference_evidence_01/`, `runs/20260909_phaded_candidate_mapping_01/`, `runs/20260909_ephaz_core_calibration_01/`, `runs/20260912_phaded_provenance_reconciliation_01/` |
| Task 1 | completed; 38-family vocabulary frozen, strict evidence remains limited | `runs/20260909_phaded_reference_evidence_01/results/reference_validation.json` |
| Task 2 | completed; 9 trained profiles and 37 reference-only records | `runs/20260909_phaded_reference_evidence_01/results/profile_build_01/profile_audit.json` |
| Task 3 | completed; 109,087 candidates reconciled | `runs/20260909_phaded_candidate_mapping_01/results/mapping_reconciliation.json` |
| Task 4 | completed candidate-only; Pfam/InterPro/structure/phylogeny integration and 19-record manual structure/source review are complete | `runs/20260912_phaded_structure_evidence_review_01/`, `runs/20260912_phaded_candidate_architecture_review_01/`, `runs/20260913_phaded_19_candidate_manual_review_01/` |
| Task 5 | completed; all three subtype decisions are `reference_only_insufficient_panel` | `runs/20260909_ephaz_core_calibration_01/results/calibration/` |
| Task 6 | completed candidate-only; candidate ranking, summaries, trees, 19-record manual structure/source review, Figure 4 repair, and experimental validation plan are complete; wet-lab execution remains future work | `runs/20260911_phaded_priority_review_03/`, `runs/20260911_phaded_competition_phylo_review_01/`, `runs/20260913_progress_report_figure4_repair_01/`, `runs/20260913_phaded_experimental_validation_plan_01/`, `runs/20260913_phaded_19_candidate_manual_review_01/` |
| Task 7 | not started by design; registry/rescan authorization gate remains closed | no formal-change plan, run or deploy exists |

The original plan was written before execution and its unchecked boxes are not
by themselves evidence that the corresponding run was not performed.

**Architecture:** Keep the accepted `run-13`, `20260905_ephaz_iphaz_layering_01`, and `20260906_overall_pha_phb_archaea_audit_04` outputs frozen. First construct a versioned PhaDED evidence library, then route only the existing ePhaZ/iPhaZ candidates through the eight PhaDED superfamilies and, where the evidence supports it, the 38 homologous families. Add independent localization, motif, catalytic-residue, domain-architecture, and sequence-integrity evidence before assigning a final candidate-only label. In parallel, use accession-level experimental evidence and held-out controls to decide whether the four-seed ePhaZ model can be replaced by separately calibrated subtype models.

**Tech Stack:** Existing Python standard-library pipeline, `unittest`, HMMER 3.4, MAFFT, SignalP6 results already bound to the accepted candidate layers, independently versioned Pfam/InterPro or equivalent domain annotations when authorized, TSV/FASTA/JSON manifests, and the PhaDED framework of Knoll et al. 2009 (DOI `10.1186/1471-2105-10-89`).

## Global Constraints

- The accepted main line is the 109,087 deduplicated `ePhaZ`/`iPhaZ` protein candidates and 76,141 `genome x family` records. It represents computational homology/function-potential evidence, not PHB or PHA degradation phenotype proof.
- The formal registry `pipeline/config/formal_scan_models.tsv`, all historical `results/`, and all existing `runs/` remain read-only. Preserve `run-11`, retries, overlength exclusions, and negative/control evidence.
- A new run may be created or executed only after separate explicit authorization. It must have an audit-safe dated `run_id`, `logs/`, `inputs/`, `results/`, and `input_contract.json`; any T141 execution must use a matching dated `deploy/<run_id>/` source snapshot.
- Record source database, version or release, retrieval date, exact accession, file size, path, SHA-256, tool version, and parameter values. Use `pending` for unavailable evidence; never infer or invent it.
- Do not use SignalP, taxonomy, neighborhood markers, a motif, a domain hit, or a tree as experimental proof of PHB degradation. `OTHER` in SignalP is an unassigned localization prediction, not proof of cytosolic localization.
- `ArchPhaZ_patatin`, `ArchPhaZ_hydrolase`, `OH`, `PhaJ`, `BdhA`, `PhaC`, and phasin remain outside this ePhaZ/iPhaZ PhaDED main line. They may be context fields only where their separate evidence rules permit.
- Every new script starts with a failing focused test. Before accepting a new candidate-only analysis, run the focused tests, the applicable test suite, `compileall`, and `git diff --check`. Commit and push require separate explicit authorization.

## Scientific Decisions Fixed by This Plan

1. **Use all eight PhaDED superfamilies as the first-level architecture, not all 38 families as mandatory new GTDB scan models.** The eight classes are: intracellular nPHASCL without a lipase box; intracellular nPHASCL with a lipase box; periplasmic/native-PHA; intracellular nPHAMCL; extracellular dPHASCL type 1; extracellular dPHASCL type 2; extracellular native-SCL/PhaZ7-like; and extracellular dPHAMCL.
2. **Use the 38 homologous families as a second-level reference vocabulary.** A candidate may receive one best family, an alternative family, and a score gap, or remain `ambiguous`/`unassigned_PhaDED_like`. An absent family in the GTDB candidates is a result, not a failure. The project will not create 38 global HMM scan models unless later reference coverage and calibration justify a specific model.
3. **PhaDED mapping does not wait for the four-seed ePhaZ repair.** The current candidates already exist and can be routed/annotated as candidate-only sequences. The repair is required before redefining `ePhaZ_curated_core`, claiming broad subtype classification performance, or launching any replacement/full GTDB scan.
4. **Do not enlarge one mixed ePhaZ HMM.** The current four seeds (`B2NHN2`, `O05527`, `P12625`, `Q51718`) include both SCL/PHB-associated and MCL-PHA-associated space. The repair must split experimentally supported reference sets by PhaDED-compatible subtype and substrate/localization evidence before any model is built.
5. **PhaDED is a sequence/architecture classification scaffold, not an assay substitute.** Substrate designations in the reference table come from named experimental assays; GTDB assignments remain `candidate_only` even when their architecture is consistent with a reference class.

## Planned Files and Responsibilities

| Path | Responsibility |
|---|---|
| `docs/superpowers/plans/2026-09-09-phaded-ephaz-iphaz-next-stage.md` | This approved execution plan; it creates no run and changes no model. |
| `runs/20260909_phaded_reference_evidence_01/` | Planned, authorization-gated frozen PhaDED reference/evidence library. |
| `runs/20260909_phaded_candidate_mapping_01/` | Planned, authorization-gated candidate-only PhaDED mapping of existing ePhaZ/iPhaZ records. |
| `runs/20260909_ephaz_core_calibration_01/` | Planned, authorization-gated independent calibration of the four-seed core and subtype models. |
| `pipeline/scripts/validate_phaded_reference_ledger.py` | Fail-closed validation of reference accessions, evidence fields, ontology values, and source provenance. |
| `pipeline/scripts/build_phaded_profiles.py` | Build explicitly scoped superfamily/family reference profiles; it must not edit the formal scan registry. |
| `pipeline/scripts/annotate_phaded_features.py` | Compute reproducible sequence, motif, localization, domain, architecture, and integrity evidence without assigning a phenotype. |
| `pipeline/scripts/assign_phaded_candidates.py` | Merge per-candidate evidence and assign first/second PhaDED labels with a reason and uncertainty state. |
| `pipeline/scripts/calibrate_ephaz_subtype_models.py` | Evaluate frozen, subtype-specific experimental panels and held-out positives/controls; it must not create a replacement model silently. |
| `pipeline/tests/test_validate_phaded_reference_ledger.py` | Unit tests for reference-ledger contracts. |
| `pipeline/tests/test_build_phaded_profiles.py` | Unit tests for model scope, sequence leakage, and profile manifest completeness. |
| `pipeline/tests/test_annotate_phaded_features.py` | Unit tests for feature extraction and explicit `pending`/`no_record` states. |
| `pipeline/tests/test_assign_phaded_candidates.py` | Unit tests for deterministic evidence integration and ambiguous calls. |
| `pipeline/tests/test_calibrate_ephaz_subtype_models.py` | Unit tests for held-out evaluation, separate denominators, and rejection of mixed classes. |

## Phase Order and Decision Gates

| Phase | Deliverable | May run without another GTDB scan? | Gate to continue |
|---|---|---:|---|
| 0. Freeze authority | Manifest-bound input inventory and terminology | Yes, after a new-run authorization | Accepted inputs exactly match existing accepted outputs. |
| 1. Reference evidence | PhaDED ontology plus accession-level assay ledger | Yes | Every reference is source-bound; unsupported fields are `pending`. |
| 2. Candidate mapping | Eight-superfamily routing and 38-family best/second labels | Yes | All assignments have score, alternative, and evidence-status fields. |
| 3. Architecture review | SignalP, motif, triad, domain/SBD/linker/lid, integrity evidence | Yes | No missing data is converted into a negative; ambiguous candidates remain retained. |
| 4. Four-seed repair | Subtype-specific calibrated models or a documented failed calibration | Yes | Held-out positives and negative/challenge panels pass predeclared criteria. |
| 5. Phylogeny and interpretation | Representative trees, taxonomic summaries, priority list | Yes | Only mapped, quality-controlled representatives enter a tree. |
| 6. Any registry replacement/full rescan | New registry decision and formal scan plan | No | Separate explicit authorization after Phase 4 review. |
| 7. Experimental validation | Candidate-specific substrate/localization activity evidence | Not applicable | Experimental evidence is retained separately from computational labels. |

### Task 0: Freeze the current authority and create a new-run contract only after approval

**Files:**
- Read: `AGENTS.md`
- Read: `docs/T141_20260909_project_handoff.md`
- Read: `runs/20260905_ephaz_iphaz_layering_01/results/layers_retry_02/protein_layers.tsv`
- Read: `runs/20260906_overall_pha_phb_archaea_audit_04/results/overall_pha_phb.tsv`
- Planned create: `runs/20260909_phaded_reference_evidence_01/input_contract.json`
- Planned create: `runs/20260909_phaded_candidate_mapping_01/input_contract.json`
- Planned create: `runs/20260909_ephaz_core_calibration_01/input_contract.json`

**Consumes:** Accepted frozen layer and overall-audit tables, their manifests, the exact HMM assets used by the accepted analysis, and the local source snapshot.

**Produces:** Three input contracts with input paths, sizes, SHA-256 values, source snapshot hash, environment/tool records, and `formal_registry_modified=false` / `formal_scan_started=false` for all candidate-only work.

  - [x] **Step 1: Obtain explicit authorization to create the three dated run directories and to perform the bounded candidate-only analyses.**

  This plan by itself does not authorize run creation, local recomputation, server execution, installation, formal scan changes, commits, pushes, or deletions.

  - [x] **Step 2: Before copying any input, compare source hashes against the accepted manifests.**

  Required accepted source identities include the layer tables in `layers_retry_02`, `overall_pha_phb.tsv`, and the HMM/profile evidence recorded in their respective input contracts. Stop on a mismatch; create no candidate result from an unbound replacement input.

  - [x] **Step 3: Create each run with `pipeline/scripts/run_context.py` and write an input contract before executing analysis code.**

  The candidate mapping and calibration runs must declare `gtdb.taxonomy`, `gtdb.metadata`, and `gtdb.tree` as `pending` when they are not inputs. They are candidate-only analyses, not a hidden rescanning workflow.

  - [x] **Step 4: Create a source snapshot only in the matching dated `deploy/<run_id>/` before any T141 command.**

  The deploy manifest must bind the exact scripts, tests, reference ledger, candidate inputs, HMMs, and executable versions. Never execute a T141 project-root script.

**Acceptance:** All three contracts are valid JSON, contain no fabricated hash, and point only to input files that exist and are non-empty. No pre-existing run or registry file changes.

### Task 1: Build the PhaDED accession-level reference and evidence ledger

**Files:**
- Create: `pipeline/scripts/validate_phaded_reference_ledger.py`
- Create: `pipeline/tests/test_validate_phaded_reference_ledger.py`
- Planned create: `runs/20260909_phaded_reference_evidence_01/inputs/phaded_reference_ledger.tsv`
- Planned create: `runs/20260909_phaded_reference_evidence_01/inputs/phaded_family_definitions.tsv`
- Planned create: `runs/20260909_phaded_reference_evidence_01/inputs/phaded_reference.faa`
- Planned create: `runs/20260909_phaded_reference_evidence_01/results/reference_validation.json`
- Read: `research/europepmc/comp_classification.md`
- Read: `research/europepmc/comp_ft_knoll_phaded.txt`
- Read: `pipeline/seeds/ephaz_curated_evidence.tsv`

**Consumes:** PhaDED/Knoll reference records, exact public protein sequences, primary-paper or reviewed-database evidence, existing ePhaZ curation, and accession-level positive/negative/control evidence.

**Produces:** A frozen reference library whose family IDs are copied from a recorded PhaDED source, not invented from protein annotations. The ledger contains exactly the eight allowed superfamily labels and the source-defined 38 family IDs, but does not force every family to have a GTDB hit.

**Required ledger columns:**

```text
reference_id,accession,sequence_sha256,sequence_source,source_database,
source_database_version,retrieval_date,organism,taxonomy_id,phaded_superfamily,
phaded_family_id,reported_localization,substrate_class,substrate_detail,
experimental_assay,experimental_result,evidence_status,positive_negative_control,
catalytic_residues,lipase_box_or_ahsmg,oxyanion_hole_evidence,
domain_architecture,primary_doi,pmid,pmcid,notes
```

Allowed `evidence_status` values are `experimental_positive`, `experimental_negative`, `challenge_control`, `annotation_only`, and `pending_review`. Allowed assignment values are defined in the same immutable `phaded_family_definitions.tsv`; `pending_review` records cannot train a strict model.

  - [x] **Step 1: Write the failing ledger tests before the validator.**

```python
class PhaDEDLedgerTests(unittest.TestCase):
    def test_rejects_reference_without_accession_or_primary_identifier(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            ledger = write_ledger(root, accession="", primary_doi="")
            with self.assertRaisesRegex(ValueError, "accession"):
                module.validate_ledger(ledger, write_family_definitions(root), write_reference_fasta(root))

    def test_rejects_annotation_only_record_from_strict_training(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            ledger = write_ledger(root, evidence_status="annotation_only", strict_training="yes")
            with self.assertRaisesRegex(ValueError, "strict-training"):
                module.validate_ledger(ledger, write_family_definitions(root), write_reference_fasta(root))
```

  - [x] **Step 2: Run the focused test and confirm failure because `validate_phaded_reference_ledger.py` is absent.**

  Run: `python -m unittest pipeline.tests.test_validate_phaded_reference_ledger -v`

  Expected: import or implementation failure.

  - [x] **Step 3: Implement fail-closed ledger validation.**

  The validator must check required headers, unique accession and sequence SHA-256 pairs, exact FASTA-to-ledger one-to-one correspondence, valid eight-superfamily/family-definition linkage, evidence identifiers for all experimental calls, and the ban on `annotation_only`/`pending_review` training records. It writes counts by evidence status and family, duplicate reports, source file hashes, and all validation failures to `reference_validation.json`.

  - [x] **Step 4: Freeze the reference collection protocol before retrieving any additional sequence.**

  A positive reference requires a stable accession connected to a named experimental assay, substrate, and result. A record is placed in `challenge_control` when it is a related but non-target substrate/localization/fold case, including iPhaZ-like and generic lipase/esterase confounders. Do not promote annotation-only entries to positives.

  - [x] **Step 5: Run the focused test until it passes and review the ledger manually.**

  Run: `python -m unittest pipeline.tests.test_validate_phaded_reference_ledger -v`

  Expected: all ledger-validity, duplicate, family-linkage, and FASTA-binding tests pass.

**Acceptance:** The exact PhaDED family vocabulary is provenance-bound. Each strict-training sequence is experimentally supported and source-resolved; unknown substrate, localization, residue, or architecture fields remain `pending`.

### Task 2: Build only PhaDED reference profiles needed for candidate routing

**Files:**
- Create: `pipeline/scripts/build_phaded_profiles.py`
- Create: `pipeline/tests/test_build_phaded_profiles.py`
- Planned create: `runs/20260909_phaded_reference_evidence_01/results/profiles/`
- Planned create: `runs/20260909_phaded_reference_evidence_01/results/profile_manifest.tsv`
- Read: `pipeline/scripts/ephaz_bridge_loo_calibration.py`
- Read: `pipeline/scripts/run_context.py`

**Consumes:** The validated reference FASTA and ledger from Task 1.

**Produces:** Separate HMMER reference profiles for superfamilies that have sufficient independent experimental sequences, plus family profiles only where the PhaDED family has a valid and non-leaking reference set. Every profile carries a training accession list, alignment hash, HMM hash, tool version, and explicit `model_status`.

**Profile rules:**

- A model never mixes `e-dPHASCL`, `e-dPHAMCL`, `i-nPHASCL`, `i-nPHAMCL`, periplasmic, and PhaZ7-like records merely to increase seed count.
- A family with fewer than three independent experimentally supported accessions is `reference_only`; it may support best-reference similarity and manual review but does not create a fitted family HMM.
- A superfamily with insufficient, confounded, or contradictory support is `reference_only` and routes no candidate automatically.
- The four-seed current ePhaZ HMM is retained unchanged as historical candidate evidence. It is neither overwritten nor used as the PhaDED family truth set.

  - [x] **Step 1: Write failing tests for no cross-superfamily training leakage and profile manifest completeness.**

```python
class PhaDEDProfileTests(unittest.TestCase):
    def test_rejects_mixed_superfamily_training_group(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            ledger = write_validated_ledger(root, ["e_dPHASCL_type1", "i_nPHASCL_no_lipase_box"])
            with self.assertRaisesRegex(ValueError, "mixed superfamily"):
                module.build_profiles(ledger, write_reference_fasta(root), root / "profiles")

    def test_reference_only_family_has_no_hmm_but_is_recorded(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            manifest = module.build_profiles(write_singleton_family_ledger(root), write_reference_fasta(root), root / "profiles")
            self.assertEqual(manifest["family_status"], "reference_only")
            self.assertEqual(manifest["hmm_path"], "")
```

  - [x] **Step 2: Run the focused test and confirm it fails before implementation.**

  Run: `python -m unittest pipeline.tests.test_build_phaded_profiles -v`

  Expected: import or implementation failure.

  - [x] **Step 3: Implement deterministic profile creation.**

  Use a documented exact accession order, produce a separate MSA per allowed profile, run `hmmbuild`, and create `profile_manifest.tsv` with `profile_id`, `phaded_superfamily`, `phaded_family_id`, `training_accessions`, `training_count`, `model_status`, `alignment_sha256`, `hmm_sha256`, `mafft_version`, and `hmmer_version`. Reject duplicate sequence hashes across training/test partitions unless the record is explicitly excluded from one partition.

  - [x] **Step 4: Run focused tests and record every profile that remains `reference_only`.**

  Run: `python -m unittest pipeline.tests.test_build_phaded_profiles -v`

  Expected: all tests pass; `reference_only` is a valid outcome, not an exception to be patched around.

**Acceptance:** No profile spans a PhaDED superfamily boundary. No training record is hidden in a held-out panel. The formal registry remains byte-identical to its pre-run hash.

### Task 3: Route existing candidates through PhaDED and preserve ambiguity

**Files:**
- Create: `pipeline/scripts/assign_phaded_candidates.py`
- Create: `pipeline/tests/test_assign_phaded_candidates.py`
- Planned create: `runs/20260909_phaded_candidate_mapping_01/inputs/candidate_layers.tsv`
- Planned create: `runs/20260909_phaded_candidate_mapping_01/inputs/candidate_sequences.faa`
- Planned create: `runs/20260909_phaded_candidate_mapping_01/inputs/profile_manifest.tsv`
- Planned create: `runs/20260909_phaded_candidate_mapping_01/results/phaded_superfamily_assignments.tsv`
- Planned create: `runs/20260909_phaded_candidate_mapping_01/results/phaded_family_assignments.tsv`
- Planned create: `runs/20260909_phaded_candidate_mapping_01/results/mapping_summary.json`

**Consumes:** Exactly the deduplicated ePhaZ/iPhaZ candidate input, the validated reference profiles, and model-search tables against the candidate-only FASTA. The scope is the existing main-line candidates; no sequence outside that set enters this analysis.

**Produces:** For each candidate: `phaded_superfamily_best`, `phaded_superfamily_second`, `superfamily_score_gap`, `phaded_family_best`, `phaded_family_second`, `family_score_gap`, `assignment_status`, and the evidence model/threshold that produced each value.

**Allowed `assignment_status` values:** `assigned`, `ambiguous_superfamily`, `ambiguous_family`, `reference_only_nearest`, `unassigned_PhaDED_like`, `not_tested_missing_sequence`, and `excluded_input_integrity`.

  - [x] **Step 1: Write failing tests for deterministic winner selection, ties, and missing evidence.**

```python
class PhaDEDAssignmentTests(unittest.TestCase):
    def test_tied_superfamily_scores_remain_ambiguous(self):
        result = module.assign_candidate(
            accession="GCF_000001|gene_1",
            scores={"e_dPHASCL_type1": 102.0, "e_dPHASCL_type2": 102.0},
            profile_status={"e_dPHASCL_type1": "calibrated", "e_dPHASCL_type2": "calibrated"},
        )
        self.assertEqual(result["assignment_status"], "ambiguous_superfamily")
        self.assertEqual(result["phaded_superfamily_best"], "")

    def test_reference_only_family_is_not_presented_as_calibrated(self):
        result = module.assign_candidate("GCF_000001|gene_2", {"family_17": 88.0}, {"family_17": "reference_only"})
        self.assertEqual(result["assignment_status"], "reference_only_nearest")
```

  - [x] **Step 2: Run the focused test and confirm failure before implementation.**

  Run: `python -m unittest pipeline.tests.test_assign_phaded_candidates -v`

  Expected: import or implementation failure.

  - [x] **Step 3: Implement the candidate-only mapper.**

  The mapper must retain the current fields `accession`, `genome`, `family`, `layer`, `signal_type`, and `pha_subtype_evidence`; preserve the original HMM evidence separately from PhaDED scores; record best and second profile scores; and use a score-gap rule that is declared in the run contract before inspecting results. A tie or insufficient gap is never force-assigned.

  - [x] **Step 4: First analyze the highest-information strata before any bulk interpretation.**

  The required review order is: 18 `ePhaZ_competition_review` records; ePhaZ curated-core secreted records; ePhaZ curated-core nonsecreted records; iPhaZ tier1 records; then tier2 review layers. Broad discovery remains outside the current main count and is not silently added.

  - [x] **Step 5: Run focused tests and perform input/output reconciliation.**

  Run: `python -m unittest pipeline.tests.test_assign_phaded_candidates -v`

  Expected: every input accession occurs once in a primary assignment table, ambiguous/unassigned outcomes are retained, and no row reports a phenotype.

**Acceptance:** The output can answer “what PhaDED architecture is most consistent with this candidate?” while always retaining the alternative, evidence class, and uncertainty. It cannot answer “does this organism degrade PHB?”

### Task 4: Add independent sequence and architecture evidence without converting absence to a negative

**Files:**
- Create: `pipeline/scripts/annotate_phaded_features.py`
- Create: `pipeline/tests/test_annotate_phaded_features.py`
- Planned create: `runs/20260909_phaded_candidate_mapping_01/results/phaded_feature_annotation.tsv`
- Planned create: `runs/20260909_phaded_candidate_mapping_01/results/phaded_assignment_evidence.tsv`
- Planned create: `runs/20260909_phaded_candidate_mapping_01/results/feature_annotation_manifest.json`
- Read: `pipeline/scripts/review_ephaz_ambiguous_structure.py`
- Read: `runs/20260905_ephaz_iphaz_layering_01/results/layers_retry_02/protein_layers.tsv`

**Consumes:** Candidate sequences, accepted SignalP6 calls already bound to the ePhaZ curated-core layer, PhaDED mapping results, and independently versioned domain annotation when it is available and registered in the contract.

**Produces:** A feature table with exact residues/coordinates and evidence status. Required fields are:

```text
accession,sequence_integrity,signalp_evidence,signalp_class,
lipase_box_state,lipase_box_coordinates,ahsmg_state,ahsmg_coordinates,
catalytic_ser_cys_state,catalytic_ser_cys_coordinates,his_state,his_coordinates,
asp_state,asp_coordinates,oxyanion_hole_state,oxyanion_hole_coordinates,
sbd_state,sbd_coordinates,linker_state,linker_coordinates,lid_state,lid_coordinates,
pfam_interpro_state,architecture_consistency,feature_evidence_status
```

**Interpretation rules:**

- `GxSxG` with a suitable hydrophobic residue, `AHSMG`, and the proposed catalytic residues are support fields; motif occurrence alone is not a family or activity call.
- A catalytic triad is `supported` only when its residues/positions are consistent with the candidate's superfamily alignment or reference model. Otherwise use `partial`, `conflicting`, or `pending`.
- Classical extracellular dPHASCL architecture is signal peptide + catalytic domain + optional linker + C-terminal SBD. dPHAMCL is assessed separately because its N terminus can be the binding region and it lacks the classical SBD/linker organization.
- `i-nPHASCL` without a lipase box uses a Cys-His-Asp route; i-nPHAMCL requires explicit lid evidence before a lid-supported label. PhaZ7-like candidates require explicit `AHSMG`/architecture review rather than a failed `GxSxG` search.
- SignalP `OTHER`, unavailable domain scans, absent annotation, sequence truncation, and no neighborhood record are recorded as `pending`, `no_record`, or `possible_N_truncation` as applicable; none becomes a biological negative.

  - [x] **Step 1: Write failing tests for motif exceptions and missing external evidence.**

```python
class PhaDEDFeatureTests(unittest.TestCase):
    def test_phaz7_like_ahsmg_is_not_failed_for_lacking_gx_sx_g(self):
        features = module.annotate_sequence("candidate_1", "M" + "A" * 80 + "AHSMG" + "A" * 180, "e_nPHASCL_phaz7")
        self.assertEqual(features["ahsmg_state"], "present")
        self.assertEqual(features["lipase_box_state"], "not_expected_for_superfamily")

    def test_unavailable_domain_database_stays_pending(self):
        features = module.annotate_sequence("candidate_2", "M" + "A" * 240, "i_nPHAMCL", domain_records=None)
        self.assertEqual(features["pfam_interpro_state"], "pending")
        self.assertEqual(features["lid_state"], "pending")
```

  - [x] **Step 2: Run the focused test and confirm failure before implementation.**

  Run: `python -m unittest pipeline.tests.test_annotate_phaded_features -v`

  Expected: import or implementation failure.

  - [x] **Step 3: Implement feature extraction with coordinate-preserving output.**

  Write motif/triad/lipase-box findings with sequence coordinates and a tested/not-tested status. Import existing SignalP fields without recomputing them. Parse independently registered Pfam/InterPro results only if their database version, command, checksums, and candidate FASTA hash appear in the input contract.

  - [x] **Step 4: Integrate feature evidence with the mapper without overwriting the mapping table.**

  `phaded_assignment_evidence.tsv` must retain raw profile calls and add `architecture_consistency` as `consistent`, `partial`, `conflicting`, or `pending`. A conflicting motif/domain result downgrades an automatic label to review; it does not delete the candidate.

  - [x] **Step 5: Run focused tests and inspect the priority strata manually.**

  Run: `python -m unittest pipeline.tests.test_annotate_phaded_features -v`

  Expected: the PhaZ7 and Cys-type exceptions are handled explicitly; unavailable domain evidence remains non-negative.

**Acceptance:** Each assignment is reproducible to exact coordinates and evidence sources. The analysis distinguishes architecture-consistent candidates from unresolved candidates without claiming activity.

### Task 5: Repair the four-seed ePhaZ model as a separate calibration study

**Files:**
- Create: `pipeline/scripts/calibrate_ephaz_subtype_models.py`
- Create: `pipeline/tests/test_calibrate_ephaz_subtype_models.py`
- Planned create: `runs/20260909_ephaz_core_calibration_01/inputs/ephaz_seed_audit.tsv`
- Planned create: `runs/20260909_ephaz_core_calibration_01/inputs/subtype_training_panels.tsv`
- Planned create: `runs/20260909_ephaz_core_calibration_01/inputs/heldout_positive.faa`
- Planned create: `runs/20260909_ephaz_core_calibration_01/inputs/formal_negative.faa`
- Planned create: `runs/20260909_ephaz_core_calibration_01/inputs/challenge_controls.faa`
- Planned create: `runs/20260909_ephaz_core_calibration_01/results/leave_one_out_metrics.tsv`
- Planned create: `runs/20260909_ephaz_core_calibration_01/results/subtype_model_decision.json`
- Read: `pipeline/scripts/ephaz_bridge_loo_calibration.py`
- Read: `pipeline/scripts/calibrate_ephaz_external_panel.py`
- Read: `pipeline/scripts/prepare_ephaz_bridge_run.py`

**Consumes:** A corrected experimental evidence ledger for `B2NHN2`, `O05527`, `P12625`, and `Q51718`; optional bridge candidates already curated separately; independent positives; formal negatives; and related challenge controls. The current model is evaluated as a historical baseline, never overwritten.

**Produces:** A decision for each candidate subtype: `calibrated_candidate_model`, `reference_only_insufficient_panel`, or `rejected_mixed_or_nonseparable`. The outcome is a model-governance result, not a GTDB rescan.

**Required seed audit:**

| Accession | Required treatment |
|---|---|
| `B2NHN2` | Verify exact assay, substrate, localization, sequence completeness, PhaDED route, and primary reference. |
| `O05527` | Verify exact assay, substrate, localization, sequence completeness, PhaDED route, and primary reference. |
| `P12625` | Verify exact assay, substrate, localization, sequence completeness, PhaDED route, and primary reference. |
| `Q51718` | Retain as an MCL-PHA-associated control/reference only in the matching PhaDED-compatible group; it must not be treated as PHB/SCL positive by default. |

**Calibration rules:**

- Create separate panels for e-dPHASCL type 1, e-dPHASCL type 2, e-nPHASCL/PhaZ7-like, and e-dPHAMCL only where the evidence ledger supports those categories.
- Keep iPhaZ-like, MCL-PHA non-PHB, annotation-only near neighbors, generic lipase/esterase, and fragments in separately reported negative/challenge panels. Do not mix their denominators.
- Use leave-one-out testing for every trainable positive reference. A sequence or sequence-identical cluster cannot occur in both training and held-out partitions.
- Freeze score/e-value/coverage decision rules before viewing the held-out result. A stricter coverage threshold that merely removes true positives is a failed trade-off, not a validation success.
- A subtype can proceed only after it has enough independently documented experimental positives to support a held-out test, represents at least three genera, recovers the predeclared held-out positives at the frozen threshold, and has no unexplained formal-negative hit. Otherwise it remains `reference_only_insufficient_panel`.

  - [x] **Step 1: Write failing tests that prohibit mixed SCL/MCL training and score-tuned evaluation.**

```python
class EphaZSubtypeCalibrationTests(unittest.TestCase):
    def test_rejects_q51718_from_phb_scl_training_panel(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            panels = write_panel_manifest(root, subtype="e_dPHASCL_type1", accessions=["B2NHN2", "Q51718"])
            with self.assertRaisesRegex(ValueError, "substrate-class mismatch"):
                module.calibrate_subtype_models(panels, root / "out")

    def test_rejects_heldout_accession_present_in_training(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            panels = write_leaking_panel_manifest(root, accession="B2NHN2")
            with self.assertRaisesRegex(ValueError, "train-heldout leakage"):
                module.calibrate_subtype_models(panels, root / "out")
```

  - [x] **Step 2: Run the focused test and confirm failure before implementation.**

  Run: `python -m unittest pipeline.tests.test_calibrate_ephaz_subtype_models -v`

  Expected: import or implementation failure.

  - [x] **Step 3: Implement subtype-specific calibration and separate denominators.**

  Reuse the audited parsing approach of `ephaz_bridge_loo_calibration.py` and `calibrate_ephaz_external_panel.py`, but output per-subtype train/held-out membership, per-panel TP/FN/FP/TN or `challenge_detected`, all thresholds tested, model/alignment hashes, and the predeclared decision rule. The script must fail if a model group crosses substrate/localization/PhaDED boundaries or if a panel accession appears twice.

  - [x] **Step 4: Decide the four-seed model's status from frozen evidence.**

  The possible outcomes are: retain it only as historical broad candidate evidence; replace it with one or more calibrated subtype candidate models in a future, separately authorized registry change; or document that current panels cannot support replacement. None of these outcomes permits silent editing of `ePhaZ_curated_core`.

  - [x] **Step 5: Run focused tests and compare baseline versus subtype models in a written decision record.**

  Run: `python -m unittest pipeline.tests.test_calibrate_ephaz_subtype_models -v`

  Expected: mixed panels, leakage, missing source evidence, and undeclared threshold changes fail closed.

**Acceptance:** The four-seed limitation is explicitly resolved as either a calibrated subtype-model program or a documented insufficient-evidence result. It is not concealed by adding more heterogeneous seeds.

### Task 6: Produce conservative downstream biology only after mapping and calibration

**Files:**
- Planned create: `runs/20260909_phaded_candidate_mapping_01/results/priority_candidate_set.tsv`
- Planned create: `runs/20260909_phaded_candidate_mapping_01/results/taxonomy_summary.tsv`
- Planned create: `runs/20260909_phaded_candidate_mapping_01/results/representative_selection.tsv`
- Planned create: `runs/20260909_phaded_candidate_mapping_01/results/interpretation_boundary.md`
- Planned create: `docs/T141_20260909_phaded_candidate_mapping_01_status.md`
- Planned modify after results are accepted: a new figure-only run and its own status record; never overwrite `20260907_progress_report_01`.

**Consumes:** PhaDED mapping, architecture evidence, current candidate layers, taxonomy, and only quality-controlled representative sequences.

**Produces:** A ranked candidate set for manual review/experiments, family-aware taxonomy summaries, and a small representative phylogenetic input set. It does not output a list of experimentally validated degraders.

  - [x] **Step 1: Define the priority-candidate rule before calculating ranks.**

  Rank `assigned` and architecture-`consistent` candidates above `partial`; retain `ambiguous` and `unassigned_PhaDED_like` in separate review lists. Apply sequence-integrity status, profile score gap, supported localization, and reference novelty as explicit columns rather than hiding them in an opaque score.

  - [x] **Step 2: Construct family-aware summaries with correct denominators.**

  Report proteins, `genome x PhaDED_superfamily`, and `genome x PhaDED_family` separately. Never merge these values with the run-13 strict four-family union (38,741 genomes), the historic Scheme A count (44,814), or the current 76,141 `genome x main-family` records.

  - [x] **Step 3: Select a bounded representative set before any tree.**

  Per PhaDED class, select validated references, all competition cases, high-confidence taxonomically diverse representatives, near-threshold/ambiguous candidates, and known negatives/challenges. Deduplicate exact sequences and record selection reason; do not build a tree from all 109,087 candidates as the first phylogenetic analysis.

  - [x] **Step 4: Make figure repairs in a new visualization-only run after the data table is frozen.**

  Correct Figure 4 terminology to `genome-family assignments` and explicitly label `Palsa-465` as a GTDB placeholder or replace it with a formally named genus under a documented selection rule. Captions must state candidate-only interpretation and taxonomy/HMM inclusion bias.

  - [x] **Step 5: Prepare experimental validation as a separate evidence layer.**

  For a small, diverse priority set, define cloned-gene/protein identity, signal-peptide/mature-protein construct, polymer substrate composition and physical state, assay conditions, heat-inactivated/vector controls, product measurement, and replicate criteria. Only a candidate's own experimental record can change an experimental-evidence field; it does not retroactively validate homologs.

**Acceptance:** All downstream tables state their unit and candidate-only boundary. Any tree or wet-lab result is traceable to a selected candidate and does not overwrite the computational source label.

### Task 7: Gate any future registry change or GTDB-wide rescan

**Files:**
- Planned create only after separate authorization: `docs/superpowers/plans/2026-09-09-phaded-formal-registry-change.md`
- Planned create only after separate authorization: `runs/20260909_phaded_formal_scan_01/`
- Planned create only after separate authorization: `deploy/20260909_phaded_formal_scan_01/`
- Read before authorization: `pipeline/config/formal_scan_models.tsv`
- Read before authorization: accepted PhaDED mapping and ePhaZ calibration manifests.

**Consumes:** A written decision that one or more subtype models is calibrated, model/HMM hashes, input panels, threshold rationale, control metrics, and the accepted candidate-only mapping report.

**Produces:** Nothing until a user explicitly authorizes the registry change and formal scan. A plan is drafted first; scan execution is not implicit.

  - [x] **Step 1: Do not request or create a formal rescan solely because a PhaDED label exists.**

  PhaDED routing of current candidates is sufficient for the next analytical stage. A full scan is warranted only when a calibrated, provenance-bound model would materially change candidate discovery.

  - [ ] **Step 2: Before a proposed change, compare three alternatives in the decision record.**

  The alternatives are: retain the current registry and publish candidate-only PhaDED mapping; add a narrow calibrated subtype model while retaining historical models; or replace a historical model with a versioned successor. The decision record must describe changes to counts and comparability before any execution.

  - [ ] **Step 3: Require a new authorization after the decision record is reviewed.**

  The authorization must explicitly cover registry modification, new dated deploy, formal scan, server execution, and any commit/push. Without it, Phase 7 remains planned and no file in `pipeline/config/` is changed.

**Acceptance:** The project never confuses an evidence-library improvement with permission to rescan 199,923 genomes.

## Verification Commands for Future Implementation

Run these only in an authorized implementation turn after the corresponding files exist:

```powershell
python -m unittest pipeline.tests.test_validate_phaded_reference_ledger -v
python -m unittest pipeline.tests.test_build_phaded_profiles -v
python -m unittest pipeline.tests.test_assign_phaded_candidates -v
python -m unittest pipeline.tests.test_annotate_phaded_features -v
python -m unittest pipeline.tests.test_calibrate_ephaz_subtype_models -v
python -m unittest discover -s pipeline/tests -v
python -m compileall pipeline/scripts
git diff --check
Get-FileHash -Algorithm SHA256 pipeline/config/formal_scan_models.tsv
```

The final registry SHA-256 must equal the pre-run value unless Phase 7 has separate written authorization and an accepted replacement manifest.

## Completion Definition

The next-stage computational project is complete when the accepted ePhaZ/iPhaZ candidate set has a frozen PhaDED reference provenance table, first/second superfamily and family evidence where supported, architecture/feature evidence with explicit uncertainty, a separate four-seed calibration decision, and correct candidate-only reporting. It is not complete merely because every candidate has a forced family name, nor because a new HMM has been constructed.

## Execution Handoff

This plan deliberately separates two authorization boundaries:

1. **Candidate-only implementation:** new reference/mapping/calibration runs, bounded local or dated-deploy T141 analysis, and new scripts/tests. This requires explicit authorization to create and execute those runs.
2. **Formal registry/full scan:** any edit to `formal_scan_models.tsv`, deployment for a GTDB-wide scan, or publication/commit/push. This requires a second, separate authorization after calibration review.

Recommended execution order is Tasks 0-4 for the PhaDED candidate-only mapping, then Task 5 for the four-seed decision, then Task 6. Task 7 remains blocked by design until its evidence and authorization gates are satisfied.
