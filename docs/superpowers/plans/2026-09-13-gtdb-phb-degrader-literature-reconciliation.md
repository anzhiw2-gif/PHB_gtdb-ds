# GTDB PHB Degrader Literature Reconciliation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `subagent-driven-development` (recommended) or an execution-plan workflow to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Using the frozen GTDB candidate universe and existing public records, identify genomes with direct literature-supported PHB degradation evidence and separate them from genomes with only sequence, domain, structure, or phylogenetic degradation potential.

**Architecture:** Keep all existing PhaDED/ePhaZ/iPhaZ runs immutable. Build a new accession-level literature and database evidence ledger, bridge exact organism/protein/genome identifiers to GTDB records, then produce genome-level evidence tiers and a conservative report. The primary analysis reuses the existing 109,087 candidate proteins and 76,141 genome-by-family records; a new broad HMM scan is not required for this question.

**Tech Stack:** Existing TSV/FASTA/JSON manifests, Python standard-library scripts and `unittest`, GTDB taxonomy/metadata already bound to accepted runs, Europe PMC/PubMed/NCBI/UniProt/GTDB public records when queried, accession and SHA-256 provenance manifests, and existing PhaDED Pfam/InterPro/structure/phylogeny outputs.

## Global Constraints

- The target claim is split into `direct_literature_supported`, `genome_supported_candidate`, and `sequence_only_candidate`; only the first tier may be described as a PHB degrader, and only when an exact strain/genome or accession-bound record supports the claim.
- Homology, HMM, Pfam, InterPro, SignalP, neighborhood, structure, and phylogeny evidence are supporting evidence or functional potential, not phenotype proof.
- Do not call `unassigned_PhaDED_like`, missing Foldseek, missing taxonomy, or missing literature a biological negative.
- Preserve `pipeline/config/formal_scan_models.tsv`, historical `results/`, historical `runs/`, HMMs, training panels, and all previous contracts.
- Every new run uses `runs/<dated_run_id>/` with `logs/`, `inputs/`, `results/`, and `input_contract.json`; source database/version/release, retrieval date, accession, path, size, SHA-256, query parameters, and tool versions are recorded.
- The no-wet-lab scope is explicit. New outputs use `not_planned_by_scope` or `not_executed` for experimental fields; they must not imply that missing wet-lab data is negative evidence.
- No registry edit, new formal deploy, GTDB-wide rescan, commit, push, installation, or deletion is part of this plan.

## Evidence Model

| Tier | Meaning | Allowed wording |
|---|---|---|
| A | Exact strain/genome/protein accession is linked to a primary or reviewed experimental record showing PHB degradation or PHB depolymerase activity | `direct_literature_supported` / `reported PHB degrader` |
| B | GTDB genome has an accession-bound homolog and coherent PhaDED/architecture/localization evidence, but no exact genome-specific phenotype record | `genome_supported_candidate` / `PHB-degradation potential` |
| C | Candidate is supported only by sequence/domain/structure/tree evidence or unresolved mapping | `sequence_only_candidate` |
| H | Conflicting architecture, partial source, low structure, or unresolved accession mapping | `hold_unresolved` |

The report must never collapse B/C/H into A.

## Current Authority

The starting evidence is:

- `runs/20260909_phaded_candidate_mapping_01/`: 109,087 candidate proteins.
- Existing family summary: 76,141 `genome x family` records.
- `runs/20260912_phaded_provenance_reconciliation_01/`: 14 retained dPHASCL1-like candidates.
- `runs/20260913_phaded_19_candidate_manual_review_01/`: 19-record manual decision table, including 3 architecture holds and 2 source/structure holds.
- `runs/20260909_ephaz_core_calibration_01/`: subtype models remain `reference_only_insufficient_panel`; no replacement model is authorized.

## Task 0: Freeze the no-wet-lab question and create a new run

**Files:**
- Create: `runs/20260913_gtdb_phb_degrader_literature_reconciliation_01/input_contract.json`
- Create: `runs/20260913_gtdb_phb_degrader_literature_reconciliation_01/logs/`
- Create: `runs/20260913_gtdb_phb_degrader_literature_reconciliation_01/inputs/`
- Create: `runs/20260913_gtdb_phb_degrader_literature_reconciliation_01/results/`

- [x] Copy only accepted manifests and candidate summaries into the new run.
- [x] Declare `formal_registry_modified=false`, `formal_scan_started=false`, `experimental_results_generated=false`, and `wet_lab_scope="not_planned_by_scope"`.
- [x] Bind the current source commit, GTDB release/path metadata, candidate table hashes, and the unchanged registry hash.

**Gate:** Stop if any accepted input hash differs from its source manifest.

## Task 1: Build the accession-level PHB literature ledger

**Files:**
- Create: `runs/20260913_gtdb_phb_degrader_literature_reconciliation_01/inputs/phb_literature_ledger.tsv`
- Create: `pipeline/scripts/validate_phb_literature_ledger.py`
- Create: `pipeline/tests/test_validate_phb_literature_ledger.py`
- Create: `runs/20260913_gtdb_phb_degrader_literature_reconciliation_01/results/literature_validation.json`

Required columns:

```text
record_id,organism_name,strain_name,protein_accession,gene_accession,genome_accession,
ncbi_taxid,gtdb_taxonomy,source_database,source_database_version,retrieval_date,
primary_doi,pmid,pmcid,experimental_substrate,experimental_result,
evidence_type,evidence_status,exact_identifier_scope,sequence_sha256,notes
```

- [x] Include only records with a stable identifier and a traceable paper/database source.
- [x] Separate PHB, other PHA, generic polyester, lipase/esterase, and annotation-only records.
- [x] Mark records with no direct PHB assay as `challenge_control`, `annotation_only`, or `pending_review`; do not promote them to Tier A.
- [x] Validate duplicate accessions, missing citation fields, unsupported evidence statuses, and sequence/source mismatches.

**Acceptance:** Every Tier A row has exact identifier scope, a source citation, substrate/result wording, and an auditable retrieval record.

## Task 2: Bridge literature records to GTDB genomes

**Files:**
- Create: `runs/20260913_gtdb_phb_degrader_literature_reconciliation_01/results/gtdb_bridge.tsv`
- Create: `runs/20260913_gtdb_phb_degrader_literature_reconciliation_01/results/gtdb_bridge_manifest.json`
- Create: `pipeline/scripts/bridge_phb_records_to_gtdb.py`
- Create: `pipeline/tests/test_bridge_phb_records_to_gtdb.py`

- [x] Prefer exact genome assembly/accession matches.
- [x] Allow exact protein-to-genome mapping only when the protein accession and GTDB genome are both bound by a source record.
- [x] Keep strain-name or species-name matches in a separate `taxonomy_only_match` class; taxonomy transfer alone cannot create Tier A evidence.
- [x] Record one-to-many and many-to-one mappings explicitly and retain unresolved records.
- [x] Bind GTDB release, taxonomy file, metadata file, retrieval date, and mapping method for every accepted bridge.

**Acceptance:** No Tier A label is assigned from a name-only or genus-only match.

## Task 3: Reuse existing GTDB candidate evidence at genome level

**Files:**
- Create: `runs/20260913_gtdb_phb_degrader_literature_reconciliation_01/results/genome_evidence_ledger.tsv`
- Create: `runs/20260913_gtdb_phb_degrader_literature_reconciliation_01/results/genome_evidence_summary.json`
- Create: `pipeline/scripts/assemble_gtdb_phb_evidence.py`
- Create: `pipeline/tests/test_assemble_gtdb_phb_evidence.py`

For each genome, join:

- exact literature bridge status;
- PhaDED superfamily/family counts;
- ePhaZ/iPhaZ layer and candidate counts;
- Pfam/InterPro architecture state;
- SignalP/localization state without treating `OTHER` as negative;
- structure/Foldseek/phylogeny evidence;
- sequence integrity and unresolved/hold flags;
- genome-level denominator and source run IDs.

- [x] Keep protein counts, unique genome counts, `genome x superfamily`, and `genome x family` as separate units.
- [x] Assign Tier A/B/C/H using a predeclared rule, not an opaque score.
- [x] A genome with multiple candidates is counted once in genome-level totals but retains all protein-level rows.
- [x] Preserve `pending`, `not_tested`, and `not_planned_by_scope` states.

**Acceptance:** The output can answer “which GTDB genomes have direct reported PHB degradation evidence?” separately from “which genomes contain computational candidates?”

## Task 4: Manually adjudicate the bounded high-information set

**Files:**
- Create: `runs/20260913_gtdb_phb_degrader_literature_reconciliation_01/results/high_information_adjudication.tsv`
- Create: `runs/20260913_gtdb_phb_degrader_literature_reconciliation_01/results/high_information_adjudication_summary.json`

- [x] Start with the 14 retained candidates in `experimental_priority_14.tsv`.
- [x] Keep the 3 dPHASCL2 architecture conflicts and 2 source/structure holds in Tier H unless a direct literature bridge resolves the claim.
- [x] Do not upgrade a candidate solely because Foldseek, pLDDT, MRCA support, or InterPro is strong.
- [x] Record whether each candidate has an exact literature bridge, a genome-supported candidate status, or sequence-only status.

**Acceptance:** Every high-information record has a written reason for its tier and an explicit non-phenotype boundary.

## Task 5: Produce the GTDB PHB-degrader report

**Files:**
- Create: `docs/T141_20260913_gtdb_phb_degrader_literature_reconciliation_01_status.md`
- Create: `docs/T141_20260913_gtdb_phb_degrader_catalog_report.md`
- Create: `runs/20260913_gtdb_phb_degrader_literature_reconciliation_01/results/figure_source_genome_tiers.tsv`

The report must contain:

1. A Tier A table of directly reported PHB degraders with exact accession and citation.
2. A Tier B table of GTDB genome-supported candidates with the supporting computational evidence and no phenotype claim.
3. A Tier C/H table of sequence-only and unresolved candidates.
4. Separate denominators for protein, genome, genome-by-superfamily, and genome-by-family units.
5. A missing-evidence table showing unresolved taxonomy, incomplete source records, absent structures, and literature gaps.
6. A methods/provenance appendix with database versions, retrieval dates, hashes, run IDs, and all classification rules.

## Task 6: Verification and stopping rule

- [x] Add failing tests before each new script, then run focused tests.
- [x] Run `python -m unittest discover -s pipeline/tests -v`.
- [x] Run `python -m compileall -q pipeline`.
- [x] Run `git diff --check`.
- [x] Verify that `pipeline/config/formal_scan_models.tsv` retains SHA-256 `8a1c05085f882daad9fddb9cf7616def74438d3a7bf55f219e0b2046adc5b7a3`.
- [x] Stop after the catalog and report. Do not create a formal registry-change plan or GTDB-wide rescan from a candidate-only result.

## Expected Outcome

The project will not claim that every GTDB sequence candidate is a PHB-degrading bacterium. It will produce a defensible catalog with:

- exact literature-supported PHB degraders where the public evidence truly supports that claim;
- GTDB genomes with computational PHB-degradation potential;
- unresolved and conflicting records retained visibly;
- no wet-lab dependency and no silent model or registry change.
