# T141 PhaDED experimental validation plan

Date: 2026-09-13  
Run: `runs/20260913_phaded_experimental_validation_plan_01/`  
Status: `completed_candidate_only_experimental_plan`

## Scope

This dated run prepares a candidate-specific wet-lab validation template for
the 14 retained records from the provenance reconciliation run. It is a plan
and evidence-bound data-entry surface, not an experiment execution.

## Inputs and outputs

- Input copy: `inputs/candidate_review_14.tsv`, copied from the verified
  14-candidate reconciliation table without changing candidate IDs or review
  decisions.
- Output: `results/experimental_validation_template.tsv` with one row per
  candidate and fields for construct identity, accession binding, mature
  protein/signal-peptide strategy, expression and purification, polymer
  substrate, assay conditions, controls, product measurement, replicates,
  acceptance criteria, and experimental evidence status.
- The run contract binds the input and output sizes/SHA-256 and the current
  local source commit. Unavailable GTDB/source/environment metadata is marked
  `pending` rather than inferred.

## Candidate order

The template preserves the evidence-gated order:

1. Five `P1` `retain_priority_candidate_only` records.
2. Four `P2_foldseek` `retain_with_fold_support_candidate_only` records.
3. Five `P2_pending_structure` `retain_pending_structure_comparison` records.

Foldseek values are recorded only as structural homology support. The five
records without a reference-bound comparison remain pending; that state is
not a biological negative and does not block documenting them in the template.

## Pending policy and boundaries

All experimental fields that have not been executed or authorized are
`pending`, including protein/gene accession resolution, construct design,
expression host, mature-protein and signal-peptide strategy, purification,
substrate composition/state, assay conditions, controls, measurement method,
replicates, acceptance criteria, and experimental evidence status. No assay
result, phenotype label, training-positive status, registry edit, formal scan,
or server execution is implied.

Only an accession-bound experimental record for the corresponding candidate
may later change its `experimental_evidence_status`. Such a change cannot
retroactively validate homologs, alter the PhaDED mapping, or promote a
PhaC-vs-PhaZ challenge record into an ePhaZ positive panel.

## Verification

The template was parsed as tab-delimited data and checked for 14 unique
candidate IDs. The contract records the verified SHA-256 values. No historical
contract, registry, HMM, training panel, GTDB scan, or historical run was
modified.
