# T141 PhaDED provenance reconciliation and 14-candidate review

Date: 2026-09-12
Run: `runs/20260912_phaded_provenance_reconciliation_01/`
Status: `completed_candidate_only_reconciliation`

## Scope

This run reconciles two stale execution descriptions without editing their
historical contracts:

- `20260909_phaded_reference_evidence_01` remains labelled
  `planned_candidate_only`, while its verified reference-validation and
  profile-audit outputs are consumed as completed candidate-only evidence.
- `20260911_phaded_competition_phylo_review_01` retains its original
  `tree:not_run` and `deploy_dir:pending_not_created` fields, while its
  completed phylogeny manifest and four panel trees are consumed as
  candidate-only evidence. The discrepancy is recorded, not rewritten.

The merged evidence is limited to Pfam architecture, secondary InterPro
evidence, ColabFold confidence/PAE review tiers, competition phylogeny, and the
accession-bound Foldseek comparison. No registry, historical HMM, training
panel, or formal scan was changed.

## Review result for the retained 14

| Group | Count | Decision |
|---|---:|---|
| P1 high architecture/phylogeny | 5 | `retain_priority_candidate_only` |
| P2 with reference-bound Foldseek | 4 | `retain_with_fold_support_candidate_only` |
| P2 without reference-bound Foldseek | 5 | `retain_pending_structure_comparison` |

All 14 are extracellular dPHASCL1-like candidate records with explicit
PHA-related InterPro evidence and partial Pfam architecture. The four
Foldseek comparisons have moderate alignment TM-scores (0.6214-0.6974), but
their intermediate-confidence peripheral regions limit interpretation. The
other five retain `pending` structural-comparison fields; missing comparison
evidence is not a biological negative.

The detailed, accession-level table is
`runs/20260912_phaded_provenance_reconciliation_01/results/candidate_review_14.tsv`.
The machine-readable summary and discrepancy record are
`results/candidate_review_summary.json` and
`results/provenance_reconciliation.json` in the same run.

## Recommended review order

Start with the five P1 records, then the four P2 records with Foldseek support
and the four highest-burden intermediate structures. Bind additional
reference structures before comparing the remaining five P2 records. Any
structural or phylogenetic similarity remains homology/architecture evidence;
none of the 14 is a validated PHB/PHA degradation phenotype or a training
positive.

## Verification

`python -m unittest discover -s pipeline/tests -v` passed with 208 tests and
one skipped test before this documentation-only addition. The reconciliation
JSON and TSV were checked for valid JSON, 14 unique candidate IDs, and the
decision counts above. Final repository verification is rerun after this
addition.
