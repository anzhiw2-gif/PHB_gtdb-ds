# T141 PhaDED 19-candidate manual review

Date: 2026-09-13  
Run: `runs/20260913_phaded_19_candidate_manual_review_01/`  
Status: `completed_candidate_only_manual_review`

## Scope

This run consolidates the existing accession-bound audits for all 19
high-information candidates. It does not rerun ColabFold, Foldseek, Pfam,
InterPro, phylogeny, or source retrieval. It does not infer missing values.

## Decisions

| Group | Count | Decision |
|---|---:|---|
| dPHASCL2 architecture conflict | 3 | `hold_architecture_conflict` |
| source/gene-model or low-structure hold | 2 | `hold_gene_model_or_structure` |
| retained bounded comparison candidates | 14 | retain candidate-only; split into 5 P1, 4 Foldseek-supported P2, and 5 pending-structure P2 |

The three dPHASCL2 records have reliable compact PF00756 cores, but lack the
expected dPHASCL2 reference fingerprint. High pLDDT, SignalP class, or tree
support does not override that architecture conflict.

The JALYXY record is source-bound as a C-terminally partial CDS (the unified
table records `partial_C_terminal_source_record`) and remains on hold. The
202-aa `NZ_CP083952.1_5059` ORF is complete but low-confidence and
must remain separate from the independent 384-aa `..._5067` ORF.

The nine intermediate dPHASCL1-like structures are all represented in the
complete audit. Four have additional priority-burden decisions; the remaining
five remain intermediate structure records with no automatic promotion.

## Experimental priority

`results/experimental_priority_14.tsv` gives the bounded order for the 14
retained records: five P1 records first, then four records with accession-bound
Foldseek support, then five records awaiting reference-bound structure
comparison. It is a review order, not an experimental result or an instruction
to execute wet-lab work.

## Boundaries

All 19 decisions remain `candidate_only`. No record is promoted to an ePhaZ
positive or a training accession. `unassigned_PhaDED_like` and pending
structure comparison remain non-negative evidence states. The run does not
modify historical contracts, registry files, HMMs, GTDB scan outputs, or
server state.

## Verification

The two output tables were parsed as tab-delimited data: the manual review has
19 unique candidate IDs and the experimental priority table has 14 unique
candidate IDs. Input and output hashes are recorded in the run contract.
