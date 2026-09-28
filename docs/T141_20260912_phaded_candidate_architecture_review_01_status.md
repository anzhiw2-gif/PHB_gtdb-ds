# T141 PhaDED candidate architecture review

Date: 2026-09-12  
Run: `runs/20260912_phaded_candidate_architecture_review_01/`  
Status: `completed_candidate_only`

## Inputs

The review consumes the frozen 19-candidate structure/evidence table from
`runs/20260912_phaded_structure_evidence_review_01/results/phaded_structure_evidence_review.tsv`.
The table already binds ColabFold pLDDT/PAE summary, Pfam coordinates and
architecture state, secondary InterPro evidence, motif fields, and competition
phylogeny.

## Decisions

| Decision | Count | Meaning |
|---|---:|---|
| `hold_architecture_conflict` | 3 | All dPHASCL2-like candidates: InterPro is PHA-related, but Pfam shows only PF00756 and conflicts with the type-2 reference fingerprint. |
| `hold_gene_model_or_structure` | 2 | One possible N-terminal truncation and one low-confidence structure require source/gene-model review. |
| `retain_candidate_comparison` | 14 | No hard architecture or integrity hold; retain for bounded structure/domain comparison only. |

The authoritative outputs are
`results/candidate_architecture_review.tsv`,
`results/candidate_architecture_review_summary.tsv`, and
`results/candidate_architecture_review_manifest.json`.

## Important candidate-specific findings

- `GCA_963818845.1|CAWFTO010000003.1_420`,
  `GCF_000759345.1|NZ_JPZL01000003.1_113`, and
  `GCF_050269715.1|NZ_JBNNIW010000004.1_132` remain architecture-conflict
  holds despite pLDDT 85.32-92.24 and tree support 98-100.
- `GCF_006442545.2|NZ_CP083952.1_5059` remains a structure hold (pLDDT
  64.09; 49.0% residues below 70; longest low run 85; PF10503 span 2-100).
  The same genome contains a separate 384-aa ePhaZ hit (`..._5067`), so the
  202-aa record should be checked for gene-model truncation or incorrect ORF
  selection before further prediction.
- `GCA_030946875.1|JALYXY010000192.1_6` has pLDDT 88.04 but is marked
  `possible_N_truncation`; its Pfam hits begin near residue 6, so upstream ORF
  context must be verified.

## Next step

1. Retrieve the exact source contig/annotation for the two hold records and
   verify ORF boundaries and upstream/downstream alternatives.
2. Inspect the three dPHASCL2 PDBs against PF00756/InterPro coordinates,
   separating signal peptide/flexible N termini from the mature domain.
3. Review the 14 retained candidates, starting with the five P1 records and
   then the four longest low-confidence intermediate structures.
4. Bind Foldseek/TM-align only if structural comparison is required; any result
   remains architecture/homology evidence, not phenotype validation.

No candidate is promoted to `ePhaZ_curated_core`, no positive training set is
changed, and no formal GTDB rescan is authorized by this run.

