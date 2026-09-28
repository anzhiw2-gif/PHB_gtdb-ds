# T141 PhaDED structure and evidence integration review

Date: 2026-09-12  
Run: `runs/20260912_phaded_structure_evidence_review_01/`  
Status: `completed_candidate_only`

## Scope

This run joins the completed ColabFold candidate summary, the corrected v2
sequence-feature annotation, the reference-derived Pfam architecture evidence,
the secondary InterPro evidence, and the bounded competition phylogeny decisions.  It also parses the B-factor
(`pLDDT`) values in each candidate's best-ranked PDB to report low-confidence
fractions, longest low-confidence runs, and N/C-terminal 20-residue means.

The authoritative table is
`results/phaded_structure_evidence_review.tsv`; the input/output hashes are in
`results/structure_evidence_review_manifest.json` and
`input_contract.json`.

## Evidence-gated result

| Review tier | Count | Interpretation |
|---|---:|---|
| `high_confidence_candidate_only` | 5 | High mean pLDDT and no Pfam conflict in the bound evidence; retain as candidate-only comparison set. |
| `architecture_conflict_manual_review` | 3 | dPHASCL2-like phylogenetic placement and high pLDDT, but only PF00756 was observed against the type-2 fingerprint; inspect domain coordinates and structure before any conclusion. |
| `intermediate_structure_manual_review` | 9 | Phylogenetic placement is retained, but the fold contains substantial low-confidence segments; inspect PDB/PAE and domain boundaries. |
| `low_confidence_or_truncation_hold` | 2 | One sequence is a possible N-terminal truncation and one has low structural confidence; do not shortlist yet. |

All 19 sequences retain `candidate_only` phenotype claims.  No sequence was
added to `ePhaZ_curated_core`, any positive training set, or a formal scan
model.

All 19 records have `interpro_secondary_state=explicit_pha_related_domain`.
This supports a PHA-related annotation layer, but it does not override the
three Pfam architecture conflicts or establish a catalytic mechanism.
The dPHASCL1 records carry the `IPR010126;IPR029058;IPR050955` signature set;
the dPHASCL2 records carry `IPR000801;IPR029058` (with `IPR050955` in two of
the three).  The InterPro table is secondary evidence and remains bound to its
recorded run/version.

## Priority observations

1. The three explicit dPHASCL2-like candidates
   (`GCA_963818845.1|CAWFTO010000003.1_420`,
   `GCF_000759345.1|NZ_JPZL01000003.1_113`, and
   `GCF_050269715.1|NZ_JBNNIW010000004.1_132`) have high pLDDT (85.32,
   90.70, and 92.24) and strong tree support (98-100), but all are flagged by
   the reference-derived Pfam architecture as `conflicting` with PF00756.
   This is a domain-architecture review trigger, not evidence of PhaC or
   evidence of activity.
2. `GCA_963818845.1|CAWFTO010000003.1_420` has an N-terminal 20-residue mean
   pLDDT of 26.46 while its C-terminal mean is 86.87.  Treat the N terminus as
   a possible signal peptide/flexible region and compare the mature-domain
   boundary before interpreting the fold.
3. `GCF_006442545.2|NZ_CP083952.1_5059` is the hold-outlier: mean pLDDT 64.09,
   49.0% residues below 70, longest low run 85 residues, and C-terminal
   20-residue mean 37.36.  Re-check the gene model/contig and re-predict only
   after sequence integrity is confirmed.
4. `GCA_030946875.1|JALYXY010000192.1_6` has high pLDDT (88.04) but is marked
   `possible_N_truncation` in the corrected annotation; verify the upstream
   coding region before using it as a structural reference.
5. The nine intermediate candidates should be ordered by PDB low-confidence
   burden.  The largest burdens are `GCA_031366635..._14` (40.1%, 113-residue
   longest run), `GCF_009760915..._13` (35.5%, 107),
   `GCF_012927045..._126` (35.5%, 105), and
   `GCA_051861205..._261` (30.7%, 105).

The merged `merged_priority` field is stricter than the original PhaDED
assignment: 3 records are P0 architecture holds, 5 are P1 high-confidence
dPHASCL1 candidates, 9 are P2 intermediate candidates, 1 is a P2 truncation
hold, and 1 is a P3 low-structure hold.  These labels are triage only and
never curate a positive training sequence.

The current P1 set is `GCA_038324515.1|JBBQIB010000018.1_12`,
`GCF_000018545.1|NC_012586.1_319`,
`GCF_000372525.1|NZ_AQUR01000076.1_30`,
`GCF_001651855.1|NZ_LNQC01000049.1_134`, and
`GCF_002944405.1|NZ_CP024310.1_365`.  Their merged Pfam spans cover about
47-56% of each cleaned sequence and all have tree support at least 90.

## Next authorized analysis step

Perform a bounded, candidate-only structure/domain review in this order:

1. Inspect the three PF00756 architecture conflicts against their Pfam
   coordinates and predicted mature-domain boundaries.
2. Verify the possible N-terminal truncation and the low-confidence hold-outlier
   against the source contig/gene model; do not use either as a positive control.
3. Review the four highest low-confidence-burden intermediate structures with
   PAE and domain boundaries.
4. If a structure-comparison tool is needed, bind its exact version and dated
   deploy first.  Foldseek/TM-align results would remain homology/architecture
   evidence only.
5. Only after this review decide which candidates merit experimental expression
   or biochemical testing.  Do not repair the four-seed ePhaZ model or launch a
   GTDB-wide rescan based on these structures alone.
