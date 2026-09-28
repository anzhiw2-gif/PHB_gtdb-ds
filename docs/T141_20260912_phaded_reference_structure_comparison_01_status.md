# PhaDED reference-bound structure comparison status

Run: `runs/20260912_phaded_reference_structure_comparison_01/`  
Date: 2026-09-12  
Status: `completed_candidate_only`

Seven exact candidate/reference pairs have been prepared. All seven candidate PDBs are the rank-001 ColabFold outputs from `20260912_phaded_structure_prediction_01`; all seven reference PDBs are accession-bound AlphaFold DB files from `20260912_phaded_reference_structure_binding_01`. Pair-level paths, byte sizes, and SHA-256 values are recorded in `inputs/comparison_pairs.tsv`.

Foldseek was installed as a dated user-level binary at `${PHB_REMOTE_ROOT}/tools/foldseek_20260912/` and executed from the dated deploy bundle with 4 threads. The executable version identifier is `463739e0014a1549a527de589102cde98f802f37`; executable and distribution hashes are recorded in `results/tool_manifest.json`. TM-align remains unavailable and was not needed for this bounded Foldseek comparison.

All seven comparisons completed. The three dPHASCL2 pairs have alignment TM-scores 0.825-0.860, while the four dPHASCL1 pairs have alignment TM-scores 0.621-0.697. These values support fold-level structural similarity, but are not a calibrated family classifier and do not override Pfam architecture conflicts or missing reference fingerprints.

The next step is to merge these scores with the existing Pfam/InterPro, structure-confidence, and phylogeny evidence. The three dPHASCL2 candidates remain `P0_architecture_conflict_hold`; the four dPHASCL1 candidates remain candidate-only (P1/P2 according to the prior confidence review). Do not alter curated HMMs, positive sets, or phenotype labels.
