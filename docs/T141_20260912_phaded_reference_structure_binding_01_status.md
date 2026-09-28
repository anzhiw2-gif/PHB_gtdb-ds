# PhaDED reference structure binding status

Run: `runs/20260912_phaded_reference_structure_binding_01/`  
Date: 2026-09-12  
Status: `partial_binding_complete; structural_comparison_pending`

## Scope

The 19 nearest-reference controls from `20260912_phaded_structure_review_01` were checked for accession-level UniProt mappings, AlphaFold DB models, and exact-UniProt RCSB PDB matches. The result is recorded in `results/reference_structure_binding.tsv` with file sizes and SHA-256 values for every downloaded artifact.

## Bound structures

Seven dPHASCL1 reference accessions map to UniProt and have AlphaFold DB models: `ABA81473.1 -> Q3IVG1`, `ABD70293.1 -> Q21VB0`, `ABE64568.1 -> Q1QGS9`, `ABM95604.1 -> A2SJ65`, `ABR61820.1 -> A6UDU0`, `EDT08140.1 -> B1G6S0`, and `EDQ04457.1 -> A0ABM9X4N4` (the last has no AlphaFold model). The five modeled entries have API `latestVersion=6`; their global confidence values range from 76.38 to 86.00.

The dPHASCL2 reference `O05527` maps to UniProt `O05527` (Delftia acidovorans, PhaZCac) and has `AF-O05527-F1`, AlphaFold DB latest version 6, model date 2025-08-01, and global metric 88.75. `EAU98424.1` remains a DED_hfam_71 GenBank-only reference without a unique UniProt/AlphaFold mapping.

## Pending references

`EAU98424.1`, `AAQ87263.1`, and `BAB49576.1` have no unique UniProt mapping in the current REST query. `EDQ04457.1` maps to `A0ABM9X4N4`, but AlphaFold DB reports no model. These records remain `pending`; the historical T141 ESMFold files are not accession-bound and are excluded.

RCSB exact-UniProt queries returned no experimental entries for the mapped accessions. This is a negative search result, not evidence that no homologous experimental structure exists.

## Decision and next step

The bound AlphaFold structures are sufficient to create a new, explicitly reference-bound structural comparison run for the 3 dPHASCL2 conflict candidates and the 4 priority dPHASCL1 intermediate candidates. The comparison must report model version, coordinate-file SHA-256, coverage/confidence, and candidate-only interpretation. It must not upgrade a candidate to a curated positive or phenotype call. Unresolved references stay pending and are not replaced by guessed structures.
