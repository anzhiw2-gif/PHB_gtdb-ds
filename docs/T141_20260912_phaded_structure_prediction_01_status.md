# T141 PhaDED structure prediction status

Date: 2026-09-12  
Run: `runs/20260912_phaded_structure_prediction_01/`  
Status: `completed_candidate_only`

## Software and hardware

- Predictor: ColabFold `1.6.1`, AlphaFold implementation `alphafold-colabfold 2.3.13`.
- JAX `0.5.3` / JAX CUDA 12 plugin `0.5.3`; TensorFlow CPU `2.21.0`.
- T141 GPU: 2 x NVIDIA GeForce RTX 4090, 24,564 MiB each; driver `595.84`.
- Execution used GPU 0 only and set server CPU thread variables to 40 or fewer.
- GPU smoke test passed on one 396-aa candidate after the original terminal-stop input issue was corrected.

## Inputs and cleaning

- Frozen candidate input: `inputs/structure_candidates.faa`, 19 sequences, SHA-256 `fcc921c23537eb2ac8aed47ee9ce4c58a11b0cbefaa65b4584d68946d7cee50e`.
- Frozen reference input: `inputs/nearest_reference_controls.faa`, 19 sequences, SHA-256 `cb20afba7c9d1a853d58360f13dd9b61d251998ed66b32d272b2594560e423ff`.
- Original candidate sequences contained terminal `*` markers. The original files were preserved; predictor inputs remove only terminal `*`, reject internal stops and invalid residues, and record lengths in `*.cleaning.tsv`.
- Structure prediction parameters: AlphaFold2-PTM, `mmseqs2_uniref_env`, 5 models, 3 recycles, 1 seed, `save_all=true`, no Amber relaxation.

## Current result

- Candidates completed: 19/19.
- Candidate PDB files: 95 (five per candidate); per-sequence `log.txt` files contain pLDDT/pTM and ranking.
- Candidate logs show no `Error`, `Exception`, `Failed`, or `Traceback` for completed records.
- The first foreground process was stopped after Query 12 began so the remaining candidate jobs could be resumed deliberately; no generated structures were deleted. `resume_candidate_prediction.sh` resumed from the same dated deploy and skipped existing `.done.txt` records.
- Reference controls: 0/19; not started and intentionally marked `pending_not_started`.

The run remains candidate-only. Predicted structures are evidence for fold, domain boundaries, and catalytic-geometry review; they do not establish PHB/PHA degradation activity, change PhaDED assignments, repair the four-seed ePhaZ model, or authorize a formal rescan.

## First confidence triage

`results/structure_confidence_summary.tsv` contains one row per candidate. Using the best-ranked model mean pLDDT only as a structural-quality indicator, 9 candidates are `high_structure_confidence` (>=80), 9 are `intermediate_structure_confidence` (70-<80), and 1 is `low_structure_confidence` (<70). The next review should combine these values with PhaDED family placement, domain architecture and catalytic-residue evidence; pLDDT/pTM alone cannot distinguish PhaZ from PhaC or prove activity.
