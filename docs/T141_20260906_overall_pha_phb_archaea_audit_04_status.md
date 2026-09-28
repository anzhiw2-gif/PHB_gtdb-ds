# T141 20260906 overall PHA, PHB focus, and archaeal audit

## Status

- Run: `20260906_overall_pha_phb_archaea_audit_04`
- Status: `completed_candidate_only`
- Archaeal inputs were copied from the dated server deploy `20260904_archaea_phazh1_candidate_01` via `<SERVER_USER>@<SERVER_HOST>`.
- Formal GTDB scan was not started and `formal_scan_models.tsv` was not modified (SHA-256 `8a1c05085f882daad9fddb9cf7616def74438d3a7bf55f219e0b2046adc5b7a3`).

## Overall PHA and conservative PHB focus

`overall_pha_phb.tsv` contains 109,087 deduplicated `ePhaZ`/`iPhaZ` candidate proteins. The original layer is retained in `overall_pha_layer`; `pha_subtype_evidence` preserves evidence boundaries:

| Subtype evidence | Proteins |
|---|---:|
| `ePhaZ_core_nonsecreted_evidence` | 16,818 |
| `ePhaZ_SCL_secreted_evidence` | 21,856 |
| `ePhaZ_tier2_review_evidence` | 31,990 |
| `ePhaZ_competition_review_evidence` | 18 |
| `iPhaZ_intracellular_evidence` | 32,926 |
| `iPhaZ_tier2_review_evidence` | 5,479 |

All 109,087 remain `PHB_candidate_review`: the reference dual-profile panel has no accession-level match to these GTDB protein identifiers, so no PHB-specific upgrade is justified. This is a candidate classification, not an experimentally verified phenotype.

## Archaeal PhaZh1-like full audit

All 2,307 candidate rows are retained; the requested 197 `PhaZh1_like_high` rows represent 173 genomes.

| Layer | Records | Treatment |
|---|---:|---|
| `PhaZh1_like_high` | 197 | review required; candidate-only |
| `PhaZh1_like_review` | 7 | review required |
| `archaeal_patatin_exploratory` | 2,027 | exploratory only |
| `non_PhaZh1_patatin` | 76 | excluded |

All 2,307 have matching FASTA and taxonomy rows. Across all rows, 2,118 have an N-terminal `GxSxG` motif in the first 80 residues and 189 do not. In the 197 high rows, 188 are motif-present and 9 motif-absent; neighborhood markers are `BdhA` (192), `BdhA;PhaJ` (4), or `PhaJ` (1). Independent domain annotation remains pending for every row.

## Outputs

- Contract: `runs/20260906_overall_pha_phb_archaea_audit_04/input_contract.json`
- PHA/PHB table: `runs/20260906_overall_pha_phb_archaea_audit_04/results/overall_pha_phb.tsv`
- Archaeal audit: `runs/20260906_overall_pha_phb_archaea_audit_04/results/archaeal_phazh1_full_audit.tsv`
- Manifest: `runs/20260906_overall_pha_phb_archaea_audit_04/results/analysis_manifest.json`

Previous `audit_01` to `audit_03` directories are preserved as intermediate/error evidence; this run is the accepted output after fixing label compatibility, review-layer preservation, and `nonsecreted` matching.
