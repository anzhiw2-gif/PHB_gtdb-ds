# T141 20260906 overall PHA, PHB focus, and archaeal audit status

## Status

- Run: `20260906_overall_pha_phb_archaea_audit_03`
- Status: `completed_candidate_only`
- Archaeal source: `<SERVER_USER>@<SERVER_HOST>:${PHB_REMOTE_ROOT}/PHB_gtdb-ds/deploy/20260904_archaea_phazh1_candidate_01/`
- Formal GTDB scan: not started.
- `pipeline/config/formal_scan_models.tsv`: unchanged, SHA-256 `8a1c05085f882daad9fddb9cf7616def74438d3a7bf55f219e0b2046adc5b7a3`.

## Overall PHA evidence layers

The accepted output contains 109,087 deduplicated bacterial `ePhaZ`/`iPhaZ` candidate proteins. The original layer is preserved in `overall_pha_layer`; `pha_subtype_evidence` records only the combined HMM/SignalP/tier evidence.

| Evidence layer | Proteins | Interpretation |
|---|---:|---|
| `ePhaZ_SCL_secreted_evidence` | 38,674 | Curated-core extracellular candidate with secreted evidence |
| `ePhaZ_tier2_review_evidence` | 31,990 | ePhaZ review layer, not core evidence |
| `ePhaZ_competition_review_evidence` | 18 | Ambiguous competing evidence; retain for review |
| `iPhaZ_intracellular_evidence` | 32,926 | Intracellular candidate evidence |
| `iPhaZ_tier2_review_evidence` | 5,479 | Intracellular review layer, not core evidence |

## PHB focus

All 109,087 records are currently `PHB_candidate_review`. The frozen dual-profile panel uses reference accessions and does not independently validate any GTDB hit accession; therefore none are upgraded to PHB-specific evidence. Historical `MCL_like` is now handled as a competing MCL signal where accession mapping is available. This is intentionally conservative and does not make a phenotype claim.

## Full archaeal PhaZh1-like audit

All 2,307 candidate records were retained for traceability. The requested 197 `PhaZh1_like_high` records are fully represented, from 173 genomes.

| Archaeal layer | Records | Audit treatment |
|---|---:|---|
| `PhaZh1_like_high` | 197 | `review_required`, not a phenotype call |
| `PhaZh1_like_review` | 7 | `review_required` |
| `archaeal_patatin_exploratory` | 2,027 | exploratory only |
| `non_PhaZh1_patatin` | 76 | excluded from PhaZh1-like layer |

For the 197 high records: 195 are `Halobacteriota; Halobacteria`, 2 are `Thermoproteota; Nitrososphaeria`; all have matched FASTA and taxonomy; 188 contain an N-terminal `GxSxG` motif in the first 80 residues and 9 do not; neighborhood markers are `BdhA` (192), `BdhA;PhaJ` (4), or `PhaJ` (1). The independent domain annotation is pending for every record. These checks support prioritization only and do not establish archaeal PHB depolymerase activity.

## Provenance and outputs

- Input contract: `runs/20260906_overall_pha_phb_archaea_audit_03/input_contract.json`
- Overall PHA and PHB-focus table: `runs/20260906_overall_pha_phb_archaea_audit_03/results/overall_pha_phb.tsv`
- Archaeal full audit: `runs/20260906_overall_pha_phb_archaea_audit_03/results/archaeal_phazh1_full_audit.tsv`
- Result manifest: `runs/20260906_overall_pha_phb_archaea_audit_03/results/analysis_manifest.json`
- Earlier `audit_01` and `audit_02` runs are retained as intermediate evidence; `audit_03` is the accepted output because it preserves review strata and records the motif audit field.
