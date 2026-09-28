# T141 PhaDED candidate mapping status

Date: 2026-09-11  
Run: `runs/20260909_phaded_candidate_mapping_01/`  
Scope: candidate-only ePhaZ/iPhaZ mapping; no formal registry change and no GTDB-wide rescan.

## Completed

- Audited PhaDED reference profile output: 46 profile records, including 38 family records and 8 superfamily records; 9 trained HMMs and 37 `reference_only` records. All trained alignment/HMM artifacts exist and match `profile_manifest.tsv` SHA-256 values.
- Bound 109,087 unique candidate accessions from `runs/20260905_ephaz_iphaz_layering_01/results/layers_retry_02/protein_layers.tsv` to a deduplicated `candidate_union.faa`; no candidate sequence was missing.
- Ran the 9 trained profiles sequentially on T141 with HMMER 3.4 and `--cpu 40`. The raw `tblout` files, `phaded_profile_scores.tsv`, and `hmmsearch_provenance.json` are retained.
- Declared the score rule before mapping: `min_score_gap=1.0 bits`; ties and smaller gaps remain ambiguous.
- Reconciled every candidate accession exactly once into `phaded_assignment.tsv` and preserved the original layer fields.
- Added feature evidence in `phaded_feature_annotation.tsv` and a separate merged table `phaded_assignment_evidence.tsv`.

## Mapping counts

| Assignment status | Count |
| --- | ---: |
| `assigned` | 67,749 |
| `ambiguous_family` | 9,664 |
| `ambiguous_superfamily` | 42 |
| `unassigned_PhaDED_like` | 31,632 |

The 31,632 `unassigned_PhaDED_like` records mean that no score was reported by the current 9-profile trained set at the declared reporting setting. They are not biological negatives. In particular, 28,573 iPhaZ tier-1 records are in this category because the current trained set contains no calibrated i-nPHASCL or periplasmic/PhaZ7 profile.

## Evidence boundary

The current feature run has no independently bound Pfam/InterPro or structure-domain input, so all 109,087 feature rows remain `feature_evidence_status=pending` and `architecture_consistency=pending`. Motif observations are retained as support fields; they are not activity calls. All PhaDED labels remain candidate homology/architecture evidence and do not establish PHB degradation phenotype.

## Next authorized analytical steps

1. Bind an independently versioned domain/structure annotation input and rerun feature annotation without overwriting the mapping table.
2. Review the priority strata in the declared order, retaining ambiguous and unassigned records.
3. Execute the separate four-seed ePhaZ calibration study; do not edit `ePhaZ_curated_core` or `formal_scan_models.tsv` during this run.
4. Only after those evidence gates are accepted, prepare a separate decision record for any future registry change or GTDB-wide rescan.

