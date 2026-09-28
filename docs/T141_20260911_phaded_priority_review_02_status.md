# T141 PhaDED priority review and InterProScan status

Date: 2026-09-11  
Runs: `runs/20260911_phaded_priority_review_01/`, `runs/20260911_phaded_interpro_priority_02/`

## Completed

- Fixed duplicate `pfam_interpro_state` and `architecture_consistency` headers in the derived feature table. The original feature and Pfam evidence tables were not modified.
- Built the ordered manual review queue at `runs/20260911_phaded_priority_review_01/results/interpro_priority_v2/priority_review_queue.tsv`.
- Queue order is: all `conflicting`, all `ePhaZ_competition_review`, all curated-core records, then all `iPhaZ_tier1` records.
- Queue size is 71,671 unique accessions. The difference from the raw layer total is intentional accession de-duplication: a conflicting record keeps the earlier `conflicting` stage.
- Built the bounded InterProScan set at `runs/20260911_phaded_priority_review_01/results/interpro_priority_v2/interpro_priority_candidates.tsv`.
- InterProScan set contains all 83 conflicting and all 18 competition records, plus the top `priority_score` representative for each `layer + PhaDED superfamily + family` stratum, capped at 10 per stratum. Total: 244 priority candidates.

## InterProScan execution

- InterProScan: `5.78-109.0`, `TSV`, `--iprlookup`, server CPU `40`.
- The first run (`20260911_phaded_interpro_priority_01`) is retained as failure evidence: InterProScan rejected terminal `*` / non-IUPAC residues.
- Retry run submitted 242 normalized sequences. Normalization removed terminal `*` only; two records containing non-IUPAC `X` were excluded from the tool input and remain in the manual queue:
  `GCA_000245075.1|JH600030.1_43` and `GCA_000715855.1|KL571366.1_27`.
- Output: `runs/20260911_phaded_interpro_priority_02/results/interpro_priority.tsv` (242 accessions, 1,610 annotation rows), SHA-256 `6C9284B011911458DEC77F9527EA725D38AF2A55A43F0DD01E03BDF248E08BAC`.
- Execution manifest: `runs/20260911_phaded_interpro_priority_02/results/interpro_run_manifest.json`.

## Secondary evidence classification

The derived table is `runs/20260911_phaded_interpro_priority_02/results/interpro_secondary_evidence.tsv`.

| State | Count |
| --- | ---: |
| `explicit_pha_related_domain` | 129 |
| `generic_hydrolase_or_esterase` | 79 |
| `pha_synthase_like_alternative` | 34 |
| `not_submitted_tool_input_limit` | 2 |

For the first manual block, the 83 `conflicting` records split into 4 explicit extracellular PHA-related records, 67 generic hydrolase/esterase records, and 12 records with competing PHA synthase/PhaC-like evidence. The latter must not be called depolymerases. The 18 competition records all have explicit PHA/PHB-related InterPro/Pfam/NCBIfam evidence in this bounded confirmation.

## Automated conflicting review

The accession-level dossier is `runs/20260911_phaded_interpro_priority_02/results/conflicting_review_dossier.tsv`.

| Provisional label | Count | Interpretation |
| --- | ---: | --- |
| `A_retain_high_priority_architecture` | 4 | dPHASCL type 2 candidates with explicit PHA-related evidence; retain for structure/experimental shortlist, but do not erase the family conflict |
| `B_cross_fold_architecture_review` | 42 | dPHASCL type 1 HMM plus cross-fingerprint hydrolase domain; plausible homolog/fold cross-hit, not a resolved type call |
| `B_family_unresolved_architecture_review` | 17 | dPHASCL type 2-like hydrolase architecture, but family boundary remains unresolved |
| `E_hold_synthase_like_structure` | 12 | nPHAMCL conflict with PF06850 and PhaC/PHA-synthase-like annotations; require structure/phylogeny before any depolymerase interpretation |
| `E_hold_structure_assignment_review` | 8 | motif/domain non-co-localization or very low assignment gap; hold for structure/phylogeny |

The sequence-integrity correction is important: 82/83 records are `terminal_stop_only` (a normal FASTA translation terminator), and one is `valid`. The old `invalid_character` field must not be used as a biological downgrade.

## Interpretation boundary

InterProScan is a secondary domain-evidence layer. It does not overwrite PhaDED assignment or Pfam architecture state. `explicit_pha_related_domain` means database-domain support, not demonstrated PHB/PHA degradation. The two excluded sequences are tool-input limitations, not biological negatives. All current records remain candidate-only until sequence-specific experimental validation.

## Next actions

1. Review the 4 `A_retain_high_priority_architecture` records first.
2. Review the 12 `E_hold_synthase_like_structure` records as a PhaC-versus-PhaZ competition set; do not merge them with depolymerase candidates.
3. Review the 8 `E_hold_structure_assignment_review` records, then the 59 B-labeled cross-fold/family-boundary candidates.
4. Review the 18 competition records as a separate ePhaZ boundary set; do not merge them into curated-core counts.
5. Only after this review, decide whether any candidate warrants structure prediction, phylogenetic placement, or experimental validation. Do not modify `formal_scan_models.tsv` or start a GTDB-wide rescan from these results.
