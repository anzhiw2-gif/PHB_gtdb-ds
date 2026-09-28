# OH / ArchPhaZ cross-phylum audit status

Runs: `20260903_oh_arch_crossphylum_audit_03`  
Parent tier run: `20260902_formal_scan13_tier_processing_02`

## Recovery and provenance

- Two OH family FASTA files had been truncated to zero bytes by an erroneous read-only inspection command. The original zero-byte states were preserved under `damage_evidence/20260903_oh_fasta_truncation/` with SHA-256 `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`.
- `OH_validated.faa` was restored from the same run's `data/screen/tiers/OH_tier1.faa`; source and restored file SHA-256 both equal `5d9af7099312172d812481372df4b66f22c16ca9d5ad869c72d9cdbb2e94a25b`.
- `OH.faa` remains unrecovered because tier1/tier2 outputs cannot reconstruct the pre-validation input exactly. It remains zero-byte evidence and must not be used as a validated sequence source.

## Cross-phylum audit

Run `20260903_oh_arch_crossphylum_audit_03` completed as a deterministic read-only sample (seed `20260903`, maximum 5 genomes per family x coverage tier x phylum stratum): 15,915 candidate genomes represented, 559 sampled records, and 0 missing sampled sequences.

Full tier1 genome counts by domain and coverage tier:

| family | tier | Archaea | Bacteria |
|---|---:|---:|---:|
| OH | `>=0.8` | 0 | 3,061 |
| OH | `0.6-0.8` | 0 | 385 |
| ArchPhaZ_hydrolase | `>=0.8` | 0 | 499 |
| ArchPhaZ_hydrolase | `0.6-0.8` | 9 | 2,704 |
| ArchPhaZ_hydrolase | `<0.6` | 270 | 8,987 |

The result does not support calling the ArchPhaZ_hydrolase high-coverage layer an archaeal expansion: all 499 `>=0.8` genomes are bacterial. Archaeal hits are concentrated in lower-coverage layers and remain exploratory/review candidates. OH candidates in this tier table are bacterial; no archaeal OH tier1 genomes were observed.

## Other major depolymerase families

The same tier-processing run contains non-empty validated and tier1 FASTA outputs for `ePhaZ`, `iPhaZ`, and `ArchPhaZ_hydrolase`:

| family | validated FASTA bytes | tier1 FASTA bytes | tier1 genomes |
|---|---:|---:|---:|
| ePhaZ | 7,498,733 | 2,417,069 | 5,080 |
| iPhaZ | 21,629,044 | 14,542,673 | 25,564 |
| ArchPhaZ_hydrolase | 21,348,637 | 6,235,597 | 12,469 |

No empty-file or missing-output problem was found for these three families. Their HMM/tier evidence remains candidate homology/function potential, not experimental PHB degradation confirmation.

## Decision

- Keep OH formal rule `E<=1e-5, min_cov=0.6`; report `>=0.8` as high-confidence only.
- Keep ArchPhaZ_hydrolase split into `>=0.8` high-confidence, `0.6-0.8` review, and `<0.6` exploratory. Do not describe the 12,469 tier1 genomes as archaeal PHB depolymerases.
- Do not rebuild or rescan ePhaZ/iPhaZ/ArchPhaZ_hydrolase based on this audit. The next useful work is targeted manual/functional review of high-coverage bacterial ArchPhaZ candidates and selected archaeal lower-coverage candidates.
