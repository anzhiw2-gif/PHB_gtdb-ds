# T141 formal frozen scan 13 downstream status

Date: 2026-09-02
Run: `20260901_formal_frozen_scan_13`
Raw evidence: `/home/data/haoyu/PHB_gtdb-ds/runs/20260901_formal_frozen_scan_13/results/`
Derived outputs: `results/downstream_20260902/`

Downstream script: `deploy/20260901_formal_frozen_scan_13/scripts/summarize_formal_scan13.py`
SHA-256: `7961d031cef478da17ec462db6e9ec443829f08318a61c5aeb5bbb71aa9b7168`

## Processing completed

The downstream script accepted rows using the run-13 registry (`E<=1e-5`; `OH min_cov=0.6`), mapped proteins to genome accessions, retained model-level and core-family-level genome sets, and wrote:

- `model_summary.tsv`
- `genome_family.tsv`
- `genome_union_summary.tsv`
- `core_cooccurrence.tsv`
- `phylum_family.tsv`

Raw HMMER outputs were not modified.

## Accepted model counts

| model | accepted hits | unique genomes |
|---|---:|---:|
| `ePhaZ_curated_core` | 19,010 | 14,006 |
| `ePhaZ_broad_discovery` | 520,217 | 134,571 |
| `iPhaZ` | 63,035 | 38,504 |
| `OH` | 6,923 | 6,470 |
| `ArchPhaZ_patatin` | 306,216 | 134,741 |
| `ArchPhaZ_hydrolase` | 248,379 | 86,281 |
| `BdhA` | 4,243,143 | 192,640 |
| `PhaJ` | 1,226,439 | 141,285 |
| `PhaC` | 107,531 | 46,210 |
| `phasin` | 7 | 7 |

The raw accepted set contains 6,740,900 rows. OH decreases from 9,220 raw hits to 6,923 accepted rows after the registered coverage gate.

## Core-family union and co-occurrence

For the four core families (`ePhaZ` combining curated and broad, `iPhaZ`, `OH`, and `ArchPhaZ_hydrolase`), the registry-threshold union is **147,690 genomes**. The largest genome-level combinations are:

| family set | genomes |
|---|---:|
| `ArchPhaZ_hydrolase+ePhaZ` | 52,756 |
| `ePhaZ` | 44,694 |
| `ArchPhaZ_hydrolase+ePhaZ+iPhaZ` | 18,318 |
| `ePhaZ+iPhaZ` | 12,768 |
| `ArchPhaZ_hydrolase` | 9,238 |

## Interpretation boundary

These are **registry-threshold HMM hit counts**, not final tier1 or phenotype-confirmed counts. They are not directly comparable to the earlier Scheme A result of 44,814 genomes, which used additional validation/tier processing and a different input snapshot. `BdhA`, `PhaJ`, `PhaC`, `phasin`, and patatin are contextual or broad-fold families and must not be added to the core depolymerase union. HMM, domain, SignalP, neighborhood, and tree evidence indicate candidate homology or function potential; they do not prove PHB/MCL-PHA degradation.

## Next action

Review this derived table set, then decide whether a separate validation/tier-rescore run is scientifically needed. Do not overwrite the historical 44,814 result and do not register MCL candidate profiles without independent validation.
