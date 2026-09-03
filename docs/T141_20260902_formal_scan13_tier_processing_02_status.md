# T141 run 13 tier processing status

Run: `20260902_formal_scan13_tier_processing_02`
Parent: `20260901_formal_frozen_scan_13`
Deploy: `${PHB_REMOTE_ROOT}/PHB_gtdb-ds/deploy/20260902_formal_scan13_tier_processing_02/`

## Completion

- Status: `completed`
- Input contract: `verified`
- Extraction: 308,862 strict-core proteins prepared; 100 shard extraction jobs completed with GNU Parallel (`-j 20`)
- Validation: 308,862 rows reviewed by `08_validate.py`
- HMM rescore: `08c_tier_rescore.py --cpu 60`
- Raw run 13 outputs were read-only inputs and were not modified.

## Strict tier results

| family | tier1 sequences | tier1 genomes | tier2 sequences |
|---|---:|---:|---:|
| ePhaZ (curated core only) | 5,646 | 5,080 | 10,202 |
| iPhaZ | 32,226 | 25,564 | 41,439 |
| OH | 3,570 | 3,446 | 3,585 |
| ArchPhaZ_hydrolase | 14,571 | 12,469 | 37,752 |

The strict four-family union is **38,741 genomes**, with 6,578 genomes containing at least two of these families. Taxonomic and ecological tables are in `runs/20260902_formal_scan13_tier_processing_02/results/tables/`.

## Broad and contextual layers

`ePhaZ_broad_discovery` remains separate at 520,217 registry-threshold proteins and is not promoted to strict tier1. `BdhA`, `PhaJ`, `PhaC`, `phasin`, and patatin remain contextual or broad-fold families and are excluded from the strict depolymerase union.

## Comparability warning

The 38,741-genome strict union is not directly comparable to the historical Scheme A value of 44,814. The current run uses the frozen split registry, isolates broad ePhaZ, and applies the current sequence-feature validation path. Both values must be retained with their input, model, and processing definitions.

## Evidence boundary

Tier1 is a high-confidence computational candidate layer, not experimental phenotype confirmation. HMM, domain, sequence-feature, SignalP, neighborhood, and tree evidence do not by themselves prove PHB or MCL-PHA degradation.
