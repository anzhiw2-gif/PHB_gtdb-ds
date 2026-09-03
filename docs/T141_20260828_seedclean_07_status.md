# T141 20260828 seed-clean run status

- Run: `/home/data/haoyu/PHB_gtdb-ds/runs/20260828_seedclean_07/`
- Deploy: `/home/data/haoyu/PHB_gtdb-ds/deploy/20260828_seedclean_07/`
- Scope: recovered three missing `ArchPhaZ_patatin` accessions, removed 18 known non-depolymerase controls, rebuilt nine HMMs, and calibrated `ePhaZ`, `iPhaZ`, and `OH`.

## Evidence

- Original server `data/seeds/v2/ArchPhaZ_patatin.faa`: 2 bytes, SHA-256 `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`; preserved and not overwritten.
- Forensic c90 source: 44,277 bytes, SHA-256 `ad7bfc8510cfc7cdde7a01c63473189051e15173070f4f3c566e9546fadb6d1f`; lacked `A0AAV3SSL3`, `A0ABD5M1Y4`, and `A0ABV4II51`.
- UniProt recovery endpoints returned HTTP 200 on 2026-08-28; normalized FASTA records were added only to the new run.
- Cleaned output checksum verification returned no errors.
- Cleaned family record counts: `ePhaZ` 4458; `iPhaZ` 152; `OH` 708; `BdhA` 5890; `ArchPhaZ_patatin` 106; `ArchPhaZ_hydrolase` 12; `PhaJ` 1071; `PhaC` 47; `phasin` 13898; negative 18.
- All 18 excluded accessions are absent from training family FASTA and retained only in `negative.faa`.

## HMM outputs

All nine HMMs passed HMMER 3.4 `hmmstat`. Full SHA-256 records are in the run directory.

## Calibration decision

- `iPhaZ`: 17/17 positives, 0 false positives at tested thresholds.
- `OH`: 27/28 positives, 0 false positives at `1e-15`.
- `ePhaZ`: 4/10 positives, 0 false positives under default HMMER filters at `E=1e-2`. A read-only `--max -E 100` diagnostic recovered all six additional positives as weak/partial hits; two had `E` below `1e-2`, while four remained weaker than that threshold. This indicates both acceleration-filter sensitivity and incomplete profile coverage, not biological absence.

## Gate

`planned_not_run`: 268G full-proteome HMMER scan and downstream candidate analysis remain paused pending `ePhaZ` seed/model review. Existing historical results remain untouched.
