# OH / ArchPhaZ hydrolase / auxiliary-family audit

Run: `20260902_oh_arch_aux_audit_01`
Parent: `20260901_formal_frozen_scan_13`
Scope: read-only, candidate-only; no production registry or historical result was changed.

## Decisions

### OH

Keep the registered `E<=1e-5, min_cov=0.6` rule. The previous calibration showed nylon-hydrolase controls at coverage `0.31-0.53`, while true OH controls were approximately `0.99`; the `0.6` gate is therefore evidence-based. Run 13 strict tier1 contains 3,570 sequences from 3,446 genomes. An optional high-confidence report can use the existing raw coverage `>=0.8` subset (3,339 sequences, 3,061 genomes), while the `0.6-0.8` band remains reviewable rather than deleted.

The increase from the historical 1,410 genomes to 3,446 should be treated as a cohort/model-processing difference and checked with cross-phylum sampling. It is not evidence of a biological expansion by itself.

### ArchPhaZ_hydrolase

Do not report all 12,469 strict tier1 genomes as high-confidence archaeal PHB depolymerase candidates. The raw accepted HMM hits are dominated by partial matches: 203,809/248,379 have coverage `<0.4`, and the median raw coverage is `0.2029`. Among the current strict tier1 sequences, only 3,431/14,571 meet coverage `>=0.6`, and only 508 meet `>=0.8`; the corresponding genome counts are 3,166 and 497.

Recommended candidate-only presentation:

- `high_confidence`: existing sequence validation + `E<=1e-20` + raw HMM coverage `>=0.8` (497 genomes);
- `review`: coverage `0.6-0.8` (the remainder of the 3,166-genome `>=0.6` set);
- `exploratory`: coverage `<0.6`, retained for traceability but excluded from the hydrolase high-confidence count.

This is a specificity/fragmentation audit, not proof that lower-coverage homologs are biologically inactive.

## Auxiliary families

- `PhaJ`: 141,285 genomes (70.7% of the 199,923-genome reference cohort). It is a useful pathway-context marker, especially for archaeal mobilization, but not a direct depolymerase call.
- `BdhA`: 192,640 genomes (96.4%). Its near-universal prevalence makes it non-specific in this scan; keep it as background metabolism only.
- `PhaC`: 46,210 genomes (23.1%). It can support PHA-island context, but its presence alone indicates synthesis potential, not degradation. It should be used as a co-occurrence feature, not a positive label.
- `phasin`: only 7 genomes were accepted despite a 7,671-sequence HMM seed profile. This is a model/sensitivity warning, not evidence that phasin is biologically absent. Keep it optional and do not use it as a required negative or positive gate without a dedicated control panel.

In the raw genome sets, `OH` co-occurs with PhaJ in 98.79% and PhaC in 67.51%; `ArchPhaZ_hydrolase` co-occurs with PhaJ in 84.27% and PhaC in 29.69%. In the strict tier1 sets, the corresponding PhaJ/PhaC rates are 99.36%/57.20% for OH and 94.31%/35.75% for ArchPhaZ hydrolase. These are contextual associations, not independent validation of enzymatic activity.

## Statistical interpretation

Counts are repeated detections on overlapping genomes and use changed model/filter definitions. Do not apply an independent-sample chi-square test to claim improvement. Report absolute changes, accession overlap, coverage strata, and phylum-stratified rates. All HMM, domain, sequence-feature, and co-occurrence evidence remains candidate homology/function potential, not experimentally verified PHB degradation.
