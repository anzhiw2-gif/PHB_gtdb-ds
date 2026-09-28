# PhaDED hold-source retrieval and CDS audit

Date: 2026-09-12  
Run: `runs/20260912_phaded_hold_source_retrieval_01/`  
Status: `completed_candidate_only_source_bound`

## Sources

Authoritative records were retrieved from NCBI nuccore and NCBI Datasets; the
JALYXY contig was cross-checked against the ENA Browser API. Exact endpoints,
accession versions, file sizes and SHA-256 values are in
`results/source_manifest.json`.

- `GCA_030946875.1`, `ASM3094687v1`, Contig, NCBI Datasets release 2023-08-29;
  PGAP 6.1 annotation.
- `GCF_006442545.2`, `ASM644254v2`, Complete Genome, NCBI Datasets release
  2021-11-01; RefSeq annotation `GCF_006442545.2-RS_2025_02_14`.
- `JALYXY010000192.1` ENA FASTA matched the NCBI nuccore FASTA sequence
  exactly (4659 bp target record among the WGS set).

## CDS findings

### `GCA_030946875.1|JALYXY010000192.1_6`

Exact translation match: `M3Z31_14030`, protein `MDQ6620787.1`, product
`PHB depolymerase family esterase`, CDS `4030..4659` on the minus strand.
The GenBank location is `[4029:>4659](-)`, so the CDS is C-terminally partial;
the earlier `possible_N_truncation` label is incorrect for this source. The
five upstream CDSs end at 577, 1437, 2079, 3136 and 3847 on the same strand,
with no upstream continuation of this CDS. Keep a gene-model/structure hold,
but reclassify the issue as `partial_C_terminal_source_record`.

### `GCF_006442545.2|NZ_CP083952.1_5059`

Exact translation match: `FJ970_RS33950`, protein `WP_415752051.1`, product
`PHB depolymerase family esterase`, CDS `5164677..5165285` on the plus strand;
the location is complete. This is a genuine 202-aa annotated ORF, not a
translation artifact. Its low structure confidence and narrow PF10503 span
remain unresolved biological/model evidence.

### `GCF_006442545.2|NZ_CP083952.1_5067`

Exact independent match: `FJ970_RS25565`, protein `WP_140762596.1`, product
`extracellular catalytic domain type 1 short-chain-length polyhydroxyalkanoate
depolymerase`, CDS `5172176..5173327` on the minus strand. It is 383 aa and
is separated from `5059` by 6890 bp with no overlap. Do not merge the ORFs.

## Decision

The source-provenance gap is closed, but both holds remain candidate-only:

- `JALYXY..._6`: retain hold; correct reason to C-terminally partial source
  record. Do not invent an upstream CDS or reconstruct the sequence.
- `NZ_CP083952.1_5059`: retain low-structure/model hold; it is independent of
  `5067`. A new structure prediction is optional only after explicit approval,
  using the verified CDS translation.

No curated model, positive training set, PhaDED family assignment, or formal
GTDB scan was changed. Homology, domain, structure and annotation evidence do
not establish PHB degradation phenotype.

